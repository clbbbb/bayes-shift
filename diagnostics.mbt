///|
pub struct TraceDiagnostics {
  steps : Int
  mean_change_probability : Double
  maximum_change_probability : Double
  maximum_change_step : Int
  shortest_map_run_length : Int
  longest_map_run_length : Int
  mean_posterior_mean : Double
}

///|
pub fn diagnose_trace(
  updates : Array[Update],
) -> Result[TraceDiagnostics, DetectError] {
  if updates.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  let mut total_probability = 0.0
  let mut total_mean = 0.0
  let mut maximum_probability = updates[0].change_probability
  let mut maximum_step = updates[0].step
  let mut shortest_run = updates[0].most_likely_run_length
  let mut longest_run = updates[0].most_likely_run_length
  for update in updates {
    total_probability = total_probability + update.change_probability
    total_mean = total_mean + update.posterior_mean
    if update.change_probability > maximum_probability {
      maximum_probability = update.change_probability
      maximum_step = update.step
    }
    if update.most_likely_run_length < shortest_run {
      shortest_run = update.most_likely_run_length
    }
    if update.most_likely_run_length > longest_run {
      longest_run = update.most_likely_run_length
    }
  }
  Ok({
    steps: updates.length(),
    mean_change_probability: total_probability / updates.length().to_double(),
    maximum_change_probability: maximum_probability,
    maximum_change_step: maximum_step,
    shortest_map_run_length: shortest_run,
    longest_map_run_length: longest_run,
    mean_posterior_mean: total_mean / updates.length().to_double(),
  })
}

///|
pub fn TraceDiagnostics::is_spiky(
  self : TraceDiagnostics,
  threshold : Double,
) -> Result[Bool, DetectError] {
  if !finite(threshold) || threshold < 0.0 || threshold > 1.0 {
    Err(InvalidHazard)
  } else {
    Ok(self.maximum_change_probability >= threshold)
  }
}

///|
pub fn TraceDiagnostics::probability_ratio(self : TraceDiagnostics) -> Double {
  if self.mean_change_probability == 0.0 {
    0.0
  } else {
    self.maximum_change_probability / self.mean_change_probability
  }
}

///|
pub struct ProbabilityHistogram {
  edges : Array[Double]
  counts : Array[Int]
}

///|
pub fn ProbabilityHistogram::edges(
  self : ProbabilityHistogram,
) -> Array[Double] {
  self.edges.copy()
}

///|
pub fn ProbabilityHistogram::counts(self : ProbabilityHistogram) -> Array[Int] {
  self.counts.copy()
}

///|
pub fn ProbabilityHistogram::total(self : ProbabilityHistogram) -> Int {
  let mut total = 0
  for count in self.counts {
    total = total + count
  }
  total
}

///|
pub fn probability_histogram(
  updates : Array[Update],
  bins : Int,
) -> Result[ProbabilityHistogram, DetectError] {
  if updates.length() == 0 || bins <= 0 {
    return Err(InvalidRunLengthLimit)
  }
  let edges = []
  let counts = Array::make(bins, 0)
  for index in 0..<=bins {
    edges.push(index.to_double() / bins.to_double())
  }
  for update in updates {
    let mut bucket = (update.change_probability * bins.to_double()).to_int()
    if bucket >= bins {
      bucket = bins - 1
    }
    if bucket < 0 {
      bucket = 0
    }
    counts[bucket] = counts[bucket] + 1
  }
  Ok({ edges, counts, })
}

///|
pub fn format_probability_histogram(
  histogram : ProbabilityHistogram,
  width : Int,
) -> Result[String, DetectError] {
  if width <= 0 {
    return Err(InvalidRunLengthLimit)
  }
  let mut maximum = 0
  for count in histogram.counts {
    if count > maximum {
      maximum = count
    }
  }
  if maximum == 0 {
    return Ok("")
  }
  let mut text = ""
  for index in 0..<histogram.counts.length() {
    let filled = histogram.counts[index] * width / maximum
    let mut bar = ""
    for position in 0..<width {
      if position < filled {
        bar = bar + "#"
      } else {
        bar = bar + "."
      }
    }
    text = text +
      histogram.edges[index].to_string() +
      "," +
      histogram.edges[index + 1].to_string() +
      "," +
      histogram.counts[index].to_string() +
      "," +
      bar
    if index + 1 < histogram.counts.length() {
      text = text + "\n"
    }
  }
  Ok(text)
}

///|
pub struct DetectedSegment {
  start_step : Int
  end_step : Int
  size : Int
  mean_observation : Double
  mean_change_probability : Double
}

///|
fn make_segment(
  updates : Array[Update],
  start : Int,
  finish : Int,
) -> DetectedSegment {
  let mut total_observation = 0.0
  let mut total_probability = 0.0
  for index in start..<finish {
    total_observation = total_observation + updates[index].observation
    total_probability = total_probability + updates[index].change_probability
  }
  let size = finish - start
  {
    start_step: updates[start].step,
    end_step: updates[finish - 1].step,
    size,
    mean_observation: total_observation / size.to_double(),
    mean_change_probability: total_probability / size.to_double(),
  }
}

///|
pub fn detected_segments(
  updates : Array[Update],
  boundaries : Array[Boundary],
) -> Result[Array[DetectedSegment], DetectError] {
  if updates.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  let result = []
  let mut start = 0
  for boundary in boundaries {
    if boundary.step <= start + 1 || boundary.step > updates.length() + 1 {
      return Err(InvalidRunLengthLimit)
    }
    result.push(make_segment(updates, start, boundary.step - 1))
    start = boundary.step - 1
  }
  if start < updates.length() {
    result.push(make_segment(updates, start, updates.length()))
  }
  Ok(result)
}

///|
pub fn segment_boundaries_from_updates(
  updates : Array[Update],
  threshold : Double,
  cooldown : Int,
) -> Result[Array[Boundary], DetectError] {
  if cooldown < 0 {
    return Err(InvalidRunLengthLimit)
  }
  let candidates = match predicted_boundaries(updates, threshold) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let result = []
  let mut last = -cooldown - 1
  for candidate in candidates {
    if candidate.step - last > cooldown {
      result.push(candidate)
      last = candidate.step
    }
  }
  Ok(result)
}

///|
pub struct BoundaryDelay {
  truth_step : Int
  predicted_step : Int?
  delay : Int?
}

///|
pub fn boundary_delays(
  truth : Array[Boundary],
  predicted : Array[Boundary],
  tolerance : Int,
) -> Result[Array[BoundaryDelay], DetectError] {
  if tolerance < 0 {
    return Err(InvalidRunLengthLimit)
  }
  let used = Array::make(predicted.length(), false)
  let result = []
  for expected in truth {
    let mut match_index : Int? = None
    let mut best_distance = tolerance + 1
    for index in 0..<predicted.length() {
      let difference = predicted[index].step - expected.step
      let distance = if difference < 0 { -difference } else { difference }
      if !used[index] && distance <= tolerance && distance < best_distance {
        match_index = Some(index)
        best_distance = distance
      }
    }
    match match_index {
      None =>
        result.push({
          truth_step: expected.step,
          predicted_step: None,
          delay: None,
        })
      Some(index) => {
        used[index] = true
        result.push({
          truth_step: expected.step,
          predicted_step: Some(predicted[index].step),
          delay: Some(predicted[index].step - expected.step),
        })
      }
    }
  }
  Ok(result)
}

///|
pub fn mean_absolute_delay(delays : Array[BoundaryDelay]) -> Int? {
  let mut count = 0
  let mut total = 0
  for item in delays {
    match item.delay {
      None => ()
      Some(value) => {
        if value < 0 {
          total = total - value
        } else {
          total = total + value
        }
        count = count + 1
      }
    }
  }
  if count == 0 {
    None
  } else {
    Some(total / count)
  }
}

///|
pub fn format_diagnostics(value : TraceDiagnostics) -> String {
  "steps=" +
  value.steps.to_string() +
  ",mean_change_probability=" +
  value.mean_change_probability.to_string() +
  ",maximum_change_probability=" +
  value.maximum_change_probability.to_string() +
  ",maximum_change_step=" +
  value.maximum_change_step.to_string() +
  ",shortest_map_run_length=" +
  value.shortest_map_run_length.to_string() +
  ",longest_map_run_length=" +
  value.longest_map_run_length.to_string()
}
