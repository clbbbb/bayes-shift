///|
/// Aggregate diagnostics over a replay trace without changing detector state.
pub struct TraceSummary {
  points : Int
  event_count : Int
  mean_change_probability : Double
  maximum_change_probability : Double
  maximum_step : Int
  final_run_length : Int
}

///|
pub fn TraceSummary::points(self : TraceSummary) -> Int {
  self.points
}

///|
pub fn TraceSummary::event_count(self : TraceSummary) -> Int {
  self.event_count
}

///|
pub fn TraceSummary::maximum_change_probability(self : TraceSummary) -> Double {
  self.maximum_change_probability
}

///|
pub fn TraceSummary::maximum_step(self : TraceSummary) -> Int {
  self.maximum_step
}

///|
/// Summarize a trace under a caller-supplied probability threshold.
pub fn summarize(
  updates : Array[Update],
  threshold : Double,
) -> Result[TraceSummary, DetectError] {
  if !finite(threshold) || threshold <= 0.0 || threshold >= 1.0 {
    return Err(InvalidHazard)
  }
  if updates.length() == 0 {
    return Ok({
      points: 0,
      event_count: 0,
      mean_change_probability: 0.0,
      maximum_change_probability: 0.0,
      maximum_step: 0,
      final_run_length: 0,
    })
  }
  let mut total = 0.0
  let mut event_count = 0
  let mut maximum = updates[0].change_probability
  let mut maximum_step = updates[0].step
  for update in updates {
    total = total + update.change_probability
    if update.change_probability >= threshold {
      event_count = event_count + 1
    }
    if update.change_probability > maximum {
      maximum = update.change_probability
      maximum_step = update.step
    }
  }
  Ok({
    points: updates.length(),
    event_count,
    mean_change_probability: total / updates.length().to_double(),
    maximum_change_probability: maximum,
    maximum_step,
    final_run_length: updates[updates.length() - 1].most_likely_run_length,
  })
}

///|
/// Return the first update whose probability reaches a threshold.
pub fn first_crossing(updates : Array[Update], threshold : Double) -> Update? {
  for update in updates {
    if update.change_probability >= threshold {
      return Some(update)
    }
  }
  None
}
