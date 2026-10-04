///|
/// One threshold evaluation point, suitable for CSV output or model selection.
pub struct ThresholdPoint {
  threshold : Double
  score : Score
  alerts : Int
}

///|
pub fn ThresholdPoint::precision(self : ThresholdPoint) -> Double {
  self.score.precision()
}

///|
pub fn ThresholdPoint::recall(self : ThresholdPoint) -> Double {
  self.score.recall()
}

///|
pub fn ThresholdPoint::f1(self : ThresholdPoint) -> Double {
  self.score.f1()
}

///|
pub fn ThresholdPoint::false_positive_rate(
  self : ThresholdPoint,
  known_negatives : Int,
) -> Double {
  if known_negatives <= 0 {
    0.0
  } else {
    self.score.false_positives.to_double() / known_negatives.to_double()
  }
}

///|
pub fn threshold_grid(
  start : Double,
  finish : Double,
  count : Int,
) -> Result[Array[Double], DetectError] {
  if !finite(start) ||
    !finite(finish) ||
    start <= 0.0 ||
    finish >= 1.0 ||
    start > finish ||
    count <= 0 {
    return Err(InvalidHazard)
  }
  if count == 1 {
    return Ok([start])
  }
  let result = []
  for index in 0..<count {
    result.push(
      start + (finish - start) * index.to_double() / (count - 1).to_double(),
    )
  }
  Ok(result)
}

///|
pub fn threshold_points(
  updates : Array[Update],
  truth : Array[Boundary],
  thresholds : Array[Double],
  tolerance : Int,
) -> Result[Array[ThresholdPoint], DetectError] {
  if thresholds.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  let result = []
  for threshold in thresholds {
    let predictions = match predicted_boundaries(updates, threshold) {
      Ok(value) => value
      Err(error) => return Err(error)
    }
    let score = match score_boundaries(truth, predictions, tolerance) {
      Ok(value) => value
      Err(error) => return Err(error)
    }
    result.push({ threshold, score, alerts: predictions.length(), })
  }
  Ok(result)
}

///|
pub fn best_f1_point(
  points : Array[ThresholdPoint],
) -> Result[ThresholdPoint, DetectError] {
  if points.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  let mut best = points[0]
  for point in points {
    if point.f1() > best.f1() ||
      (point.f1() == best.f1() && point.threshold < best.threshold) {
      best = point
    }
  }
  Ok(best)
}

///|
pub fn point_at_precision(
  points : Array[ThresholdPoint],
  minimum_precision : Double,
) -> Result[ThresholdPoint, DetectError] {
  if points.length() == 0 ||
    !finite(minimum_precision) ||
    minimum_precision < 0.0 ||
    minimum_precision > 1.0 {
    return Err(InvalidHazard)
  }
  let mut found : ThresholdPoint? = None
  for point in points {
    if point.precision() >= minimum_precision {
      match found {
        None => found = Some(point)
        Some(current) =>
          if point.recall() > current.recall() ||
            (
              point.recall() == current.recall() &&
              point.threshold < current.threshold
            ) {
            found = Some(point)
          }
      }
    }
  }
  match found {
    Some(point) => Ok(point)
    None => Err(InvalidRunLengthLimit)
  }
}

///|
pub struct CurvePoint {
  threshold : Double
  recall : Double
  precision : Double
  false_positive_rate : Double
}

///|
pub fn precision_recall_curve(
  points : Array[ThresholdPoint],
) -> Array[CurvePoint] {
  let result = []
  for point in points {
    result.push({
      threshold: point.threshold,
      recall: point.recall(),
      precision: point.precision(),
      false_positive_rate: 0.0,
    })
  }
  result.sort_by(fn(left, right) {
    if left.recall < right.recall {
      -1
    } else if left.recall > right.recall {
      1
    } else if left.threshold < right.threshold {
      -1
    } else if left.threshold > right.threshold {
      1
    } else {
      0
    }
  })
  result
}

///|
pub fn roc_curve(
  points : Array[ThresholdPoint],
  negatives : Int,
) -> Result[Array[CurvePoint], DetectError] {
  if negatives <= 0 {
    return Err(InvalidRunLengthLimit)
  }
  let result = []
  for point in points {
    result.push({
      threshold: point.threshold,
      recall: point.recall(),
      precision: point.precision(),
      false_positive_rate: point.false_positive_rate(negatives),
    })
  }
  result.sort_by(fn(left, right) {
    if left.false_positive_rate < right.false_positive_rate {
      -1
    } else if left.false_positive_rate > right.false_positive_rate {
      1
    } else {
      0
    }
  })
  Ok(result)
}

///|
pub fn trapezoid_area(
  horizontal : Array[Double],
  vertical : Array[Double],
) -> Result[Double, DetectError] {
  if horizontal.length() < 2 || horizontal.length() != vertical.length() {
    return Err(InvalidRunLengthLimit)
  }
  let mut total = 0.0
  for index in 1..<horizontal.length() {
    if horizontal[index] < horizontal[index - 1] {
      return Err(InvalidRunLengthLimit)
    }
    let width = horizontal[index] - horizontal[index - 1]
    total = total + width * (vertical[index] + vertical[index - 1]) / 2.0
  }
  Ok(total)
}

///|
pub fn area_under_roc(curve : Array[CurvePoint]) -> Result[Double, DetectError] {
  if curve.length() < 2 {
    return Err(InvalidRunLengthLimit)
  }
  let horizontal = []
  let vertical = []
  for point in curve {
    horizontal.push(point.false_positive_rate)
    vertical.push(point.recall)
  }
  trapezoid_area(horizontal, vertical)
}

///|
pub struct HazardPoint {
  hazard : Double
  best_threshold : Double
  best_score : Score
  maximum_probability : Double
}

///|
pub fn HazardPoint::f1(self : HazardPoint) -> Double {
  self.best_score.f1()
}

///|
pub fn evaluate_hazards(
  scenario : NumericScenario,
  model : GaussianModel,
  hazards : Array[Double],
  max_run_length : Int,
  thresholds : Array[Double],
  tolerance : Int,
) -> Result[Array[HazardPoint], DetectError] {
  if hazards.length() == 0 || thresholds.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  let result = []
  for hazard in hazards {
    let config = match Config::new(hazard, max_run_length, model) {
      Ok(value) => value
      Err(error) => return Err(error)
    }
    let detector = Detector::new(config)
    let updates = match detector.replay(scenario.values) {
      Ok(value) => value
      Err(error) => return Err(error)
    }
    let points = match
      threshold_points(updates, scenario.boundaries, thresholds, tolerance) {
      Ok(value) => value
      Err(error) => return Err(error)
    }
    let best = match best_f1_point(points) {
      Ok(value) => value
      Err(error) => return Err(error)
    }
    let mut maximum = 0.0
    for update in updates {
      if update.change_probability > maximum {
        maximum = update.change_probability
      }
    }
    result.push({
      hazard,
      best_threshold: best.threshold,
      best_score: best.score,
      maximum_probability: maximum,
    })
  }
  result.sort_by(fn(left, right) {
    if left.f1() > right.f1() {
      -1
    } else if left.f1() < right.f1() {
      1
    } else if left.hazard < right.hazard {
      -1
    } else {
      1
    }
  })
  Ok(result)
}

///|
pub fn posterior_quantile(
  probabilities : Array[Double],
  quantile : Double,
) -> Result[Int, DetectError] {
  if probabilities.length() == 0 ||
    !finite(quantile) ||
    quantile < 0.0 ||
    quantile > 1.0 {
    return Err(InvalidHazard)
  }
  let mut total = 0.0
  for probability in probabilities {
    if probability < 0.0 || !finite(probability) {
      return Err(InvalidVariance)
    }
    total = total + probability
  }
  if total <= 0.0 {
    return Err(InvalidVariance)
  }
  let mut cumulative = 0.0
  for index in 0..<probabilities.length() {
    cumulative = cumulative + probabilities[index] / total
    if cumulative >= quantile {
      return Ok(index)
    }
  }
  Ok(probabilities.length() - 1)
}

///|
pub fn mean_change_probability(
  updates : Array[Update],
) -> Result[Double, DetectError] {
  if updates.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  let mut total = 0.0
  for update in updates {
    total = total + update.change_probability
  }
  Ok(total / updates.length().to_double())
}
