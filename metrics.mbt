///|
/// A labeled boundary position for deterministic detector evaluation.
pub struct Boundary {
  step : Int
}

///|
pub fn Boundary::new(step : Int) -> Result[Boundary, DetectError] {
  if step <= 0 {
    return Err(InvalidRunLengthLimit)
  }
  Ok({ step, })
}

///|
/// Confusion-matrix counts for boundary matching within a tolerance window.
pub struct Score {
  true_positives : Int
  false_positives : Int
  false_negatives : Int
  mean_delay : Double
}

///|
pub fn Score::precision(self : Score) -> Double {
  let denominator = self.true_positives + self.false_positives
  if denominator == 0 {
    0.0
  } else {
    self.true_positives.to_double() / denominator.to_double()
  }
}

///|
pub fn Score::recall(self : Score) -> Double {
  let denominator = self.true_positives + self.false_negatives
  if denominator == 0 {
    0.0
  } else {
    self.true_positives.to_double() / denominator.to_double()
  }
}

///|
pub fn Score::f1(self : Score) -> Double {
  let precision = self.precision()
  let recall = self.recall()
  if precision + recall == 0.0 {
    0.0
  } else {
    2.0 * precision * recall / (precision + recall)
  }
}

///|
/// Extract candidate boundary positions from a posterior trace.
pub fn predicted_boundaries(
  updates : Array[Update],
  threshold : Double,
) -> Result[Array[Boundary], DetectError] {
  if !finite(threshold) || threshold <= 0.0 || threshold >= 1.0 {
    return Err(InvalidHazard)
  }
  let result : Array[Boundary] = []
  for update in updates {
    if update.change_probability >= threshold {
      result.push({ step: update.step, })
    }
  }
  Ok(result)
}

///|
/// Match each truth boundary to at most one closest predicted boundary.
pub fn score_boundaries(
  truth : Array[Boundary],
  predicted : Array[Boundary],
  tolerance : Int,
) -> Result[Score, DetectError] {
  if tolerance < 0 {
    return Err(InvalidRunLengthLimit)
  }
  let used = Array::make(predicted.length(), false)
  let mut true_positives = 0
  let mut total_delay = 0.0
  for expected in truth {
    let mut matched = -1
    let mut best_distance = tolerance + 1
    for index in 0..<predicted.length() {
      let difference = predicted[index].step - expected.step
      let distance = if difference < 0 { -difference } else { difference }
      if !used[index] && distance <= tolerance && distance < best_distance {
        matched = index
        best_distance = distance
      }
    }
    if matched >= 0 {
      used[matched] = true
      true_positives = true_positives + 1
      total_delay = total_delay +
        (predicted[matched].step - expected.step).to_double()
    }
  }
  let false_positives = predicted.length() - true_positives
  let false_negatives = truth.length() - true_positives
  Ok({
    true_positives,
    false_positives,
    false_negatives,
    mean_delay: if true_positives == 0 {
      0.0
    } else {
      total_delay / true_positives.to_double()
    },
  })
}

///|
/// Average run-length posterior entropy in bits for a trace of posterior rows.
pub fn posterior_entropy(
  probabilities : Array[Double],
) -> Result[Double, DetectError] {
  if probabilities.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  let mut total = 0.0
  let mut entropy = 0.0
  for probability in probabilities {
    if !finite(probability) || probability < 0.0 {
      return Err(NonFiniteObservation)
    }
    total = total + probability
    if probability > 0.0 {
      entropy = entropy - probability * @math.ln(probability) / @math.ln(2.0)
    }
  }
  if total < 0.999999 || total > 1.000001 {
    return Err(InvalidVariance)
  }
  Ok(entropy)
}
