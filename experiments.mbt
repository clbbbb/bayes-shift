///|
/// A named numeric fixture and its expected segment boundaries.
pub struct NumericScenario {
  name : String
  values : Array[Double]
  boundaries : Array[Boundary]
}

///|
pub fn NumericScenario::mean_shift() -> NumericScenario {
  {
    name: "mean-shift",
    values: piecewise_constant([30, 30], [0.0, 5.0], 0.2, 73).unwrap(),
    boundaries: [Boundary::new(31).unwrap()],
  }
}

///|
pub fn NumericScenario::no_change() -> NumericScenario {
  {
    name: "no-change",
    values: piecewise_constant([60], [0.0], 0.2, 73).unwrap(),
    boundaries: [],
  }
}

///|
pub fn NumericScenario::two_shifts() -> NumericScenario {
  {
    name: "two-shifts",
    values: piecewise_constant([20, 20, 20], [0.0, 4.0, -3.0], 0.15, 11).unwrap(),
    boundaries: [Boundary::new(21).unwrap(), Boundary::new(41).unwrap()],
  }
}

///|
/// A fully specified evaluation result for one detector configuration.
pub struct NumericExperiment {
  scenario : String
  hazard : Double
  threshold : Double
  score : Score
  trace : TraceSummary
}

///|
pub fn NumericExperiment::f1(self : NumericExperiment) -> Double {
  self.score.f1()
}

///|
pub fn NumericExperiment::true_positives(self : NumericExperiment) -> Int {
  self.score.true_positives
}

///|
pub fn NumericExperiment::false_positives(self : NumericExperiment) -> Int {
  self.score.false_positives
}

///|
/// Evaluate Gaussian BOCPD against a deterministic labeled scenario.
pub fn evaluate_numeric(
  scenario : NumericScenario,
  model : GaussianModel,
  hazard : Double,
  max_run_length : Int,
  threshold : Double,
  tolerance : Int,
) -> Result[NumericExperiment, DetectError] {
  let config = match Config::new(hazard, max_run_length, model) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let detector = Detector::new(config)
  let updates = match detector.replay(scenario.values) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let predictions = match predicted_boundaries(updates, threshold) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let score = match
    score_boundaries(scenario.boundaries, predictions, tolerance) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let trace = match summarize(updates, threshold) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  Ok({ scenario: scenario.name, hazard, threshold, score, trace, })
}

///|
/// Pick the best threshold by F1; ties prefer the lower threshold.
pub fn sweep_thresholds(
  scenario : NumericScenario,
  model : GaussianModel,
  hazard : Double,
  max_run_length : Int,
  thresholds : Array[Double],
  tolerance : Int,
) -> Result[Array[NumericExperiment], DetectError] {
  if thresholds.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  let result : Array[NumericExperiment] = []
  for threshold in thresholds {
    match
      evaluate_numeric(
        scenario, model, hazard, max_run_length, threshold, tolerance,
      ) {
      Ok(experiment) => result.push(experiment)
      Err(error) => return Err(error)
    }
  }
  result.sort_by(fn(left, right) {
    if left.f1() > right.f1() {
      -1
    } else if left.f1() < right.f1() {
      1
    } else if left.threshold < right.threshold {
      -1
    } else {
      1
    }
  })
  Ok(result)
}

///|
/// Count stream fixture for the Poisson-Gamma model.
pub struct CountScenario {
  name : String
  values : Array[Int]
  boundaries : Array[Boundary]
}

///|
pub fn CountScenario::rate_shift() -> CountScenario {
  {
    name: "rate-shift",
    values: [1, 0, 2, 1, 1, 0, 1, 24, 28, 31, 27, 30],
    boundaries: [Boundary::new(8).unwrap()],
  }
}

///|
pub struct CountExperiment {
  scenario : String
  threshold : Double
  score : Score
  maximum_probability : Double
}

///|
/// Evaluate a count-rate detector with the same matching rule as numeric BOCPD.
pub fn evaluate_count(
  scenario : CountScenario,
  model : PoissonModel,
  hazard : Double,
  max_run_length : Int,
  threshold : Double,
  tolerance : Int,
) -> Result[CountExperiment, DetectError] {
  let config = match PoissonConfig::new(hazard, max_run_length, model) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let detector = PoissonDetector::new(config)
  let updates = match detector.replay(scenario.values) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let predictions : Array[Boundary] = []
  let mut maximum = 0.0
  for update in updates {
    if update.change_probability > maximum {
      maximum = update.change_probability
    }
    if update.change_probability >= threshold {
      predictions.push(Boundary::new(update.step).unwrap())
    }
  }
  let score = match
    score_boundaries(scenario.boundaries, predictions, tolerance) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  Ok({
    scenario: scenario.name,
    threshold,
    score,
    maximum_probability: maximum,
  })
}
