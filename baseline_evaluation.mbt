///|
/// Comparable result for a baseline detector on one labeled numeric trace.
pub struct BaselineExperiment {
  detector : String
  scenario : String
  threshold : Double
  score : Score
  alerts : Int
}

///|
pub fn BaselineExperiment::f1(self : BaselineExperiment) -> Double {
  self.score.f1()
}

///|
pub fn BaselineExperiment::precision(self : BaselineExperiment) -> Double {
  self.score.precision()
}

///|
pub fn BaselineExperiment::recall(self : BaselineExperiment) -> Double {
  self.score.recall()
}

///|
pub fn ewma_boundaries(updates : Array[EwmaUpdate]) -> Array[Boundary] {
  let result = []
  for update in updates {
    if update.drift {
      result.push(Boundary::new(update.step).unwrap())
    }
  }
  result
}

///|
pub fn cusum_boundaries(updates : Array[CusumUpdate]) -> Array[Boundary] {
  let result = []
  for update in updates {
    if update.drift {
      result.push(Boundary::new(update.step).unwrap())
    }
  }
  result
}

///|
pub fn evaluate_ewma(
  scenario : NumericScenario,
  config : EwmaConfig,
  tolerance : Int,
) -> Result[BaselineExperiment, DetectError] {
  let detector = EwmaDetector::new(config)
  let updates = match detector.replay(scenario.values) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let boundaries = ewma_boundaries(updates)
  let score = match
    score_boundaries(scenario.boundaries, boundaries, tolerance) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  Ok({
    detector: "ewma",
    scenario: scenario.name,
    threshold: config.z_threshold,
    score,
    alerts: boundaries.length(),
  })
}

///|
pub fn evaluate_cusum(
  scenario : NumericScenario,
  config : CusumConfig,
  tolerance : Int,
) -> Result[BaselineExperiment, DetectError] {
  let detector = CusumDetector::new(config)
  let updates = []
  for value in scenario.values {
    let update = match detector.update(value) {
      Ok(result) => result
      Err(error) => return Err(error)
    }
    updates.push(update)
  }
  let boundaries = cusum_boundaries(updates)
  let score = match
    score_boundaries(scenario.boundaries, boundaries, tolerance) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  Ok({
    detector: "cusum",
    scenario: scenario.name,
    threshold: config.threshold,
    score,
    alerts: boundaries.length(),
  })
}

///|
pub fn sweep_ewma(
  scenario : NumericScenario,
  alphas : Array[Double],
  z_thresholds : Array[Double],
  warmup : Int,
  tolerance : Int,
) -> Result[Array[BaselineExperiment], DetectError] {
  if alphas.length() == 0 || z_thresholds.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  let result = []
  for alpha in alphas {
    for threshold in z_thresholds {
      let config = match EwmaConfig::new(alpha, threshold, warmup) {
        Ok(value) => value
        Err(error) => return Err(error)
      }
      let experiment = match evaluate_ewma(scenario, config, tolerance) {
        Ok(value) => value
        Err(error) => return Err(error)
      }
      result.push(experiment)
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
pub fn sweep_cusum(
  scenario : NumericScenario,
  reference_mean : Double,
  slacks : Array[Double],
  thresholds : Array[Double],
  tolerance : Int,
) -> Result[Array[BaselineExperiment], DetectError] {
  if slacks.length() == 0 || thresholds.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  let result = []
  for slack in slacks {
    for threshold in thresholds {
      let config = match CusumConfig::new(reference_mean, slack, threshold) {
        Ok(value) => value
        Err(error) => return Err(error)
      }
      let experiment = match evaluate_cusum(scenario, config, tolerance) {
        Ok(value) => value
        Err(error) => return Err(error)
      }
      result.push(experiment)
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
pub enum DetectorKind {
  Bocpd
  Ewma
  Cusum
}

///|
pub fn DetectorKind::label(self : DetectorKind) -> String {
  match self {
    Bocpd => "bocpd"
    Ewma => "ewma"
    Cusum => "cusum"
  }
}

///|
pub struct ComparisonRow {
  detector : DetectorKind
  score : Score
  alert_count : Int
}

///|
pub fn ComparisonRow::f1(self : ComparisonRow) -> Double {
  self.score.f1()
}

///|
pub fn rank_comparisons(rows : Array[ComparisonRow]) -> Array[ComparisonRow] {
  rows.sort_by(fn(left, right) {
    if left.f1() > right.f1() {
      -1
    } else if left.f1() < right.f1() {
      1
    } else if left.alert_count < right.alert_count {
      -1
    } else if left.alert_count > right.alert_count {
      1
    } else {
      0
    }
  })
  rows
}

///|
pub fn compare_default_detectors(
  scenario : NumericScenario,
  model : GaussianModel,
  tolerance : Int,
) -> Result[Array[ComparisonRow], DetectError] {
  let bocpd = match
    evaluate_numeric(scenario, model, 0.05, 64, 0.8, tolerance) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let ewma = match
    evaluate_ewma(scenario, EwmaConfig::new(0.15, 2.0, 5).unwrap(), tolerance) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let cusum = match
    evaluate_cusum(
      scenario,
      CusumConfig::new(0.0, 0.3, 3.0).unwrap(),
      tolerance,
    ) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let bocpd_updates = match
    Detector::new(Config::new(0.05, 64, model).unwrap()).replay(scenario.values) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let bocpd_alerts = match predicted_boundaries(bocpd_updates, 0.8) {
    Ok(value) => value.length()
    Err(error) => return Err(error)
  }
  Ok(
    rank_comparisons([
      { detector: Bocpd, score: bocpd.score, alert_count: bocpd_alerts, },
      { detector: Ewma, score: ewma.score, alert_count: ewma.alerts, },
      { detector: Cusum, score: cusum.score, alert_count: cusum.alerts, },
    ]),
  )
}

///|
pub fn baseline_experiment_csv_header() -> String {
  "detector,scenario,threshold,precision,recall,f1,alerts"
}

///|
pub fn format_baseline_experiment(experiment : BaselineExperiment) -> String {
  experiment.detector +
  "," +
  experiment.scenario +
  "," +
  experiment.threshold.to_string() +
  "," +
  experiment.precision().to_string() +
  "," +
  experiment.recall().to_string() +
  "," +
  experiment.f1().to_string() +
  "," +
  experiment.alerts.to_string()
}
