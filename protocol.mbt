///|
/// Immutable plan for one reproducible Gaussian BOCPD evaluation.
pub struct ExperimentPlan {
  name : String
  hazard : Double
  max_run_length : Int
  threshold : Double
  tolerance : Int
  seed : Int
}

///|
pub fn ExperimentPlan::new(
  name : String,
  hazard : Double,
  max_run_length : Int,
  threshold : Double,
  tolerance : Int,
  seed : Int,
) -> Result[ExperimentPlan, DetectError] {
  if name.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  if !finite(hazard) || hazard <= 0.0 || hazard >= 1.0 {
    return Err(InvalidHazard)
  }
  if max_run_length < 2 || tolerance < 0 || seed <= 0 {
    return Err(InvalidRunLengthLimit)
  }
  if !finite(threshold) || threshold <= 0.0 || threshold >= 1.0 {
    return Err(InvalidHazard)
  }
  Ok({ name, hazard, max_run_length, threshold, tolerance, seed, })
}

///|
pub fn ExperimentPlan::describe(self : ExperimentPlan) -> String {
  "name=" +
  self.name +
  ",hazard=" +
  self.hazard.to_string() +
  ",max_run_length=" +
  self.max_run_length.to_string() +
  ",threshold=" +
  self.threshold.to_string() +
  ",tolerance=" +
  self.tolerance.to_string() +
  ",seed=" +
  self.seed.to_string()
}

///|
pub fn default_experiment_plan() -> ExperimentPlan {
  ExperimentPlan::new("mean-shift-default", 0.05, 64, 0.8, 3, 73).unwrap()
}

///|
pub fn conservative_experiment_plan() -> ExperimentPlan {
  ExperimentPlan::new("mean-shift-conservative", 0.02, 128, 0.95, 3, 73).unwrap()
}

///|
pub fn sensitive_experiment_plan() -> ExperimentPlan {
  ExperimentPlan::new("mean-shift-sensitive", 0.1, 64, 0.6, 3, 73).unwrap()
}

///|
pub struct PlanResult {
  plan : ExperimentPlan
  experiment : NumericExperiment
  diagnostics : TraceDiagnostics
  predictions : Array[Boundary]
}

///|
pub fn PlanResult::summary(self : PlanResult) -> String {
  self.plan.name +
  ",f1=" +
  self.experiment.f1().to_string() +
  ",precision=" +
  self.experiment.score.precision().to_string() +
  ",recall=" +
  self.experiment.score.recall().to_string() +
  ",maximum_probability=" +
  self.diagnostics.maximum_change_probability.to_string()
}

///|
pub fn execute_plan(
  plan : ExperimentPlan,
  scenario : NumericScenario,
  model : GaussianModel,
) -> Result[PlanResult, DetectError] {
  let config = match Config::new(plan.hazard, plan.max_run_length, model) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let detector = Detector::new(config)
  let updates = match detector.replay(scenario.values) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let predictions = match predicted_boundaries(updates, plan.threshold) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let score = match
    score_boundaries(scenario.boundaries, predictions, plan.tolerance) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let trace = match summarize(updates, plan.threshold) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let diagnostics = match diagnose_trace(updates) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  Ok({
    plan,
    experiment: {
      scenario: scenario.name,
      hazard: plan.hazard,
      threshold: plan.threshold,
      score,
      trace,
    },
    diagnostics,
    predictions,
  })
}

///|
pub fn execute_default_plans(
  scenario : NumericScenario,
  model : GaussianModel,
) -> Result[Array[PlanResult], DetectError] {
  let plans = [
    default_experiment_plan(),
    conservative_experiment_plan(),
    sensitive_experiment_plan(),
  ]
  let result = []
  for plan in plans {
    let value = match execute_plan(plan, scenario, model) {
      Ok(item) => item
      Err(error) => return Err(error)
    }
    result.push(value)
  }
  result.sort_by(fn(left, right) {
    if left.experiment.f1() > right.experiment.f1() {
      -1
    } else if left.experiment.f1() < right.experiment.f1() {
      1
    } else if left.plan.threshold < right.plan.threshold {
      -1
    } else {
      1
    }
  })
  Ok(result)
}

///|
pub struct ReproducibilityRecord {
  project : String
  scenario : String
  observation_count : Int
  boundary_count : Int
  plan : String
}

///|
pub fn reproducibility_record(
  scenario : NumericScenario,
  plan : ExperimentPlan,
) -> ReproducibilityRecord {
  {
    project: "BayesShift",
    scenario: scenario.name,
    observation_count: scenario.values.length(),
    boundary_count: scenario.boundaries.length(),
    plan: plan.describe(),
  }
}

///|
pub fn ReproducibilityRecord::to_csv(self : ReproducibilityRecord) -> String {
  self.project +
  "," +
  self.scenario +
  "," +
  self.observation_count.to_string() +
  "," +
  self.boundary_count.to_string() +
  "," +
  self.plan
}

///|
pub fn reproducibility_csv_header() -> String {
  "project,scenario,observation_count,boundary_count,plan"
}

///|
pub fn validate_scenario(
  scenario : NumericScenario,
) -> Result[Unit, DetectError] {
  if scenario.values.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  let mut previous = 0
  for boundary in scenario.boundaries {
    if boundary.step <= previous || boundary.step > scenario.values.length() {
      return Err(InvalidRunLengthLimit)
    }
    previous = boundary.step
  }
  for value in scenario.values {
    if !finite(value) {
      return Err(NonFiniteObservation)
    }
  }
  Ok(())
}

///|
pub fn plan_grid(
  hazards : Array[Double],
  run_lengths : Array[Int],
  thresholds : Array[Double],
  tolerance : Int,
  seed : Int,
) -> Result[Array[ExperimentPlan], DetectError] {
  if hazards.length() == 0 ||
    run_lengths.length() == 0 ||
    thresholds.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  let result = []
  let mut sequence = 0
  for hazard in hazards {
    for run_length in run_lengths {
      for threshold in thresholds {
        sequence = sequence + 1
        let name = "grid-" + sequence.to_string()
        let plan = match
          ExperimentPlan::new(
            name,
            hazard,
            run_length,
            threshold,
            tolerance,
            seed + sequence,
          ) {
          Ok(value) => value
          Err(error) => return Err(error)
        }
        result.push(plan)
      }
    }
  }
  Ok(result)
}

///|
pub fn execute_grid(
  plans : Array[ExperimentPlan],
  scenario : NumericScenario,
  model : GaussianModel,
) -> Result[Array[PlanResult], DetectError] {
  let result = []
  for plan in plans {
    let evaluation = match execute_plan(plan, scenario, model) {
      Ok(value) => value
      Err(error) => return Err(error)
    }
    result.push(evaluation)
  }
  result.sort_by(fn(left, right) {
    if left.experiment.f1() > right.experiment.f1() {
      -1
    } else if left.experiment.f1() < right.experiment.f1() {
      1
    } else {
      0
    }
  })
  Ok(result)
}

///|
pub fn plan_result_csv_header() -> String {
  "plan,scenario,hazard,threshold,precision,recall,f1,max_probability,max_step"
}

///|
pub fn format_plan_result(result : PlanResult) -> String {
  result.plan.name +
  "," +
  result.experiment.scenario +
  "," +
  result.plan.hazard.to_string() +
  "," +
  result.plan.threshold.to_string() +
  "," +
  result.experiment.score.precision().to_string() +
  "," +
  result.experiment.score.recall().to_string() +
  "," +
  result.experiment.f1().to_string() +
  "," +
  result.diagnostics.maximum_change_probability.to_string() +
  "," +
  result.diagnostics.maximum_change_step.to_string()
}
