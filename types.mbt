///|
pub enum DetectError {
  InvalidHazard
  InvalidVariance
  InvalidPriorStrength
  InvalidRunLengthLimit
  NonFiniteObservation
}

///|
pub struct GaussianModel {
  prior_mean : Double
  prior_strength : Double
  observation_variance : Double
}

///|
pub fn GaussianModel::new(
  prior_mean : Double,
  prior_strength : Double,
  observation_variance : Double,
) -> Result[GaussianModel, DetectError] {
  if !finite(prior_mean) {
    return Err(NonFiniteObservation)
  }
  if !finite(prior_strength) || prior_strength <= 0.0 {
    return Err(InvalidPriorStrength)
  }
  if !finite(observation_variance) || observation_variance <= 0.0 {
    return Err(InvalidVariance)
  }
  Ok({ prior_mean, prior_strength, observation_variance, })
}

///|
pub struct Config {
  hazard : Double
  hazard_log : Double
  survival_log : Double
  max_run_length : Int
  model : GaussianModel
}

///|
pub fn Config::new(
  hazard : Double,
  max_run_length : Int,
  model : GaussianModel,
) -> Result[Config, DetectError] {
  if !finite(hazard) || hazard <= 0.0 || hazard >= 1.0 {
    return Err(InvalidHazard)
  }
  if max_run_length < 2 {
    return Err(InvalidRunLengthLimit)
  }
  Ok({
    hazard,
    hazard_log: @math.ln(hazard),
    survival_log: @math.ln(1.0 - hazard),
    max_run_length,
    model,
  })
}

///|
pub struct Update {
  step : Int
  observation : Double
  change_probability : Double
  most_likely_run_length : Int
  posterior_mean : Double
  retained_run_lengths : Int
}

///|
pub struct Snapshot {
  step : Int
  log_posterior : Array[Double]
  counts : Array[Int]
  sums : Array[Double]
}

///|
pub fn Update::step(self : Update) -> Int {
  self.step
}

///|
pub fn Update::change_probability(self : Update) -> Double {
  self.change_probability
}

///|
pub fn Update::most_likely_run_length(self : Update) -> Int {
  self.most_likely_run_length
}
