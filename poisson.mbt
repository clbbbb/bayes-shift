///|
/// Conjugate Gamma prior for a Poisson segment rate.
pub struct PoissonModel {
  prior_shape : Double
  prior_rate : Double
}

///|
pub fn PoissonModel::new(
  prior_shape : Double,
  prior_rate : Double,
) -> Result[PoissonModel, DetectError] {
  if !finite(prior_shape) || prior_shape <= 0.0 {
    return Err(InvalidPriorStrength)
  }
  if !finite(prior_rate) || prior_rate <= 0.0 {
    return Err(InvalidVariance)
  }
  Ok({ prior_shape, prior_rate, })
}

///|
/// Bounded BOCPD configuration for non-negative integer observations.
pub struct PoissonConfig {
  hazard : Double
  max_run_length : Int
  model : PoissonModel
}

///|
pub fn PoissonConfig::new(
  hazard : Double,
  max_run_length : Int,
  model : PoissonModel,
) -> Result[PoissonConfig, DetectError] {
  if !finite(hazard) || hazard <= 0.0 || hazard >= 1.0 {
    return Err(InvalidHazard)
  }
  if max_run_length < 2 {
    return Err(InvalidRunLengthLimit)
  }
  Ok({ hazard, max_run_length, model, })
}

///|
/// Posterior result for one count observation.
pub struct CountUpdate {
  step : Int
  count : Int
  change_probability : Double
  most_likely_run_length : Int
  posterior_rate : Double
}

///|
pub fn CountUpdate::step(self : CountUpdate) -> Int {
  self.step
}

///|
pub fn CountUpdate::change_probability(self : CountUpdate) -> Double {
  self.change_probability
}

///|
pub fn CountUpdate::most_likely_run_length(self : CountUpdate) -> Int {
  self.most_likely_run_length
}

///|
pub fn CountUpdate::posterior_rate(self : CountUpdate) -> Double {
  self.posterior_rate
}

///|
/// Separate detector avoids making numeric and count-model callers pay for an enum dispatch.
pub struct PoissonDetector {
  config : PoissonConfig
  mut step : Int
  mut log_posterior : Array[Double]
  mut totals : Array[Int]
  mut observations : Array[Int]
}

///|
pub fn PoissonDetector::new(config : PoissonConfig) -> PoissonDetector {
  { config, step: 0, log_posterior: [0.0], totals: [0], observations: [0], }
}

///|
pub fn PoissonDetector::steps(self : PoissonDetector) -> Int {
  self.step
}

///|
pub fn PoissonDetector::run_length_posterior(
  self : PoissonDetector,
) -> Array[Double] {
  self.log_posterior.map(fn(value) { @math.exp(value) })
}

///|
/// Gamma-Poisson posterior predictive probability in log form.
fn PoissonDetector::predictive(
  self : PoissonDetector,
  total : Int,
  observations : Int,
  value : Int,
) -> Double {
  let model = self.config.model
  let shape = model.prior_shape + total.to_double()
  let rate = model.prior_rate + observations.to_double()
  let next_rate = rate + 1.0
  let numerator = match log_gamma(shape + value.to_double()) {
    Ok(result) => result
    Err(_) => negative_infinity
  }
  let denominator = match log_gamma(shape) {
    Ok(result) => result
    Err(_) => negative_infinity
  }
  let factorial = match log_factorial(value) {
    Ok(result) => result
    Err(_) => negative_infinity
  }
  numerator -
  denominator -
  factorial +
  shape * @math.ln(rate / next_rate) +
  value.to_double() * @math.ln(1.0 / next_rate)
}

///|
/// Update with a non-negative count. The change branch always predicts from the prior.
pub fn PoissonDetector::update(
  self : PoissonDetector,
  value : Int,
) -> Result[CountUpdate, DetectError] {
  if value < 0 {
    return Err(InvalidRunLengthLimit)
  }
  let old_size = self.log_posterior.length()
  let new_size = if old_size < self.config.max_run_length {
    old_size + 1
  } else {
    old_size
  }
  let next_log = Array::make(new_size, negative_infinity)
  let next_totals = Array::make(new_size, 0)
  let next_observations = Array::make(new_size, 0)
  let reset_predictive = self.predictive(0, 0, value)
  let mut changed = negative_infinity
  for run in 0..<old_size {
    changed = log_add(
      changed,
      self.log_posterior[run] + @math.ln(self.config.hazard) + reset_predictive,
    )
    let grown = run + 1
    if grown < new_size {
      let predictive = self.predictive(
        self.totals[run],
        self.observations[run],
        value,
      )
      next_log[grown] = self.log_posterior[run] +
        @math.ln(1.0 - self.config.hazard) +
        predictive
      next_totals[grown] = self.totals[run] + value
      next_observations[grown] = self.observations[run] + 1
    }
  }
  next_log[0] = changed
  next_totals[0] = value
  next_observations[0] = 1
  normalize_log_weights(next_log)
  let mut best = 0
  for index in 1..<new_size {
    if next_log[index] > next_log[best] {
      best = index
    }
  }
  let model = self.config.model
  let rate = (model.prior_shape + next_totals[best].to_double()) /
    (model.prior_rate + next_observations[best].to_double())
  self.step = self.step + 1
  self.log_posterior = next_log
  self.totals = next_totals
  self.observations = next_observations
  Ok({
    step: self.step,
    count: value,
    change_probability: @math.exp(self.log_posterior[0]),
    most_likely_run_length: best,
    posterior_rate: rate,
  })
}

///|
pub fn PoissonDetector::replay(
  self : PoissonDetector,
  values : Array[Int],
) -> Result[Array[CountUpdate], DetectError] {
  let result : Array[CountUpdate] = []
  for value in values {
    match self.update(value) {
      Ok(update) => result.push(update)
      Err(error) => return Err(error)
    }
  }
  Ok(result)
}
