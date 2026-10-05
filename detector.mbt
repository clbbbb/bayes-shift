///|
pub struct Detector {
  config : Config
  mut step : Int
  mut log_posterior : Array[Double]
  mut counts : Array[Int]
  mut sums : Array[Double]
}

///|
pub fn Detector::new(config : Config) -> Detector {
  { config, step: 0, log_posterior: [0.0], counts: [0], sums: [0.0], }
}

///|
pub fn Detector::steps(self : Detector) -> Int {
  self.step
}

///|
pub fn Detector::run_length_posterior(self : Detector) -> Array[Double] {
  self.log_posterior.map(fn(value) { @math.exp(value) })
}

///|
pub fn Detector::snapshot(self : Detector) -> Snapshot {
  {
    step: self.step,
    log_posterior: self.log_posterior.copy(),
    counts: self.counts.copy(),
    sums: self.sums.copy(),
  }
}

///|
pub fn Detector::restore(self : Detector, snapshot : Snapshot) -> Unit {
  self.step = snapshot.step
  self.log_posterior = snapshot.log_posterior.copy()
  self.counts = snapshot.counts.copy()
  self.sums = snapshot.sums.copy()
}

///|
fn Detector::predictive(
  self : Detector,
  count : Int,
  sum : Double,
  value : Double,
) -> Double {
  let model = self.config.model
  let effective = model.prior_strength + count.to_double()
  let mean = (model.prior_strength * model.prior_mean + sum) / effective
  normal_log_density(
    value,
    mean,
    model.observation_variance * (1.0 + 1.0 / effective),
  )
}

///|
pub fn Detector::update(
  self : Detector,
  observation : Double,
) -> Result[Update, DetectError] {
  if !finite(observation) {
    return Err(NonFiniteObservation)
  }
  let old_size = self.log_posterior.length()
  let new_size = if old_size < self.config.max_run_length {
    old_size + 1
  } else {
    old_size
  }
  let next_log = Array::make(new_size, negative_infinity)
  let next_counts = Array::make(new_size, 0)
  let next_sums = Array::make(new_size, 0.0)
  let reset_predictive = self.predictive(0, 0.0, observation)
  // The previous run-length posterior is normalized, so reset mass is
  // hazard times prior predictive without a per-state log-sum-exp reduction.
  let changed = self.config.hazard_log + reset_predictive
  for run in 0..<old_size {
    let grown = run + 1
    if grown < new_size {
      let predictive = self.predictive(
        self.counts[run],
        self.sums[run],
        observation,
      )
      next_log[grown] = self.log_posterior[run] +
        self.config.survival_log +
        predictive
      next_counts[grown] = self.counts[run] + 1
      next_sums[grown] = self.sums[run] + observation
    }
  }
  next_log[0] = changed
  next_counts[0] = 1
  next_sums[0] = observation
  let best = normalize_log_weights(next_log)
  let model = self.config.model
  let mean = (model.prior_strength * model.prior_mean + next_sums[best]) /
    (model.prior_strength + next_counts[best].to_double())
  self.step = self.step + 1
  self.log_posterior = next_log
  self.counts = next_counts
  self.sums = next_sums
  Ok({
    step: self.step,
    observation,
    change_probability: @math.exp(self.log_posterior[0]),
    most_likely_run_length: best,
    posterior_mean: mean,
    retained_run_lengths: new_size,
  })
}

///|
pub fn Detector::replay(
  self : Detector,
  observations : Array[Double],
) -> Result[Array[Update], DetectError] {
  let updates : Array[Update] = []
  for observation in observations {
    match self.update(observation) {
      Ok(value) => updates.push(value)
      Err(error) => return Err(error)
    }
  }
  Ok(updates)
}
