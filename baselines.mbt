///|
/// Exponentially weighted moving-average baseline for online drift comparison.
pub struct EwmaConfig {
  alpha : Double
  z_threshold : Double
  warmup : Int
}

///|
pub fn EwmaConfig::new(
  alpha : Double,
  z_threshold : Double,
  warmup : Int,
) -> Result[EwmaConfig, DetectError] {
  if !finite(alpha) || alpha <= 0.0 || alpha >= 1.0 {
    return Err(InvalidHazard)
  }
  if !finite(z_threshold) || z_threshold <= 0.0 || warmup < 2 {
    return Err(InvalidRunLengthLimit)
  }
  Ok({ alpha, z_threshold, warmup, })
}

///|
pub struct EwmaUpdate {
  step : Int
  value : Double
  mean : Double
  variance : Double
  z : Double
  drift : Bool
}

///|
pub struct EwmaDetector {
  config : EwmaConfig
  mut step : Int
  mut mean : Double
  mut variance : Double
  mut initialized : Bool
}

///|
pub fn EwmaDetector::new(config : EwmaConfig) -> EwmaDetector {
  { config, step: 0, mean: 0.0, variance: 0.0, initialized: false, }
}

///|
/// Update the EWMA mean and variance before comparing the new point to its prior state.
pub fn EwmaDetector::update(
  self : EwmaDetector,
  value : Double,
) -> Result[EwmaUpdate, DetectError] {
  if !finite(value) {
    return Err(NonFiniteObservation)
  }
  if !self.initialized {
    self.initialized = true
    self.step = 1
    self.mean = value
    return Ok({
      step: 1,
      value,
      mean: value,
      variance: 0.0,
      z: 0.0,
      drift: false,
    })
  }
  let prior_mean = self.mean
  let delta = value - prior_mean
  let alpha = self.config.alpha
  self.mean = prior_mean + alpha * delta
  self.variance = (1.0 - alpha) * (self.variance + alpha * delta * delta)
  self.step = self.step + 1
  let scale = if self.variance <= 0.0 {
    0.0
  } else {
    @math.exp(@math.ln(self.variance) * 0.5)
  }
  let z = if scale == 0.0 { 0.0 } else { delta / scale }
  Ok({
    step: self.step,
    value,
    mean: self.mean,
    variance: self.variance,
    z,
    drift: self.step >= self.config.warmup &&
    (z > self.config.z_threshold || z < -self.config.z_threshold),
  })
}

///|
pub fn EwmaDetector::replay(
  self : EwmaDetector,
  values : Array[Double],
) -> Result[Array[EwmaUpdate], DetectError] {
  let result : Array[EwmaUpdate] = []
  for value in values {
    match self.update(value) {
      Ok(update) => result.push(update)
      Err(error) => return Err(error)
    }
  }
  Ok(result)
}

///|
/// CUSUM is a second deterministic baseline for a known reference mean.
pub struct CusumConfig {
  reference_mean : Double
  slack : Double
  threshold : Double
}

///|
pub fn CusumConfig::new(
  reference_mean : Double,
  slack : Double,
  threshold : Double,
) -> Result[CusumConfig, DetectError] {
  if !finite(reference_mean) ||
    !finite(slack) ||
    !finite(threshold) ||
    slack < 0.0 ||
    threshold <= 0.0 {
    return Err(InvalidVariance)
  }
  Ok({ reference_mean, slack, threshold, })
}

///|
pub struct CusumDetector {
  config : CusumConfig
  mut positive : Double
  mut negative : Double
  mut step : Int
}

///|
pub struct CusumUpdate {
  step : Int
  positive : Double
  negative : Double
  drift : Bool
}

///|
pub fn CusumDetector::new(config : CusumConfig) -> CusumDetector {
  { config, positive: 0.0, negative: 0.0, step: 0, }
}

///|
pub fn CusumDetector::update(
  self : CusumDetector,
  value : Double,
) -> Result[CusumUpdate, DetectError] {
  if !finite(value) {
    return Err(NonFiniteObservation)
  }
  let residual = value - self.config.reference_mean
  let positive = self.positive + residual - self.config.slack
  let negative = self.negative - residual - self.config.slack
  self.positive = if positive > 0.0 { positive } else { 0.0 }
  self.negative = if negative > 0.0 { negative } else { 0.0 }
  self.step = self.step + 1
  Ok({
    step: self.step,
    positive: self.positive,
    negative: self.negative,
    drift: self.positive >= self.config.threshold ||
    self.negative >= self.config.threshold,
  })
}
