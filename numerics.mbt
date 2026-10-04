///|
/// Numerically stable scalar utilities shared by alternative observation models.
pub fn log_sum_exp(values : Array[Double]) -> Result[Double, DetectError] {
  if values.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  let mut maximum = negative_infinity
  for value in values {
    if !finite(value) {
      return Err(NonFiniteObservation)
    }
    if value > maximum {
      maximum = value
    }
  }
  let mut total = 0.0
  for value in values {
    total = total + @math.exp(value - maximum)
  }
  Ok(maximum + @math.ln(total))
}

///|
pub fn log1mexp(value : Double) -> Result[Double, DetectError] {
  if !finite(value) || value >= 0.0 {
    return Err(InvalidHazard)
  }
  if value < -0.6931471805599453 {
    Ok(@math.ln(1.0 - @math.exp(value)))
  } else {
    Ok(@math.ln(-value) + @math.ln(1.0 - value / 2.0 + value * value / 6.0))
  }
}

///|
/// Lanczos approximation for log Gamma on positive arguments.
pub fn log_gamma(value : Double) -> Result[Double, DetectError] {
  if !finite(value) || value <= 0.0 {
    return Err(InvalidVariance)
  }
  if value < 0.5 {
    let reflected = match log_gamma(1.0 - value) {
      Ok(result) => result
      Err(error) => return Err(error)
    }
    return Ok(
      @math.ln(3.141592653589793) -
      @math.ln(@math.sin(3.141592653589793 * value)) -
      reflected,
    )
  }
  let coefficients = [
    0.9999999999998099, 676.5203681218851, -1259.1392167224028, 771.3234287776531,
    -176.6150291621406, 12.507343278686905, -0.13857109526572012, 0.000009984369578019572,
    0.00000015056327351493116,
  ]
  let shifted = value - 1.0
  let mut accumulator = coefficients[0]
  for index in 1..<coefficients.length() {
    accumulator = accumulator +
      coefficients[index] / (shifted + index.to_double())
  }
  let t = shifted + 7.5
  Ok(
    0.9189385332046727 +
    (shifted + 0.5) * @math.ln(t) -
    t +
    @math.ln(accumulator),
  )
}

///|
pub fn log_factorial(value : Int) -> Result[Double, DetectError] {
  if value < 0 {
    return Err(InvalidRunLengthLimit)
  }
  log_gamma(value.to_double() + 1.0)
}

///|
pub fn mean(values : Array[Double]) -> Result[Double, DetectError] {
  if values.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  let mut total = 0.0
  for value in values {
    if !finite(value) {
      return Err(NonFiniteObservation)
    }
    total = total + value
  }
  Ok(total / values.length().to_double())
}

///|
pub fn sample_variance(values : Array[Double]) -> Result[Double, DetectError] {
  if values.length() < 2 {
    return Err(InvalidRunLengthLimit)
  }
  let average = match mean(values) {
    Ok(result) => result
    Err(error) => return Err(error)
  }
  let mut total = 0.0
  for value in values {
    let difference = value - average
    total = total + difference * difference
  }
  Ok(total / (values.length() - 1).to_double())
}

///|
pub fn z_score(
  value : Double,
  center : Double,
  variance : Double,
) -> Result[Double, DetectError] {
  if !finite(value) || !finite(center) {
    return Err(NonFiniteObservation)
  }
  if !finite(variance) || variance <= 0.0 {
    return Err(InvalidVariance)
  }
  Ok((value - center) / @math.exp(@math.ln(variance) * 0.5))
}
