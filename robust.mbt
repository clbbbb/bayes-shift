///|
/// Copy values into ascending order after validating all observations.
pub fn sorted_values(
  values : Array[Double],
) -> Result[Array[Double], DetectError] {
  if values.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  let result = values.copy()
  for value in result {
    if !finite(value) {
      return Err(NonFiniteObservation)
    }
  }
  result.sort_by(fn(left, right) {
    if left < right {
      -1
    } else if left > right {
      1
    } else {
      0
    }
  })
  Ok(result)
}

///|
pub fn quantile(
  values : Array[Double],
  probability : Double,
) -> Result[Double, DetectError] {
  if !finite(probability) || probability < 0.0 || probability > 1.0 {
    return Err(InvalidHazard)
  }
  let sorted = match sorted_values(values) {
    Ok(result) => result
    Err(error) => return Err(error)
  }
  if sorted.length() == 1 {
    return Ok(sorted[0])
  }
  let position = probability * (sorted.length() - 1).to_double()
  let lower = position.to_int()
  let upper = if lower + 1 < sorted.length() { lower + 1 } else { lower }
  let fraction = position - lower.to_double()
  Ok(sorted[lower] + (sorted[upper] - sorted[lower]) * fraction)
}

///|
pub fn median(values : Array[Double]) -> Result[Double, DetectError] {
  quantile(values, 0.5)
}

///|
pub fn interquartile_range(
  values : Array[Double],
) -> Result[Double, DetectError] {
  let lower = match quantile(values, 0.25) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let upper = match quantile(values, 0.75) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  Ok(upper - lower)
}

///|
pub fn median_absolute_deviation(
  values : Array[Double],
) -> Result[Double, DetectError] {
  let center = match median(values) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let deviations = []
  for value in values {
    let difference = value - center
    deviations.push(if difference < 0.0 { -difference } else { difference })
  }
  median(deviations)
}

///|
pub fn robust_scale(values : Array[Double]) -> Result[Double, DetectError] {
  let deviation = match median_absolute_deviation(values) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  Ok(deviation * 1.4826)
}

///|
pub fn robust_z_scores(
  values : Array[Double],
) -> Result[Array[Double], DetectError] {
  let center = match median(values) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let scale = match robust_scale(values) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let result = []
  if scale == 0.0 {
    for _ in values {
      result.push(0.0)
    }
    return Ok(result)
  }
  for value in values {
    result.push((value - center) / scale)
  }
  Ok(result)
}

///|
pub fn winsorize(
  values : Array[Double],
  lower_probability : Double,
  upper_probability : Double,
) -> Result[Array[Double], DetectError] {
  if lower_probability < 0.0 ||
    upper_probability > 1.0 ||
    lower_probability > upper_probability {
    return Err(InvalidHazard)
  }
  let lower = match quantile(values, lower_probability) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let upper = match quantile(values, upper_probability) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let result = []
  for value in values {
    if value < lower {
      result.push(lower)
    } else if value > upper {
      result.push(upper)
    } else {
      result.push(value)
    }
  }
  Ok(result)
}

///|
pub fn trimmed_mean(
  values : Array[Double],
  trim_fraction : Double,
) -> Result[Double, DetectError] {
  if !finite(trim_fraction) || trim_fraction < 0.0 || trim_fraction >= 0.5 {
    return Err(InvalidHazard)
  }
  let sorted = match sorted_values(values) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let trim = (sorted.length().to_double() * trim_fraction).to_int()
  if trim * 2 >= sorted.length() {
    return Err(InvalidRunLengthLimit)
  }
  let mut total = 0.0
  for index in trim..<(sorted.length() - trim) {
    total = total + sorted[index]
  }
  Ok(total / (sorted.length() - trim * 2).to_double())
}

///|
pub fn rolling_median(
  values : Array[Double],
  width : Int,
) -> Result[Array[Double], DetectError] {
  if width <= 0 || width > values.length() {
    return Err(InvalidRunLengthLimit)
  }
  let result = []
  for start in 0..<(values.length() - width + 1) {
    let window = []
    for offset in 0..<width {
      window.push(values[start + offset])
    }
    let value = match median(window) {
      Ok(result) => result
      Err(error) => return Err(error)
    }
    result.push(value)
  }
  Ok(result)
}

///|
pub fn huber_weight(
  residual : Double,
  cutoff : Double,
) -> Result[Double, DetectError] {
  if !finite(residual) || !finite(cutoff) || cutoff <= 0.0 {
    return Err(InvalidVariance)
  }
  let magnitude = if residual < 0.0 { -residual } else { residual }
  if magnitude <= cutoff {
    Ok(1.0)
  } else {
    Ok(cutoff / magnitude)
  }
}

///|
pub fn huber_loss(
  residual : Double,
  cutoff : Double,
) -> Result[Double, DetectError] {
  if !finite(residual) || !finite(cutoff) || cutoff <= 0.0 {
    return Err(InvalidVariance)
  }
  let magnitude = if residual < 0.0 { -residual } else { residual }
  if magnitude <= cutoff {
    Ok(0.5 * residual * residual)
  } else {
    Ok(cutoff * (magnitude - 0.5 * cutoff))
  }
}

///|
pub fn theil_sen_slope(values : Array[Double]) -> Result[Double, DetectError] {
  if values.length() < 2 {
    return Err(InvalidRunLengthLimit)
  }
  let slopes = []
  for left in 0..<(values.length() - 1) {
    for right in (left + 1)..<values.length() {
      slopes.push((values[right] - values[left]) / (right - left).to_double())
    }
  }
  median(slopes)
}

///|
pub fn detrend_linear(
  values : Array[Double],
) -> Result[Array[Double], DetectError] {
  let slope = match theil_sen_slope(values) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let intercept = match median(values) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let center = (values.length() - 1).to_double() / 2.0
  let result = []
  for index in 0..<values.length() {
    result.push(
      values[index] - (intercept + slope * (index.to_double() - center)),
    )
  }
  Ok(result)
}
