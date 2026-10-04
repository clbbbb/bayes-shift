///|
/// Description of a constant segment for generated numeric traces.
pub struct SegmentSpec {
  length : Int
  level : Double
}

///|
pub fn SegmentSpec::new(
  length : Int,
  level : Double,
) -> Result[SegmentSpec, DetectError] {
  if length <= 0 {
    return Err(InvalidRunLengthLimit)
  }
  if !finite(level) {
    return Err(NonFiniteObservation)
  }
  Ok({ length, level, })
}

///|
pub fn segment_boundaries(specs : Array[SegmentSpec]) -> Array[Boundary] {
  let result = []
  let mut accumulated = 0
  for index in 0..<specs.length() {
    accumulated = accumulated + specs[index].length
    if index + 1 < specs.length() {
      result.push(Boundary::new(accumulated + 1).unwrap())
    }
  }
  result
}

///|
pub fn uniform_noise(
  length : Int,
  lower : Double,
  upper : Double,
  seed : Int,
) -> Result[Array[Double], DetectError] {
  if length <= 0 || !finite(lower) || !finite(upper) || lower >= upper {
    return Err(InvalidVariance)
  }
  let rng = ParkMiller::new(seed).unwrap()
  let result = []
  for _ in 0..<length {
    result.push(lower + (upper - lower) * rng.next_unit())
  }
  Ok(result)
}

///|
/// Box-Muller normal samples, deterministic for a fixed seed.
pub fn normal_noise(
  length : Int,
  center : Double,
  scale : Double,
  seed : Int,
) -> Result[Array[Double], DetectError] {
  if length <= 0 || !finite(center) || !finite(scale) || scale < 0.0 {
    return Err(InvalidVariance)
  }
  let rng = ParkMiller::new(seed).unwrap()
  let result = []
  while result.length() < length {
    let mut uniform = rng.next_unit()
    if uniform < 0.000000001 {
      uniform = 0.000000001
    }
    let radius = @math.exp(-@math.ln(uniform) * 0.5 * 2.0)
    let angle = 6.283185307179586 * rng.next_unit()
    result.push(center + scale * radius * @math.cos(angle))
    if result.length() < length {
      result.push(center + scale * radius * @math.sin(angle))
    }
  }
  Ok(result)
}

///|
pub fn gaussian_segments(
  specs : Array[SegmentSpec],
  scale : Double,
  seed : Int,
) -> Result[NumericSeries, DetectError] {
  if specs.length() == 0 || scale < 0.0 || !finite(scale) {
    return Err(InvalidVariance)
  }
  let values = []
  let mut offset = 0
  for spec in specs {
    let noise = match
      normal_noise(spec.length, spec.level, scale, seed + offset + 1) {
      Ok(result) => result
      Err(error) => return Err(error)
    }
    for value in noise {
      values.push(value)
    }
    offset = offset + spec.length
  }
  NumericSeries::new(values)
}

///|
pub fn linear_ramp(
  length : Int,
  start : Double,
  finish : Double,
) -> Result[NumericSeries, DetectError] {
  if length < 2 || !finite(start) || !finite(finish) {
    return Err(InvalidRunLengthLimit)
  }
  let values = []
  for index in 0..<length {
    values.push(
      start + (finish - start) * index.to_double() / (length - 1).to_double(),
    )
  }
  NumericSeries::new(values)
}

///|
pub fn sinusoid(
  length : Int,
  center : Double,
  amplitude : Double,
  period : Double,
) -> Result[NumericSeries, DetectError] {
  if length <= 0 ||
    !finite(center) ||
    !finite(amplitude) ||
    !finite(period) ||
    period <= 0.0 {
    return Err(InvalidRunLengthLimit)
  }
  let values = []
  for index in 0..<length {
    values.push(
      center +
      amplitude * @math.sin(6.283185307179586 * index.to_double() / period),
    )
  }
  NumericSeries::new(values)
}

///|
pub fn add_series(
  left : NumericSeries,
  right : NumericSeries,
) -> Result[NumericSeries, DetectError] {
  if left.length() != right.length() {
    return Err(InvalidRunLengthLimit)
  }
  let values = []
  let left_values = left.values()
  let right_values = right.values()
  for index in 0..<left.length() {
    values.push(left_values[index] + right_values[index])
  }
  NumericSeries::new(values)
}

///|
pub fn add_normal_noise(
  series : NumericSeries,
  scale : Double,
  seed : Int,
) -> Result[NumericSeries, DetectError] {
  let noise = match normal_noise(series.length(), 0.0, scale, seed) {
    Ok(result) => result
    Err(error) => return Err(error)
  }
  add_series(series, NumericSeries::new(noise).unwrap())
}

///|
pub fn inject_spikes(
  series : NumericSeries,
  positions : Array[Int],
  amount : Double,
) -> Result[NumericSeries, DetectError] {
  if !finite(amount) {
    return Err(NonFiniteObservation)
  }
  let values = series.values()
  for position in positions {
    if position < 0 || position >= values.length() {
      return Err(InvalidRunLengthLimit)
    }
    values[position] = values[position] + amount
  }
  NumericSeries::new(values)
}

///|
pub fn clip_series(
  series : NumericSeries,
  lower : Double,
  upper : Double,
) -> Result[NumericSeries, DetectError] {
  if !finite(lower) || !finite(upper) || lower > upper {
    return Err(InvalidVariance)
  }
  let values = []
  for value in series.values() {
    if value < lower {
      values.push(lower)
    } else if value > upper {
      values.push(upper)
    } else {
      values.push(value)
    }
  }
  NumericSeries::new(values)
}

///|
pub fn poisson_sample(
  rate : Double,
  rng : ParkMiller,
) -> Result[Int, DetectError] {
  if !finite(rate) || rate < 0.0 || rate > 100.0 {
    return Err(InvalidVariance)
  }
  let cutoff = @math.exp(-rate)
  let mut product = 1.0
  let mut count = 0
  while product > cutoff {
    count = count + 1
    product = product * rng.next_unit()
  }
  Ok(count - 1)
}

///|
pub fn poisson_segments(
  lengths : Array[Int],
  rates : Array[Double],
  seed : Int,
) -> Result[Array[Int], DetectError] {
  if lengths.length() == 0 || lengths.length() != rates.length() {
    return Err(InvalidRunLengthLimit)
  }
  let rng = ParkMiller::new(seed).unwrap()
  let values = []
  for index in 0..<lengths.length() {
    if lengths[index] <= 0 || rates[index] < 0.0 || !finite(rates[index]) {
      return Err(InvalidVariance)
    }
    for _ in 0..<lengths[index] {
      let sample = match poisson_sample(rates[index], rng) {
        Ok(value) => value
        Err(error) => return Err(error)
      }
      values.push(sample)
    }
  }
  Ok(values)
}

///|
pub fn contamination_mask(
  length : Int,
  probability : Double,
  seed : Int,
) -> Result[Array[Bool], DetectError] {
  if length <= 0 ||
    !finite(probability) ||
    probability < 0.0 ||
    probability > 1.0 {
    return Err(InvalidHazard)
  }
  let rng = ParkMiller::new(seed).unwrap()
  let result = []
  for _ in 0..<length {
    result.push(rng.next_unit() < probability)
  }
  Ok(result)
}
