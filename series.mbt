///|
/// An ordered, validated numeric sequence used by offline evaluators.
pub struct NumericSeries {
  values : Array[Double]
}

///|
pub fn NumericSeries::new(
  values : Array[Double],
) -> Result[NumericSeries, DetectError] {
  if values.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  for value in values {
    if !finite(value) {
      return Err(NonFiniteObservation)
    }
  }
  Ok({ values: values.copy(), })
}

///|
pub fn NumericSeries::values(self : NumericSeries) -> Array[Double] {
  self.values.copy()
}

///|
pub fn NumericSeries::length(self : NumericSeries) -> Int {
  self.values.length()
}

///|
pub fn NumericSeries::at(self : NumericSeries, index : Int) -> Double? {
  self.values.get(index)
}

///|
pub fn NumericSeries::first(self : NumericSeries) -> Double {
  self.values[0]
}

///|
pub fn NumericSeries::last(self : NumericSeries) -> Double {
  self.values[self.values.length() - 1]
}

///|
pub fn NumericSeries::difference(self : NumericSeries) -> Array[Double] {
  let result = []
  for index in 1..<self.values.length() {
    result.push(self.values[index] - self.values[index - 1])
  }
  result
}

///|
pub fn NumericSeries::cumulative_sum(self : NumericSeries) -> Array[Double] {
  let result = []
  let mut total = 0.0
  for value in self.values {
    total = total + value
    result.push(total)
  }
  result
}

///|
pub fn NumericSeries::centered(
  self : NumericSeries,
) -> Result[Array[Double], DetectError] {
  let average = match mean(self.values) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let result = []
  for value in self.values {
    result.push(value - average)
  }
  Ok(result)
}

///|
pub fn NumericSeries::standardized(
  self : NumericSeries,
) -> Result[Array[Double], DetectError] {
  let average = match mean(self.values) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let variance = match sample_variance(self.values) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  if variance == 0.0 {
    return Ok(self.centered().unwrap())
  }
  let scale = @math.exp(@math.ln(variance) * 0.5)
  let result = []
  for value in self.values {
    result.push((value - average) / scale)
  }
  Ok(result)
}

///|
pub fn NumericSeries::moving_means(
  self : NumericSeries,
  width : Int,
) -> Result[Array[Double], DetectError] {
  if width <= 0 || width > self.values.length() {
    return Err(InvalidRunLengthLimit)
  }
  let result = []
  let mut sum = 0.0
  for index in 0..<self.values.length() {
    sum = sum + self.values[index]
    if index >= width {
      sum = sum - self.values[index - width]
    }
    if index + 1 >= width {
      result.push(sum / width.to_double())
    }
  }
  Ok(result)
}

///|
pub fn NumericSeries::moving_ranges(
  self : NumericSeries,
  width : Int,
) -> Result[Array[Double], DetectError] {
  if width <= 0 || width > self.values.length() {
    return Err(InvalidRunLengthLimit)
  }
  let result = []
  for start in 0..<(self.values.length() - width + 1) {
    let mut smallest = self.values[start]
    let mut largest = self.values[start]
    for offset in 1..<width {
      let value = self.values[start + offset]
      if value < smallest {
        smallest = value
      }
      if value > largest {
        largest = value
      }
    }
    result.push(largest - smallest)
  }
  Ok(result)
}

///|
pub fn NumericSeries::subseries(
  self : NumericSeries,
  start : Int,
  finish : Int,
) -> Result[NumericSeries, DetectError] {
  if start < 0 || finish <= start || finish > self.values.length() {
    return Err(InvalidRunLengthLimit)
  }
  let result = []
  for index in start..<finish {
    result.push(self.values[index])
  }
  NumericSeries::new(result)
}

///|
pub fn NumericSeries::split(
  self : NumericSeries,
  width : Int,
) -> Result[Array[NumericSeries], DetectError] {
  if width <= 0 {
    return Err(InvalidRunLengthLimit)
  }
  let result = []
  let mut start = 0
  while start < self.values.length() {
    let mut finish = start + width
    if finish > self.values.length() {
      finish = self.values.length()
    }
    let part = match self.subseries(start, finish) {
      Ok(value) => value
      Err(error) => return Err(error)
    }
    result.push(part)
    start = finish
  }
  Ok(result)
}

///|
pub fn NumericSeries::rank_indices(self : NumericSeries) -> Array[Int] {
  let indices = []
  for index in 0..<self.values.length() {
    indices.push(index)
  }
  indices.sort_by(fn(left, right) {
    if self.values[left] < self.values[right] {
      -1
    } else if self.values[left] > self.values[right] {
      1
    } else {
      left - right
    }
  })
  indices
}

///|
pub fn concatenate_series(
  series : Array[NumericSeries],
) -> Result[NumericSeries, DetectError] {
  let result = []
  for item in series {
    for value in item.values {
      result.push(value)
    }
  }
  NumericSeries::new(result)
}

///|
pub struct TimedValue {
  time : Int
  value : Double
}

///|
pub fn TimedValue::new(
  time : Int,
  value : Double,
) -> Result[TimedValue, DetectError] {
  if time < 0 {
    return Err(InvalidRunLengthLimit)
  }
  if !finite(value) {
    return Err(NonFiniteObservation)
  }
  Ok({ time, value, })
}

///|
pub struct TimedSeries {
  points : Array[TimedValue]
}

///|
pub fn TimedSeries::new(
  points : Array[TimedValue],
) -> Result[TimedSeries, DetectError] {
  if points.length() == 0 {
    return Err(InvalidRunLengthLimit)
  }
  let mut previous = -1
  for point in points {
    if point.time <= previous {
      return Err(InvalidRunLengthLimit)
    }
    previous = point.time
  }
  Ok({ points: points.copy(), })
}

///|
pub fn TimedSeries::points(self : TimedSeries) -> Array[TimedValue] {
  self.points.copy()
}

///|
pub fn TimedSeries::values(self : TimedSeries) -> Array[Double] {
  let values = []
  for point in self.points {
    values.push(point.value)
  }
  values
}

///|
pub fn TimedSeries::intervals(self : TimedSeries) -> Array[Int] {
  let result = []
  for index in 1..<self.points.length() {
    result.push(self.points[index].time - self.points[index - 1].time)
  }
  result
}

///|
pub fn TimedSeries::is_regular(self : TimedSeries) -> Bool {
  if self.points.length() < 3 {
    return true
  }
  let expected = self.points[1].time - self.points[0].time
  for index in 2..<self.points.length() {
    if self.points[index].time - self.points[index - 1].time != expected {
      return false
    }
  }
  true
}

///|
pub fn TimedSeries::numeric(self : TimedSeries) -> NumericSeries {
  NumericSeries::new(self.values()).unwrap()
}
