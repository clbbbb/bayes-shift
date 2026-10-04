///|
/// Welford running statistics for a single numeric stream.
pub struct OnlineStats {
  mut count : Int
  mut average : Double
  mut sum_of_squares : Double
  mut minimum_value : Double?
  mut maximum_value : Double?
}

///|
pub fn OnlineStats::new() -> OnlineStats {
  {
    count: 0,
    average: 0.0,
    sum_of_squares: 0.0,
    minimum_value: None,
    maximum_value: None,
  }
}

///|
pub fn OnlineStats::count(self : OnlineStats) -> Int {
  self.count
}

///|
pub fn OnlineStats::mean(self : OnlineStats) -> Double? {
  if self.count == 0 {
    None
  } else {
    Some(self.average)
  }
}

///|
pub fn OnlineStats::population_variance(self : OnlineStats) -> Double? {
  if self.count == 0 {
    None
  } else {
    Some(self.sum_of_squares / self.count.to_double())
  }
}

///|
pub fn OnlineStats::sample_variance(self : OnlineStats) -> Double? {
  if self.count < 2 {
    None
  } else {
    Some(self.sum_of_squares / (self.count - 1).to_double())
  }
}

///|
pub fn OnlineStats::minimum(self : OnlineStats) -> Double? {
  self.minimum_value
}

///|
pub fn OnlineStats::maximum(self : OnlineStats) -> Double? {
  self.maximum_value
}

///|
pub fn OnlineStats::push(
  self : OnlineStats,
  value : Double,
) -> Result[Unit, DetectError] {
  if !finite(value) {
    return Err(NonFiniteObservation)
  }
  self.count = self.count + 1
  let delta = value - self.average
  self.average = self.average + delta / self.count.to_double()
  self.sum_of_squares = self.sum_of_squares + delta * (value - self.average)
  match self.minimum_value {
    None => self.minimum_value = Some(value)
    Some(current) => if value < current { self.minimum_value = Some(value) }
  }
  match self.maximum_value {
    None => self.maximum_value = Some(value)
    Some(current) => if value > current { self.maximum_value = Some(value) }
  }
  Ok(())
}

///|
pub fn OnlineStats::push_all(
  self : OnlineStats,
  values : Array[Double],
) -> Result[Unit, DetectError] {
  for value in values {
    match self.push(value) {
      Ok(_) => ()
      Err(error) => return Err(error)
    }
  }
  Ok(())
}

///|
pub fn OnlineStats::z_score(
  self : OnlineStats,
  value : Double,
) -> Result[Double, DetectError] {
  if !finite(value) {
    return Err(NonFiniteObservation)
  }
  let variance = match self.sample_variance() {
    Some(result) => result
    None => return Err(InvalidVariance)
  }
  if variance <= 0.0 {
    return Err(InvalidVariance)
  }
  Ok((value - self.average) / @math.exp(@math.ln(variance) * 0.5))
}

///|
pub fn OnlineStats::merge(
  self : OnlineStats,
  other : OnlineStats,
) -> OnlineStats {
  if self.count == 0 {
    return other
  }
  if other.count == 0 {
    return self
  }
  let total = self.count + other.count
  let difference = other.average - self.average
  let correction = difference *
    difference *
    self.count.to_double() *
    other.count.to_double() /
    total.to_double()
  let minimum_value = match (self.minimum_value, other.minimum_value) {
    (Some(left), Some(right)) =>
      if left < right {
        Some(left)
      } else {
        Some(right)
      }
    (Some(value), None) => Some(value)
    (None, Some(value)) => Some(value)
    (None, None) => None
  }
  let maximum_value = match (self.maximum_value, other.maximum_value) {
    (Some(left), Some(right)) =>
      if left > right {
        Some(left)
      } else {
        Some(right)
      }
    (Some(value), None) => Some(value)
    (None, Some(value)) => Some(value)
    (None, None) => None
  }
  {
    count: total,
    average: self.average +
    difference * other.count.to_double() / total.to_double(),
    sum_of_squares: self.sum_of_squares + other.sum_of_squares + correction,
    minimum_value,
    maximum_value,
  }
}

///|
pub fn OnlineStats::reset(self : OnlineStats) -> Unit {
  self.count = 0
  self.average = 0.0
  self.sum_of_squares = 0.0
  self.minimum_value = None
  self.maximum_value = None
}

///|
pub fn OnlineStats::copy(self : OnlineStats) -> OnlineStats {
  {
    count: self.count,
    average: self.average,
    sum_of_squares: self.sum_of_squares,
    minimum_value: self.minimum_value,
    maximum_value: self.maximum_value,
  }
}

///|
pub fn stats_of(values : Array[Double]) -> Result[OnlineStats, DetectError] {
  let stats = OnlineStats::new()
  match stats.push_all(values) {
    Ok(_) => Ok(stats)
    Err(error) => Err(error)
  }
}

///|
pub struct RunningMoments {
  mut count : Int
  mut mean_value : Double
  mut m2 : Double
  mut m3 : Double
  mut m4 : Double
}

///|
pub fn RunningMoments::new() -> RunningMoments {
  { count: 0, mean_value: 0.0, m2: 0.0, m3: 0.0, m4: 0.0, }
}

///|
pub fn RunningMoments::push(
  self : RunningMoments,
  value : Double,
) -> Result[Unit, DetectError] {
  if !finite(value) {
    return Err(NonFiniteObservation)
  }
  let old_count = self.count.to_double()
  self.count = self.count + 1
  let count = self.count.to_double()
  let delta = value - self.mean_value
  let delta_over_count = delta / count
  let term = delta * delta_over_count * old_count
  self.m4 = self.m4 +
    term *
    delta_over_count *
    delta_over_count *
    (count * count - 3.0 * count + 3.0) +
    6.0 * delta_over_count * delta_over_count * self.m2 -
    4.0 * delta_over_count * self.m3
  self.m3 = self.m3 +
    term * delta_over_count * (count - 2.0) -
    3.0 * delta_over_count * self.m2
  self.m2 = self.m2 + term
  self.mean_value = self.mean_value + delta_over_count
  Ok(())
}

///|
pub fn RunningMoments::skewness(self : RunningMoments) -> Double? {
  if self.count < 3 || self.m2 == 0.0 {
    None
  } else {
    Some(
      @math.exp(@math.ln(self.count.to_double()) * 0.5) *
      self.m3 /
      @math.exp(@math.ln(self.m2) * 1.5),
    )
  }
}

///|
pub fn RunningMoments::excess_kurtosis(self : RunningMoments) -> Double? {
  if self.count < 4 || self.m2 == 0.0 {
    None
  } else {
    Some(self.count.to_double() * self.m4 / (self.m2 * self.m2) - 3.0)
  }
}
