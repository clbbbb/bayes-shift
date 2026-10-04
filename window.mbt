///|
/// Fixed-capacity numeric window with O(1) mean updates and explicit quality state.
pub struct Window {
  capacity : Int
  mut values : Array[Double]
  mut sum : Double
  mut sum_squares : Double
}

///|
pub fn Window::new(capacity : Int) -> Result[Window, DetectError] {
  if capacity < 1 {
    return Err(InvalidRunLengthLimit)
  }
  Ok({ capacity, values: [], sum: 0.0, sum_squares: 0.0, })
}

///|
pub fn Window::length(self : Window) -> Int {
  self.values.length()
}

///|
pub fn Window::is_full(self : Window) -> Bool {
  self.values.length() == self.capacity
}

///|
pub fn Window::values(self : Window) -> Array[Double] {
  self.values.copy()
}

///|
/// Add one finite value and evict the oldest value after reaching capacity.
pub fn Window::push(self : Window, value : Double) -> Result[Unit, DetectError] {
  if !finite(value) {
    return Err(NonFiniteObservation)
  }
  if self.values.length() == self.capacity {
    let evicted = self.values.remove(0)
    self.sum = self.sum - evicted
    self.sum_squares = self.sum_squares - evicted * evicted
  }
  self.values.push(value)
  self.sum = self.sum + value
  self.sum_squares = self.sum_squares + value * value
  Ok(())
}

///|
pub fn Window::mean(self : Window) -> Double? {
  if self.values.length() == 0 {
    None
  } else {
    Some(self.sum / self.values.length().to_double())
  }
}

///|
/// Population variance; useful for short online windows where a sample correction is not wanted.
pub fn Window::variance(self : Window) -> Double? {
  if self.values.length() == 0 {
    return None
  }
  let size = self.values.length().to_double()
  let variance = self.sum_squares / size - self.sum / size * (self.sum / size)
  Some(if variance > 0.0 { variance } else { 0.0 })
}

///|
pub fn Window::minimum(self : Window) -> Double? {
  if self.values.length() == 0 {
    return None
  }
  let mut result = self.values[0]
  for value in self.values {
    if value < result {
      result = value
    }
  }
  Some(result)
}

///|
pub fn Window::maximum(self : Window) -> Double? {
  if self.values.length() == 0 {
    return None
  }
  let mut result = self.values[0]
  for value in self.values {
    if value > result {
      result = value
    }
  }
  Some(result)
}

///|
pub fn Window::clear(self : Window) -> Unit {
  self.values = []
  self.sum = 0.0
  self.sum_squares = 0.0
}

///|
/// Online quality facts for a numeric sequence.
pub struct Quality {
  values : Int
  non_finite : Int
  minimum : Double?
  maximum : Double?
  monotonic_increasing : Bool
}

///|
pub fn inspect_values(values : Array[Double]) -> Quality {
  let mut non_finite = 0
  let mut minimum : Double? = None
  let mut maximum : Double? = None
  let mut previous : Double? = None
  let mut monotonic = true
  for value in values {
    if !finite(value) {
      non_finite = non_finite + 1
    } else {
      match minimum {
        None => minimum = Some(value)
        Some(current) => if value < current { minimum = Some(value) }
      }
      match maximum {
        None => maximum = Some(value)
        Some(current) => if value > current { maximum = Some(value) }
      }
      match previous {
        None => ()
        Some(current) => if value < current { monotonic = false }
      }
      previous = Some(value)
    }
  }
  {
    values: values.length(),
    non_finite,
    minimum,
    maximum,
    monotonic_increasing: monotonic,
  }
}
