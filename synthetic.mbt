///|
/// A small deterministic PRNG used only for reproducible fixtures and examples.
pub struct ParkMiller {
  mut state : Int
}

///|
pub fn ParkMiller::new(seed : Int) -> Result[ParkMiller, DetectError] {
  if seed <= 0 || seed >= 2147483647 {
    return Err(InvalidRunLengthLimit)
  }
  Ok({ state: seed, })
}

///|
pub fn ParkMiller::next_unit(self : ParkMiller) -> Double {
  self.state = self.state * 48271 % 2147483647
  self.state.to_double() / 2147483647.0
}

///|
/// Generate a piecewise-constant series with deterministic uniform jitter.
pub fn piecewise_constant(
  lengths : Array[Int],
  means : Array[Double],
  half_width : Double,
  seed : Int,
) -> Result[Array[Double], DetectError] {
  if lengths.length() == 0 ||
    lengths.length() != means.length() ||
    !finite(half_width) ||
    half_width < 0.0 {
    return Err(InvalidRunLengthLimit)
  }
  let rng = ParkMiller::new(seed).unwrap()
  let result : Array[Double] = []
  for index in 0..<lengths.length() {
    if lengths[index] <= 0 || !finite(means[index]) {
      return Err(InvalidRunLengthLimit)
    }
    for _ in 0..<lengths[index] {
      result.push(means[index] + (rng.next_unit() * 2.0 - 1.0) * half_width)
    }
  }
  Ok(result)
}
