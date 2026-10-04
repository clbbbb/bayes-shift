///|
fn finite(value : Double) -> Bool {
  value == value && value != 1.0 / 0.0 && value != -1.0 / 0.0
}

///|
let negative_infinity : Double = -1.0 / 0.0

///|
fn log_add(left : Double, right : Double) -> Double {
  if left == negative_infinity {
    return right
  }
  if right == negative_infinity {
    return left
  }
  let high = if left > right { left } else { right }
  let low = if left > right { right } else { left }
  high + @math.ln(1.0 + @math.exp(low - high))
}

///|
fn normalize_log_weights(weights : Array[Double]) -> Unit {
  let mut maximum = negative_infinity
  for value in weights {
    if value > maximum {
      maximum = value
    }
  }
  let mut total = 0.0
  for value in weights {
    if value != negative_infinity {
      total = total + @math.exp(value - maximum)
    }
  }
  let evidence = maximum + @math.ln(total)
  for index in 0..<weights.length() {
    weights[index] = weights[index] - evidence
  }
}

///|
fn normal_log_density(
  value : Double,
  mean : Double,
  variance : Double,
) -> Double {
  let difference = value - mean
  -0.5 *
  (
    @math.ln(2.0 * 3.141592653589793 * variance) +
    difference * difference / variance
  )
}
