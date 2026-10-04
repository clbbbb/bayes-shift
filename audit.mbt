///|
/// Input audit for a numeric stream before it is handed to a detector.
pub struct SeriesAudit {
  length : Int
  finite_values : Int
  non_finite_values : Int
  repeated_adjacent_values : Int
  increasing_steps : Int
  decreasing_steps : Int
  zero_steps : Int
  minimum : Double?
  maximum : Double?
}

///|
pub fn audit_series(values : Array[Double]) -> SeriesAudit {
  let mut finite_values = 0
  let mut non_finite_values = 0
  let mut repeated = 0
  let mut increasing = 0
  let mut decreasing = 0
  let mut zero = 0
  let mut minimum : Double? = None
  let mut maximum : Double? = None
  let mut previous : Double? = None
  for value in values {
    if !finite(value) {
      non_finite_values = non_finite_values + 1
    } else {
      finite_values = finite_values + 1
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
        Some(current) =>
          if value > current {
            increasing = increasing + 1
          } else if value < current {
            decreasing = decreasing + 1
          } else {
            repeated = repeated + 1
            zero = zero + 1
          }
      }
      previous = Some(value)
    }
  }
  {
    length: values.length(),
    finite_values,
    non_finite_values,
    repeated_adjacent_values: repeated,
    increasing_steps: increasing,
    decreasing_steps: decreasing,
    zero_steps: zero,
    minimum,
    maximum,
  }
}

///|
pub fn SeriesAudit::is_usable(self : SeriesAudit) -> Bool {
  self.length > 0 && self.finite_values == self.length
}

///|
pub fn SeriesAudit::is_constant(self : SeriesAudit) -> Bool {
  self.finite_values > 0 &&
  self.increasing_steps == 0 &&
  self.decreasing_steps == 0
}

///|
pub fn SeriesAudit::monotonic_direction(self : SeriesAudit) -> String {
  if self.increasing_steps > 0 && self.decreasing_steps == 0 {
    "increasing"
  } else if self.decreasing_steps > 0 && self.increasing_steps == 0 {
    "decreasing"
  } else if self.increasing_steps == 0 && self.decreasing_steps == 0 {
    "constant"
  } else {
    "mixed"
  }
}

///|
pub fn SeriesAudit::range(self : SeriesAudit) -> Double? {
  match (self.minimum, self.maximum) {
    (Some(lower), Some(upper)) => Some(upper - lower)
    _ => None
  }
}

///|
pub fn SeriesAudit::summary(self : SeriesAudit) -> String {
  "length=" +
  self.length.to_string() +
  ",finite=" +
  self.finite_values.to_string() +
  ",non_finite=" +
  self.non_finite_values.to_string() +
  ",repeated_adjacent=" +
  self.repeated_adjacent_values.to_string() +
  ",direction=" +
  self.monotonic_direction()
}

///|
pub fn validate_audit(audit : SeriesAudit) -> Result[Unit, DetectError] {
  if audit.length == 0 {
    return Err(InvalidRunLengthLimit)
  }
  if audit.non_finite_values > 0 {
    return Err(NonFiniteObservation)
  }
  Ok(())
}

///|
pub fn audit_difference(
  values : Array[Double],
) -> Result[SeriesAudit, DetectError] {
  let series = match NumericSeries::new(values) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  Ok(audit_series(series.difference()))
}

///|
pub struct RunLengthAudit {
  observed_steps : Int
  retained_states : Int
  truncation_ratio : Double
}

///|
pub fn audit_run_length_capacity(
  updates : Array[Update],
  max_run_length : Int,
) -> Result[RunLengthAudit, DetectError] {
  if updates.length() == 0 || max_run_length < 2 {
    return Err(InvalidRunLengthLimit)
  }
  let mut retained = 0
  for update in updates {
    retained = retained + update.retained_run_lengths
  }
  Ok({
    observed_steps: updates.length(),
    retained_states: retained,
    truncation_ratio: retained.to_double() /
    (updates.length() * max_run_length).to_double(),
  })
}

///|
pub fn SeriesAudit::quality_score(self : SeriesAudit) -> Int {
  if self.length == 0 {
    return 0
  }
  let mut score = 100
  score = score - self.non_finite_values * 100 / self.length
  if self.is_constant() {
    score = score - 20
  }
  if self.repeated_adjacent_values * 2 > self.length {
    score = score - 10
  }
  if score < 0 {
    0
  } else {
    score
  }
}

///|
pub fn SeriesAudit::recommendation(self : SeriesAudit) -> String {
  if self.length == 0 {
    "provide at least one observation"
  } else if self.non_finite_values > 0 {
    "remove or impute non-finite observations before detection"
  } else if self.is_constant() {
    "constant series is valid but cannot demonstrate a change point"
  } else if self.quality_score() < 80 {
    "inspect repeated or irregular observations before calibration"
  } else {
    "series is suitable for detector calibration"
  }
}
