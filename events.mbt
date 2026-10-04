///|
/// A threshold policy turns posterior summaries into application-level events.
pub struct EventPolicy {
  probability_threshold : Double
  minimum_run_length : Int
  cooldown : Int
}

///|
/// Build a policy without embedding alert transport or side effects in the library.
pub fn EventPolicy::new(
  probability_threshold : Double,
  minimum_run_length : Int,
  cooldown : Int,
) -> Result[EventPolicy, DetectError] {
  if !finite(probability_threshold) ||
    probability_threshold <= 0.0 ||
    probability_threshold >= 1.0 {
    return Err(InvalidHazard)
  }
  if minimum_run_length < 0 || cooldown < 0 {
    return Err(InvalidRunLengthLimit)
  }
  Ok({ probability_threshold, minimum_run_length, cooldown, })
}

///|
/// An event contains the evidence needed by a caller to make an alert decision.
pub struct ChangeEvent {
  step : Int
  observation : Double
  probability : Double
  prior_run_length : Int
}

///|
pub fn ChangeEvent::step(self : ChangeEvent) -> Int {
  self.step
}

///|
pub fn ChangeEvent::probability(self : ChangeEvent) -> Double {
  self.probability
}

///|
/// Stateful debounce helper. It never mutates the detector itself.
pub struct EventGate {
  policy : EventPolicy
  mut last_emitted_step : Int
}

///|
pub fn EventGate::new(policy : EventPolicy) -> EventGate {
  { policy, last_emitted_step: -1000000000, }
}

///|
/// Emit at most one event per cooldown interval after a sufficiently long segment.
pub fn EventGate::observe(self : EventGate, update : Update) -> ChangeEvent? {
  let since_last = update.step - self.last_emitted_step
  if update.change_probability >= self.policy.probability_threshold &&
    update.most_likely_run_length <= self.policy.minimum_run_length &&
    since_last > self.policy.cooldown {
    self.last_emitted_step = update.step
    Some({
      step: update.step,
      observation: update.observation,
      probability: update.change_probability,
      prior_run_length: update.most_likely_run_length,
    })
  } else {
    None
  }
}

///|
/// Reset only the debounce state, retaining threshold configuration.
pub fn EventGate::reset(self : EventGate) -> Unit {
  self.last_emitted_step = -1000000000
}
