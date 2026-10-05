///|
/// Direct BOCPD recurrence used only as a small correctness oracle. Unlike the
/// production detector, it explicitly sums every reset branch.
fn reference_gaussian_posterior(
  observations : Array[Double],
  config : Config,
) -> Array[Double] {
  let mut log_posterior : Array[Double] = [0.0]
  let mut counts : Array[Int] = [0]
  let mut sums : Array[Double] = [0.0]
  for observation in observations {
    let old_size = log_posterior.length()
    let next_log = Array::make(old_size + 1, negative_infinity)
    let next_counts = Array::make(old_size + 1, 0)
    let next_sums = Array::make(old_size + 1, 0.0)
    let reset_terms : Array[Double] = []
    for run in 0..<old_size {
      let model = config.model
      let effective = model.prior_strength + counts[run].to_double()
      let mean = (model.prior_strength * model.prior_mean + sums[run]) /
        effective
      let predictive = normal_log_density(
        observation,
        mean,
        model.observation_variance * (1.0 + 1.0 / effective),
      )
      let prior_predictive = normal_log_density(
        observation,
        model.prior_mean,
        model.observation_variance * (1.0 + 1.0 / model.prior_strength),
      )
      reset_terms.push(
        log_posterior[run] + config.hazard_log + prior_predictive,
      )
      let grown = run + 1
      next_log[grown] = log_posterior[run] + config.survival_log + predictive
      next_counts[grown] = counts[run] + 1
      next_sums[grown] = sums[run] + observation
    }
    next_log[0] = log_sum_exp(reset_terms).unwrap()
    next_counts[0] = 1
    next_sums[0] = observation
    let _ = normalize_log_weights(next_log)
    log_posterior = next_log
    counts = next_counts
    sums = next_sums
  }
  log_posterior.map(fn(value) { @math.exp(value) })
}

///|
fn reference_poisson_predictive(
  total : Int,
  observations : Int,
  value : Int,
  value_log_factorial : Double,
  model : PoissonModel,
) -> Double {
  let shape = model.prior_shape + total.to_double()
  let rate = model.prior_rate + observations.to_double()
  let next_rate = rate + 1.0
  log_gamma(shape + value.to_double()).unwrap() -
  log_gamma(shape).unwrap() -
  value_log_factorial +
  shape * @math.ln(rate / next_rate) +
  value.to_double() * @math.ln(1.0 / next_rate)
}

///|
/// Direct Gamma-Poisson BOCPD recurrence used only as a correctness oracle.
fn reference_poisson_posterior(
  values : Array[Int],
  config : PoissonConfig,
) -> Array[Double] {
  let mut log_posterior : Array[Double] = [0.0]
  let mut totals : Array[Int] = [0]
  let mut observations : Array[Int] = [0]
  for value in values {
    let old_size = log_posterior.length()
    let next_log = Array::make(old_size + 1, negative_infinity)
    let next_totals = Array::make(old_size + 1, 0)
    let next_observations = Array::make(old_size + 1, 0)
    let reset_terms : Array[Double] = []
    let value_log_factorial = log_factorial(value).unwrap()
    let prior_predictive = reference_poisson_predictive(
      0,
      0,
      value,
      value_log_factorial,
      config.model,
    )
    for run in 0..<old_size {
      reset_terms.push(
        log_posterior[run] + config.hazard_log + prior_predictive,
      )
      let grown = run + 1
      next_log[grown] = log_posterior[run] +
        config.survival_log +
        reference_poisson_predictive(
          totals[run],
          observations[run],
          value,
          value_log_factorial,
          config.model,
        )
      next_totals[grown] = totals[run] + value
      next_observations[grown] = observations[run] + 1
    }
    next_log[0] = log_sum_exp(reset_terms).unwrap()
    next_totals[0] = value
    next_observations[0] = 1
    let _ = normalize_log_weights(next_log)
    log_posterior = next_log
    totals = next_totals
    observations = next_observations
  }
  log_posterior.map(fn(value) { @math.exp(value) })
}

///|
test "configuration rejects unsafe values" {
  let model = GaussianModel::new(0.0, 1.0, 1.0).unwrap()
  assert_true(Config::new(0.0, 32, model) is Err(InvalidHazard))
  assert_true(Config::new(0.1, 1, model) is Err(InvalidRunLengthLimit))
}

///|
test "log1mexp stays accurate near zero and across its branch point" {
  let near_zero = log1mexp(-0.000000000001).unwrap()
  let moderate = log1mexp(-0.5).unwrap()
  assert_true(near_zero > -27.631021116 && near_zero < -27.631021115)
  assert_true(moderate > -0.932752130 && moderate < -0.932752129)
}

///|
test "posterior remains normalized" {
  let model = GaussianModel::new(0.0, 1.0, 1.0).unwrap()
  let detector = Detector::new(Config::new(0.05, 32, model).unwrap())
  let _ = detector.replay([0.1, -0.2, 0.0, 0.3]).unwrap()
  let mut total = 0.0
  for probability in detector.run_length_posterior() {
    total = total + probability
  }
  assert_true(total > 0.999999999 && total < 1.000000001)
}

///|
test "bounded posterior remains normalized after reaching its limit" {
  let model = GaussianModel::new(0.0, 1.0, 1.0).unwrap()
  let detector = Detector::new(Config::new(0.05, 4, model).unwrap())
  let _ = detector.replay([0.0, 0.1, -0.1, 0.0, 0.2, -0.2, 0.1]).unwrap()
  assert_eq(detector.run_length_posterior().length(), 4)
  let mut total = 0.0
  for probability in detector.run_length_posterior() {
    total = total + probability
  }
  assert_true(total > 0.999999999 && total < 1.000000001)
}

///|
test "Gaussian detector rejects non-finite observations without advancing" {
  let model = GaussianModel::new(0.0, 1.0, 1.0).unwrap()
  let detector = Detector::new(Config::new(0.05, 8, model).unwrap())
  assert_true(detector.update(1.0 / 0.0) is Err(NonFiniteObservation))
  assert_eq(detector.steps(), 0)
}

///|
test "first Gaussian update preserves the configured hazard posterior" {
  let detector = Detector::new(
    Config::new(0.17, 8, GaussianModel::new(2.0, 3.0, 0.5).unwrap()).unwrap(),
  )
  let update = detector.update(-4.0).unwrap()
  assert_true(
    update.change_probability > 0.169999999999 &&
    update.change_probability < 0.170000000001,
  )
}

///|
test "Gaussian detector matches the direct BOCPD recurrence" {
  let config = Config::new(
    0.13,
    16,
    GaussianModel::new(0.5, 2.0, 0.75).unwrap(),
  ).unwrap()
  let observations = [-0.3, 0.2, 0.1, 2.4, 2.1]
  let expected = reference_gaussian_posterior(observations, config)
  let detector = Detector::new(config)
  let _ = detector.replay(observations).unwrap()
  let actual = detector.run_length_posterior()
  assert_eq(actual.length(), expected.length())
  for index in 0..<actual.length() {
    let difference = actual[index] - expected[index]
    assert_true(difference > -0.000000000001 && difference < 0.000000000001)
  }
}

///|
test "a large shift favors a short run" {
  let model = GaussianModel::new(0.0, 1.0, 0.05).unwrap()
  let detector = Detector::new(Config::new(0.1, 64, model).unwrap())
  let _ = detector.replay([0.0, 0.1, -0.1, 0.0]).unwrap()
  let update = detector.update(8.0).unwrap()
  assert_true(update.most_likely_run_length <= 1)
}

///|
test "snapshot restores an exact continuation" {
  let model = GaussianModel::new(0.0, 1.0, 1.0).unwrap()
  let detector = Detector::new(Config::new(0.05, 32, model).unwrap())
  let _ = detector.replay([0.2, 0.3, -0.1]).unwrap()
  let snapshot = detector.snapshot()
  let first = detector.update(0.4).unwrap()
  detector.restore(snapshot)
  let second = detector.update(0.4).unwrap()
  assert_eq(first.most_likely_run_length, second.most_likely_run_length)
}

///|
test "event gate emits a single high-confidence change" {
  let model = GaussianModel::new(0.0, 1.0, 0.05).unwrap()
  let detector = Detector::new(Config::new(0.1, 64, model).unwrap())
  let gate = EventGate::new(EventPolicy::new(0.9, 1, 3).unwrap())
  let _ = detector.replay([0.0, 0.1, -0.1]).unwrap()
  assert_true(gate.observe(detector.update(8.0).unwrap()) is Some(_))
}

///|
test "synthetic fixture is deterministic and summary finds its boundary" {
  let values = piecewise_constant([5, 5], [0.0, 7.0], 0.1, 42).unwrap()
  let repeat = piecewise_constant([5, 5], [0.0, 7.0], 0.1, 42).unwrap()
  assert_eq(values, repeat)
  let model = GaussianModel::new(0.0, 1.0, 0.05).unwrap()
  let detector = Detector::new(Config::new(0.05, 32, model).unwrap())
  let summary = summarize(detector.replay(values).unwrap(), 0.8).unwrap()
  assert_true(summary.maximum_change_probability > 0.8)
  assert_true(summary.maximum_step >= 6)
}

///|
test "boundary scoring uses each prediction at most once" {
  let truth = [Boundary::new(10).unwrap(), Boundary::new(20).unwrap()]
  let predicted = [
    Boundary::new(9).unwrap(),
    Boundary::new(11).unwrap(),
    Boundary::new(25).unwrap(),
  ]
  let score = score_boundaries(truth, predicted, 2).unwrap()
  assert_eq(score.true_positives, 1)
  assert_eq(score.false_positives, 2)
  assert_eq(score.false_negatives, 1)
}

///|
test "poisson detector recognizes a count-rate jump" {
  let model = PoissonModel::new(1.0, 1.0).unwrap()
  let detector = PoissonDetector::new(
    PoissonConfig::new(0.05, 64, model).unwrap(),
  )
  let _ = detector.replay([1, 0, 2, 1, 1]).unwrap()
  let update = detector.update(30).unwrap()
  assert_true(update.change_probability > 0.8)
  assert_true(update.most_likely_run_length <= 1)
}

///|
test "Poisson detector keeps a normalized bounded posterior" {
  let detector = PoissonDetector::new(
    PoissonConfig::new(0.05, 4, PoissonModel::new(1.0, 1.0).unwrap()).unwrap(),
  )
  assert_true(detector.update(-1) is Err(InvalidRunLengthLimit))
  let _ = detector.replay([0, 1, 2, 0, 3, 1]).unwrap()
  assert_eq(detector.steps(), 6)
  assert_eq(detector.run_length_posterior().length(), 4)
  let mut total = 0.0
  for probability in detector.run_length_posterior() {
    total = total + probability
  }
  assert_true(total > 0.999999999 && total < 1.000000001)
}

///|
test "Poisson detector matches the direct BOCPD recurrence" {
  let config = PoissonConfig::new(
    0.19,
    16,
    PoissonModel::new(1.5, 0.75).unwrap(),
  ).unwrap()
  let values = [1, 0, 2, 5, 4]
  let expected = reference_poisson_posterior(values, config)
  let detector = PoissonDetector::new(config)
  let _ = detector.replay(values).unwrap()
  let actual = detector.run_length_posterior()
  assert_eq(actual.length(), expected.length())
  for index in 0..<actual.length() {
    let difference = actual[index] - expected[index]
    assert_true(difference > -0.000000000001 && difference < 0.000000000001)
  }
}

///|
test "first Poisson update matches the conjugate posterior rate" {
  let detector = PoissonDetector::new(
    PoissonConfig::new(0.2, 8, PoissonModel::new(2.0, 1.0).unwrap()).unwrap(),
  )
  let update = detector.update(3).unwrap()
  assert_true(
    update.change_probability > 0.199999999999 &&
    update.change_probability < 0.200000000001,
  )
  assert_eq(update.posterior_rate, 2.5)
}

///|
test "numeric scenario sweep scores the known mean shift" {
  let model = GaussianModel::new(0.0, 1.0, 0.05).unwrap()
  let results = sweep_thresholds(
    NumericScenario::mean_shift(),
    model,
    0.05,
    64,
    [0.5, 0.8, 0.95],
    3,
  ).unwrap()
  assert_true(results[0].f1() > 0.5)
}

///|
test "count experiment locates its known rate boundary" {
  let result = evaluate_count(
    CountScenario::rate_shift(),
    PoissonModel::new(1.0, 1.0).unwrap(),
    0.05,
    32,
    0.8,
    2,
  ).unwrap()
  assert_true(result.score.true_positives == 1)
}

///|
test "ewma baseline detects a sustained change" {
  let detector = EwmaDetector::new(EwmaConfig::new(0.1, 2.0, 3).unwrap())
  let updates = detector.replay([0.0, 0.1, -0.1, 0.0, 10.0]).unwrap()
  assert_true(updates[4].drift)
}

///|
test "cusum accumulates positive residuals" {
  let detector = CusumDetector::new(CusumConfig::new(0.0, 0.2, 2.0).unwrap())
  let _ = detector.update(1.0).unwrap()
  let update = detector.update(1.5).unwrap()
  assert_true(update.drift)
}

///|
test "report helpers produce stable CSV and charts" {
  assert_eq(
    update_csv_header(),
    "step,observation,change_probability,most_likely_run_length,posterior_mean,retained_run_lengths",
  )
  assert_eq(probability_bar(0.5, 4).unwrap(), "##..")
  assert_eq(posterior_chart([0.5, 0.0], 2).unwrap(), "0 #.\n1 ..")
}

///|
test "window evicts oldest values and maintains statistics" {
  let window = Window::new(3).unwrap()
  let _ = window.push(1.0)
  let _ = window.push(2.0)
  let _ = window.push(3.0)
  let _ = window.push(4.0)
  let _ = window.push(5.0)
  assert_eq(window.values(), [3.0, 4.0, 5.0])
  assert_eq(window.mean(), Some(4.0))
  assert_eq(window.minimum(), Some(3.0))
  assert_eq(window.maximum(), Some(5.0))
}

///|
test "numeric series supports deterministic rolling summaries" {
  let series = NumericSeries::new([1.0, 2.0, 3.0, 4.0]).unwrap()
  assert_eq(series.difference(), [1.0, 1.0, 1.0])
  assert_eq(series.moving_means(2).unwrap(), [1.5, 2.5, 3.5])
  assert_eq(series.subseries(1, 3).unwrap().values(), [2.0, 3.0])
}

///|
test "robust statistics resist a single large spike" {
  let values = [1.0, 1.0, 1.0, 1.0, 100.0]
  assert_eq(median(values).unwrap(), 1.0)
  assert_eq(median_absolute_deviation(values).unwrap(), 0.0)
  assert_true(trimmed_mean(values, 0.2).unwrap() < 2.0)
}

///|
test "threshold calibration picks a usable mean-shift point" {
  let scenario = NumericScenario::mean_shift()
  let detector = Detector::new(
    Config::new(0.05, 64, GaussianModel::new(0.0, 1.0, 0.05).unwrap()).unwrap(),
  )
  let updates = detector.replay(scenario.values).unwrap()
  let points = threshold_points(
    updates,
    scenario.boundaries,
    [0.5, 0.8, 0.95],
    3,
  ).unwrap()
  assert_true(best_f1_point(points).unwrap().f1() > 0.5)
}

///|
test "protocol output records diagnostics and audit catches constants" {
  let scenario = NumericScenario::mean_shift()
  let result = execute_plan(
    default_experiment_plan(),
    scenario,
    GaussianModel::new(0.0, 1.0, 0.05).unwrap(),
  ).unwrap()
  assert_true(result.diagnostics.maximum_change_probability > 0.8)
  let audit = audit_series([2.0, 2.0, 2.0])
  assert_true(audit.is_constant())
  assert_eq(audit.quality_score(), 70)
}
