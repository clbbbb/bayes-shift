///|
/// Stable text formatting for a numeric BOCPD update.
pub fn format_update(update : Update) -> String {
  "step=" +
  update.step.to_string() +
  " value=" +
  update.observation.to_string() +
  " change_probability=" +
  update.change_probability.to_string() +
  " run_length=" +
  update.most_likely_run_length.to_string() +
  " posterior_mean=" +
  update.posterior_mean.to_string()
}

///|
/// Stable CSV header for numeric update rows.
pub fn update_csv_header() -> String {
  "step,observation,change_probability,most_likely_run_length,posterior_mean,retained_run_lengths"
}

///|
pub fn update_csv_row(update : Update) -> String {
  update.step.to_string() +
  "," +
  update.observation.to_string() +
  "," +
  update.change_probability.to_string() +
  "," +
  update.most_likely_run_length.to_string() +
  "," +
  update.posterior_mean.to_string() +
  "," +
  update.retained_run_lengths.to_string()
}

///|
pub fn format_count_update(update : CountUpdate) -> String {
  "step=" +
  update.step.to_string() +
  " count=" +
  update.count.to_string() +
  " change_probability=" +
  update.change_probability.to_string() +
  " run_length=" +
  update.most_likely_run_length.to_string() +
  " posterior_rate=" +
  update.posterior_rate.to_string()
}

///|
pub fn count_update_csv_header() -> String {
  "step,count,change_probability,most_likely_run_length,posterior_rate"
}

///|
pub fn count_update_csv_row(update : CountUpdate) -> String {
  update.step.to_string() +
  "," +
  update.count.to_string() +
  "," +
  update.change_probability.to_string() +
  "," +
  update.most_likely_run_length.to_string() +
  "," +
  update.posterior_rate.to_string()
}

///|
/// One line intended for reproducible threshold-sweep output.
pub fn experiment_csv_header() -> String {
  "scenario,hazard,threshold,true_positives,false_positives,false_negatives,precision,recall,f1,mean_delay,maximum_probability,maximum_step"
}

///|
pub fn experiment_csv_row(experiment : NumericExperiment) -> String {
  experiment.scenario +
  "," +
  experiment.hazard.to_string() +
  "," +
  experiment.threshold.to_string() +
  "," +
  experiment.score.true_positives.to_string() +
  "," +
  experiment.score.false_positives.to_string() +
  "," +
  experiment.score.false_negatives.to_string() +
  "," +
  experiment.score.precision().to_string() +
  "," +
  experiment.score.recall().to_string() +
  "," +
  experiment.score.f1().to_string() +
  "," +
  experiment.score.mean_delay.to_string() +
  "," +
  experiment.trace.maximum_change_probability.to_string() +
  "," +
  experiment.trace.maximum_step.to_string()
}

///|
pub fn format_score(score : Score) -> String {
  "tp=" +
  score.true_positives.to_string() +
  " fp=" +
  score.false_positives.to_string() +
  " fn=" +
  score.false_negatives.to_string() +
  " precision=" +
  score.precision().to_string() +
  " recall=" +
  score.recall().to_string() +
  " f1=" +
  score.f1().to_string() +
  " mean_delay=" +
  score.mean_delay.to_string()
}

///|
/// Bounded ASCII bar for console output; NaN is not rendered silently.
pub fn probability_bar(
  probability : Double,
  width : Int,
) -> Result[String, DetectError] {
  if !finite(probability) || probability < 0.0 || probability > 1.0 || width < 1 {
    return Err(InvalidVariance)
  }
  let filled = (probability * width.to_double()).to_int()
  let mut output = ""
  for index in 0..<width {
    if index < filled {
      output = output + "#"
    } else {
      output = output + "."
    }
  }
  Ok(output)
}

///|
/// Render a small posterior chart in run-length order for debugging and demos.
pub fn posterior_chart(
  probabilities : Array[Double],
  width : Int,
) -> Result[String, DetectError] {
  if probabilities.length() == 0 || width < 1 {
    return Err(InvalidRunLengthLimit)
  }
  let mut output = ""
  for index in 0..<probabilities.length() {
    let bar = match probability_bar(probabilities[index], width) {
      Ok(value) => value
      Err(error) => return Err(error)
    }
    output = output + index.to_string() + " " + bar
    if index + 1 < probabilities.length() {
      output = output + "\n"
    }
  }
  Ok(output)
}
