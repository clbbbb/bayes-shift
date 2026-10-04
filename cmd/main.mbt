///|
fn main {
  let model = @shift.GaussianModel::new(0.0, 1.0, 0.25).unwrap()
  let detector = @shift.Detector::new(
    @shift.Config::new(0.04, 64, model).unwrap(),
  )
  for
    update in detector
    .replay([0.0, 0.2, -0.1, 0.1, 0.0, 4.8, 5.1, 4.9, 5.0])
    .unwrap() {
    println(
      "step=\{update.step()} change=\{update.change_probability()} run=\{update.most_likely_run_length()}",
    )
  }
}
