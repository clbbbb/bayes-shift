///|
fn main {
  let model = @shift.GaussianModel::new(0.0, 1.0, 0.05).unwrap()
  let results = @shift.execute_default_plans(
    @shift.NumericScenario::mean_shift(),
    model,
  ).unwrap()
  println(@shift.plan_result_csv_header())
  for result in results {
    println(@shift.format_plan_result(result))
  }
}
