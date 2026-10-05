///|
/// Fixed workloads make benchmark output comparable across revisions.
fn gaussian_benchmark_values() -> Array[Double] {
  let values = []
  for index in 0..<512 {
    let level = if index < 256 { 0.0 } else { 4.0 }
    values.push(level + (index % 7).to_double() * 0.03)
  }
  values
}

///|
fn poisson_benchmark_values() -> Array[Int] {
  let values = []
  for index in 0..<512 {
    if index < 256 {
      values.push(index % 3)
    } else {
      values.push(20 + index % 5)
    }
  }
  values
}

///|
test (bench : @bench.T) {
  let gaussian_values = gaussian_benchmark_values()
  let poisson_values = poisson_benchmark_values()
  let gaussian_model = GaussianModel::new(0.0, 1.0, 0.25).unwrap()
  let poisson_model = PoissonModel::new(1.0, 1.0).unwrap()
  bench.bench(name="gaussian_bocpd_512x64", fn() {
    let detector = Detector::new(Config::new(0.04, 64, gaussian_model).unwrap())
    bench.keep(detector.replay(gaussian_values).unwrap().length())
  })
  bench.bench(name="poisson_bocpd_512x64", fn() {
    let detector = PoissonDetector::new(
      PoissonConfig::new(0.04, 64, poisson_model).unwrap(),
    )
    bench.keep(detector.replay(poisson_values).unwrap().length())
  })
}
