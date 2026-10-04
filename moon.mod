// Learn more about moon.mod configuration:
// https://docs.moonbitlang.com/en/latest/toolchain/moon/module.html
//
// To add a dependency, run this command in your terminal:
//   moon add moonbitlang/x
//
// Or manually declare it in `import`, for example:
// import {
//   "moonbitlang/x@0.4.6",
// }

name = "clbbbb/bayes-shift"

version = "0.1.0"

readme = "README.mbt.md"

repository = "https://github.com/clbbbb/bayes-shift"

license = "Apache-2.0"

keywords = [ "bayesian", "changepoint", "streaming", "timeseries" ]

preferred_target = "wasm"

description = "Bayesian online change point detection for streaming numeric data"
