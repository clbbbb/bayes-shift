# BayesShift

BayesShift 是一个 MoonBit 的在线变点检测库。它在每次收到新观测值时维护“当前段已经持续多久”的后验分布，并返回该位置开始新段的概率。适合把连续数值或事件计数流接入告警、监控、实验分析和离线评估程序。

目前提供两套共轭观测模型：

- 高斯未知均值模型：用于延迟、温度、传感器读数、评分等连续数值；
- Gamma-Poisson 模型：用于请求量、错误数、事件数等非负整数计数。

## 安装与验证

项目要求 MoonBit 0.10.14 或更高版本。

~~~sh
moon fmt
moon info
moon check --deny-warn
moon test --deny-warn
moon run cmd
~~~

cmd 会运行一段均值突变的可执行示例，并打印每一步的变点概率和最大后验运行长度。
持续集成还会在 Linux 上运行 Wasm、WasmGC、JavaScript 和 Native 全目标测试。

## 性能基准

仓库提供固定的 512 点工作负载：前半段稳定、后半段发生均值或计数率变化，运行长度上限为 64。可使用下面的命令在自己的机器上复测：

~~~sh
moon bench --release --deny-warn
~~~

基准同时运行高斯和 Gamma-Poisson 检测器。它避免把随机数据生成或输出格式化混入计时，衡量的是完整在线重放路径。CI 会执行该基准以确保性能工作负载持续可运行。

## 最小使用示例

~~~moonbit nocheck
let model = @shift.GaussianModel::new(0.0, 1.0, 0.25).unwrap()
let config = @shift.Config::new(0.04, 128, model).unwrap()
let detector = @shift.Detector::new(config)

for update in detector.replay([0.0, 0.2, -0.1, 0.0, 5.1, 4.9]).unwrap() {
  println(update.change_probability())
}
~~~

hazard 是每一步发生新段的先验概率。较小的值会减少告警，较大的值会提高对短段变化的灵敏度。max_run_length 给运行长度后验设上限，限制在线状态量。

## 主要接口

| 模块 | 内容 |
| --- | --- |
| Detector / GaussianModel | 连续值 BOCPD、快照恢复、批量重放 |
| PoissonDetector / PoissonModel | 计数流 BOCPD |
| EventGate | 概率阈值、最小段长与冷却窗口组成的事件门控 |
| NumericSeries / TimedSeries | 输入验证、滚动统计、时间序列规则性检查 |
| threshold_points / evaluate_hazards | 阈值和 hazard 网格评估 |
| evaluate_ewma / evaluate_cusum | 与 EWMA、CUSUM 基线对照 |
| diagnose_trace / audit_series | 后验轨迹诊断、输入质量检查 |
| gaussian_segments / poisson_segments | 固定种子的合成实验数据 |

## 可复现实验

NumericScenario::mean_shift()、two_shifts() 和 no_change() 是内置确定性场景。ExperimentPlan 固定 hazard、阈值、容忍窗口和随机种子；execute_plan 返回评分、预测边界与轨迹诊断。

~~~moonbit nocheck
let result = @shift.execute_plan(
  @shift.default_experiment_plan(),
  @shift.NumericScenario::mean_shift(),
  @shift.GaussianModel::new(0.0, 1.0, 0.05).unwrap(),
).unwrap()
println(result.summary())
~~~

边界评分会在给定容忍窗口内一对一匹配预测和真实变点，避免同一真实变点被多个告警重复计为真阳性。

## 工程边界

本项目关注单变量、顺序到达的变点检测与可复现实验。它不负责持久化时序数据、解析 PCAP、绘图服务或外部告警投递。应用层可以通过 Update.change_probability 和 EventGate 接入自己的存储、界面或通知系统。

## 方法与参考

高斯和 Poisson 检测器实现 Adams 与 MacKay 在 *Bayesian Online Changepoint Detection*（2007）描述的运行长度递推，采用对数空间归一化以避免长序列下溢。本仓库以 MoonBit 重新实现，未复制其他项目代码。

- Adams, R. P. and MacKay, D. J. C. (2007), *Bayesian Online Changepoint Detection*.
- 许可证：Apache-2.0，详见 [LICENSE](LICENSE)。

完整递推、观测模型、数值稳定策略、截断近似和性能测量方式见 [算法说明](docs/ALGORITHM.md)。
