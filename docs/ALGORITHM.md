# BayesShift 算法说明

BayesShift 实现单变量 Bayesian Online Change Point Detection（BOCPD）。它在每个时刻维护运行长度 r_t 的后验；r_t 表示当前观测段从最近一次变点开始已经包含的观测数。应用不必等待完整历史数据，在接收到一个新观测后即可得到变点概率和最可能的运行长度。

实现依据 Adams 与 MacKay 的 [Bayesian Online Changepoint Detection](https://arxiv.org/abs/0710.3742)。论文是算法参考资料；本仓库没有移植或复制任何现有代码实现。

## 递推

设常数 hazard 为 H，观测为 x_t，前一步运行长度后验为 P(r_(t-1) | x_(1:t-1))。对于每一个被保留的运行长度，BayesShift 计算其观测模型给出的预测密度 p(x_t | r_(t-1))。

增长分支为：

~~~text
P(r_t = r + 1, x_(1:t))
  = (1 - H) p(x_t | r) P(r_(t-1) = r | x_(1:t-1))
~~~

变点分支从先验段模型开始：

~~~text
P(r_t = 0, x_(1:t))
  = H p(x_t | prior) Sum_r P(r_(t-1) = r | x_(1:t-1))
~~~

每一步都会把这些未归一化权重除以证据 P(x_t | x_(1:t-1))。由于前一步后验已经归一化，第二个求和恒为 1。实现直接使用 log(H) + log p(x_t | prior)，避免对每个候选运行长度重复一次 log-sum-exp。

## 高斯未知均值模型

GaussianModel 假定观测方差已知，段均值未知：

~~~text
mu ~ Normal(mu0, sigma^2 / kappa0)
x | mu ~ Normal(mu, sigma^2)
~~~

给定一个包含 n 个观测、观测和为 s 的候选段，预测均值为：

~~~text
(kappa0 mu0 + s) / (kappa0 + n)
~~~

预测方差为：

~~~text
sigma^2 (1 + 1 / (kappa0 + n))
~~~

因此每个候选段只需保存 n 和 s。更新不需要保留所有历史观测。

## Gamma-Poisson 计数模型

PoissonModel 用 Gamma 先验描述每段的非负整数事件率：

~~~text
lambda ~ Gamma(shape, rate)
x | lambda ~ Poisson(lambda)
~~~

对于候选段的累计计数 k 和观测数 n，后验参数是：

~~~text
shape' = shape + k
rate'  = rate + n
~~~

实现以对数 Gamma 函数计算相应的负二项预测概率。单次更新只计算一次当前观测的 log(x!)，所有候选运行长度复用该值。

## 数值稳定性

后验权重以对数形式保存。每轮更新先减去最大对数权重，再对指数权重求和并归一化。这样可避免长序列、极低密度或高计数输入下的浮点下溢。

输入构造函数拒绝无效 hazard、非正方差、非正先验强度和负计数。Detector.update 还会拒绝 NaN 与无穷连续观测。

## 运行长度上限

max_run_length 控制每轮保留的状态数。完整 BOCPD 的状态数随时间增长；BayesShift 保留从 0 到 max_run_length - 1 的状态，在达到上限后丢弃超出上限的增长分支，并重新归一化剩余后验。

这使每次更新的空间复杂度保持为 O(max_run_length)，时间复杂度为 O(min(t, max_run_length))。它是明确的截断近似：若业务需要识别远长于上限的稳定段，应提高 max_run_length 并用项目提供的基准检查成本。

## 复现与性能

benchmark.mbt 固定使用 512 个观测、前后各一个稳定段、运行长度上限 64 的高斯与 Poisson 工作负载：

~~~sh
moon bench --release --deny-warn
~~~

工作负载不把随机生成、文件读写或输出格式化放入计时区间。基准用于比较同一机器上的修订版本；不同硬件和后端的绝对时间不应直接比较。

核心回归测试验证：

- 高斯后验归一化、变点响应和快照恢复；
- 高斯与 Gamma-Poisson 生产递推和显式分支求和参考递推的完整后验对照；
- 高斯达到运行长度上限后的归一化；
- Poisson 计数率跳变、负计数拒绝和有界后验；
- 合成场景上的边界评分、阈值校准与基线对照。
