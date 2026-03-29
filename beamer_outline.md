# Beamer 演讲大纲

基于评分标准（满分 100 分），15 分钟演讲。

---

## Slide 1 — 封面

- 标题：Dynamic Predictive Model for Growth of *Salmonella* Enteritidis in Egg Yolk
- 副标题：Gompertz Model with Modified Ratkowsky Secondary Model
- 课程：Modeling Methods in Biosystems Engineering, Spring 2026
- 小组成员

## Slide 2 — 背景与目标（5 分）

- 肠炎沙门氏菌（SE）是全球主要食源性致病菌，鸡蛋是主要传播载体
- SE 穿透蛋壳后在富铁蛋黄中快速增殖
- USDA-FSIS 要求鸡蛋产后 12 小时内冷却至 ≤ 7.2°C
- 问题：冷却、储存、运输过程中温度持续波动
- 目标：建立正弦变温条件下（3–43°C，24 h 周期）SE 在蛋黄中的动态 Gompertz 生长模型
- 因变量：log₁₀(N) (CFU/mL)；自变量：时间 (hr)
- 7 个参数（5 个估计，2 个固定）
- 参考文献：Gumudavelli et al. (2007)

## Slide 3 — 模型描述

- 一级模型：Gompertz 解析式（恒温）+ 微分形式（变温）
- 二级模型：修正 Ratkowsky 方程
- 动态模型：一级 + 二级整合，`ode45` 数值求解
- 展示公式

## Slide 4 — 数据描述

- 菌落数据：12 个时间点（0–17 hr），20 个观测值（0–2 hr 单次，3–17 hr 双重复）
- 温度数据：47 个点，正弦波 ~7–43°C，通过 `interp1` 插值

## Slide 5 — 正向问题 + SSC（10 分）

- 图 1：初始猜测值的 Ypred vs 观测数据（含温度曲线）→ `fig01_forward_guess.png`
- 图 2：初始猜测的缩放灵敏度系数（SSC）→ `fig02_ssc_initial.png`
- 分析：
  - 哪些参数可估计，为什么
  - 参数估计精度排序
  - Tmin/Tmax 的 SSC 很小 → 固定

## Slide 6 — OLS 统计结果（18 分）

- 参数估计表：估计值、标准误、相对误差、95% 置信区间
- 相关矩阵
- RMSE、Pseudo-R²
- 结果解读

## Slide 7 — 拟合曲线 + CB/PB（12 分）

- 图 3：观测值（点）、预测值（线）、渐近置信带（CB）、预测带（PB）→ `fig03_fit_CB_PB.png`
- 右轴温度曲线
- CB = 均值响应的不确定性；PB = 新观测值的不确定性
- 使用 `nlpredci` 计算

## Slide 8 — 残差分析（12 分）

- 图 4：残差散点图 → `fig04_residual_scatter.png`
- 图 5：残差直方图 → `fig05_residual_histogram.png`
- 五项标准统计假设（逐项报告通过/未通过）：
  1. 模型正确（目视检查）
  2. 误差随机（Durbin-Watson 统计量）
  3. 方差齐性（残差图目视）
  4. 误差不相关（Durbin-Watson）
  5. 正态分布（Shapiro-Wilk 检验）

## Slide 9 — 最终 SSC（7 分）

- 图 6：用估计参数重新计算的 SSC → `fig06_ssc_final.png`
- 与初始 SSC 对比
- 确认参数在解处仍可辨识
- 参数精度排序

## Slide 10 — 最优实验设计（8 分）

- 图 7：Delta 准则 det(X'X) vs 最后测量时间 → `fig07_delta_criterion.png`
- 图 8：Cii 曲线 → `fig08_Cii_curves.png`
- 解读：最优实验时长、关键测量时间窗口
- 如果没有明显最优点，解释原因

## Slide 11 — Bootstrap（10 分）

- 方法：残差重采样，1000 次迭代
- 各参数的 Bootstrap 95% CI（与渐近 CI 对比）
- 图 9：Bootstrap CB 和 PB → `fig09_bootstrap_CB_PB.png`
- 对比 Bootstrap 与渐近带宽（哪个更宽/更窄）

## Slide 12 — 总结与结论

- 模型拟合良好（RMSE < 0.3 log₁₀ CFU/mL）
- 5 个估计参数均可辨识
- 五项统计假设已评估
- 最优实验设计已分析
- Bootstrap CI 与渐近 CI 一致
- 实际意义：可预测鸡蛋冷却/储存/运输过程中的 SE 风险
- 支持 USDA-FSIS 7.2°C 储存建议

## Slide 13 — 参考文献

- Gumudavelli et al. (2007), J. Food Sci. 72:M254–62
- Baranyi & Roberts (1994), Int J Food Microbiol 23:277–94
- Zwietering et al. (1991), Appl Environ Microbiol 57:1094–101
- Beck & Arnold (1977), Parameter Estimation in Engineering and Science

---

## 演讲注意事项（演讲质量 10 分）

- 图表：大标记、粗线条、大字号坐标轴
- 每页 slide 文字尽量少，图为主
- 少用图例，能标注就标注
- 所有图必须让观众看清楚
- 语速适中，重点强调
- 所有组员都要上台讲
- 总时长不超过 15 分钟
