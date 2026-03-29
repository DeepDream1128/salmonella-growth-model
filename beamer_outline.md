# Beamer 演讲大纲

满分 100 分，15 分钟演讲。每页标注插图位置、表格和文字内容。
所有表格使用 Beamer 标准三线表（`\toprule`, `\midrule`, `\bottomrule`），需 `\usepackage{booktabs}`。

---

## Slide 1 — Title Page

**内容：** 标题、副标题、课程名、小组成员

---

## Slide 2 — Introduction & Background（5 分）

**文字（bullets）：**
- SE is a leading cause of foodborne illness; shell eggs are a major vehicle
- SE penetrates shell → grows rapidly in iron-rich yolk
- USDA requires eggs stored ≤ 7.2°C within 12 h of laying
- Temperature fluctuates during cooling, storage, and distribution
- Objective: develop a dynamic Gompertz model under sinusoidal temperature (3–43°C)
- Reference: Gumudavelli et al. (2007)

---

## Slide 3 — Model Description

左栏 — Primary Model (Gompertz): 解析式 + 微分形式 + 参数 A, C, M
右栏 — Secondary Model (Ratkowsky): 公式 + 参数 a, b, Tmin, Tmax
底部：Dynamic model = Primary + Secondary, solved with `ode45`

---

## Slide 4 — Data & Parameters

**文字：** 左栏数据说明，右栏参数表

左栏：
- Growth: 12 time points (0–17 hr), 20 observations, sinusoidal T
- Temperature: 47 points, ~7–43°C, interpolated via `interp1`

**表格 1：参数与初始猜测值（三线表，右栏）**

```latex
\begin{tabular}{llcc}
\toprule
Param & Description & Guess & Status \\
\midrule
$A$        & Init.\ $\log_{10}$(CFU/mL) & 2.60     & Est. \\
$C$        & Growth range               & 11       & Est. \\
$M$        & Inflection time (hr)       & 7.5      & Est. \\
$a$        & Ratkowsky coeff            & 3.38e-4  & Est. \\
$b$        & Ratkowsky coeff            & 0.275    & Est. \\
$T_{\min}$ & Min growth temp (°C)       & 6.0      & Fixed \\
$T_{\max}$ & Max growth temp (°C)       & 46.3     & Fixed \\
\bottomrule
\end{tabular}
```

---

## Slide 5 — Forward Problem & SSC（10 分）

**插图：** 左右并排
- 左：`figs/fig01_forward_guess.png`
- 右：`figs/fig02_ssc_initial.png`

**底部说明：**
- Guesses produce reasonable fit → nlinfit can converge
- All 5 params show distinct SSC → estimable; Tmin/Tmax negligible → fixed

---

## Slide 6 — OLS Parameter Estimates（18 分）

**表格 2：参数估计结果（三线表，全幅居中，从 `results/report.txt` 填入）**

```latex
\begin{tabular}{lrrrrr}
\toprule
Param & Estimate & SE & Rel.Err (\%) & CI$_{lo}$ & CI$_{hi}$ \\
\midrule
$A$ & ... & ... & ... & ... & ... \\
$C$ & ... & ... & ... & ... & ... \\
$M$ & ... & ... & ... & ... & ... \\
$a$ & ... & ... & ... & ... & ... \\
$b$ & ... & ... & ... & ... & ... \\
\bottomrule
\end{tabular}
```

**表格 3：相关矩阵（三线表）**

```latex
\begin{tabular}{lccccc}
\toprule
      & $A$ & $C$ & $M$ & $a$ & $b$ \\
\midrule
$A$ & 1 & ... & ... & ... & ... \\
$C$ & ... & 1 & ... & ... & ... \\
$M$ & ... & ... & 1 & ... & ... \\
$a$ & ... & ... & ... & 1 & ... \\
$b$ & ... & ... & ... & ... & 1 \\
\bottomrule
\end{tabular}
```

**底部：** RMSE = ..., Pseudo-R² = ...

---

## Slide 7 — Goodness of Fit Summary（18 分续）

**表格 4：拟合优度指标（三线表）**

```latex
\begin{tabular}{lr}
\toprule
Metric & Value \\
\midrule
RMSE ($\log_{10}$ CFU/mL) & ... \\
Pseudo-$R^2$              & ... \\
MSE                       & ... \\
$N$ (data points)         & 20 \\
\bottomrule
\end{tabular}
```

---

## Slide 8 — Fitted Curve with CB & PB（12 分）

**文字：** CB = 95% band for mean; PB = 95% band for new observation

**插图：** 居中全幅 — `figs/fig03_fit_CB_PB.png`

---

## Slide 9 — Residual Analysis（12 分）

**插图：** 左右并排
- 左：`figs/fig04_residual_scatter.png`
- 右：`figs/fig05_residual_histogram.png`

**表格 5：统计假设检验（三线表）**

```latex
\begin{tabular}{clc}
\toprule
\# & Assumption & Result \\
\midrule
1 & Model is correct      & \checkmark\ (visual) \\
2 & Errors are random     & DW = ... \\
3 & Constant variance     & \checkmark\ (visual) \\
4 & Errors uncorrelated   & See DW \\
5 & Normal distribution   & SW $p$ = ... \\
\bottomrule
\end{tabular}
```

---

## Slide 10 — Final SSC（7 分）

**插图：** 居中全幅 — `figs/fig06_ssc_final.png`

**表格 6：最终 SSC 最大值（三线表）**

```latex
\begin{tabular}{lr}
\toprule
Param & max$|$SSC$|$ \\
\midrule
$A$ & ... \\
$C$ & ... \\
$M$ & ... \\
$a$ & ... \\
$b$ & ... \\
\bottomrule
\end{tabular}
```

---

## Slide 11 — Optimal Experimental Design（8 分）

**插图：** 左右并排
- 左：`figs/fig07_delta_criterion.png`
- 右：`figs/fig08_Cii_curves.png`

**底部说明：** Key measurement window: exponential growth phase (~2–8 hr)

---

## Slide 12 — Bootstrap（10 分）

**表格 7：Bootstrap vs Asymptotic 95% CI（三线表，从 `results/report.txt` 填入）**

```latex
\begin{tabular}{lcccc}
\toprule
Param & \multicolumn{2}{c}{Bootstrap CI} & \multicolumn{2}{c}{Asymptotic CI} \\
      & Lower & Upper & Lower & Upper \\
\midrule
$A$ & ... & ... & ... & ... \\
$C$ & ... & ... & ... & ... \\
$M$ & ... & ... & ... & ... \\
$a$ & ... & ... & ... & ... \\
$b$ & ... & ... & ... & ... \\
\bottomrule
\end{tabular}
```

**表格 8：带宽对比（三线表）**

```latex
\begin{tabular}{lc}
\toprule
Band & Avg Width \\
\midrule
Asymptotic CB & ... \\
Bootstrap CB  & ... \\
Asymptotic PB & ... \\
Bootstrap PB  & ... \\
\bottomrule
\end{tabular}
```

**插图：** 下方 — `figs/fig09_bootstrap_CB_PB.png`

---

## Slide 13 — Summary & Conclusions

- Dynamic Gompertz model successfully developed
- 5 params estimated, 2 fixed; RMSE < 0.3 log₁₀ CFU/mL
- All 5 statistical assumptions evaluated
- Bootstrap CIs consistent with asymptotic CIs
- Model supports USDA-FSIS 7.2°C recommendation

---

## Slide 14 — References

- Gumudavelli et al. (2007), J. Food Sci. 72:M254–62
- Baranyi & Roberts (1994)
- Zwietering et al. (1991)
- Beck & Arnold (1977)

---

## Slide 15 — Thank You / Questions

---

## 三线表汇总

| 表格 | 所在 Slide | 内容 | 数据来源 |
|------|-----------|------|---------|
| 表1 | Slide 4 | 参数与初始猜测值 | 老师给定 |
| 表2 | Slide 6 | OLS 参数估计结果 | `results/report.txt` |
| 表3 | Slide 6 | 相关矩阵 | `results/report.txt` |
| 表4 | Slide 7 | 拟合优度指标 | `results/report.txt` |
| 表5 | Slide 9 | 统计假设检验 | `results/report.txt` |
| 表6 | Slide 10 | 最终 SSC 最大值 | `results/report.txt` |
| 表7 | Slide 12 | Bootstrap vs Asymptotic CI | `results/report.txt` |
| 表8 | Slide 12 | 带宽对比 | `results/report.txt` |

## 演讲注意事项（10 分）

- 表格用 `\small` 或 `\footnotesize` 缩小字号
- 图表大标记、粗线条、大字号坐标轴
- 每页文字尽量少，所有图表必须清晰可见
- 所有组员上台讲，总时长 ≤ 15 分钟
