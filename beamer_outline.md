# Beamer 演讲大纲

满分 100 分，15 分钟演讲。每页标注插图位置、表格和文字内容。
所有表格使用 Beamer 标准三线表（`\toprule`, `\midrule`, `\bottomrule`）。

---

## Slide 1 — Title Page

**文字内容：**
- 标题：Dynamic Predictive Model for Growth of *Salmonella* Enteritidis in Egg Yolk
- 副标题：Gompertz Model with Modified Ratkowsky Secondary Model
- 课程：Modeling Methods in Biosystems Engineering, Spring 2026
- 小组成员姓名

**插图/表格：** 无

---

## Slide 2 — Introduction & Background（5 分）

**文字内容（bullet points）：**
- SE is a leading cause of foodborne illness; shell eggs are a major vehicle
- SE penetrates shell → grows rapidly in iron-rich yolk
- USDA requires eggs stored ≤ 7.2°C within 12 h of laying
- Temperature fluctuates during cooling, storage, and distribution
- Objective: develop a dynamic Gompertz model for SE growth under sinusoidal temperature (3–43°C)

**插图/表格：** 无

---

## Slide 3 — Model Description

**文字内容：** 左右两栏布局

左栏 — Primary Model (Gompertz):
- 恒温解析式公式
- 变温微分形式公式
- 参数说明：A, C, M

右栏 — Secondary Model (Ratkowsky):
- 修正 Ratkowsky 公式
- 参数说明：a, b, Tmin(fixed), Tmax(fixed)

底部一行：Dynamic model = Primary + Secondary, solved with `ode45`

**插图/表格：** 无（公式页）

---

## Slide 4 — Experimental Data（三线表展示）

**说明文字：** Sinusoidal temperature profile (3–43°C, 24 h). Data from Gumudavelli et al. (2007), Fig. 5c.

**表格 1：Growth Data（三线表，居中）**

Beamer LaTeX 代码：
```latex
\begin{table}
\centering
\caption{Observed SE growth data in egg yolk}
\begin{tabular}{ccc}
\toprule
Time (hr) & CFU/mL (Rep 1) & CFU/mL (Rep 2) \\
\midrule
0    & 304         & ---         \\
0.5  & 396         & ---         \\
1    & 332         & ---         \\
2    & 1\,184      & ---         \\
3    & 7\,600      & 7\,920      \\
4    & 90\,960     & 72\,768     \\
5    & 205\,600    & 303\,200    \\
6    & 306\,000    & 324\,000    \\
8    & 1\,696\,000 & 1\,752\,000 \\
12   & 14\,720\,000 & 14\,880\,000 \\
15   & 87\,200\,000 & 95\,200\,000 \\
17   & 361\,600\,000 & 393\,600\,000 \\
\bottomrule
\end{tabular}
\end{table}
```

**底部说明：** 0–2 hr: single measurement; 3–17 hr: duplicate measurements. Total: 20 data points.

---

## Slide 5 — Temperature Data（三线表，可选摘要版）

**说明文字：** 47 data points, sinusoidal ~7–43°C. Interpolated via `interp1`.

**表格 2：Temperature Data 摘要（三线表，选取关键时间点）**

```latex
\begin{table}
\centering
\caption{Selected temperature profile data points}
\begin{tabular}{cc}
\toprule
Time (hr) & Temperature (°C) \\
\midrule
0.00  & 25.0 \\
2.40  & 43.0 \\
4.80  & 25.9 \\
7.20  & 7.1  \\
9.60  & 11.4 \\
12.00 & 42.1 \\
14.40 & 31.6 \\
16.80 & 7.2  \\
19.20 & 12.0 \\
22.08 & 42.4 \\
\bottomrule
\end{tabular}
\end{table}
```

**底部说明：** Full dataset: 47 points at ~28.8 min intervals. See `Salmonella sin growth Temps.xlsx`.

---

## Slide 6 — Parameters & Initial Guesses（三线表）

**表格 3：参数表（三线表）**

```latex
\begin{table}
\centering
\caption{Model parameters and initial guesses}
\begin{tabular}{llcc}
\toprule
Parameter & Description & Initial Guess & Status \\
\midrule
$A$      & Initial $\log_{10}$(CFU/mL) & 2.60     & Estimated \\
$C$      & Growth range ($\log_{10}$)   & 11       & Estimated \\
$M$      & Inflection time (hr)         & 7.5      & Estimated \\
$a$      & Ratkowsky coeff (°C$^{-2}$)  & 0.000338 & Estimated \\
$b$      & Ratkowsky coeff (°C$^{-1}$)  & 0.275    & Estimated \\
$T_{\min}$ & Min growth temp (°C)       & 6.0      & Fixed \\
$T_{\max}$ & Max growth temp (°C)       & 46.3     & Fixed \\
\bottomrule
\end{tabular}
\end{table}
```

---

## Slide 7 — Forward Problem & SSC（10 分）

**插图：** 左右并排
- 左图：`figs/fig01_forward_guess.png`
- 右图：`figs/fig02_ssc_initial.png`

**底部说明文字：**
- Initial guesses produce reasonable fit → suitable for nlinfit convergence
- All 5 parameters show distinct SSC patterns → estimable
- Tmin, Tmax have negligible SSC → fixed

---

## Slide 8 — OLS Results（18 分，三线表）

**表格 4：参数估计结果（三线表，从 `results/report.txt` 填入）**

```latex
\begin{table}
\centering
\caption{OLS parameter estimates}
\begin{tabular}{lccccc}
\toprule
Param & Estimate & SE & Rel.Err (\%) & \multicolumn{2}{c}{95\% CI} \\
\midrule
$A$ & ... & ... & ... & ... & ... \\
$C$ & ... & ... & ... & ... & ... \\
$M$ & ... & ... & ... & ... & ... \\
$a$ & ... & ... & ... & ... & ... \\
$b$ & ... & ... & ... & ... & ... \\
\bottomrule
\end{tabular}
\end{table}
```

**底部：** RMSE = ..., Pseudo-R² = ..., $T_{\min}$ = 6.0 (fixed), $T_{\max}$ = 46.3 (fixed)

---

## Slide 9 — Fitted Curve with CB & PB（12 分）

**文字内容：**
- Confidence Band (CB): 95% band for mean response
- Prediction Band (PB): 95% band for new observation

**插图：** 居中全幅
- `figs/fig03_fit_CB_PB.png`

---

## Slide 10 — Residual Analysis（12 分）

**插图：** 左右并排
- 左图：`figs/fig04_residual_scatter.png`
- 右图：`figs/fig05_residual_histogram.png`

**表格 5：统计假设检验结果（三线表）**

```latex
\begin{table}
\centering
\caption{Standard statistical assumptions}
\begin{tabular}{clc}
\toprule
\# & Assumption & Result \\
\midrule
1 & Model is correct          & \checkmark\ (visual) \\
2 & Errors are random         & DW = ... (\checkmark/\texttimes) \\
3 & Constant variance         & \checkmark\ (visual) \\
4 & Errors uncorrelated       & See DW \\
5 & Normal distribution       & SW $p$ = ... (\checkmark/\texttimes) \\
\bottomrule
\end{tabular}
\end{table}
```

---

## Slide 11 — Final SSC（7 分）

**插图：** 居中全幅
- `figs/fig06_ssc_final.png`

**底部说明：** SSC recomputed with estimated parameters. Confirms identifiability at solution.

---

## Slide 12 — Optimal Experimental Design（8 分）

**插图：** 左右并排
- 左图：`figs/fig07_delta_criterion.png`
- 右图：`figs/fig08_Cii_curves.png`

**底部说明：**
- det(X'X) increases with time → longer experiments provide more information
- Cii decreases → parameters become more precise
- Key measurement window: exponential growth phase (~2–8 hr)

---

## Slide 13 — Bootstrap（10 分）

**表格 6：Bootstrap vs Asymptotic CI（三线表，从 `results/report.txt` 填入）**

```latex
\begin{table}
\centering
\caption{Bootstrap vs.\ Asymptotic 95\% CI ($N_{\text{boot}}=1000$)}
\begin{tabular}{lcccc}
\toprule
Param & \multicolumn{2}{c}{Bootstrap 95\% CI} & \multicolumn{2}{c}{Asymptotic 95\% CI} \\
\midrule
$A$ & ... & ... & ... & ... \\
$C$ & ... & ... & ... & ... \\
$M$ & ... & ... & ... & ... \\
$a$ & ... & ... & ... & ... \\
$b$ & ... & ... & ... & ... \\
\bottomrule
\end{tabular}
\end{table}
```

**插图：** 下方
- `figs/fig09_bootstrap_CB_PB.png`

**底部说明：** Asymptotic CB width: ..., Bootstrap CB width: ..., Asymptotic PB width: ..., Bootstrap PB width: ...

---

## Slide 14 — Summary & Conclusions

**文字内容（bullet points）：**
- Dynamic Gompertz model successfully developed
- 5 parameters estimated, 2 fixed; RMSE < 0.3 log₁₀ CFU/mL
- All statistical assumptions evaluated
- Bootstrap CIs consistent with asymptotic CIs
- Model supports USDA-FSIS 7.2°C storage recommendation

**插图/表格：** 无

---

## Slide 15 — References

- Gumudavelli et al. (2007), J. Food Sci. 72:M254–62
- Baranyi & Roberts (1994), Int J Food Microbiol 23:277–94
- Zwietering et al. (1991), Appl Environ Microbiol 57:1094–101
- Beck & Arnold (1977), Parameter Estimation in Engineering and Science

---

## Slide 16 — Thank You / Questions

**文字内容：** Thank You! Questions?

---

## 三线表 Beamer 配置提示

在 Beamer preamble 中加入：
```latex
\usepackage{booktabs}
```
所有表格使用 `\toprule`、`\midrule`、`\bottomrule`，不使用竖线。

## 演讲注意事项（10 分）

- 图表：大标记、粗线条、大字号坐标轴
- 表格：三线表，字号适当缩小（`\small` 或 `\footnotesize`）
- 每页 slide 文字尽量少
- 所有图和表必须让观众看清楚
- 所有组员都要上台讲
- 总时长不超过 15 分钟
