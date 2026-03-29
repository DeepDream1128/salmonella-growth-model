# Beamer 演讲大纲

满分 100 分，15 分钟演讲。每页标注插图位置和文字内容。

---

## Slide 1 — Title Page

**文字内容：**
- 标题：Dynamic Predictive Model for Growth of *Salmonella* Enteritidis in Egg Yolk
- 副标题：Gompertz Model with Modified Ratkowsky Secondary Model
- 课程：Modeling Methods in Biosystems Engineering, Spring 2026
- 小组成员姓名

**插图：** 无

---

## Slide 2 — Introduction & Background（5 分）

**文字内容（bullet points）：**
- SE is a leading cause of foodborne illness; shell eggs are a major vehicle
- SE penetrates shell → grows rapidly in iron-rich yolk
- USDA requires eggs stored ≤ 7.2°C within 12 h of laying
- Temperature fluctuates during cooling, storage, and distribution
- Objective: develop a dynamic Gompertz model for SE growth under sinusoidal temperature (3–43°C)

**插图：** 无（纯文字背景介绍）

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

**插图：** 无（公式页）

---

## Slide 4 — Data Description

**文字内容：** 左右两栏

左栏 — Growth Data:
- 12 time points (0–17 hr), 20 observations
- 0–2 hr: single measurement; 3–17 hr: duplicates

右栏 — Temperature Data:
- 47 points, sinusoidal ~7–43°C, 24 h cycle
- Interpolated via `interp1`

**插图：** 无（或可选插入原始数据散点图）

---

## Slide 5 — Parameters & Initial Guesses

**文字内容：** 参数表格

| Parameter | Description | Initial Guess | Status |
|-----------|-------------|---------------|--------|
| A | Initial log₁₀(CFU/mL) | 2.60 | Estimated |
| C | Growth range | 11 | Estimated |
| M | Inflection time (hr) | 7.5 | Estimated |
| a | Ratkowsky coeff | 0.000338 | Estimated |
| b | Ratkowsky coeff | 0.275 | Estimated |
| Tmin | Min growth temp | 6°C | Fixed |
| Tmax | Max growth temp | 46.3°C | Fixed |

**插图：** 无

---

## Slide 6 — Forward Problem（10 分，含 SSC）

**文字内容：**
- Left: Forward prediction with initial guesses
- Right: Scaled Sensitivity Coefficients

**插图：** 左右并排
- 左图：`figs/fig01_forward_guess.png`
- 右图：`figs/fig02_ssc_initial.png`

**底部说明文字：**
- Initial guesses produce reasonable fit → suitable for nlinfit convergence
- All 5 parameters show distinct SSC patterns → estimable
- Tmin, Tmax have negligible SSC → fixed

---

## Slide 7 — OLS Results（18 分）

**文字内容：** 参数估计结果表格（从 `results/report.txt` 填入）

| Param | Estimate | SE | Rel.Err% | 95% CI |
|-------|----------|----|----------|--------|
| A | ... | ... | ...% | [..., ...] |
| C | ... | ... | ...% | [..., ...] |
| M | ... | ... | ...% | [..., ...] |
| a | ... | ... | ...% | [..., ...] |
| b | ... | ... | ...% | [..., ...] |

底部：RMSE = ..., Pseudo-R² = ...

**插图：** 无（数据表格页）

---

## Slide 8 — Fitted Curve with CB & PB（12 分）

**文字内容：**
- Confidence Band (CB): 95% band for mean response
- Prediction Band (PB): 95% band for new observation

**插图：** 居中全幅
- `figs/fig03_fit_CB_PB.png`

---

## Slide 9 — Residual Analysis（12 分）

**文字内容：** 底部列出五项统计假设通过/未通过

**插图：** 左右并排
- 左图：`figs/fig04_residual_scatter.png`
- 右图：`figs/fig05_residual_histogram.png`

**底部说明文字：**
1. Model correct: ✓ (visual)
2. Errors random: DW = ... (✓/✗)
3. Constant variance: ✓ (visual)
4. Errors uncorrelated: see DW
5. Normal distribution: SW p = ... (✓/✗)

---

## Slide 10 — Final SSC（7 分）

**文字内容：**
- SSC recomputed with estimated parameters
- Confirms identifiability at solution

**插图：** 居中全幅
- `figs/fig06_ssc_final.png`

**底部说明文字：**
- Parameter accuracy ranking: (from report.txt)

---

## Slide 11 — Optimal Experimental Design（8 分）

**文字内容：**
- Delta criterion: det(X'X) measures information content
- Cii: diagonal of (X'X)⁻¹, proportional to parameter variance

**插图：** 左右并排
- 左图：`figs/fig07_delta_criterion.png`
- 右图：`figs/fig08_Cii_curves.png`

**底部说明文字：**
- Interpretation of optimal experiment duration

---

## Slide 12 — Bootstrap（10 分）

**文字内容：**
- Method: residual resampling, N = 1000
- Bootstrap vs Asymptotic CI comparison table (from report.txt)
- Band width comparison

**插图：** 右侧或下方
- `figs/fig09_bootstrap_CB_PB.png`

---

## Slide 13 — Summary & Conclusions

**文字内容（bullet points）：**
- Dynamic Gompertz model successfully developed
- 5 parameters estimated, 2 fixed; RMSE < 0.3 log₁₀ CFU/mL
- All statistical assumptions evaluated
- Bootstrap CIs consistent with asymptotic CIs
- Model supports USDA-FSIS 7.2°C storage recommendation

**插图：** 无

---

## Slide 14 — References

- Gumudavelli et al. (2007), J. Food Sci. 72:M254–62
- Baranyi & Roberts (1994), Int J Food Microbiol 23:277–94
- Zwietering et al. (1991), Appl Environ Microbiol 57:1094–101
- Beck & Arnold (1977), Parameter Estimation in Engineering and Science

**插图：** 无

---

## Slide 15 — Thank You / Questions

**文字内容：** Thank You! Questions?

**插图：** 无

---

## 演讲注意事项（10 分）

- 图表：大标记、粗线条、大字号坐标轴
- 每页 slide 文字尽量少
- 所有图必须让观众看清楚
- 所有组员都要上台讲
- 总时长不超过 15 分钟
