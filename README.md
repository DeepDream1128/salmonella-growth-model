# Salmonella Growth Model — Final Project

Dynamic Gompertz Model for predicting *Salmonella* Enteritidis growth in egg yolk under sinusoidal temperature profiles.

## Files

- `global_example.m` — Complete MATLAB code (based on instructor template)
- `Salmonella sin growth.xlsx` — Experimental growth data (CFU/mL vs time)
- `Salmonella sin growth Temps.xlsx` — Temperature profile data (sinusoidal)
- `final_project.m` — Standalone version (alternative)

## Usage

1. Open MATLAB, `cd` to this repo
2. Run `global_example.m`
3. Console outputs all statistics; figures auto-generated

## Model

- **Primary**: Gompertz differential form — `dy/dt = −μ(T)·e·(y−A)/C·ln((y−A)/C)`
- **Secondary**: Modified Ratkowsky — `μ(T) = a·(T−Tmin)²·(1−exp(b·(T−Tmax)))`
- **Solver**: `ode45` + `interp1` for temperature interpolation
- **Estimation**: `nlinfit` (OLS, Levenberg-Marquardt)

## Parameters

| Parameter | Description | Initial Guess | Status |
|-----------|-------------|---------------|--------|
| A | Initial log₁₀(CFU/mL) | log₁₀(400) ≈ 2.60 | Estimated |
| C | Growth range (log₁₀) | 11 | Estimated |
| M | Inflection time (hr) | 7.5 | Estimated |
| a | Ratkowsky coeff (°C⁻²) | 0.000338 | Estimated |
| b | Ratkowsky coeff (°C⁻¹) | 0.275 | Estimated |
| Tmin | Min growth temp (°C) | 6 | Fixed |
| Tmax | Max growth temp (°C) | 46.3 | Fixed |

---

## Beamer Presentation Outline

Based on the scoring rubric (100 points total).

### Slide 1 — Title
- Title: Dynamic Predictive Model for Growth of *Salmonella* Enteritidis in Egg Yolk
- Subtitle: Gompertz Model with Modified Ratkowsky Secondary Model
- Course: Modeling Methods in Biosystems Engineering, Spring 2026

### Slide 2 — Introduction & Background (5 pts)
- SE is a leading cause of foodborne illness; shell eggs are a major vehicle
- SE penetrates shell → grows rapidly in iron-rich yolk
- USDA-FSIS requires eggs stored ≤ 7.2°C within 12 h of laying
- Problem: temperature fluctuates during cooling, storage, distribution
- Objective: develop and validate a dynamic Gompertz model for SE growth under sinusoidal temperature (3–43°C, 24 h)
- Dependent variable: log₁₀(N) (CFU/mL); Independent variable: time (hr)
- 7 parameters total (5 estimated, 2 fixed)
- Reference: Gumudavelli et al. (2007), J. Food Sci. 72:M254–62

### Slide 3 — Model Description
- Primary model: Gompertz analytical form (constant T) and differential form (varying T)
- Secondary model: Modified Ratkowsky equation
- Dynamic model: integrate primary + secondary, solve with `ode45`
- Show equations

### Slide 4 — Data Description
- Growth data: 12 time points (0–17 hr), 20 total observations (singles at 0–2 hr, duplicates at 3–17 hr)
- Temperature data: 47 points, sinusoidal ~7–43°C, interpolated via `interp1`

### Slide 5 — Forward Problem: Ypred with Guesses + SSC (10 pts)
- Plot 1: Ypred (initial guesses) vs observed data, with temperature on right axis
- Plot 2: Scaled Sensitivity Coefficients (initial guesses)
  - Which parameters can be estimated and why
  - Ranking of estimability: which will be most accurate
  - Tmin/Tmax have small SSC → fixed

### Slide 6 — OLS Statistics (18 pts)
- Parameter estimates table: value, SE, relative error, 95% CI
- Correlation matrix
- RMSE, Pseudo-R²
- Interpretation of results

### Slide 7 — Fitted Curve with CB and PB (12 pts)
- Plot: Yobs (dots), Ypred (line), asymptotic Confidence Band, Prediction Band
- Temperature on right axis
- CB = uncertainty in mean response; PB = uncertainty for new observation
- Computed via `nlpredci`

### Slide 8 — Residual Analysis (12 pts)
- Residual scatter plot (residuals vs time)
- Residual histogram
- Five standard statistical assumptions (pass/fail):
  1. Model is correct (visual)
  2. Errors are random (Durbin-Watson)
  3. Constant variance (visual)
  4. Errors uncorrelated (Durbin-Watson)
  5. Normally distributed (Shapiro-Wilk test)

### Slide 9 — Final SSC (7 pts)
- SSC recomputed with estimated parameters
- Compare with initial SSC
- Confirm parameter identifiability at solution
- Ranking of parameter accuracy

### Slide 10 — Optimal Experimental Design (8 pts)
- Delta criterion: det(X'X) vs last measurement time
- Cii curves: C₁₁, C₂₂, ..., C₅₅ vs time
- Interpretation: optimal experiment duration, key measurement windows
- If no clear optimum, explain why

### Slide 11 — Bootstrap (10 pts)
- Method: residual resampling, 1000 iterations
- Bootstrap 95% CI for each parameter (compare with asymptotic CI)
- Bootstrap CB and PB plot
- Compare bootstrap vs asymptotic band widths (wider or narrower?)

### Slide 12 — Summary & Conclusions
- Model fits well (RMSE < 0.3 log₁₀ CFU/mL)
- All 5 estimated parameters identifiable
- Statistical assumptions evaluated
- Practical significance: predict SE risk during egg cooling/storage/distribution
- Supports USDA-FSIS 7.2°C recommendation

### Slide 13 — References
- Gumudavelli et al. (2007)
- Baranyi & Roberts (1994)
- Zwietering et al. (1991)
- Beck & Arnold (1977)

### Presentation Tips (10 pts for quality)
- Large marker size, thick lines, large fonts on axes
- Minimize legends, minimize text per slide
- Every plot must be clearly visible to audience
- Speak clearly, proper volume, emphasis
- All team members present
- Max 15 minutes

## Reference

Gumudavelli V, Subbiah J, Thippareddi H, Velugoti PR, Froning G. 2007. Dynamic predictive model for growth of *Salmonella* Enteritidis in egg yolk. *J. Food Sci.* 72:M254–62.
