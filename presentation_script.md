# Presentation Script (English)

~13 minutes total. Matches beamer_outline.md (15 slides, 8 three-line tables for results).

---

## Slide 1 — Title Page (~15 sec)

> Good morning everyone. Our presentation is on developing a dynamic predictive model for the growth of *Salmonella* Enteritidis in egg yolk, using the Gompertz model with the modified Ratkowsky secondary model.

---

## Slide 2 — Introduction & Background (~1.5 min)

> *Salmonella* Enteritidis is one of the leading causes of foodborne illness worldwide. Shell eggs are a major vehicle — SE penetrates the eggshell and grows rapidly in the iron-rich yolk.
>
> The USDA requires eggs be cooled to 7.2 degrees within 12 hours of laying. However, temperature fluctuates during cooling, storage, and distribution.
>
> Our objective is to develop a dynamic model that predicts SE growth under continuously changing temperatures — specifically a sinusoidal profile cycling between 3 and 43 degrees over 24 hours. This work is based on Gumudavelli et al., 2007.

---

## Slide 3 — Model Description (~1.5 min)

> The primary model is the Gompertz equation. At constant temperature it has an analytical solution; for varying temperature we use the differential form shown here. The parameters are A, the initial concentration; C, the growth range; and M, the inflection time.
>
> The secondary model is the modified Ratkowsky equation, which describes how the growth rate mu depends on temperature. Parameters a and b are regression coefficients; T-min and T-max are theoretical growth limits, fixed at 6 and 46.3 degrees.
>
> The dynamic model integrates both: at each time step, we compute mu from the current temperature and feed it into the Gompertz ODE, solved with ode45.

---

## Slide 4 — Data & Parameters (~45 sec)

> On the left, our data summary: 20 growth observations across 12 time points from 0 to 17 hours, and 47 temperature points following a sinusoidal wave. Temperature is interpolated via interp1.
>
> On the right is the parameter table. Five parameters are estimated, two are fixed from the literature. The initial guesses were provided by the instructor.

---

## Slide 5 — Forward Problem & SSC (~1.5 min)

> On the left, the forward prediction with initial guesses. The blue curve follows the general trend of the observed data, confirming our starting values are reasonable.
>
> On the right, the Scaled Sensitivity Coefficients. All five parameters show significant, distinct patterns — they are estimable. T-min and T-max showed negligible SSC, which is why we fixed them.

---

## Slide 6 — OLS Parameter Estimates (~1 min)

> Here are our estimation results. This table shows the estimated values, standard errors, relative errors, and 95% confidence intervals for each parameter. Below is the correlation matrix — we check for high correlations that might indicate identifiability issues.
>
> The RMSE and pseudo-R-squared are shown at the bottom. An RMSE below 0.3 log CFU per mL is within the precision of standard plating methods, indicating a good fit.

---

## Slide 7 — Goodness of Fit Summary (~30 sec)

> This table summarizes our key fit metrics: RMSE, pseudo-R-squared, MSE, and the number of data points. These confirm the model describes the data well.

---

## Slide 8 — Fitted Curve with CB & PB (~1 min)

> This figure shows the fitted curve with asymptotic confidence and prediction bands. The darker band is the 95% confidence band for the mean response. The lighter band is the prediction band, which is wider because it also accounts for measurement error. All observed points fall within the prediction band.

---

## Slide 9 — Residual Analysis (~1.5 min)

> On the left, the residual scatter plot — residuals appear randomly distributed around zero with roughly constant spread. On the right, the histogram is approximately bell-shaped.
>
> The table below lists the five standard statistical assumptions. The Durbin-Watson statistic tests for randomness and autocorrelation. The Shapiro-Wilk test checks normality. Results are filled in from our MATLAB output.

---

## Slide 10 — Final SSC (~45 sec)

> We recomputed the SSC with estimated parameters. The patterns are consistent with the initial SSC, confirming all five parameters remain identifiable at the solution. The table shows the maximum absolute SSC for each parameter.

---

## Slide 11 — Optimal Experimental Design (~1 min)

> On the left, the delta criterion — determinant of X-transpose-X — increases with time, meaning longer experiments provide more information. On the right, the C-ii curves decrease, meaning parameters become more precisely estimated. The key measurement window is the exponential growth phase, roughly 2 to 8 hours.

---

## Slide 12 — Bootstrap (~1 min)

> We performed residual resampling bootstrap with 1000 iterations. The first table compares bootstrap and asymptotic 95% confidence intervals — they are generally consistent. The second table compares average band widths. The figure below shows the bootstrap CB and PB.

---

## Slide 13 — Summary & Conclusions (~45 sec)

> To summarize: we successfully developed a dynamic Gompertz model for SE growth in egg yolk. Five parameters were estimated with good precision. All statistical assumptions were evaluated. Bootstrap confirmed the robustness of our estimates. The model supports the USDA recommendation to store eggs at or below 7.2 degrees.

---

## Slide 14 — References (~5 sec)

> Here are our references.

---

## Slide 15 — Thank You (~5 sec)

> Thank you for your attention. We're happy to take questions.
