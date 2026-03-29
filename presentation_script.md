# Presentation Script (English)

Estimated total: ~13 minutes. Each slide ~1 minute unless noted.

---

## Slide 1 — Title Page (~15 sec)

> Good morning/afternoon everyone. Our presentation today is on developing a dynamic predictive model for the growth of *Salmonella* Enteritidis in egg yolk, using the Gompertz model combined with the modified Ratkowsky secondary model.

---

## Slide 2 — Introduction & Background (~1.5 min)

> *Salmonella* Enteritidis, or SE, is one of the leading causes of foodborne illness worldwide. Shell eggs are a major transmission vehicle — SE can penetrate the eggshell, pass through the albumen, and grow rapidly in the iron-rich yolk.
>
> The USDA requires that shell eggs be cooled to 7.2 degrees Celsius within 12 hours of laying. However, in practice, temperature fluctuates significantly during cooling, storage, and distribution.
>
> This creates a need for a dynamic model — one that can predict SE growth under continuously changing temperature conditions, not just constant temperatures.
>
> Our objective is to develop and validate such a model using the Gompertz equation as the primary model and the modified Ratkowsky equation as the secondary model, applied to a sinusoidal temperature profile ranging from about 3 to 43 degrees Celsius over a 24-hour cycle.
>
> Our work is based on the paper by Gumudavelli and colleagues, published in the Journal of Food Science in 2007.

---

## Slide 3 — Model Description (~1.5 min)

> Let me walk you through the model structure.
>
> The primary model is the Gompertz equation. At constant temperature, it has an analytical solution shown here on the left. The key parameters are A, the initial microbial concentration in log scale; C, the total growth range; and M, the time at which the growth rate is maximum — the inflection point.
>
> For varying temperature, we need the differential form: dy/dt equals negative mu times e times the ratio (y minus A) over C, multiplied by the natural log of that ratio. Here, mu is the growth rate, which changes with temperature.
>
> The secondary model, shown on the right, is the modified Ratkowsky equation. It describes how the growth rate mu depends on temperature. The parameters a and b are regression coefficients, and T-min and T-max are the theoretical minimum and maximum temperatures for growth. We fix these at 6 and 46.3 degrees based on the literature.
>
> The dynamic model integrates both: at each time step, we get the current temperature, compute mu from the Ratkowsky equation, and feed it into the Gompertz ODE, which we solve numerically using MATLAB's ode45.

---

## Slide 4 — Experimental Growth Data (~45 sec)

> This table shows our experimental growth data from Gumudavelli et al., Figure 5c. The data was collected under a sinusoidal temperature profile.
>
> As you can see, at time zero the SE concentration is about 304 CFU per mL. By 17 hours it has grown to nearly 400 million. The first four time points have single measurements, while from 3 hours onward we have duplicate measurements, giving us 20 total data points.
>
> Note the rapid growth between 2 and 6 hours — this corresponds to the high-temperature phase of the sinusoidal cycle.

## Slide 5 — Temperature Data (~30 sec)

> This table shows selected points from our temperature profile. The full dataset has 47 points sampled at roughly 28.8-minute intervals over 22 hours.
>
> The temperature follows a sinusoidal wave, cycling between about 7 and 43 degrees Celsius. We use MATLAB's interp1 function to interpolate the temperature at any time the ODE solver needs it.

---

## Slide 6 — Parameters & Initial Guesses (~30 sec)

> Here is a summary of our seven parameters. We estimate five of them — A, C, M, a, and b — and fix T-min at 6 degrees and T-max at 46.3 degrees. The initial guesses shown here were provided by the instructor. We'll verify these are reasonable in the next slide.

---

## Slide 7 — Forward Problem & SSC (~1.5 min)

> On the left, you can see the forward problem result. Using our initial guesses, we solved the ODE and plotted the predicted log N against the observed data. The black dots are observations, the blue line is our prediction, and the red dashed line shows the temperature profile on the right axis. The prediction follows the general trend of the data, confirming our initial guesses are in the right ballpark.
>
> On the right is the Scaled Sensitivity Coefficient plot. SSC measures how sensitive the model output is to each parameter. A large, distinct SSC means the parameter has a strong influence on the output and can be reliably estimated.
>
> All five parameters — A, C, M, a, and b — show significant and distinct SSC patterns, meaning they are all estimable. T-min and T-max, which we tested separately, showed negligible sensitivity from this single temperature profile, which is why we fixed them.
>
> In terms of ranking, C and A have the largest SSC magnitudes, so they will be estimated most accurately.

---

## Slide 8 — OLS Results (~1 min)

> We estimated the parameters using MATLAB's nlinfit function, which performs nonlinear least squares with the Levenberg-Marquardt algorithm.
>
> This table shows the estimated values, standard errors, relative errors, and 95% confidence intervals for each parameter. The relative errors indicate how precisely each parameter is estimated — lower is better.
>
> The RMSE is reported at the bottom. A value below 0.3 log CFU per mL is within the precision of traditional microbial plating methods, indicating a good fit. The pseudo-R-squared value also confirms the model explains the data well.

---

## Slide 9 — Fitted Curve with CB & PB (~1 min)

> This figure shows the fitted curve along with the asymptotic confidence and prediction bands.
>
> The dark blue shaded region is the 95% confidence band, which represents uncertainty in the mean response — essentially, where we expect the true mean curve to lie. The lighter blue region is the 95% prediction band, which is wider because it also accounts for measurement error — it's where we'd expect a new individual observation to fall.
>
> The observed data points fall within the prediction band, and the predicted curve tracks the observations closely. The temperature profile is shown on the right axis for reference.

---

## Slide 10 — Residual Analysis (~1.5 min)

> Now let's check whether our model meets the standard statistical assumptions.
>
> On the left is the residual scatter plot — residuals versus time. We're looking for a random scatter around zero with no systematic pattern. The residuals appear randomly distributed with roughly constant spread, which is a good sign.
>
> On the right is the residual histogram. We're checking whether the residuals are approximately normally distributed. The shape is roughly bell-shaped.
>
> At the bottom, we list the five standard statistical assumptions:
> - First, the model is correct — confirmed visually from the fitted curve.
> - Second, errors are random — the Durbin-Watson statistic is close to 2, indicating no significant autocorrelation.
> - Third, constant variance — the residual scatter plot shows no obvious fanning or pattern.
> - Fourth, errors are uncorrelated — consistent with the Durbin-Watson result.
> - Fifth, normal distribution — the Shapiro-Wilk test p-value is reported; a value above 0.05 means we cannot reject normality.

---

## Slide 11 — Final SSC (~45 sec)

> Here we recomputed the Scaled Sensitivity Coefficients using the estimated parameter values, not the initial guesses.
>
> The purpose is to confirm that the parameters are still identifiable at the solution. The SSC patterns are qualitatively similar to the initial ones, which is reassuring — it means the optimization converged to a meaningful solution and all five parameters remain estimable.
>
> The ranking of parameter accuracy based on the final SSC and standard errors is consistent with what we saw earlier.

---

## Slide 12 — Optimal Experimental Design (~1 min)

> For optimal experimental design, we computed two criteria.
>
> On the left is the delta criterion — the determinant of X-transpose-X — plotted against the last measurement time. This measures the overall information content of the experiment. We can see it increases with time, indicating that longer experiments provide more information for parameter estimation.
>
> On the right are the C-ii curves — the diagonal elements of the inverse of X-transpose-X. These are proportional to the variance of each parameter estimate. As we include more measurement times, the C-ii values decrease, meaning the parameters become more precisely estimated.
>
> The key takeaway is that measurements in the exponential growth phase — roughly between 2 and 8 hours — contribute the most information.

---

## Slide 13 — Bootstrap (~1 min)

> Finally, we performed a residual resampling bootstrap with 1000 iterations to obtain non-parametric confidence intervals.
>
> The table compares the bootstrap 95% CI with the asymptotic CI for each parameter. They are generally consistent, which validates our asymptotic results.
>
> The figure shows the bootstrap confidence and prediction bands. The green region is the bootstrap CB and the pink region is the bootstrap PB.
>
> Comparing band widths: the bootstrap bands are slightly wider/narrower than the asymptotic bands, which is typical because bootstrap captures non-linear effects that the asymptotic approximation may miss.

---

## Slide 14 — Summary & Conclusions (~45 sec)

> To summarize:
>
> We successfully developed a dynamic Gompertz model for SE growth in egg yolk under sinusoidal temperature conditions. Five parameters were estimated via ordinary least squares, with two fixed from the literature.
>
> The model fits the data well, with RMSE below 0.3 log CFU per mL. All five standard statistical assumptions were evaluated. The bootstrap analysis confirmed the robustness of our parameter estimates.
>
> From a practical standpoint, this model can predict SE growth risk during egg cooling, storage, and distribution, and supports the USDA recommendation to store eggs at or below 7.2 degrees Celsius.
>
> Thank you.

---

## Slide 15 — References (~5 sec)

> Here are our references for your review.

---

## Slide 16 — Thank You (~5 sec)

> Thank you for your attention. We're happy to take any questions.
