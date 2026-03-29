# Salmonella Growth Model - Final Project

Dynamic Gompertz Model for predicting *Salmonella* Enteritidis growth in egg yolk under sinusoidal temperature profiles.

## Files

- `final_project.m` — Main MATLAB script (ODE model + SSC analysis + nlinfit parameter estimation + plotting)
- `Salmonella sin growth.xlsx` — Experimental growth data (CFU/mL vs time)
- `Salmonella sin growth Temps.xlsx` — Temperature profile data (sinusoidal)
- `global_example.m` — Instructor-provided template for global variable usage

## Usage

1. Open MATLAB, set current directory to this repo
2. Run `final_project.m`
3. Examine SSC plot to verify which parameters are estimable
4. Adjust `est_idx` / `fix_idx` if needed, re-run

## Model

- **Primary**: Gompertz differential form — `dy/dt = -μ(T)·e·(y-A)/C·ln((y-A)/C)`
- **Secondary**: Modified Ratkowsky — `μ(T) = a·(T-Tmin)²·(1-exp(b·(T-Tmax)))`
- **Numerical solver**: `ode45` with `interp1` for temperature interpolation

## Reference

Gumudavelli V, Subbiah J, Thippareddi H, Velugoti PR, Froning G. 2007. Dynamic predictive model for growth of *Salmonella* Enteritidis in egg yolk. *J. Food Sci.* 72:M254–62.
