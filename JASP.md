#### [JASP](https://jasp-stats.org/getting-started/) multiple linear regression - *13 May 2026*
In JASP top row, chose Regression, then Classical -> Linear Regression
- in Linear Regression window, put one variable in `Dependent Variable`,
	- Method `Backward` allows JASP to prune useless independent variables
	- then put independent variables in `Covariates`, *NOT* `Factors`
- in Model, confirm independent variables are in Model 0
	- set `Include intercept`
- in Statistics, set
	- Model Summary: R squared change, F change
	- Coeffient: Estimates, 95% Confidence intervals
	- Display: Model fit, Descriptives, Regression equation
- in Method Specification
	- Use p value
- in Plots
	- Residuals vs. dependent
