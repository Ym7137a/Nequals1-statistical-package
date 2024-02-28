The purpose of this package is to guide the user through the calculation of confidence intervals for Normal distributions with a sample size of n = 1 or n = 2.  This work directly references the work of Wall, Boen, and Tweedie (2001) and is based off of the following form for constructing confidence intervals where n = 1:

x ± (zeta) * abs(x).

x in this case represents the single observation, and zeta can be approximated using the following formula:

(zeta) = (2(tau)/level)) + 1

where tau is the standard normal density.

For this package, the user can calculate zeta one of two ways, by either numerically solving an equation or by closed-form approximation.

The "numerical" method essentially takes an equation and sets it to zero, and then relies on the uniroot function to solve for zeta.

The "closed-form" method solves for zeta using the (2(tau)/level)) + 1 equation.

After solving for zeta, the user can then use the following function to construct a confidence interval:

calculate_ci <- function(x, a = 0, level = 0.95).

These are the default parameters and their values.  In this case, x is the single observation for n = 1, a is an arbitrary value the user can adjust the equation by, and the level is the confidence level used to construct the interval.  This function creates a confidence interval using the following boundaries:

ci_lower <- (x - a) - zeta * abs(x - a)
ci_upper <- (x - a) + zeta * abs(x - a)

The package also allows for the user to construct a confidence interval for when n = 2.  When switching to a two-observation distribution, the function calculates using a new function that more closely resembles a normal confidence interval calculation.  The user should use the following syntax for this operation:

calculate_ci <- function(x = c(x1, x2), a = 0, level = 0.95).

The value for a should be set to zero, since there is no shifting of the value.

The package can also produce a p-value against a null value of zero for hypothesis construction, using the following function:

calculate_ci_p <- function(x, a, null = 0, level = 0.90).

Lastly, the package can produce a summary and basic plot of the confidence interval by using summary() and plot() functions against the ci_result class.

Any questions or suggestions for improvement (and they are always up for consideration), please contact Kevin Norris at kn7295a@american.edu, Yashvi Malviya at ym7137a@american.edu, or Mia Monintja at am2366a@american.edu.  Thank you.



