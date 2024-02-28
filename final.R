#' CALCULATE CONFIDENCE INTERVALS FOR MEANS
#' WITH SAMPLE SIZE N = 1 OR N = 2 USING ZETA AND OTHER METHODS
#'
#'
#' Calculate Zeta Numerically
#'
#' Numerically solves the equation to calculate zeta.
#'
#' @param level The confidence level.
#' @return The calculated zeta.
#'
calculate_zeta_num <- function(level) {
  # Define the equation to solve: (2 * tau / zeta) + 1 - level = 0
  equation <- function(zeta, level) {
    (2 * dnorm(0) / zeta) + 1 - level
  }

  # Use numerical solver to find the root of the equation
  result <- uniroot(equation, interval = c(0.001, 10), level = level, extendInt = "yes")
  zeta <- result$root

  return(zeta)
}

result <- calculate_zeta_num(level = 0.95)
print(result)

#' Calculate Zeta Using Closed-Form Approximation
#'
#' Calculates zeta using a closed-form approximation.
#'
#' @param level The confidence level.
#' @return The calculated zeta.
#'
calculate_zeta_cf <- function(level) {
  # Implement the closed-form approximation
  tau <- dnorm(0)
  zeta <- (2 * tau / level) + 1

  return(zeta)
}

result <- calculate_zeta_cf(level = 0.95)
print(result)

#' Define new class
#'
#'
#'
#'
setClass("ci_result", representation(
  confidence_level = "numeric",
  ci_lower = "numeric",
  ci_upper = "numeric"
))

#' Calculate Confidence Interval for Normal Distribution (Sample Size = 1 or 2)
#'
#' This function calculates the confidence interval for a Normal distribution
#' when the sample size is either 1 or 2, based on the efforts of Wall, Boen, and Tweedie (2001).
#'
#' @param x The sample observation(s) in a Normal Distribution.
#' @param a A constant (default is 0) the user can include (optionally) to shift x arbitrarily.
#' @param level The confidence level (default is 0.95).
#' @return The confidence interval.
#' @export
#'
calculate_ci <- function(x, a = 0, level = 0.95, method = "numerical") {
    if (!(method %in% c("numerical", "closed-form"))) {
      stop("Invalid method. Choose 'numerical' or 'closed-form'")
    }

    # Calculate zeta based on the specified method
    zeta <- ifelse(method == "numerical", calculate_zeta_num(level),
                   calculate_zeta_cf(level))


  # if x is a single observation or a vector of observations
  if (length(x) == 1) {
    # For sample size 1
    tau <- dnorm(0)

    zeta <- (2 * tau / level) + 1
    ci_lower <- (x - a) - zeta * abs(x - a)
    ci_upper <- (x - a) + zeta * abs(x - a)
  }
  else if (length(x) == 2) {

    # For sample size 2, calculating the mean and standard deviation
    mean_x <- mean(x)
    sd_x <- sd(x)

    # t-distribution for sample size 2
    t_quantile <- qt((1 + level) / 2, df = 1)

    ci_lower <- mean_x - t_quantile * sd_x / sqrt(2)
    ci_upper <- mean_x + t_quantile * sd_x / sqrt(2)
  } else {
    stop("Sample size must be either 1 or 2.")
  }

  # Return result
  return(new("ci_result", confidence_level = level, ci_lower = ci_lower, ci_upper = ci_upper))
}

result <- calculate_ci (x = 4.6, a = 0, level = 0.95, method = "closed-form")
print(result)

#' Calculate CI and P-Value for Normal Distribution (n = 1)
#'
#' calculates the confidence interval and p-value
#' for a normal distribution with one observation
#'
#' @param x This is our only data point from a normal distribution.
#' @param a A constant (default is 0) the user can include (optionally) to shift x arbitrarily.
#' @param null This is the value we want to test against our hypothesis!
#' @param level Confidence level (between 0 and 1, default is 0.95).
#' @return A list with the confidence interval and p-value.
#' \code{CI = x ± zeta * |x|}, where zeta is from a standard normal distribution.
#' P-value is computed against the null hypothesis.
#'
#' @export
#'
calculate_ci_p <- function(x, a, null = 0, level = 0.90) {
  if (level < 0 || level > 1) stop(" the confidence level is between 0 and 1.")

  # To find  zeta from the equation below
  tau <- dnorm(0)  # Standard normal density
  zeta <- (2 * dnorm(0) / level) + 1

  # Obtain the confidence interval
  ci_lower <- (x - a) - zeta * abs(x - a)
  ci_upper <- (x - a) + zeta * abs(x - a)

  # Find the p-value against our null hypothesis
  p_value <- 2 * pnorm(q = abs(null - x), lower.tail = FALSE)

  # return
  return(list(ci_lower = ci_lower, ci_upper = ci_upper, p_value = p_value))
}

result <- calculate_ci_p(x = 3.8, a = 0, null = 0, level = 0.90)
print(result)

#' ALTERNATIVE METHODS TO FINDING CONFIDENCE INTERVALS
#' BASED ON PORTNOY, DISGUPTA (2019)
#'
#' Calculate Confidence Interval I1
#'
#' Calculates the confidence interval I1 based on a single observation X.
#'
#' @param x The single observation.
#' @param c The constant parameter.
#' @return The confidence interval I1.
#'
calculate_ci_I1 <- function(x, c) {
  ci_lower <- x - c * abs(x)
  ci_upper <- x + c * abs(x)
  return(c(ci_lower, ci_upper))
}

#' Coverage Probability for I1
#'
#' Calculates the coverage probability P1 for confidence interval I1.
#'
#' @param lambda The parameter lambda.
#' @param c The constant parameter.
#' @return The coverage probability P1.
#'
calculate_cov_P1 <- function(lambda, c) {
  term1 <- pnorm((c / (c + 1)) * lambda)
  term2 <- 1 - pnorm((c / (c - 1)) * lambda)
  return(term1 + term2)
}

#' Optimal Lambda for I1
#'
#' Calculates the optimal lambda for confidence interval I1.
#'
#' @param c The constant parameter.
#' @return The optimal lambda for I1.
#'
calculate_ol_I1 <- function(c) {
  lambda_star <- ((c^2 - 1) / (sqrt(2) * c^(3/2))) * sqrt(log(c + 1 / c - 1))
  return(lambda_star)
}

# Example usage:
x <- 3.5
c_value <- 1.96 # Adjust as needed
lambda_star_I1 <- calculate_ol_I1(c_value)
ci_I1 <- calculate_ci_I1(x, c_value)
cov_P1 <- calculate_cov_P1(lambda_star_I1, c_value)

# Check coverage probabilities and intervals for I1
print(ci_I1)
print(cov_P1)

#' Calculate Confidence Interval I2
#'
#' Calculates the confidence interval I2 based on a single observation X.
#'
#' @param x The single observation.
#' @param c The constant parameter.
#' @return The confidence interval I2.
#'
calculate_ci_I2 <- function(x, c) {
  ci_lower <- -c * abs(x)
  ci_upper <- c * abs(x)
  return(c(ci_lower, ci_upper))
}

#' Coverage Probability for I2
#'
#' Calculates the coverage probability P2 for confidence interval I2.
#'
#' @param lambda The parameter lambda.
#' @param c The constant parameter.
#' @return The coverage probability P2.
#'
calculate_cov_P2 <- function(lambda, c) {
  term1 <- pnorm((c - 1) / c * lambda)
  term2 <- 1 - pnorm((c + 1) / c * lambda)
  return(term1 + term2)
}

#' Optimal Lambda for I2
#'
#' Calculates the optimal lambda for confidence interval I2.
#'
#' @param c The constant parameter.
#' @return The optimal lambda for I2.
#'
calculate_ol_I2 <- function(c) {
  lambda_star <- sqrt((c / 2) * log(c + 1 / c - 1))
  return(lambda_star)
}

# Example usage:
x <- 3.5
c_value <- 1.96 # Adjust as needed
lambda_star_I2 <- calculate_ol_I2(c_value)
ci_I2 <- calculate_ci_I1(x, c_value)
cov_P2 <- calculate_cov_P1(lambda_star_I1, c_value)

# Check coverage probabilities and intervals for I2
print(ci_I2)
print(cov_P2)

#' Plot Method
#'
#' Generates a plot based on the ci_result class.
#'
#' @param x An object of class ci_result
#' @return NULL (used for side effect of plotting)
#' @export
#'
plot.ci_result <- function(x) {
  # Results from ci_result object
  level <- x@confidence_level
  ci_lower <- x@ci_lower
  ci_upper <- x@ci_upper

  # Calculate the midpoint of the confidence interval
  midpoint <- (ci_lower + ci_upper) / 2

  # Create a plot
  plot(midpoint, type = "n", ylim = c(min(ci_lower), max(ci_upper)),
       main = "Confidence Interval",
       xlab = "Observation", ylab = "Confidence Interval")

  # Points for confidence interval values
  points(1, ci_lower, col = "blue", pch = 16)
  points(1, ci_upper, col = "blue", pch = 16)

  # Point for the Midpoint
  points(1, midpoint, col = "red", pch = 19)

  # Add labels for ci_lower, midpoint, and ci_upper
  text(1, ci_lower, labels = sprintf("Lower Value: %.2f", ci_lower), pos = 3, col = "blue")
  text(1, midpoint, labels = sprintf("Midpoint: %.2f", midpoint), pos = 3, col = "red")
  text(1, ci_upper, labels = sprintf("Upper Value: %.2f", ci_upper), pos = 1, col = "blue")
}

#' Summary Method for ci_result Class
#'
#' Generates a summary based on the ci_result class.
#'
#' @param y An object of class ci_result.
#' @return A summary of the ci_result object.
#' @export
#'
summary.ci_result <- function(y) {
  # Results from ci_result object
  level <- y@confidence_level
  ci_lower <- y@ci_lower
  ci_upper <- y@ci_upper

  # Calculate the midpoint (original x)
  midpoint <- (ci_lower + ci_upper) / 2

  cat("Summary of Confidence Interval Object \n")
  cat("Confidence Level:", level, "\n")
  cat("Lower Confidence Interval:", ci_lower, "\n")
  cat("Upper Confidence Interval:", ci_upper, "\n")
  cat("Midpoint:", midpoint, "\n")

  invisible(NULL)
}

#Example
result <- calculate_ci (x = 10, a = 0, level = 0.95)
summary(result)
plot(result)

result <- calculate_ci(x = c(5.3, 4.6), a = 0, level = 0.95)
summary(result)
plot(result)


