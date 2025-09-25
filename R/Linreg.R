#' @importFrom stats model.matrix pt
#' @import methods
#' @importFrom gridExtra grid.arrange
#' @import ggplot2

# Define the RC (Reference Class) with a distinct name to avoid conflict
LinregClass <- setRefClass("LinregClass",
                           fields = list(
                             call = "call",                # Stores the original function call (like lm())
                             formula = "formula",          # The input formula
                             data_name = "character",      # Name of the data frame
                             beta_hat = "matrix",          # Estimated regression coefficients (column matrix)
                             y_hat = "numeric",            # Fitted values
                             e = "numeric",                # Residuals
                             df = "numeric",               # Degrees of freedom
                             st_e = "numeric",             # Standardized residuals
                             sigma2_hat = "numeric",       # Residual variance (sigma squared)
                             var_beta_hat = "matrix",      # Variance-covariance matrix of coefficients
                             se_beta = "numeric",          # Standard errors of coefficients
                             t_beta = "matrix",            # t-values for coefficients (column matrix)
                             p_values = "matrix",          # p-values for coefficients (column matrix)
                             coeff_names = "character"     # Names of coefficients (e.g., "(Intercept)", predictors)
                           ),
                           methods = list(
                             # Print method to mimic lm() output
                             print = function() {
                               cat("Call:\n")
                               base::print(.self$call)
                               cat("\nCoefficients:\n")
                               base::print(as.numeric(.self$beta_hat))
                               invisible(.self)
                             },

                             # Plot method with two ggplot2 plots using grid.arrange
                             plot = function() {
                               med_resid <- stats::median(.self$e, na.rm = TRUE)
                               med_sqrt <- stats::median(sqrt(abs(.self$st_e)), na.rm = TRUE)

                               # Plot 1: Residuals vs Fitted with median line and loess smooth
                               plot1 <- ggplot(data.frame(fitted = .self$y_hat, residuals = .self$e), aes(x = fitted, y = residuals)) +
                                 geom_point() +
                                 geom_hline(yintercept = med_resid, color = "red", linetype = "dashed") +
                                 geom_smooth(method = "loess", se = FALSE, color = "blue") +
                                 ggtitle("Residuals vs Fitted") +
                                 xlab("Fitted values") +
                                 ylab("Residuals")

                               # Plot 2: Scale-Location with loess smooth
                               plot2 <- ggplot(data.frame(fitted = .self$y_hat, std_res_sqrt = sqrt(abs(.self$st_e))), aes(x = fitted, y = std_res_sqrt)) +
                                 geom_point() +
                                 geom_smooth(method = "loess", se = FALSE, color = "red") +
                                 ggtitle("Scale-Location") +
                                 xlab("Fitted values") +
                                 ylab("sqrt(|Standardized residuals|)")

                               grid.arrange(plot1, plot2, ncol = 1)
                               invisible(.self)
                             },

                             # Residuals method to return the vector of residuals
                             resid = function() {
                               return(.self$e)
                             },

                             # Predicted values method to return the fitted values
                             pred = function() {
                               return(.self$y_hat)
                             },

                             # Coefficients method to return a named vector of coefficients
                             coef = function() {
                               coeffs <- as.vector(.self$beta_hat)
                               names(coeffs) <- .self$coeff_names
                               return(coeffs)
                             },

                             # Summary method to return a table similar to lm()
                             summary = function() {
                               coeff_table <- cbind(
                                 Estimate = as.vector(.self$beta_hat),
                                 `Std. Error` = .self$se_beta,
                                 `t value` = as.vector(.self$t_beta),
                                 `Pr(>|t|)` = as.vector(.self$p_values)
                               )
                               rownames(coeff_table) <- .self$coeff_names

                               cat("Coefficients:\n")
                               base::print(coeff_table)
                               cat("\nResidual standard error:", sqrt(.self$sigma2_hat), "on", .self$df, "degrees of freedom\n")
                               invisible(.self)
                             }
                           )
)

#' Linear Regression Function
#'
#' Fits a linear regression model using OLS (Ordinary Least Squares).
#' Returns an RC (Reference Class) object of class "LinregClass".
#'
#' @param formula A formula object (e.g., y ~ x1 + x2).
#' @param data A data frame containing the variables.
#' @return An object of class "LinregClass" (RC) with regression results.
#' @export
#' @examples
#' data(iris)
#' mod <- linreg(Petal.Length ~ Species, iris)
linreg <- function(formula, data) {
  # --- Input Handling ---
  X <- model.matrix(formula, data)
  dep_var <- all.vars(formula)[1]
  y <- data[[dep_var]]

  # --- Computations: Ordinary Least Squares (OLS) ---
  XtX_inv <- solve(t(X) %*% X)
  beta_hat <- XtX_inv %*% t(X) %*% y

  # --- Common Computations (Shared Post-Estimation) ---
  y_hat <- as.numeric(X %*% beta_hat)
  e <- y - y_hat

  n <- nrow(X)
  p <- ncol(X)
  df <- n - p
  sigma2_hat <- as.numeric((t(e) %*% e) / df)
  var_beta_hat <- sigma2_hat * XtX_inv
  st_e <- e / (sqrt(sigma2_hat) * sqrt(1 - diag(X %*% XtX_inv %*% t(X))))
  se_beta <- sqrt(diag(var_beta_hat))
  t_beta <- beta_hat / se_beta
  p_values <- 2 * pt(-abs(t_beta), df)

  # --- Store Results in RC Object ---
  result <- LinregClass$new(
    call = match.call(),
    formula = formula,
    data_name = deparse(substitute(data)),
    beta_hat = beta_hat,
    y_hat = y_hat,
    e = e,
    st_e = st_e,
    df = df,
    sigma2_hat = sigma2_hat,
    var_beta_hat = var_beta_hat,
    se_beta = se_beta,
    t_beta = t_beta,
    p_values = p_values,
    coeff_names = colnames(X)
  )

  return(result)
}
