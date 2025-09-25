#' @importFrom stats model.matrix pt
#' @import methods

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
                             sigma2_hat = "numeric",       # Residual variance (sigma squared)
                             var_beta_hat = "matrix",      # Variance-covariance matrix of coefficients
                             se_beta = "numeric",          # Standard errors of coefficients
                             t_beta = "matrix",            # t-values for coefficients (column matrix)
                             p_values = "matrix",          # p-values for coefficients (column matrix)
                             coeff_names = "character"     # Names of coefficients (e.g., "(Intercept)", predictors)
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
  # Create the design matrix X (independent variables + intercept)
  X <- model.matrix(formula, data)
  # Extract the dependent variable name and vector y
  dep_var <- all.vars(formula)[1]
  y <- data[[dep_var]]

  # --- Computations: Ordinary Least Squares (OLS) ---
  # Compute (X'X)^(-1), the inverse of X transpose times X
  XtX_inv <- solve(t(X) %*% X)
  # Regression coefficients: beta_hat = (X'X)^(-1) * X'y
  beta_hat <- XtX_inv %*% t(X) %*% y

  # --- Common Computations (Shared Post-Estimation) ---
  # Fitted values: y_hat = X * beta_hat
  y_hat <- as.numeric(X %*% beta_hat)  # Flatten to numeric vector
  # Residuals: e = y - y_hat
  e <- y - y_hat
  # Number of observations (n) and parameters (p)
  n <- nrow(X)
  p <- ncol(X)
  # Degrees of freedom: df = n - p
  df <- n - p
  # Residual variance: sigma2_hat = (e'e) / df
  sigma2_hat <- as.numeric((t(e) %*% e) / df)
  # Variance-covariance of betas: var_beta_hat = sigma2_hat * (X'X)^(-1)
  var_beta_hat <- sigma2_hat * XtX_inv
  # Standard errors: se_beta = sqrt(diag(var_beta_hat))
  se_beta <- sqrt(diag(var_beta_hat))
  # t-values: t_beta = beta_hat / se_beta
  t_beta <- beta_hat / se_beta
  # p-values: 2 * pt(-|t_beta|, df) for two-tailed test
  p_values <- 2 * pt(-abs(t_beta), df)

  # --- Store Results in RC Object ---
  # Create a new instance of the LinregClass RC class
  # Populate all fields with the computed values
  result <- LinregClass$new(
    call = match.call(),           # Capture the full call for print() method
    formula = formula,             # Store formula
    data_name = deparse(substitute(data)),  # Store data name as string
    beta_hat = beta_hat,           # Coefficients matrix
    y_hat = y_hat,                 # Fitted values vector
    e = e,                         # Residuals vector
    df = df,                       # Degrees of freedom scalar
    sigma2_hat = sigma2_hat,       # Residual variance scalar
    var_beta_hat = var_beta_hat,   # Variance matrix
    se_beta = se_beta,             # Standard errors vector
    t_beta = t_beta,               # t-values matrix
    p_values = p_values,           # p-values matrix
    coeff_names = colnames(X)      # Coefficient names (from X columns)
  )

  # Return the RC object
  return(result)
}
