#' Linear Regression Function
#'
#' Fits a linear regression model using OLS (basic) or QR decomposition (advanced).
#'
#' @param formula A formula object (e.g., y ~ x1 + x2).
#' @param data A data frame containing the variables.
#' @return An object of class "linreg" with regression results (S3 list or RC object).
#' @export
#' @examples
#' data(iris)
#' mod <- linreg(Petal.Length ~ Species, iris)
linreg <- function(formula, data) {
  # Extract X matrix and y vector
  X <- model.matrix(formula, data)
  dep_var <- all.vars(formula)[1]
  y <- data[[dep_var]]  # Dependent variable

  # --- Computations: Choose OLS (basic) or QR (*) ---
  # Basic: Ordinary Least Squares
  XtX_inv <- solve(t(X) %*% X)
  beta_hat <- XtX_inv %*% t(X) %*% y

  # --- Alternative: (*) QR Decomposition (uncomment to use instead of OLS) ---
  # qr_decomp <- qr(X)
  # Q <- qr.Q(qr_decomp)
  # R <- qr.R(qr_decomp)
  # beta_hat <- backsolve(R, t(Q) %*% y)  # Efficient for upper triangular R
  # XtX_inv <- solve(t(R) %*% R)  # For var_beta_hat below

  # Common computations (same for OLS or QR)
  y_hat <- X %*% beta_hat
  e <- y - y_hat
  n <- nrow(X)
  p <- ncol(X)
  df <- n - p
  sigma2_hat <- as.numeric((t(e) %*% e) / df)
  var_beta_hat <- sigma2_hat * XtX_inv
  se_beta <- sqrt(diag(var_beta_hat))  # Standard errors (for later use in methods)
  t_beta <- beta_hat / se_beta
  p_values <- 2 * pt(-abs(t_beta), df)

  # --- Store Results: Choose S3 (basic) or RC (*) ---
  # Basic: S3 class (list)
  result <- list(
    call = match.call(),  # Store the call for print()
    formula = formula,
    data_name = deparse(substitute(data)),
    beta_hat = beta_hat,
    y_hat = y_hat,
    e = e,
    df = df,
    sigma2_hat = sigma2_hat,
    var_beta_hat = var_beta_hat,
    se_beta = se_beta,
    t_beta = t_beta,
    p_values = p_values,
    coeff_names = colnames(X)  # For naming coefficients
  )
  class(result) <- "linreg"

  # --- Alternative: (*) RC class (uncomment to use instead of S3; define class first) ---
  # # First, define the RC class at the top of this file (outside the function):
  # linreg <- setRefClass("linreg",
  #   fields = list(
  #     call = "call",
  #     formula = "formula",
  #     data_name = "character",
  #     beta_hat = "matrix",
  #     y_hat = "numeric",
  #     e = "numeric",
  #     df = "numeric",
  #     sigma2_hat = "numeric",
  #     var_beta_hat = "matrix",
  #     se_beta = "numeric",
  #     t_beta = "matrix",
  #     p_values = "matrix",
  #     coeff_names = "character"
  #   )
  # )
  # # Then, in the function (replace the S3 return):
  # result <- linreg$new(
  #   call = match.call(),
  #   formula = formula,
  #   data_name = deparse(substitute(data)),
  #   beta_hat = beta_hat,
  #   y_hat = y_hat,
  #   e = e,
  #   df = df,
  #   sigma2_hat = sigma2_hat,
  #   var_beta_hat = var_beta_hat,
  #   se_beta = se_beta,
  #   t_beta = t_beta,
  #   p_values = p_values,
  #   coeff_names = colnames(X)
  # )

  return(result)
}
