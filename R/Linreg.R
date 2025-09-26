#' @importFrom stats model.matrix pt model.frame model.response terms
#' @import methods
#' @importFrom gridExtra grid.arrange
#' @import ggplot2

linreg <- setRefClass(
  "linreg",
  fields = list(
    call          = "call",
    formula       = "formula",
    data          = "data.frame",
    data_name     = "character",
    beta_hat      = "numeric",
    y_hat         = "numeric",
    e             = "numeric",
    df            = "numeric",
    st_e          = "numeric",
    sigma2_hat    = "numeric",
    var_beta_hat  = "matrix",
    se_beta       = "numeric",
    t_beta        = "numeric",
    p_values      = "numeric",
    coeff_names   = "character"
  ),
  methods = list(

    initialize = function(formula, data, ...) {
      callSuper(...)

      data_sym        <- substitute(data)
      .self$data_name <- deparse(data_sym)
      .self$formula   <- formula
      .self$data      <- data
      .self$call <- as.call(list(
        as.name("linreg"),
        formula = formula,
        data    = data_sym
      ))


      mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
      tt <- stats::terms(mf)
      X  <- stats::model.matrix(tt, mf)
      y  <- as.numeric(stats::model.response(mf))

      n <- nrow(X); p <- ncol(X)
      XtX     <- crossprod(X)
      XtX_inv <- solve(XtX)
      Xty     <- crossprod(X, y)

      ## Estimates
      beta_hat_ <- as.numeric(solve(XtX, Xty))
      names(beta_hat_) <- colnames(X)

      y_hat_   <- as.numeric(X %*% beta_hat_)
      e_       <- y - y_hat_
      df_      <- n - p
      sigma2_  <- as.numeric(crossprod(e_) / df_)
      vcov_    <- sigma2_ * XtX_inv
      se_      <- sqrt(diag(vcov_))
      t_       <- beta_hat_ / se_
      p_       <- 2 * stats::pt(abs(t_), df = df_, lower.tail = FALSE)

      ## standardized residuals
      Hdiag <- rowSums((X %*% XtX_inv) * X)
      den   <- sqrt(pmax(1e-12, 1 - Hdiag))
      st_e_ <- e_ / (sqrt(sigma2_) * den)

      ## Store into fields
      .self$beta_hat     <- beta_hat_
      .self$y_hat        <- y_hat_
      .self$e            <- e_
      .self$df           <- df_
      .self$sigma2_hat   <- sigma2_
      .self$var_beta_hat <- vcov_
      .self$se_beta      <- se_
      .self$t_beta       <- t_
      .self$p_values     <- p_
      .self$st_e         <- st_e_
      .self$coeff_names  <- colnames(X)

      invisible(.self)
    },

    print = function() {
      cat("Call:\n"); base::print(.self$call)
      cat("\nCoefficients:\n"); base::print(.self$beta_hat)
      invisible(.self)
    },

    resid = function() .self$e,
    pred  = function() .self$y_hat,
    coef  = function() stats::setNames(.self$beta_hat, .self$coeff_names),

    summary = function() {
      stars <- function(p)
        ifelse(p < 0.001,"***",
               ifelse(p < 0.01,"**",
                      ifelse(p < 0.05,"*",
                             ifelse(p < 0.1,"."," "))))

      tbl <- data.frame(
        Estimate     = .self$beta_hat,
        `Std. Error` = .self$se_beta,
        `t value`    = .self$t_beta,
        `Pr(>|t|)`   = .self$p_values,
        check.names  = FALSE
      )
      tbl$` ` <- stars(tbl$`Pr(>|t|)`)
      rownames(tbl) <- .self$coeff_names

      cat("Coefficients:\n"); base::print(tbl)
      cat(sprintf(
        "\nResidual standard error: %.6f on %d degrees of freedom\n",
        sqrt(.self$sigma2_hat), as.integer(.self$df)
      ))
      invisible(tbl)
    },

    plot = function() {
      med_resid <- stats::median(.self$e, na.rm = TRUE)
      plot1 <- ggplot2::ggplot(
        data.frame(fitted = .self$y_hat, residuals = .self$e),
        ggplot2::aes(x = fitted, y = residuals)
      ) +
        ggplot2::geom_point() +
        ggplot2::geom_hline(yintercept = med_resid, color = "red", linetype = "dashed") +
        ggplot2::geom_smooth(method = "loess", se = FALSE, color = "blue", formula = y ~ x) +
        ggplot2::ggtitle("Residuals vs Fitted") +
        ggplot2::xlab("Fitted values") +
        ggplot2::ylab("Residuals")

      plot2 <- ggplot2::ggplot(
        data.frame(fitted = .self$y_hat, std_res_sqrt = sqrt(abs(.self$st_e))),
        ggplot2::aes(x = fitted, y = std_res_sqrt)
      ) +
        ggplot2::geom_point() +
        ggplot2::geom_smooth(method = "loess", se = FALSE, color = "red", formula = y ~ x) +
        ggplot2::ggtitle("Scale-Location") +
        ggplot2::xlab("Fitted values") +
        ggplot2::ylab("sqrt(|Standardized residuals|)")

      gridExtra::grid.arrange(plot1, plot2, ncol = 1)
      invisible(.self)
    }
  )
)

#' Fit a linear regression using OLS (wrapper around the RC class)
#'
#' This convenience function constructs and returns an object of class **"linreg"**
#' by calling its generator with `$new()`.
#'
#' @param formula A formula like \code{y ~ x1 + x2}.
#' @param data A \code{data.frame} containing the variables.
#' @return An object of class \code{"linreg"} (RC) with fields and methods.
#' @examples
#' data(iris)
#' m <- linreg_fit(Petal.Length ~ Sepal.Width + Sepal.Length, iris)
#' m$summary()
#' @export
linreg_fit <- function(formula, data) {
  linreg$new(formula = formula, data = data)
}

