#' Residualize one feature on the random effects (two-step modelling)
#'
#' Fits a null model containing only the random intercepts, \code{feature ~ (1 | re)},
#' and returns its cell-level residuals, that is, the feature with the
#' random-effect structure (and its mean) removed. If the model cannot be fitted
#' the feature is simply mean-centred.
#'
#' @param targetExpr Numeric vector with the values of the feature across cells.
#' @param metadata Data.frame with cell metadata. Row names must be cell IDs.
#' @param randomEffects Character vector of \code{metadata} columns used
#'  as random intercepts.
#' @param targetOmicType Numeric flag: \code{0} = Gaussian (\code{lmer}),
#'   \code{1} = binary (\code{glmer}).
#' @param REML Logical. Restricted maximum likelihood (Gaussian models only).
#' @param verbose Logical. If \code{TRUE}, prints optimization progress.
#'
#' @return Numeric vector of cell-level residuals. The attribute \code{fallback}
#'   is \code{TRUE} when the model failed and the mean-centred values were
#'   returned instead.
#'
#' @importFrom stats residuals as.formula
#' @keywords internal
fitTwoStep <- function(targetExpr,
                            metadata,
                            randomEffects,
                            targetOmicType,
                            REML,
                            verbose) {

  # ~~~~~~~~~~ 1. Create data.frame ~~~~~~~~~~ #
  model_data <- data.frame(
    TargetExpression = as.numeric(targetExpr),
    metadata[, randomEffects, drop = FALSE]
  )

  re_terms      <- paste0("(1 | ", randomEffects, ")", collapse = " + ")
  model_formula <- stats::as.formula(paste("TargetExpression ~", re_terms), env = globalenv())

  # ~~~~~~~~~~ 2. Fit the null model ~~~~~~~~~~ #
  fit_res <- runGLMM(
    formula        = model_formula,
    data           = model_data,
    targetOmicType = targetOmicType,
    REML           = REML,
    verbose        = verbose
  )

  # ~~~~~~~~~~ 3. Residuals (or centred values as fallback) ~~~~~~~~~~ #
  if (!is.null(fit_res$fit)) {
    resid <- as.numeric(stats::residuals(fit_res$fit))
    attr(resid, "fallback") <- FALSE
  } else {
    resid <- as.numeric(targetExpr) - mean(targetExpr, na.rm = TRUE)
    attr(resid, "fallback") <- TRUE
  }
  resid
}


#' Residualize every regulatory layer on the random effects (step 1)
#'
#' For each regulator of each layer fits a null mixed model with
#' \code{\link{fitTwoStep}} and replaces the regulator by its residuals.
#' This removes the sample (donor, batch) structure from the regulators before
#' the target-level models of step 2 are fitted.
#'
#' @details
#' Sparse layers are transposed once so that extracting the values of a
#' regulator is cheap. The residual matrices are dense. Note that residualizing
#' on a random intercept removes the between-sample differences of the regulator,
#' so the regulator effects estimated in step 2 are within-sample (cell-level)
#' effects.
#'
#' @param regulatoryData Named list of matrices, sparse matrices or data.frames
#'   (regulators x cells), one per regulatory layer.
#' @param metadata Data.frame with cell metadata. Row names must be cell IDs.
#' @param randomEffects Character vector of \code{metadata} columns used
#'  as random intercepts.
#' @param omicType Named numeric vector with the type of each layer (0 =
#'   continuous, 1 = binary).
#' @param REML Logical. Restricted maximum likelihood (Gaussian models only).
#' @param verbose Logical. If \code{TRUE}, prints optimization progress.
#'
#' @return A named list of dense numeric matrices of residuals with the same
#'   dimensions and dimnames as the corresponding element of
#'   \code{regulatoryData}.
#'
#' @importFrom cli cli_alert_info cli_alert_warning
#' @importFrom pbapply pblapply
#' @keywords internal
prepareTwoSteps <- function(regulatoryData,
                            metadata,
                            randomEffects,
                            omicType,
                            REML    = FALSE,
                            verbose = FALSE) {

  residuals_list <- lapply(names(regulatoryData), function(layer) {
    mat <- regulatoryData[[layer]]
    n_regs <- nrow(mat)
    cli::cli_alert_info("Layer {.field {layer}}: residualizing {n_regs} regulator(s)...")

    # Create access to one regulator across cells
    if (inherits(mat, "sparseMatrix")) {
      tmat   <- Matrix::t(mat)
      getReg <- function(i) as.numeric(tmat[, i])
    } else {
      getReg <- function(i) as.numeric(mat[i, ])
    }

    res <- pbapply::pblapply(seq_len(n_regs), function(i) {
      fitTwoStep(
        targetExpr     = getReg(i),
        metadata       = metadata,
        randomEffects  = randomEffects,
        targetOmicType = omicType[[layer]],
        REML           = REML,
        verbose        = verbose
      )
    })

    n_fallback <- sum(vapply(res, function(r) isTRUE(attr(r, "fallback")), logical(1)))
    if (n_fallback > 0) {
      cli::cli_alert_warning(
        "Layer {.field {layer}}: the null model failed for {n_fallback} regulator(s); mean-centred values were used instead."
      )
    }

    resid <- matrix(unlist(lapply(res, as.numeric), use.names = FALSE),
                    nrow = n_regs, byrow = TRUE)
    dimnames(resid) <- dimnames(mat)
    resid
  })

  names(residuals_list) <- names(regulatoryData)
  residuals_list
}
