#' Fit an individual LMM or GLMM model
#'
#' Low-level internal function that fits a Gaussian linear mixed model (LMM)
#' via \code{lme4::lmer} or a binomial generalized linear mixed model (GLMM)
#' via \code{lme4::glmer}. When the formula has no random terms it falls back to
#' \code{stats::lm} or \code{stats::glm}. It captures warnings and errors,
#' assesses convergence and checks for singularity, so a failing model never
#' interrupts a loop over many features.
#'
#' @param formula Object of class \code{formula} specifying the model.
#' @param data A \code{data.frame} containing the response and the covariates.
#' @param targetOmicType Numeric flag: \code{0} = Gaussian (\code{lmer}),
#'   \code{1} = binary (\code{glmer}).
#' @param REML Logical. Restricted maximum likelihood (Gaussian models only).
#' @param verbose Logical. If \code{TRUE}, prints optimization progress.
#'
#' @return A named list with:
#' \itemize{
#'   \item \code{fit}: The fitted model (\code{merMod}, \code{lm} or \code{glm}),
#'     or \code{NULL} if fitting failed.
#'   \item \code{converged}: Logical. \code{TRUE} when the model was fitted
#'     without errors and the optimizer reported no convergence problems.
#'   \item \code{is_singular}: Logical. \code{TRUE} if the mixed model fit is
#'     singular (variance components close to zero); \code{NA} if not fitted.
#'   \item \code{warnings}: Character string with the collapsed warnings, or
#'     \code{NULL}.
#'   \item \code{error}: Character string with the error message if fitting
#'     failed, or \code{NULL}.
#' }
#'
#' @importFrom lme4 lmer glmer glmerControl isSingular findbars
#' @importFrom stats lm glm binomial
#' @keywords internal
runGLMM <- function(formula,
                    data,
                    targetOmicType,
                    REML,
                    verbose){
  res <- list(
    fit         = NULL,
    converged   = FALSE,
    is_singular = NA,
    warnings    = NULL,
    error       = NULL
  )

  warn_list <- character()
  has_random <- !is.null(lme4::findbars(formula))

  # ~~~~~~~~~~ 1. Model fitting & error/warning handling ~~~~~~~~~~ #
  fit <- tryCatch({
    withCallingHandlers({
      if (targetOmicType == 0) {
        # ~~~~~~~~~~ 1.0: LMM Gaussian ~~~~~~~~~~ #
        # if (has_random) {
        lme4::lmer(formula = formula, data = data, REML = REML, verbose = verbose)
          # } else {
          #   stats::lm(formula = formula, data = data)
          #   }
      } else if (targetOmicType == 1) {
        # ~~~~~~~~~~ 1.1: GLMM Binomial Logit ~~~~~~~~~~ #
        # if (has_random) {
        lme4::glmer(
          formula = formula,
          data    = data,
          family  = stats::binomial(link = "logit"),
          control = lme4::glmerControl(optimizer = "bobyqa",
                                       optCtrl = list(maxfun = 2e5)),
          verbose = verbose
        )
        # } else {
        #   stats::glm(formula = formula, data = data, family = stats::binomial(link = "logit"))
        # }
      } else {
        stop("Unsupported `targetOmicType`. Use 0 for Gaussian or 1 for Binomial.")
      }
    }, warning = function(w) {
      warn_list <<- c(warn_list, conditionMessage(w))
      invokeRestart("muffleWarning")
    }, message = function(m) {
      # lme4 reports singular fits as messages; they are assessed below
      if (!verbose) invokeRestart("muffleMessage")
    })
  }, error = function(e) {
    res$error <<- conditionMessage(e)
    return(NULL)
  })

  # Return early if model fitting crashed
  if (is.null(fit)) return(res)

  # ~~~~~~~~~~ 2. Convergence and Singularity Assessment ~~~~~~~~~~ #
  res$warnings <- if (length(warn_list) > 0) paste(warn_list, collapse = " | ") else NULL

  if (inherits(fit, "merMod")) {
    conv_code <- fit@optinfo$conv$opt
    conv_msgs <- fit@optinfo$conv$lme4$messages
    # A boundary (singular) fit is not a convergence failure
    conv_msgs <- conv_msgs[!grepl("singular", conv_msgs, ignore.case = TRUE)]
    res$converged   <- (is.null(conv_code) || conv_code == 0) && length(conv_msgs) == 0
    res$is_singular <- lme4::isSingular(fit, tol = 1e-4)
  } else if (inherits(fit, "glm")) {
    res$converged   <- isTRUE(fit$converged)
    res$is_singular <- FALSE
  } else {
  res$converged   <- TRUE
  res$is_singular <- TRUE
  }

  res$fit <- fit

  return(res)

}
