#' Construct the formula and fit the mixed model of one target feature
#'
#' Prepares the cell-level data of a target feature, builds the model formula
#' with its prior regulators, the fixed effects and (optionally) the
#' regulator-by-design interactions, adds random intercepts and fits the model
#' with \code{\link{runGLMM}}. Optionally, the fitted model is replaced by a
#' lightweight list of \code{emmeans} trend objects (one per regulator).
#'
#' @details
#' With \code{interaction = TRUE} the fixed part of the formula is
#' \code{(reg_1 + ... + reg_k) * (fe_1 * ... * fe_m)}, that is, every regulator
#' is crossed with all fixed effects and their interactions. With
#' \code{interaction = FALSE} it is \code{reg_1 + ... + reg_k + fe_1 + ... + fe_m}.
#' The random part is \code{(1 | re_1) + ... + (1 | re_j)}.
#'
#' When \code{asEmmGrid = TRUE} the regulator slopes for every combination of
#' the fixed effects are estimated with \code{\link[emmeans]{emtrends}} (using
#' asymptotic degrees of freedom for mixed models, which is much faster than
#' Satterthwaite or Kenward-Roger) and stored as lists that can be turned back
#' into \code{emmGrid} objects with \code{emmeans::as.emmGrid}.
#'
#' @param target Character. Name of the target feature (row of \code{targetData}).
#' @param targetData Matrix, sparse \code{Matrix} or data.frame of target
#'   features (features x cells).
#' @param regData Numeric matrix of regulatory features (regulators x cells),
#'   typically the stacked regulatory layers.
#' @param metadata Data.frame with cell metadata. Row names must be cell IDs.
#' @param regulators Character vector of regulators of \code{target}
#'   (rownames of \code{regData}).
#' @param fixedEffects Character vector of \code{metadata} columns used as
#'  fixed tearm.
#' @param randomEffects Character vector of \code{metadata} columns used
#'  as random intercepts.
#' @param interaction Logical. If \code{TRUE}, regulators are crossed
#'   with all fixed effects and their interactions.
#' @param targetOmicType Numeric flag: \code{0} = Gaussian (\code{lmer}),
#'   \code{1} = binary (\code{glmer}).
#' @param REML Logical. Restricted maximum likelihood (Gaussian models only).
#' @param verbose Logical. If \code{TRUE}, prints optimization progress.
#' @param asEmmGrid Logical. If \code{TRUE}, the fitted model is converted into a
#'   named list (one element per regulator) of \code{emtrends} results.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{fit}: The fitted model, or a named list of \code{emtrends}
#'     lists per regulator if \code{asEmmGrid = TRUE} (an element is \code{NULL}
#'     if its estimation failed; the whole \code{fit} is \code{NULL} if the model
#'     failed).
#'   \item \code{converged}, \code{is_singular}, \code{warnings}, \code{error}:
#'     Diagnostics returned by \code{\link{runGLMM}}.
#'   \item \code{target}: Identifier of the target feature.
#'   \item \code{formula}: The model formula.
#'   \item \code{regulators}: Regulators included in the model.
#' }
#'
#' @importFrom stats as.formula
#' @importFrom emmeans emtrends
#' @keywords internal
fitTargetModel <- function(target,
                         targetData,
                         regData,
                         metadata,
                         regulators,
                         fixedEffects,
                         randomEffects,
                         interaction,
                         targetOmicType,
                         REML,
                         verbose,
                         asEmmGrid) {

  regulators <- intersect(unique(regulators), rownames(regData))
  if (length(regulators) == 0) {
    return(list(fit = NULL, converged = FALSE, is_singular = NA, warnings = NULL,
                error = "No regulators of this target are available in `regData`.",
                target = target, formula = NULL, regulators = character(0)))
  }

  # ~~~~~~~~~~ 1. Design matrix construction ~~~~~~~~~~ #
  # Extract regulatory features associated with the target
  design_vars <- unique(c(fixedEffects, randomEffects))
  target_expr <- targetData[target, ]
  reg_sub <- t(as.matrix(regData[regulators, , drop = FALSE]))

  model_data <- data.frame(
    TargetExpression = target_expr,
    reg_sub,
    metadata[, design_vars, drop = FALSE],
    check.names = FALSE
  )

  # ~~~~~~~~~~ 2. Dynamic formula generation ~~~~~~~~~~ #
  has_fixed <- length(fixedEffects) > 0

  if (interaction && has_fixed) {
    fe_terms <- paste0("(", paste(regulators, collapse = " + "), ") * (",
                       paste(fixedEffects, collapse = " * "), ")")
  } else {
    fe_terms <- paste(c(regulators, fixedEffects), collapse = " + ")
  }

  re_terms <- if (length(randomEffects) > 0) {
    paste0("(1 | ", randomEffects, ")", collapse = " + ")
  } else {
    NULL
  }

  formula_str   <- paste("TargetExpression ~", paste(c(fe_terms, re_terms), collapse = " + "))
  model_formula <- stats::as.formula(formula_str, env = globalenv())

  # ~~~~~~~~~~ 4. Model execution ~~~~~~~~~~ #
  fit_res <- runGLMM(
    formula        = model_formula,
    data           = model_data,
    targetOmicType = targetOmicType,
    REML           = REML,
    verbose        = verbose
  )

  # ~~~~~~~~~~ 5. Model to EmmGrid ~~~~~~~~~~ #
  if(asEmmGrid && !is.null(fit_res$fit)){
    fit   <- fit_res$fit
    is_me <- inherits(fit, "merMod")
    specs <- stats::as.formula(paste("~", if (has_fixed) paste(fixedEffects, collapse = " + ") else "1"))

    trends = sapply(regulators, function(reg){
      tryCatch({
        p <- if (is_me) {
          emmeans::emtrends(fit, specs = specs, var = reg, mode = "asymptotic")
        } else {
          emmeans::emtrends(fit, specs = specs, var = reg)
        }
        as.list(p, model.info.slot = TRUE)
      }, error = function(e) NULL)
    }, simplify = FALSE)

    names(trends) <- regulators
    fit_res$fit <- trends
  }

  fit_res$target     <- target
  fit_res$formula    <- model_formula
  fit_res$regulators <- regulators

  return(fit_res)
}
