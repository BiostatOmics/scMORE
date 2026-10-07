#' Compare regulator effects between conditions, cell types or other groups
#'
#' Uses the regulator slopes stored by \code{\link{fitMore}} to test, for every
#' target and regulator, whether the effect of the regulator differs between the
#' levels of one or several fixed effects (for example treatment versus control,
#' or one cell type versus another). Differences can be computed within the
#' levels of other fixed effects, and interaction contrasts (difference of
#' differences) are also available.
#'
#' @details
#' The slope of a regulator is the change in the expression of the target per
#' unit of the regulator. A contrast such as \code{"treat - ctrl"} estimates
#' \eqn{\beta_{treat} - \beta_{ctrl}}{beta_treat - beta_ctrl}, so a positive value
#' means a larger slope in \code{treat}. Note that if the slopes are negative, a
#' positive difference means a weaker repression, not a stronger effect. A
#' regulator can be significant in one group and not in another without the
#' difference between them being significant: this function tests the
#' difference itself, whereas \code{\link{testRegulation}} tests each group
#' against zero.
#'
#' The arguments define which comparison is made:
#' \itemize{
#'   \item \code{compare} gives the fixed effect(s) whose levels are compared.
#'   \item \code{by} gives the fixed effects within whose levels the comparison
#'     is made. By default, all the remaining fixed effects (nothing is
#'     averaged). For example, \code{compare = "condition"} gives the difference
#'     between conditions separately in every cell type. Fixed effects that are
#'     in neither \code{compare} nor \code{by} are averaged over (see
#'     \code{weights}); averaging can hide opposite effects.
#'   \item \code{method} gives the contrasts:
#'     \code{"pairwise"} (all pairs, \code{A - B} with levels in their order),
#'     \code{"revpairwise"} (\code{B - A}),
#'     \code{"trt.vs.ctrl"} (each level against a reference level, set with
#'     \code{ref}), \code{"trt.vs.ctrlk"}, \code{"trt.vs.ctrl1"},
#'     \code{"eff"} (each level against the mean of all) and
#'     \code{"consec"} (consecutive levels). A named list of coefficient vectors
#'     can also be given for custom contrasts, one coefficient per level of
#'     \code{compare}, for example
#'     \code{list("treat - ctrl" = c(-1, 1))}.
#'   \item \code{interaction = TRUE} computes interaction contrasts, that is,
#'     differences of differences such as
#'     \eqn{(\beta_{treat,B} - \beta_{ctrl,B}) - (\beta_{treat,T} - \beta_{ctrl,T})}{(b_treat,B - b_ctrl,B) - (b_treat,T - b_ctrl,T)}.
#'     It needs two or more fixed effects in \code{compare}, and \code{method} is
#'     then applied to each of them.
#'   \item \code{at} restricts the analysis to specific levels before comparing,
#'     for example to compare only two of three cell types.
#' }
#'
#' Each contrast is tested with a Wald test of no difference, using the
#' asymptotic (normal) approximation, which can be optimistic when the number of
#' samples is small. The model must have been fitted with \code{interaction = TRUE},
#' because otherwise the slopes do not depend on the fixed effects.
#'
#' \strong{Multiple testing.} As in \code{\link{testRegulation}}, p-values are
#' adjusted only at the level of the test, with \code{emmeans}: the family is the
#' set of contrasts of one target-regulator pair within each level of \code{by}.
#' When the family has a single contrast (for example two conditions compared
#' within each cell type) the adjusted p-value equals the unadjusted one. No
#' correction is applied across targets, regulators or layers; the unadjusted
#' p-value is returned in \code{p.value} so a correction over the whole table can
#' be applied afterwards. Confidence intervals are never adjusted.
#'
#' @param fit An object of class \code{scMoreFit} returned by \code{\link{fitMore}}.
#' @param compare Character vector with the fixed effect(s) whose levels are
#'   compared (for example \code{"condition"}). Required.
#' @param by Character vector with the fixed effects within whose levels the
#'   comparison is made. If \code{NULL} (default) all the other fixed effects
#'   are used. Use \code{character(0)} to average over all of them.
#' @param method Contrast method: one of \code{"pairwise"} (default),
#'   \code{"revpairwise"}, \code{"trt.vs.ctrl"}, \code{"trt.vs.ctrl1"},
#'   \code{"trt.vs.ctrlk"}, \code{"eff"}, \code{"del.eff"} or \code{"consec"}, or a
#'   named list of coefficient vectors. See Details.
#' @param interaction Logical. If \code{TRUE}, interaction contrasts (differences
#'   of differences) are computed between the factors in \code{compare}.
#' @param ref Reference level for \code{method = "trt.vs.ctrl"} (and its
#'   variants), given as a level name or a position. Only valid with a single
#'   fixed effect in \code{compare}.
#' @param at Optional named list restricting the analysis to some levels, for
#'   example \code{list(cell_type = c("Tcell", "Bcell"))}. The order of the
#'   levels in \code{at} does not change the direction of the contrasts, which
#'   follows the order of the levels of the factor.
#' @param weights Either \code{"equal"} (default) or \code{"proportional"}: how
#'   the fixed effects that are in neither \code{compare} nor \code{by} are
#'   averaged (same weight per level, or weight by number of cells).
#' @param alpha Numeric. Significance threshold applied to the adjusted p-value
#'   (\code{p.adj}).
#' @param adjust Method used by \code{emmeans} to adjust the p-values of each
#'   family of contrasts. One of \code{"fdr"} (default, the same as \code{"BH"}),
#'   \code{"none"}, \code{"holm"}, \code{"bonferroni"}, \code{"BH"},
#'   \code{"BY"}, \code{"hochberg"}, \code{"hommel"}, \code{"sidak"},
#'   \code{"mvt"}, \code{"tukey"} (for pairwise contrasts), \code{"dunnettx"}
#'   (for contrasts against a reference) or \code{"scheffe"}. If a method is not
#'   valid for the chosen contrasts, \code{emmeans} replaces it by a valid one.
#' @param minEffect Numeric. Minimum absolute difference of slopes for a contrast
#'   to be called significant. Defaults to 0 (no threshold).
#' @param targets,regulators,layers Optional character vectors to restrict the
#'   analysis to some targets, regulators or regulatory layers. Targets and
#'   regulators are given with the sanitized names used by the models (see
#'   \code{GlobalSummary$nameMap}).
#' @param onlySignificant Logical. If \code{TRUE}, only significant rows are
#'   returned.
#' @param useOriginalNames Logical. If \code{TRUE}, sanitized target and
#'   regulator names are translated back to the original identifiers.
#'
#' @return A data.frame, sorted by adjusted p-value, with one row per target,
#'   regulator, contrast and level of \code{by}, and the columns \code{target},
#'   \code{regulator}, \code{layer}, \code{contrast}, one column per variable in
#'   \code{by}, \code{group}, \code{estimate} (difference of slopes),
#'   \code{std.error}, \code{conf.low}, \code{conf.high}, \code{statistic},
#'   \code{p.value} (unadjusted), \code{p.adj} (adjusted within the family),
#'   \code{direction} (\code{"positive"} or \code{"negative"}, the sign of the
#'   difference) and \code{significant}. Pairs whose slopes could not be
#'   compared are listed in the attribute \code{failed}.
#'
#' @examples
#' \dontrun{
#' # Treatment versus control in every cell type
#' res <- testContrasts(fit, compare = "condition")
#'
#' # Same, written as treat - ctrl, keeping only significant regulators
#' testContrasts(fit, compare = "condition", method = "revpairwise",
#'               onlySignificant = TRUE)
#'
#' # Cell types against a reference cell type, within each condition
#' testContrasts(fit, compare = "cell_type", method = "trt.vs.ctrl", ref = "Bcell")
#'
#' # Does the effect of the treatment depend on the cell type?
#' testContrasts(fit, compare = c("condition", "cell_type"), interaction = TRUE)
#'
#' # Condition difference averaged over cell types
#' testContrasts(fit, compare = "condition", by = character(0))
#' }
#'
#' @importFrom emmeans as.emmGrid emmeans contrast test
#' @importFrom stats as.formula
#' @export
testContrasts <- function(fit,
                          compare,
                          by               = NULL,
                          method           = "pairwise",
                          interaction      = FALSE,
                          ref              = NULL,
                          at               = NULL,
                          weights          = c("equal", "proportional"),
                          alpha            = 0.05,
                          adjust           = c("fdr", "none", "holm", "bonferroni", "BH",
                                               "BY", "hochberg", "hommel", "sidak", "mvt",
                                               "tukey", "dunnettx", "scheffe"),
                          minEffect        = 0,
                          targets          = NULL,
                          regulators       = NULL,
                          layers           = NULL,
                          onlySignificant  = FALSE,
                          useOriginalNames = FALSE) {

  weights <- match.arg(weights)
  adjust  <- match.arg(adjust)

  # ~~~~~~~~~~ 1. Validate arguments ~~~~~~~~~~ #
  if (!inherits(fit, "scMoreFit")) {
    cli::cli_abort("{.arg fit} must be an object returned by {.fn fitMore}.")
  }
  fixedEffects <- fit$GlobalSummary$args$fixedEffects
  if (length(fixedEffects) == 0) {
    cli::cli_abort("The model has no fixed effects, so there is nothing to compare.")
  }
  if (isFALSE(fit$GlobalSummary$args$interaction)) {
    cli::cli_abort(c(
      "Contrasts between groups need a model fitted with {.code interaction = TRUE}.",
      "i" = "Without interactions the slopes do not depend on the fixed effects, so all differences are zero."
    ))
  }

  if (missing(compare) || is.null(compare) || length(compare) == 0) {
    cli::cli_abort(c("{.arg compare} is required.", "i" = "Fixed effects of the model: {.val {fixedEffects}}."))
  }
  bad_cmp <- setdiff(compare, fixedEffects)
  if (length(bad_cmp) > 0) {
    cli::cli_abort(c("{.arg compare} must be fixed effects of the model.",
                     "x" = "Not found: {.val {bad_cmp}}.", "i" = "Available: {.val {fixedEffects}}."))
  }

  # By default the comparison is made within every combination of the other fixed effects
  if (is.null(by)) by <- setdiff(fixedEffects, compare)
  bad_by <- setdiff(by, fixedEffects)
  if (length(bad_by) > 0) {
    cli::cli_abort(c("{.arg by} must be fixed effects of the model.",
                     "x" = "Not found: {.val {bad_by}}.", "i" = "Available: {.val {fixedEffects}}."))
  }
  if (length(intersect(by, compare)) > 0) {
    cli::cli_abort("A fixed effect cannot be in both {.arg compare} and {.arg by}.")
  }

  # Contrast method
  std_methods <- c("pairwise", "revpairwise", "trt.vs.ctrl", "trt.vs.ctrl1",
                   "trt.vs.ctrlk", "eff", "del.eff", "consec")
  if (is.character(method)) {
    if (length(method) != 1 || !method %in% std_methods) {
      cli::cli_abort("{.arg method} must be one of {.val {std_methods}} or a named list of coefficients.")
    }
  } else if (!(is.list(method) && !is.null(names(method)) && all(nzchar(names(method))))) {
    cli::cli_abort("{.arg method} must be one of {.val {std_methods}} or a named list of coefficients.")
  }

  if (isTRUE(interaction)) {
    if (length(compare) < 2) {
      cli::cli_abort("{.code interaction = TRUE} needs two or more fixed effects in {.arg compare}.")
    }
    if (!is.character(method)) {
      cli::cli_abort("{.code interaction = TRUE} needs {.arg method} to be a method name, not a list.")
    }
  }

  if (!is.null(ref)) {
    if (!is.character(method) || !startsWith(method, "trt.vs.ctrl")) {
      cli::cli_abort("{.arg ref} can only be used with {.code method = \"trt.vs.ctrl\"} (or its variants).")
    }
    if (length(compare) != 1 || isTRUE(interaction)) {
      cli::cli_abort("{.arg ref} needs a single fixed effect in {.arg compare} and {.code interaction = FALSE}.")
    }
  }

  # ~~~~~~~~~~ 2. Select models and validate against the stored grid ~~~~~~~~~~ #
  models  <- selectModels(fit, targets)
  layerOf <- layerLookup(fit)
  checkLayers(layers, fit)

  grid0 <- firstTrendGrid(models)
  checkAt(at, fixedEffects, grid0)

  for (v in compare) {
    if (is.numeric(grid0[[v]])) {
      cli::cli_abort("{.field {v}} is a continuous fixed effect and cannot be compared by levels.")
    }
    if (nLevelsAt(v, at, grid0) < 2) {
      cli::cli_abort("{.field {v}} needs at least two levels to be compared (check {.arg at}).")
    }
  }

  # Validate the reference level against the levels that will be compared
  if (!is.null(ref)) {
    lv <- unique(as.character(grid0[[compare]]))
    if (!is.null(at[[compare]])) lv <- intersect(lv, as.character(at[[compare]]))
    ok_ref <- if (is.character(ref)) ref %in% lv else (is.numeric(ref) && all(ref >= 1 & ref <= length(lv)))
    if (!ok_ref) {
      cli::cli_abort(c(
        "{.arg ref} is not one of the levels being compared.",
        "i" = "Levels of {.field {compare}}: {.val {lv}}."
      ))
    }
  }

  # Fixed effects in neither `compare` nor `by` are averaged over
  averaged <- Filter(function(v) nLevelsAt(v, at, grid0) > 1, setdiff(fixedEffects, c(compare, by)))
  if (length(averaged) > 0) {
    cli::cli_alert_info(c(
      "Slopes are averaged over {.field {averaged}} ({weights} weights). ",
      "An averaged difference can hide opposite effects in different levels; ",
      "add {.field {averaged}} to {.arg by} to compare within each level."
    ))
  }

  # emmeans needs the factors being compared and the stratifying factors in the grid
  specs <- stats::as.formula(paste("~", paste(c(compare, by), collapse = " + ")))

  # ~~~~~~~~~~ 3. Compare every target-regulator pair ~~~~~~~~~~ #
  testOne <- function(trend, reg) {
    emm <- emmeans::as.emmGrid(trend)
    if (!is.null(at)) {
      emm <- subsetTrendGrid(emm, at)
      if (is.null(emm)) stop("none of the requested levels in `at` exist in the model")
    }

    # Average over the fixed effects that are not in `compare` or `by`.
    # emmeans prints a NOTE about interactions on every call: it is reported once above
    sub <- suppressMessages(emmeans::emmeans(emm, specs = specs, weights = weights))

    # Build the contrasts. With `by`, they are computed within each level of `by`
    ct_args <- list(object = sub)
    if (isTRUE(interaction)) {
      ct_args$interaction <- rep(method, length(compare))     # one method per compared factor
    } else {
      ct_args$method <- method
    }
    if (length(by) > 0) ct_args$by <- by
    if (!is.null(ref)) {
      # A level name is converted into its position among the compared levels
      ct_args$ref <- if (is.character(ref)) match(ref, sub@levels[[compare]]) else ref
      if (anyNA(ct_args$ref)) stop("reference level not found among the compared levels")
    }
    ct <- suppressMessages(do.call(emmeans::contrast, ct_args))

    # Unadjusted differences, confidence intervals and Wald tests against 0
    s <- as.data.frame(summary(ct, infer = c(TRUE, TRUE), null = 0, adjust = "none"))

    # Adjusted p-value: family = contrasts of this pair within each level of `by`
    # (same as test(adjust = ) in emmeans). Confidence intervals are left unadjusted
    s$p.adj <- if (adjust == "none") {
      s$p.value
    } else {
      as.data.frame(emmeans::test(ct, null = 0, adjust = adjust))$p.value
    }

    # One tidy `contrast` column. Interaction contrasts have one column per factor
    # (for example condition_pairwise and cell_type_pairwise)
    est_i <- match("estimate", names(s))
    cols  <- setdiff(names(s)[seq_len(est_i - 1)], by)
    lab   <- lapply(s[cols], as.character)
    s$contrast <- if (length(cols) > 1) {
      do.call(paste, c(lapply(lab, function(x) paste0("(", x, ")")), sep = " x "))
    } else {
      lab[[1]]
    }
    # Drop the original contrast columns (the tidy `contrast` column is kept)
    s[, setdiff(names(s), setdiff(cols, "contrast")), drop = FALSE]
  }

  # Compare every pair. A failure in one pair never stops the loop: it is collected
  # and reported at the end through attr(, "failed") (see testPairs())
  pairs_out <- testPairs(models, regulators, layers, layerOf, testOne, verb = "compared")
  res       <- pairs_out$res
  failed    <- pairs_out$failed

  # ~~~~~~~~~~ 4. Tidy columns ~~~~~~~~~~ #
  res <- tidyEmmeansColumns(res)

  res$group <- if (length(by) > 0) do.call(paste, c(lapply(res[by], as.character), sep = " / ")) else "overall"

  # ~~~~~~~~~~ 5. Multiple testing correction ~~~~~~~~~~ #
  # `p.adj` was computed in testOne() within each family of contrasts. No correction
  # across targets, regulators or layers is applied for now (see testRegulation()
  # for the code left commented to enable a global correction in the future).

  res$direction   <- ifelse(res$estimate >= 0, "positive", "negative")
  res$significant <- !is.na(res$p.adj) & res$p.adj < alpha & abs(res$estimate) >= minEffect

  # ~~~~~~~~~~ 6. Order columns, filter and names ~~~~~~~~~~ #
  first      <- c("target", "regulator", "layer", "contrast", by, "group")
  stats_cols <- c("estimate", "std.error", "df", "conf.low", "conf.high",
                  "statistic", "p.value", "p.adj", "direction", "significant")
  res <- res[, c(first, intersect(stats_cols, names(res))), drop = FALSE]

  if (useOriginalNames) res <- restoreOriginalNames(res, fit)
  if (onlySignificant)  res <- res[res$significant, , drop = FALSE]
  res <- res[order(res$p.adj, na.last = TRUE), , drop = FALSE]
  rownames(res) <- NULL

  n_sig   <- sum(res$significant)
  adj_txt <- if (adjust == "none") "no p-value adjustment" else paste0(adjust, " adjustment within each family of contrasts")
  cli::cli_alert_success(
    "Tested {nrow(res)} contrast(s); {n_sig} significant (p < {alpha}, {adj_txt})."
  )
  if (!is.null(failed) && nrow(failed) > 0) {
    cli::cli_alert_warning("{nrow(failed)} regulator-target pair(s) could not be compared (see {.code attr(, \"failed\")}).")
  }

  attr(res, "failed") <- failed
  attr(res, "args")   <- list(compare = compare, by = by, method = method,
                              interaction = interaction, ref = ref, at = at,
                              weights = weights, alpha = alpha, adjust = adjust,
                              minEffect = minEffect)
  res
}
