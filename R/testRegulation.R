#' Test regulator activity per condition, cell type or any fixed effect
#'
#' Uses the regulator slopes stored by \code{\link{fitMore}} to test, for every
#' target and regulator, whether the regulator has an effect (slope different
#' from zero) in each group defined by one or several fixed effects of the model.
#' P-values can be adjusted within each target-regulator pair (see Details) and
#' the result is returned as a tidy table.
#'
#' @details
#' \code{fitMore} stores the slope of every regulator for each combination of all
#' the fixed effects (for example condition x cell type). This function works on
#' those stored slopes with \pkg{emmeans}:
#' \itemize{
#'   \item \code{by} selects the fixed effects that define the groups to report.
#'     Slopes are averaged over the fixed effects that are not in \code{by}. For
#'     example, with \code{by = "condition"} the slope of a condition is the
#'     average of its slopes across cell types. Averaging can hide opposite
#'     effects in different cell types, so use all fixed effects (the default)
#'     when the interaction is of interest.
#'   \item \code{at} restricts the analysis to specific levels before averaging.
#'     For example, \code{by = "condition"} with
#'     \code{at = list(cell_type = "Tcell")} reports the slope of each condition
#'     within T cells.
#'   \item \code{weights} sets how the omitted fixed effects are averaged:
#'     \code{"equal"} gives the same weight to each level, \code{"proportional"}
#'     weights each level by its number of cells.
#' }
#' Each test is a Wald test of \eqn{H_0: \text{slope} = 0}{H0: slope = 0} using
#' the asymptotic (normal) approximation, which can be optimistic when the number
#' of samples is small.
#'
#' \strong{Multiple testing.} P-values are adjusted only at the level of the
#' test, as in \code{emmeans::test(..., adjust = )}: the family of tests is the
#' set of groups returned for one target-regulator pair (for example the
#' condition x cell type combinations of one regulator of one gene). No
#' correction is applied across targets, regulators or layers, so a list of
#' significant regulators obtained with this function does not control the false
#' discovery rate of the whole analysis. The unadjusted p-value is always
#' returned in \code{p.value}, so a correction over the whole table can be applied
#' afterwards, for example \code{res$q <- p.adjust(res$p.value, "BH")}. Confidence
#' intervals are never adjusted. A global correction (over all tests, by group or
#' by layer) is not implemented yet: the code is left commented in the source.
#'
#' If the model was fitted with \code{interaction = FALSE}, slopes do not depend
#' on the fixed effects, so every group returns the same result.
#'
#' @param fit An object of class \code{scMoreFit} returned by \code{\link{fitMore}}.
#' @param by Character vector with the fixed effects defining the groups (for
#'   example \code{"condition"} or \code{c("condition", "cell_type")}). If
#'   \code{NULL} (default) all the fixed effects of the model are used. Use
#'   \code{character(0)} for a single overall estimate.
#' @param at Optional named list restricting the analysis to some levels, for
#'   example \code{list(cell_type = c("Tcell", "Bcell"))}. Names must be fixed
#'   effects of the model.
#' @param weights Either \code{"equal"} (default) or \code{"proportional"}. See
#'   Details.
#' @param alpha Numeric. Significance threshold applied to the adjusted p-value
#'   (\code{p.adj}).
#' @param adjust Method used by \code{emmeans} to adjust the p-values of the
#'   groups of each target-regulator pair. One of \code{"fdr"} (default, the same
#'   as \code{"BH"}), \code{"none"}, \code{"holm"}, \code{"bonferroni"},
#'   \code{"BH"}, \code{"BY"}, \code{"hochberg"}, \code{"hommel"},
#'   \code{"sidak"} or \code{"mvt"}. See Details.
#' @param minEffect Numeric. Minimum absolute slope for a regulator to be called
#'   significant. Defaults to 0 (no threshold).
#' @param targets,regulators,layers Optional character vectors to restrict the
#'   analysis to some targets, regulators or regulatory layers. Targets and
#'   regulators are given with the sanitized names used by the models (see
#'   \code{GlobalSummary$nameMap}).
#' @param onlySignificant Logical. If \code{TRUE}, only significant rows are
#'   returned.
#' @param useOriginalNames Logical. If \code{TRUE}, sanitized target and
#'   regulator names are translated back to the original identifiers using
#'   \code{GlobalSummary$nameMap}.
#'
#' @return A data.frame, sorted by adjusted p-value, with one row per target,
#'   regulator and group and the columns \code{target}, \code{regulator},
#'   \code{layer}, one column per variable in \code{by}, \code{group},
#'   \code{estimate} (slope), \code{std.error}, \code{conf.low},
#'   \code{conf.high}, \code{statistic}, \code{p.value} (unadjusted),
#'   \code{p.adj} (adjusted within the target-regulator pair),
#'   \code{direction} (\code{"activator"} or \code{"repressor"}, according to the
#'   sign of the slope) and \code{significant}. Regulator-target pairs whose
#'   slopes could not be estimated are listed in the attribute \code{failed}.
#'
#' @examples
#' \dontrun{
#' # Regulators active in each condition (averaged over cell types)
#' res <- testRegulation(fit, by = "condition")
#'
#' # Each condition within T cells only, keeping significant regulators
#' testRegulation(fit, by = "condition", at = list(cell_type = "Tcell"),
#'                onlySignificant = TRUE)
#'
#' # Holm instead of FDR within each regulator-target pair
#' testRegulation(fit, adjust = "holm")
#'
#' # Correction over the whole table, applied afterwards on the raw p-values
#' res <- testRegulation(fit, adjust = "none")
#' res$q <- p.adjust(res$p.value, method = "BH")
#' }
#'
#' @importFrom emmeans as.emmGrid emmeans
#' @importFrom stats as.formula
#' @export
testRegulation <- function(fit,
                           by               = NULL,
                           at               = NULL,
                           weights          = c("equal", "proportional"),
                           alpha            = 0.05,
                           adjust           = c("fdr", "none", "holm", "bonferroni", "BH",
                                                "BY", "hochberg", "hommel", "sidak", "mvt"),
                           # FUTURE: global correction (see section 5 below)
                           # adjustBy       = c("global", "group", "layer", "pair"),
                           minEffect        = 0,
                           targets          = NULL,
                           regulators       = NULL,
                           layers           = NULL,
                           onlySignificant  = FALSE,
                           useOriginalNames = FALSE) {

  weights <- match.arg(weights)
  adjust  <- match.arg(adjust)
  # adjustBy <- match.arg(adjustBy)   # FUTURE: see section 5

  # ~~~~~~~~~~ 1. Validate arguments ~~~~~~~~~~ #
  if (!inherits(fit, "scMoreFit")) {
    cli::cli_abort("{.arg fit} must be an object returned by {.fn fitMore}.")
  }
  fixedEffects <- fit$GlobalSummary$args$fixedEffects

  if (is.null(by)) by <- fixedEffects
  bad_by <- setdiff(by, fixedEffects)
  if (length(bad_by) > 0) {
    cli::cli_abort(c(
      "{.arg by} must be fixed effects of the model.",
      "x" = "Not found: {.val {bad_by}}.",
      "i" = "Available: {.val {fixedEffects}}."
    ))
  }
  if (isFALSE(fit$GlobalSummary$args$interaction)) {
    cli::cli_alert_warning("The model was fitted with {.code interaction = FALSE}: slopes do not depend on the fixed effects, so all groups give the same result.")
  }

  # ~~~~~~~~~~ 2. Select models and validate against the stored grid ~~~~~~~~~~ #
  models  <- selectModels(fit, targets)
  layerOf <- layerLookup(fit)
  checkLayers(layers, fit)

  # The stored grid of the first model is used to validate `at` and to know
  # which fixed effects will be averaged over
  grid0 <- firstTrendGrid(models)
  checkAt(at, fixedEffects, grid0)

  # Fixed effects not in `by` that still have several levels are averaged over
  averaged <- Filter(function(v) nLevelsAt(v, at, grid0) > 1, setdiff(fixedEffects, by))
  if (length(averaged) > 0 && !isFALSE(fit$GlobalSummary$args$interaction)) {
    cli::cli_alert_info(c(
      "Slopes are averaged over {.field {averaged}} ({weights} weights). ",
      "As the model includes interactions, an averaged slope can hide opposite effects in different levels; ",
      "add {.field {averaged}} to {.arg by} to see each level."
    ))
  }

  specs <- stats::as.formula(paste("~", if (length(by) > 0) paste(by, collapse = " + ") else "1"))

  # ~~~~~~~~~~ 3. Test every target-regulator pair ~~~~~~~~~~ #
  testOne <- function(trend, reg) {
    emm <- emmeans::as.emmGrid(trend)
    if (!is.null(at)) {
      emm <- subsetTrendGrid(emm, at)
      if (is.null(emm)) stop("none of the requested levels in `at` exist in the model")
    }
    # Average over the fixed effects that are not in `by` (see `weights`).
    # emmeans prints a NOTE about interactions on every call: it is reported once above
    sub <- suppressMessages(emmeans::emmeans(emm, specs = specs, weights = weights))

    # Slopes, standard errors, confidence intervals and Wald tests against 0.
    # Here everything is unadjusted: `p.value` is the raw p-value
    s <- as.data.frame(summary(sub, infer = c(TRUE, TRUE), null = 0, adjust = "none"))
    names(s)[names(s) == paste0(reg, ".trend")] <- "estimate"

    # Adjusted p-value. The family of tests is every row of `sub`, that is, the
    # groups of this single target-regulator pair (same as test(adjust = ) in emmeans).
    # Confidence intervals are left unadjusted on purpose
    s$p.adj <- if (adjust == "none") {
      s$p.value
    } else {
      as.data.frame(emmeans::test(sub, null = 0, adjust = adjust))$p.value
    }
    if (length(by) == 0) s <- s[, setdiff(names(s), "1"), drop = FALSE]
    s
  }

  # Test every pair. A failure in one pair never stops the loop: it is collected
  # and reported at the end through attr(, "failed") (see testPairs())
  pairs_out <- testPairs(models, regulators, layers, layerOf, testOne, verb = "tested")
  res       <- pairs_out$res
  failed    <- pairs_out$failed

  # ~~~~~~~~~~ 4. Tidy columns ~~~~~~~~~~ #
  res <- tidyEmmeansColumns(res)

  res$group <- if (length(by) > 0) do.call(paste, c(lapply(res[by], as.character), sep = " / ")) else "overall"

  # ~~~~~~~~~~ 5. Multiple testing correction ~~~~~~~~~~ #
  # `p.adj` was already computed in testOne(): emmeans adjusts the groups of each
  # target-regulator pair. No correction across targets, regulators or layers is
  # applied for now (the user can use `p.value` for that, see the documentation).
  #
  # FUTURE: global correction across the whole table. To enable it, restore the
  # `adjustBy` argument (signature, match.arg and documentation) and use:
  #
  # family <- switch(adjustBy,
  #                  global = rep(1L, nrow(res)),                # all tests together
  #                  group  = res$group,                         # one family per group
  #                  layer  = res$layer,                         # one family per layer
  #                  pair   = paste(res$target, res$regulator))  # as the current behaviour
  # res$p.adj <- stats::ave(res$p.value, family,
  #                         FUN = function(p) stats::p.adjust(p, method = adjust))
  #
  # Notes: (1) p.adjust() only knows "holm", "hochberg", "hommel", "bonferroni",
  # "BH", "BY", "fdr" and "none", so "sidak" and "mvt" would have to be excluded;
  # (2) with adjustBy = "pair" the result equals the current one; (3) remember to
  # add "importFrom stats p.adjust ave" to the roxygen header.
  # Rationale and simulation: in a global BH the proportion of false discoveries
  # stays at the nominal level, whereas adjusting within each pair does not
  # control it across thousands of pairs.

  res$direction   <- ifelse(res$estimate >= 0, "activator", "repressor")
  res$significant <- !is.na(res$p.adj) & res$p.adj < alpha & abs(res$estimate) >= minEffect

  # ~~~~~~~~~~ 6. Order columns, filter and names ~~~~~~~~~~ #
  first <- c("target", "regulator", "layer", by, "group")
  stats_cols <- c("estimate", "std.error", "df", "conf.low", "conf.high",
                  "statistic", "p.value", "p.adj", "direction", "significant")
  res <- res[, c(first, intersect(stats_cols, names(res))), drop = FALSE]

  if (useOriginalNames) res <- restoreOriginalNames(res, fit)
  if (onlySignificant) res <- res[res$significant, , drop = FALSE]
  res <- res[order(res$p.adj, na.last = TRUE), , drop = FALSE]
  rownames(res) <- NULL

  n_sig   <- sum(res$significant)
  adj_txt <- if (adjust == "none") "no p-value adjustment" else paste0(adjust, " adjustment within each regulator-target pair")
  cli::cli_alert_success(
    "Tested {nrow(res)} regulator-target-group combination(s); {n_sig} significant (p < {alpha}, {adj_txt})."
  )
  if (!is.null(failed) && nrow(failed) > 0) {
    cli::cli_alert_warning("{nrow(failed)} regulator-target pair(s) could not be tested (see {.code attr(, \"failed\")}).")
  }

  attr(res, "failed") <- failed
  attr(res, "args")   <- list(by = by, at = at, weights = weights, alpha = alpha,
                              adjust = adjust, minEffect = minEffect)
  res
}
