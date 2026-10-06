#' Test regulator activity per condition, cell type or any fixed effect
#'
#' Uses the regulator slopes stored by \code{\link{fitMore}} to test, for every
#' target and regulator, whether the regulator has an effect (slope different
#' from zero) in each group defined by one or several fixed effects of the model.
#' P-values are adjusted for multiple testing and the result is returned as a
#' tidy table.
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
#' of samples is small. The family of tests used for the multiple testing
#' correction is set with \code{adjustBy}.
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
#' @param alpha Numeric. Significance threshold applied to the adjusted p-value.
#' @param adjust Method passed to \code{\link[stats]{p.adjust}}. Defaults to
#'   \code{"BH"}.
#' @param adjustBy Family of tests for the correction: \code{"global"} (all tests
#'   together), \code{"group"} (separately within each group) or \code{"layer"}
#'   (separately within each regulatory layer).
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
#'   \code{conf.high}, \code{statistic}, \code{p.value}, \code{p.adj},
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
#' # Every condition x cell type combination, FDR within each group
#' testRegulation(fit, adjustBy = "group")
#' }
#'
#' @importFrom emmeans as.emmGrid emmeans
#' @importFrom stats p.adjust ave as.formula
#' @importFrom pbapply pblapply
#' @export
testRegulation <- function(fit,
                           by               = NULL,
                           at               = NULL,
                           weights          = c("equal", "proportional"),
                           alpha            = 0.05,
                           adjust           = "BH",
                           adjustBy         = c("global", "group", "layer"),
                           minEffect        = 0,
                           targets          = NULL,
                           regulators       = NULL,
                           layers           = NULL,
                           onlySignificant  = FALSE,
                           useOriginalNames = FALSE) {

  weights  <- match.arg(weights)
  adjustBy <- match.arg(adjustBy)

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
  if (!is.null(at)) {
    if (!is.list(at) || is.null(names(at))) {
      cli::cli_abort("{.arg at} must be a named list, for example {.code list(cell_type = \"Tcell\")}.")
    }
    bad_at <- setdiff(names(at), fixedEffects)
    if (length(bad_at) > 0) {
      cli::cli_abort("{.arg at} contains variables that are not fixed effects of the model: {.val {bad_at}}.")
    }
  }
  if (isFALSE(fit$GlobalSummary$args$interaction)) {
    cli::cli_alert_warning("The model was fitted with {.code interaction = FALSE}: slopes do not depend on the fixed effects, so all groups give the same result.")
  }

  # ~~~~~~~~~~ 2. Select models ~~~~~~~~~~ #
  models <- fit$Models
  if (!is.null(targets)) models <- models[intersect(targets, names(models))]
  if (length(models) == 0) cli::cli_abort("No targets to analyse.")

  priors  <- fit$GlobalSummary$priorAssociations
  layerOf <- unlist(lapply(names(priors), function(l) {
    regs <- unique(as.character(priors[[l]]$regulator))
    stats::setNames(rep(l, length(regs)), regs)
  }))

  if (!is.null(layers)) {
    bad_layers <- setdiff(layers, names(priors))
    if (length(bad_layers) > 0) {
      cli::cli_abort("Unknown layer(s): {.val {bad_layers}}. Available: {.val {names(priors)}}.")
    }
  }

  # Validate the requested levels against the stored grid
  if (!is.null(at)) {
    first_trend <- NULL
    for (m in models) {
      ok <- Filter(Negate(is.null), m$fit)
      if (length(ok) > 0) { first_trend <- ok[[1]]; break }
    }
    if (!is.null(first_trend)) {
      grid0 <- emmeans::as.emmGrid(first_trend)@grid
      for (v in names(at)) {
        avail <- unique(as.character(grid0[[v]]))
        miss  <- setdiff(as.character(at[[v]]), avail)
        if (length(miss) > 0) {
          cli::cli_abort(c(
            "Level(s) {.val {miss}} not found in {.field {v}}.",
            "i" = "Available levels: {.val {avail}}."
          ))
        }
      }
    }
  }

  specs <- stats::as.formula(paste("~", if (length(by) > 0) paste(by, collapse = " + ") else "1"))

  # Restrict the stored grid to the requested levels
  subsetGrid <- function(emm, at) {
    keep <- rep(TRUE, nrow(emm@grid))
    for (v in names(at)) keep <- keep & (as.character(emm@grid[[v]]) %in% as.character(at[[v]]))
    if (!any(keep)) return(NULL)
    emm[which(keep)]
  }

  # ~~~~~~~~~~ 3. Test every target-regulator pair ~~~~~~~~~~ #
  failed <- list()

  testOne <- function(trend, reg) {
    emm <- emmeans::as.emmGrid(trend)
    if (!is.null(at)) {
      emm <- subsetGrid(emm, at)
      if (is.null(emm)) stop("none of the requested levels in `at` exist in the model")
    }
    sub <- emmeans::emmeans(emm, specs = specs, weights = weights)
    s   <- as.data.frame(summary(sub, infer = c(TRUE, TRUE), null = 0, adjust = "none"))
    names(s)[names(s) == paste0(reg, ".trend")] <- "estimate"
    if (length(by) == 0) s <- s[, setdiff(names(s), "1"), drop = FALSE]
    s
  }

  res_list <- pbapply::pblapply(names(models), function(g) {
    trends <- models[[g]]$fit
    if (is.null(trends)) return(NULL)

    regs <- names(trends)
    if (!is.null(regulators)) regs <- intersect(regs, regulators)
    if (!is.null(layers))     regs <- regs[layerOf[regs] %in% layers]

    out <- lapply(regs, function(reg) {
      if (is.null(trends[[reg]])) return(list(ok = NULL, fail = c(g, reg, "slopes could not be estimated")))
      tryCatch({
        s <- testOne(trends[[reg]], reg)
        list(ok = cbind(target = g, regulator = reg, layer = unname(layerOf[reg]), s,
                        stringsAsFactors = FALSE), fail = NULL)
      }, error = function(e) list(ok = NULL, fail = c(g, reg, conditionMessage(e))))
    })
    list(ok   = do.call(rbind, lapply(out, `[[`, "ok")),
         fail = do.call(rbind, lapply(out, `[[`, "fail")))
  })

  res <- do.call(rbind, lapply(res_list, `[[`, "ok"))
  failed <- do.call(rbind, lapply(res_list, `[[`, "fail"))
  if (!is.null(failed)) {
    failed <- as.data.frame(failed, stringsAsFactors = FALSE)
    names(failed) <- c("target", "regulator", "reason")
  }

  if (is.null(res) || nrow(res) == 0) {
    cli::cli_abort("No regulator-target pairs could be tested. Check {.arg at}, {.arg regulators} and {.arg layers}.")
  }

  # ~~~~~~~~~~ 4. Tidy columns ~~~~~~~~~~ #
  rename <- c(SE = "std.error", asymp.LCL = "conf.low", lower.CL = "conf.low",
              asymp.UCL = "conf.high", upper.CL = "conf.high",
              z.ratio = "statistic", t.ratio = "statistic")
  for (old in names(rename)) names(res)[names(res) == old] <- rename[[old]]
  res$null <- NULL

  res$group <- if (length(by) > 0) do.call(paste, c(lapply(res[by], as.character), sep = " / ")) else "overall"

  # ~~~~~~~~~~ 5. Multiple testing correction ~~~~~~~~~~ #
  family <- switch(adjustBy,
                   global = rep(1L, nrow(res)),
                   group  = res$group,
                   layer  = res$layer)
  res$p.adj <- stats::ave(res$p.value, family,
                          FUN = function(p) stats::p.adjust(p, method = adjust))

  res$direction   <- ifelse(res$estimate >= 0, "activator", "repressor")
  res$significant <- !is.na(res$p.adj) & res$p.adj < alpha & abs(res$estimate) >= minEffect

  # ~~~~~~~~~~ 6. Order columns, filter and names ~~~~~~~~~~ #
  first <- c("target", "regulator", "layer", by, "group")
  stats_cols <- c("estimate", "std.error", "df", "conf.low", "conf.high",
                  "statistic", "p.value", "p.adj", "direction", "significant")
  res <- res[, c(first, intersect(stats_cols, names(res))), drop = FALSE]

  if (useOriginalNames) {
    nm  <- fit$GlobalSummary$nameMap
    inv <- function(map) stats::setNames(names(map), unname(map))
    if (length(nm$target) > 0) {
      i <- inv(nm$target); hit <- res$target %in% names(i); res$target[hit] <- i[res$target[hit]]
    }
    for (l in names(priors)) {
      if (length(nm[[l]]) == 0) next
      i <- inv(nm[[l]]); hit <- res$layer == l & res$regulator %in% names(i)
      res$regulator[hit] <- i[res$regulator[hit]]
    }
  }

  if (onlySignificant) res <- res[res$significant, , drop = FALSE]
  res <- res[order(res$p.adj, na.last = TRUE), , drop = FALSE]
  rownames(res) <- NULL

  n_sig <- sum(res$significant)
  cli::cli_alert_success(
    "Tested {nrow(res)} regulator-target-group combination(s); {n_sig} significant (adjusted p < {alpha}, {adjust}, {adjustBy})."
  )
  if (!is.null(failed) && nrow(failed) > 0) {
    cli::cli_alert_warning("{nrow(failed)} regulator-target pair(s) could not be tested (see {.code attr(, \"failed\")}).")
  }

  attr(res, "failed") <- failed
  attr(res, "args")   <- list(by = by, at = at, weights = weights, alpha = alpha,
                              adjust = adjust, adjustBy = adjustBy, minEffect = minEffect)
  res
}
