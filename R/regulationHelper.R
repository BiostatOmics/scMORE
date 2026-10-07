# Internal helpers for the functions that analyse the output of fitMore()
# (testContrasts() and future exploration functions) --------------------------

#' Map every regulator to its regulatory layer
#'
#' @param fit An object of class \code{scMoreFit}.
#'
#' @return A named character vector: names are regulators and values are the
#'   names of their layers.
#'
#' @keywords internal
layerLookup <- function(fit) {
  priors <- fit$GlobalSummary$priorAssociations
  unlist(lapply(names(priors), function(l) {
    regs <- unique(as.character(priors[[l]]$regulator))
    stats::setNames(rep(l, length(regs)), regs)
  }))
}


#' Grid of fixed-effect combinations stored in a scMoreFit
#'
#' Reads the \code{emmGrid} of the first target-regulator pair that has a
#' fitted model. All pairs share the same grid, so it is used to validate user
#' input (variables and levels) without looping over every pair.
#'
#' @param models Named list of models (\code{fit$Models}, possibly subset).
#'
#' @return A data.frame with one row per combination of fixed-effect levels.
#'
#' @importFrom emmeans as.emmGrid
#' @keywords internal
firstTrendGrid <- function(models) {
  for (m in models) {
    ok <- Filter(Negate(is.null), m$fit)
    if (length(ok) > 0) return(emmeans::as.emmGrid(ok[[1]])@grid)
  }
  cli::cli_abort("None of the selected targets has a fitted model.")
}


#' Validate a list of levels against the stored grid
#'
#' @param at Named list of levels, or \code{NULL}.
#' @param fixedEffects Character vector with the fixed effects of the model.
#' @param grid Data.frame returned by \code{\link{firstTrendGrid}}.
#'
#' @return \code{NULL}, invisibly. Aborts with an informative message if a
#'   variable is not a fixed effect or a level does not exist.
#'
#' @keywords internal
checkAt <- function(at, fixedEffects, grid) {
  if (is.null(at)) return(invisible(NULL))
  if (!is.list(at) || is.null(names(at))) {
    cli::cli_abort("{.arg at} must be a named list, for example {.code list(cell_type = \"Tcell\")}.")
  }
  bad_at <- setdiff(names(at), fixedEffects)
  if (length(bad_at) > 0) {
    cli::cli_abort("{.arg at} contains variables that are not fixed effects of the model: {.val {bad_at}}.")
  }
  for (v in names(at)) {
    avail <- unique(as.character(grid[[v]]))
    miss  <- setdiff(as.character(at[[v]]), avail)
    if (length(miss) > 0) {
      cli::cli_abort(c(
        "Level(s) {.val {miss}} not found in {.field {v}}.",
        "i" = "Available levels: {.val {avail}}."
      ))
    }
  }
  invisible(NULL)
}


#' Restrict a trend grid to the requested levels
#'
#' @param emm An \code{emmGrid} reconstructed from the stored slopes.
#' @param at Named list of levels to keep.
#'
#' @return The \code{emmGrid} restricted to the rows that match \code{at}, or
#'   \code{NULL} if no row matches.
#'
#' @keywords internal
subsetTrendGrid <- function(emm, at) {
  keep <- rep(TRUE, nrow(emm@grid))
  for (v in names(at)) {
    keep <- keep & (as.character(emm@grid[[v]]) %in% as.character(at[[v]]))
  }
  if (!any(keep)) return(NULL)
  emm[which(keep)]
}


#' Number of levels of a fixed effect after applying \code{at}
#'
#' @param v Character. Name of the fixed effect.
#' @param at Named list of levels, or \code{NULL}.
#' @param grid Data.frame returned by \code{\link{firstTrendGrid}}.
#'
#' @return An integer.
#'
#' @keywords internal
nLevelsAt <- function(v, at, grid) {
  if (!is.null(at[[v]])) length(unique(at[[v]])) else length(unique(as.character(grid[[v]])))
}


#' Translate sanitized names back to the original identifiers
#'
#' @param res Data.frame with the columns \code{target}, \code{regulator} and
#'   \code{layer}.
#' @param fit An object of class \code{scMoreFit} (uses
#'   \code{GlobalSummary$nameMap}).
#'
#' @return \code{res} with \code{target} and \code{regulator} translated.
#'
#' @keywords internal
restoreOriginalNames <- function(res, fit) {
  nm  <- fit$GlobalSummary$nameMap
  inv <- function(map) stats::setNames(names(map), unname(map))

  if (length(nm$target) > 0) {
    i   <- inv(nm$target)
    hit <- res$target %in% names(i)
    res$target[hit] <- i[res$target[hit]]
  }
  for (l in names(fit$GlobalSummary$priorAssociations)) {
    if (length(nm[[l]]) == 0) next
    i   <- inv(nm[[l]])
    hit <- res$layer == l & res$regulator %in% names(i)
    res$regulator[hit] <- i[res$regulator[hit]]
  }
  res
}


#' Validate the requested regulatory layers
#'
#' @param layers Character vector of layer names, or \code{NULL}.
#' @param fit An object of class \code{scMoreFit}.
#'
#' @return \code{NULL}, invisibly. Aborts if a layer does not exist.
#'
#' @keywords internal
checkLayers <- function(layers, fit) {
  if (is.null(layers)) return(invisible(NULL))
  available  <- names(fit$GlobalSummary$priorAssociations)
  bad_layers <- setdiff(layers, available)
  if (length(bad_layers) > 0) {
    cli::cli_abort("Unknown layer(s): {.val {bad_layers}}. Available: {.val {available}}.")
  }
  invisible(NULL)
}


#' Select the models of the requested targets
#'
#' @param fit An object of class \code{scMoreFit}.
#' @param targets Character vector of targets (sanitized names), or \code{NULL}
#'   for all of them.
#'
#' @return The sublist of \code{fit$Models} for the requested targets. Aborts if
#'   there is none.
#'
#' @keywords internal
selectModels <- function(fit, targets = NULL) {
  models <- fit$Models
  if (!is.null(targets)) models <- models[intersect(targets, names(models))]
  if (length(models) == 0) cli::cli_abort("No targets to analyse.")
  models
}


#' Give common names to the columns returned by emmeans
#'
#' \pkg{emmeans} names columns differently for mixed models (asymptotic z) and
#' for \code{lm} (t with degrees of freedom), so both sets of names are mapped to
#' common ones.
#'
#' @param res Data.frame built from \code{summary()} of \code{emmeans} objects.
#'
#' @return \code{res} with the columns \code{std.error}, \code{conf.low},
#'   \code{conf.high} and \code{statistic}, and without the \code{null} column.
#'
#' @keywords internal
tidyEmmeansColumns <- function(res) {
  rename <- c(SE = "std.error", asymp.LCL = "conf.low", lower.CL = "conf.low",
              asymp.UCL = "conf.high", upper.CL = "conf.high",
              z.ratio = "statistic", t.ratio = "statistic")
  for (old in names(rename)) names(res)[names(res) == old] <- rename[[old]]
  res$null <- NULL
  res
}


#' Apply a test function to every target-regulator pair
#'
#' Loops over the targets (with a progress bar) and, inside, over their
#' regulators. A failure in one pair never stops the loop: it is recorded
#' together with its reason and returned, so the caller can report it.
#'
#' @param models Named list of models (see \code{\link{selectModels}}).
#' @param regulators Optional character vector to restrict the regulators.
#' @param layers Optional character vector to restrict the layers.
#' @param layerOf Named vector regulator -> layer (see \code{\link{layerLookup}}).
#' @param testOne Function with arguments \code{(trend, reg)} that returns a
#'   data.frame for one pair, or fails.
#' @param verb Character. Used in the error message when no pair can be
#'   processed (for example \code{"tested"} or \code{"compared"}).
#'
#' @return A list with \code{res} (data.frame with one block of rows per pair,
#'   preceded by the columns \code{target}, \code{regulator} and \code{layer})
#'   and \code{failed} (data.frame with the pairs that failed and the reason, or
#'   \code{NULL}).
#'
#' @importFrom pbapply pblapply
#' @keywords internal
testPairs <- function(models, regulators, layers, layerOf, testOne, verb = "tested") {

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

  res    <- do.call(rbind, lapply(res_list, `[[`, "ok"))
  failed <- do.call(rbind, lapply(res_list, `[[`, "fail"))
  if (!is.null(failed)) {
    failed <- as.data.frame(failed, stringsAsFactors = FALSE)
    names(failed) <- c("target", "regulator", "reason")
  }

  if (is.null(res) || nrow(res) == 0) {
    cli::cli_abort(c(
      "No regulator-target pairs could be {verb}.",
      "i" = if (!is.null(failed) && nrow(failed) > 0) "First error: {failed$reason[1]}" else
        "Check {.arg at}, {.arg regulators} and {.arg layers}."
    ))
  }

  list(res = res, failed = failed)
}
