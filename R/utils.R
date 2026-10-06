#' Check for feature names that would break model fitting
#'
#' Regulators of all layers are stacked in a single matrix and added as columns
#' of the model data. This function aborts when (i) the same regulator name is
#' used in more than one layer or (ii) a regulator name collides with a
#' metadata column or with the reserved name \code{TargetExpression}.
#'
#' @param regulatoryData Named list of matrices, sparse matrices or data.frames
#'   (regulators x cells), one per regulatory layer.
#' @param metadata Data.frame with cell metadata. Row names must be cell IDs.
#'
#' @return \code{NULL}, invisibly. Called for its side effect of aborting.
#'
#' @keywords internal
checkNameCollisions <- function(regulatoryData, metadata) {
  all_regs <- unlist(lapply(regulatoryData, rownames), use.names = FALSE)

  dup <- unique(all_regs[duplicated(all_regs)])
  if (length(dup) > 0) {
    cli::cli_abort(c(
      "Regulator names must be unique across regulatory layers.",
      "x" = "Duplicated: {.val {utils::head(dup, 5)}}{if (length(dup) > 5) ' ...' else ''}."
    ))
  }

  reserved <- intersect(all_regs, c("TargetExpression", colnames(metadata)))
  if (length(reserved) > 0) {
    cli::cli_abort(c(
      "Regulator names collide with {.arg metadata} columns or with the reserved name {.val TargetExpression}.",
      "x" = "Colliding names: {.val {utils::head(reserved, 5)}}."
    ))
  }
  invisible(NULL)
}

# checkInputNames
# Features are sanitize
sanitize <- function(raw, label) {
  if (anyDuplicated(raw)) {
    cli::cli_abort("Duplicated feature names found in {.field {label}}.")
  }
  clean <- make.names(raw)
  if (anyDuplicated(clean)) {
    cli::cli_abort(c(
      "Sanitizing feature names of {.field {label}} produced duplicated identifiers.",
      "i" = "Please rename the features so that they remain unique after {.fn make.names}."
    ))
  }
  changed <- raw != clean
  if (any(changed)) {
    cli::cli_alert_warning(
      "Non-syntactic characters detected in {.field {label}} ({sum(changed)} feature(s)). Sanitizing identifiers..."
    )
  }
  list(clean = clean, map = stats::setNames(clean[changed], raw[changed]))
}

# fitMORE
# A single matrix / data.frame is wrapped in a list
asLayerList <- function(x) {
  if (is.data.frame(x) || is.matrix(x) || inherits(x, "Matrix")) list(x) else x
}

# filterInputData()
addExcluded <- function(tbl, features, layer, reason) {
  if (length(features) == 0) return(tbl)
  rbind(tbl, data.frame(feature = features, layer = layer,
                        exclusion_reason = reason, stringsAsFactors = FALSE))
}

# Excludes non-finite and zero-variance rows of one matrix
qualityFilter <- function(mat, layer) {
  q <- rowQuality(mat)
  excludedFeatures <<- addExcluded(excludedFeatures, rownames(mat)[q$nonFinite],
                                   layer, "Contains NA, NaN, or Inf values")
  excludedFeatures <<- addExcluded(excludedFeatures, rownames(mat)[q$zeroVar],
                                   layer, "Zero variance across cells")
  if (any(q$nonFinite) || any(q$zeroVar)) {
    cli::cli_alert_info(
      "{.field {layer}}: excluded {sum(q$nonFinite)} feature(s) with NA/NaN/Inf and {sum(q$zeroVar)} with zero variance."
    )
  }
  mat[!q$nonFinite & !q$zeroVar, , drop = FALSE]
}
