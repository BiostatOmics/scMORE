# Internal helpers shared across scMore functions -----------------------------

#' Coerce an input into a numeric feature matrix
#'
#' Converts data.frames and base matrices into a base matrix, and leaves
#' \pkg{Matrix} objects (for example \code{dgCMatrix}) untouched so that sparse
#' data are never densified unnecessarily.
#'
#' @param x A matrix, sparse \code{Matrix} or data.frame (features x cells).
#' @param arg Character. Argument name used in error messages.
#'
#' @return A matrix or \code{Matrix} object with the same dimnames as \code{x}.
#'
#' @keywords internal
asFeatureMatrix <- function(x, arg = "input") {
  if (inherits(x, "Matrix")) return(x)
  if (is.data.frame(x) || is.matrix(x)) {
    x <- as.matrix(x)
    if (!is.numeric(x) && !is.logical(x)) {
      cli::cli_abort("{.arg {arg}} must be numeric.")
    }
    return(x)
  }
  cli::cli_abort("{.arg {arg}} must be a matrix, a sparse Matrix or a data.frame.")
}


#' Row-wise data quality metrics (sparse aware)
#'
#' Flags features (rows) with non-finite values (NA, NaN, Inf) and features
#' with (numerically) zero variance. Works on dense and sparse matrices without
#' converting the latter to dense. A row is considered constant when its
#' variance is below \code{tol} times its mean square, which makes the criterion
#' independent of the measurement scale.
#'
#' @param mat Numeric matrix or sparse \code{Matrix} (features x cells).
#' @param tol Numeric. Relative tolerance used to declare zero variance.
#'
#' @return A list with the logical vectors \code{nonFinite} and \code{zeroVar}
#'   (both of length \code{nrow(mat)}) and the numeric vector \code{variance}.
#'
#' @keywords internal
rowQuality <- function(mat, tol = 1e-12) {
  n  <- ncol(mat)
  nr <- nrow(mat)

  if (inherits(mat, "sparseMatrix")) {
    tm <- methods::as(mat, "TsparseMatrix")
    x  <- if ("x" %in% methods::slotNames(tm)) tm@x else rep(1, length(tm@i))
    nonFinite <- logical(nr)
    nonFinite[unique(tm@i[!is.finite(x)]) + 1L] <- TRUE
  } else {
    nonFinite <- rowSums(!is.finite(mat)) > 0
  }

  s1 <- Matrix::rowSums(mat, na.rm = TRUE)
  s2 <- Matrix::rowSums(mat^2, na.rm = TRUE)
  variance <- (s2 - s1^2 / n) / (n - 1)
  zeroVar  <- !nonFinite & (!is.finite(variance) | variance <= tol * (s2 / n))

  list(nonFinite = nonFinite, zeroVar = zeroVar, variance = variance)
}


#' Distinct values of a matrix (sparse aware)
#'
#' Returns the distinct values of a matrix. For sparse matrices only the stored
#' entries are inspected and the implicit zero is added when the matrix has
#' empty cells, so the matrix is never densified.
#'
#' @param mat Numeric matrix or sparse \code{Matrix}.
#'
#' @return A numeric vector of unique values.
#'
#' @keywords internal
uniqueValues <- function(mat) {
  if (inherits(mat, "sparseMatrix")) {
    tm <- methods::as(mat, "TsparseMatrix")
    x  <- if ("x" %in% methods::slotNames(tm)) tm@x else rep(1, length(tm@i))
    vals <- unique(x)
    if (length(x) < prod(dim(mat))) vals <- unique(c(vals, 0))
    return(vals)
  }
  unique(as.vector(mat))
}


#' Build a lookup list of regulators per target feature
#'
#' Combines the prior association tables of all regulatory layers and groups
#' the regulators by target feature, so the regulators of a given target can be
#' retrieved instantly with \code{dict[[target]]}.
#'
#' @param priorAssociations Named list of data.frames (one per regulatory
#'   layer) or a single data.frame. Each table must contain the columns
#'   \code{regulator} and \code{gene} (as produced by \code{\link{checkInputNames}}).
#'   Additional columns are ignored.
#'
#' @return A named list. Names are target features and each element is a
#'   character vector with the unique regulators associated with that target.
#'
#' @keywords internal
preparePriorDict <- function(priorAssociations) {
  if (is.data.frame(priorAssociations)) {
    priorAssociations <- list(priorAssociations)
  }

  prior_df <- do.call(rbind, lapply(priorAssociations, function(df) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    df[, c("regulator", "gene"), drop = FALSE]
  }))
  prior_df <- unique(prior_df)
  rownames(prior_df) <- NULL

  split(as.character(prior_df$regulator), as.character(prior_df$gene))
}


#' Check for feature names that would break model fitting
#'
#' Regulators of all layers are stacked in a single matrix and added as columns
#' of the model data. This function aborts when (i) the same regulator name is
#' used in more than one layer or (ii) a regulator name collides with a
#' metadata column or with the reserved name \code{TargetExpression}.
#'
#' @param regulatoryData Named list of regulatory matrices.
#' @param metadata Data.frame with cell metadata.
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


#' Sanitize feature names
#'
#' Converts feature identifiers into syntactically valid R names with
#' \code{\link[base]{make.names}} and keeps track of the ones that changed.
#' Aborts if the original or the sanitized names contain duplicates, because
#' features could no longer be told apart.
#'
#' @param raw Character vector with the original feature names.
#' @param label Character. Name of the dataset or layer, used in messages.
#'
#' @return A list with \code{clean} (sanitized names, same order as
#'   \code{raw}) and \code{map} (named vector: original -> sanitized, only for
#'   the names that changed).
#'
#' @keywords internal
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


#' Wrap a single matrix or data.frame in a list
#'
#' A data.frame is itself a list, so \code{is.list()} cannot be used to detect
#' a single layer. Lists are returned unchanged.
#'
#' @param x A matrix, sparse \code{Matrix}, data.frame or a list of them.
#'
#' @return A list.
#'
#' @keywords internal
asLayerList <- function(x) {
  if (is.data.frame(x) || is.matrix(x) || inherits(x, "Matrix")) list(x) else x
}


#' Append excluded features to the exclusion table
#'
#' @param tbl Data.frame with columns \code{feature}, \code{layer} and
#'   \code{exclusion_reason}.
#' @param features Character vector of excluded features (can be empty).
#' @param layer Character. Layer the features belong to.
#' @param reason Character. Reason for the exclusion.
#'
#' @return \code{tbl} with the new rows appended (unchanged if \code{features}
#'   is empty).
#'
#' @keywords internal
addExcluded <- function(tbl, features, layer, reason) {
  if (length(features) == 0) return(tbl)
  rbind(tbl, data.frame(feature = features, layer = layer,
                        exclusion_reason = reason, stringsAsFactors = FALSE))
}
