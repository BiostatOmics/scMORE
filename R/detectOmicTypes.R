#' Normalize a user-supplied omic type vector
#'
#' Recycles a single value to all layers, assigns layer names to unnamed
#' vectors and matches named vectors to the layer names.
#'
#' @param omicType Numeric vector with values 0 (continuous) or 1 (binary).
#' @param layerNames Character vector with the names of the regulatory layers.
#'
#' @return A named numeric vector with one value per layer, in the order of
#'   \code{layerNames}.
#'
#' @importFrom cli cli_abort
#' @keywords internal
normalizeOmicType <- function(omicType, layerNames) {
  if (!all(omicType %in% c(0, 1))) {
    cli::cli_abort("{.arg omicType} values must be {.val 0} (continuous) or {.val 1} (binary).")
  }
  if (length(omicType) == 1 && length(layerNames) > 1) {
    omicType <- stats::setNames(rep(unname(omicType), length(layerNames)), layerNames)
  }
  if (is.null(names(omicType))) {
    if (length(omicType) != length(layerNames)) {
      cli::cli_abort("Length of {.arg omicType} ({length(omicType)}) must match the number of regulatory layers ({length(layerNames)}).")
    }
    names(omicType) <- layerNames
  } else {
    absent <- setdiff(layerNames, names(omicType))
    if (length(absent) > 0) {
      cli::cli_abort("{.arg omicType} has no value for layer(s): {.val {absent}}.")
    }
    omicType <- omicType[layerNames]
  }
  omicType
}


#' Infers (or validates) whether each regulatory layer is continuous (0) or
#' binary (1). Automatic detection labels a layer as binary when all its values
#' belong to \code{c(0, 1, NA)}. Sparse matrices are inspected without being
#' densified.
#'
#' @param regulatoryData Named list of matrices, sparse matrices or data.frames
#'   (regulators x cells), one per regulatory layer.
#' @param omicType Optional numeric vector (0 = continuous, 1 = binary). It can
#'   be a single value (applied to all layers), an unnamed vector with one value
#'   per layer, or a named vector. If \code{NULL} the types are inferred.
#'
#' @return A named numeric vector with the type (0 or 1) of each regulatory
#'   layer.
#'
#' @importFrom cli cli_alert_info cli_text cli_ul cli_li cli_end cli_h2
#' @keywords internal
detectOmicTypes <- function(regulatoryData, omicType = NULL) {

  cli::cli_text("Encoding standard: {.val 0} = Continuous/Numeric | {.val 1} = Binary")
  # ~~~~~~~~~~ 1. Automatic detection or user validation ~~~~~~~~~~ #
  if (is.null(omicType)) {
    # ~~~~~~~~~~ 1.1 Automatic detection ~~~~~~~~~~ #
    cli::cli_h2("Detecting regulatory omic types")
    omicType <- sapply(regulatoryData, function(mat) {
      u_vals <- unique(as.vector(mat))
      if (all(u_vals %in% c(0, 1, NA))) 1 else 0
    })
    names(omicType) <- names(regulatoryData)
  } else {
    # ~~~~~~~~~~ 1.2 User validation ~~~~~~~~~~ #
    omicType <- normalizeOmicType(omicType, names(regulatoryData))
  }

  # ~~~~~~~~~~ 2. Display layer classification summary ~~~~~~~~~~ #
  cli::cli_ul()
  for (name in names(omicType)) {
    type_label <- if (omicType[[name]] == 1) "Binary (1)" else "Continuous (0)"
    cli::cli_li("{.field {name}}: {type_label}")
  }
  cli::cli_end()

  cli::cli_alert_info("If this classification is incorrect, stop execution and specify {.arg omicType} manually.")

  return(omicType)
}
