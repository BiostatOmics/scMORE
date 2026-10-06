#' Validate and sanitize input identifiers across datasets
#'
#' Aligns cell identifiers across the target matrix, the metadata and every
#' regulatory layer (keeping the same cells in the same order everywhere), and
#' sanitizes non-syntactic feature names so they can be used in model formulas.
#' Prior association tables are updated with the sanitized names.
#'
#' @details
#' Feature names are converted with \code{\link[base]{make.names}} (for example
#' \code{"chr1-100-200"} becomes \code{"chr1.100.200"}). The same conversion is
#' applied to the first two columns of every prior table, so matrices and priors
#' stay consistent as long as they used the same original identifiers. The
#' function aborts if the sanitization produces duplicated names.
#'
#' @param targetData Matrix, sparse \code{Matrix} or data.frame of target
#'   features (features x cells).
#' @param regulatoryData Named list of matrices, sparse matrices or data.frames
#'   (regulators x cells), one per regulatory layer.
#' @param metadata Data.frame with cell metadata. Row names must be cell IDs.
#' @param priorAssociations Named list of data.frames with the prior
#'   regulator-target pairs of each layer (first column = regulator, second
#'   column = target feature). Any additional columns are preserved.
#'
#' @return A named list with the aligned and sanitized \code{targetData},
#'   \code{regulatoryData}, \code{metadata} and \code{priorAssociations}
#'   (whose first two columns are renamed to \code{regulator} and \code{target}),
#'   plus \code{nameMap}: a list with the original and sanitized names of the
#'   features that were renamed (\code{target} and one entry per layer).
#'
#' @importFrom cli cli_abort cli_alert_info cli_alert_warning cli_alert_success
#' @keywords internal
checkInputNames <- function(targetData,
                            regulatoryData,
                            metadata,
                            priorAssociations) {

  # ~~~~~~~~~~ 1. Coercion and early validation ~~~~~~~~~~ #
  targetData <- asFeatureMatrix(targetData, "targetData")
  for (layer in names(regulatoryData)) {
    regulatoryData[[layer]] <- asFeatureMatrix(regulatoryData[[layer]],
                                               paste0("regulatoryData$", layer))
  }
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (is.null(colnames(targetData)) || is.null(rownames(targetData))) {
    cli::cli_abort("{.arg targetData} must have valid row and column names.")
  }
  if (is.null(rownames(metadata))) {
    cli::cli_abort("Rownames of {.arg metadata} must contain valid cell IDs.")
  }
  for (layer in names(regulatoryData)) {
    if (is.null(colnames(regulatoryData[[layer]])) || is.null(rownames(regulatoryData[[layer]]))) {
      cli::cli_abort("Regulatory layer {.val {layer}} must have valid row and column names.")
    }
  }

  # ~~~~~~~~~~ 2. Aligning cell identifiers (colnames) ~~~~~~~~~~ #
  # Align cells (columns) across expression matrices and metadata
  reg_cell_lists <- lapply(regulatoryData, colnames)
  all_cell_sources <- c(list(target = colnames(targetData),
                             metadata = rownames(metadata)),
                        reg_cell_lists)

  # Intersect all cell ids
  common_cells <- Reduce(intersect, all_cell_sources)

  if (length(common_cells) == 0) {
    cli::cli_abort("No overlapping cell IDs found across {.arg targetData}, {.arg metadata}, and {.arg regulatoryData}.")
  }

  # Calculate total initial cells per input to detect if subsetting is needed
  initial_cell_counts <- sapply(all_cell_sources, length)

  if (any(initial_cell_counts > length(common_cells))) {
    cli::cli_alert_info(
      "Subsetting dataset to {.val {length(common_cells)}} intersecting cells across all inputs."
    )
    targetData     <- targetData[, common_cells, drop = FALSE]
    metadata       <- metadata[common_cells, , drop = FALSE]
    regulatoryData <- lapply(regulatoryData, function(mat) mat[, common_cells, drop = FALSE])
  }

  # ~~~~~~~~~~ 3. Sanitizing feture names (rownames) ~~~~~~~~~~ #
  nameMap <- list()
  # Target
  tg <- sanitize(rownames(targetData), "target")
  rownames(targetData) <- tg$clean
  nameMap$target <- tg$map
  # Regulators
  for (m_name in names(regulatoryData)) {
    rg <- sanitize(rownames(regulatoryData[[m_name]]), m_name)
    rownames(regulatoryData[[m_name]]) <- rg$clean
    nameMap[[m_name]] <- rg$map
  }

  # ~~~~~~~~~~ 4. Sanitizing prior associations ~~~~~~~~~~ #
  # Update priorAssociations entries with sanitized names

  sanitized_priors <- lapply(names(priorAssociations), function(layer_name) {
    assoc_table <- priorAssociations[[layer_name]]

    if (is.null(assoc_table) || nrow(assoc_table) == 0) {
      return(assoc_table)
    }

    # Convert to data.frame to ensure column indexing stability
    assoc_df <- as.data.frame(assoc_table, stringsAsFactors = FALSE)

    if (ncol(assoc_df) < 2) {
      cli::cli_abort("Layer {.field {layer_name}} in {.arg priorAssociations} must contain at least 2 columns.")
    }

    # First column = regulator, second column = target. Extra columns are kept
    colnames(assoc_df) = c("regulator", "target")
    colnames(assoc_df)[1:2] <- c("regulator", "target")
    assoc_df[[1]] <- make.names(as.character(assoc_df[[1]]))
    assoc_df[[2]] <- make.names(as.character(assoc_df[[2]]))

    return(assoc_df)
  })
  names(sanitized_priors) <- names(priorAssociations)

  cli::cli_alert_success("Identifiers successfully validated and sanitized.")
  # Return data
  return(list(targetData        = targetData,
              regulatoryData    = regulatoryData,
              metadata          = metadata,
              priorAssociations = sanitized_priors)
  )
}
