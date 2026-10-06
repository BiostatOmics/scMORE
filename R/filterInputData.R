#' Filter data quality and prior associations
#'
#' Evaluates numeric validity (NA, NaN, Inf and zero variance) of the target
#' features and of every regulatory layer, intersects the prior associations
#' with the features that survived, drops layers without usable priors, removes
#' targets without regulators and regulators without targets, and records every
#' excluded feature together with the reason.
#'
#' @details
#' The quality metrics are computed with sparse-aware routines, so sparse
#' inputs (for example ATAC counts stored as \code{dgCMatrix}) are not
#' converted to dense matrices. Duplicated prior pairs are removed.
#'
#' @param targetData Matrix, sparse \code{Matrix} or data.frame of target
#'   features (features x cells).
#' @param regulatoryData Named list of matrices, sparse matrices or data.frames
#'   (regulators x cells), one per regulatory layer.
#' @param priorAssociations Named list of data.frames with the prior
#'   regulator-target pairs of each layer (first column = regulator, second
#'   column = target feature). Any additional columns are preserved.
#'
#' @return A list with the filtered \code{targetData}, \code{regulatoryData} and
#'   \code{priorAssociations}, and the data.frame \code{excludedFeatures}
#'   (columns \code{feature}, \code{layer} and \code{exclusion_reason}).
#'
#' @importFrom cli cli_abort cli_alert_info cli_alert_warning cli_alert_success
#' @keywords internal
filterInputData <- function(targetData,
                            regulatoryData,
                            priorAssociations) {

  excludedFeatures <- data.frame(
    feature = character(0),
    layer = character(0),
    exclusion_reason = character(0),
    stringsAsFactors = FALSE
  )

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

  # ~~~~~~~~~~ 1. Zero variance and quality check on targets ~~~~~~~~~~ #
  targetData <- qualityFilter(targetData, "target")
  if (nrow(targetData) == 0) {
    cli::cli_abort("No valid target features remaining after filtering out NA values and zero-variance features.")
  }

  # ~~~~~~~~~~ 2. Zero variance and quality check on regulatory layers ~~~~~~~~~~ #
  cleaned_regulatoryData <- lapply(names(regulatoryData), function(m_name) {
    qualityFilter(regulatoryData[[m_name]], m_name)
  })
  names(cleaned_regulatoryData) <- names(regulatoryData)

  # ~~~~~~~~~~ 3. Validate layers and prior associations ~~~~~~~~~~ #
  valid_layers <- character(0)
  filtered_priors <- list()

  for (m_name in names(cleaned_regulatoryData)) {
    assoc <- priorAssociations[[m_name]]

    if (is.null(assoc) || nrow(assoc) == 0) {
      cli::cli_alert_warning("Layer {.field {m_name}}: {.arg priorAssociations} is NULL or empty. Excluding this regulatory omic.")
      next
    }

    if (ncol(assoc) < 2) {
      cli::cli_alert_warning("Layer {.field {m_name}}: {.arg priorAssociations} has fewer than 2 columns. Excluding layer.")
      next
    }

    # Intersect prior pairs with available valid features (col 1 = regulator, col 2 = target)
    reg_features <- rownames(cleaned_regulatoryData[[m_name]])
    target_features <- rownames(targetData)

    keep_idx <- (assoc[[1]] %in% reg_features) & (assoc[[2]] %in% target_features)
    valid_assoc <- unique(assoc[keep_idx, , drop = FALSE])
    rownames(valid_assoc) <- NULL

    if (nrow(valid_assoc) == 0) {
      cli::cli_alert_warning("Layer {.field {m_name}}: No matching target-regulator pairs found in expression matrices. Excluding layer.")
      next
    }

    filtered_priors[[m_name]] <- valid_assoc
    valid_layers <- c(valid_layers, m_name)
  }

  if (length(valid_layers) == 0) {
    cli::cli_abort("No regulatory layers remaining with valid prior associations and matching expression data.")
  }

  regulatoryData <- cleaned_regulatoryData[valid_layers]
  priorAssociations <- filtered_priors[valid_layers]

  # ~~~~~~~~~~ 4. Subset active targets and regulators ~~~~~~~~~~ #
  # Filter target features without any active prior association
  active_targets_in_priors <- unique(unlist(lapply(priorAssociations, function(df) df[[2]])))
  targets_without_regulators <- setdiff(rownames(targetData), active_targets_in_priors)

  if (length(targets_without_regulators) > 0) {
    cli::cli_alert_info("{length(targets_without_regulators)} target feature(s) have no prior regulatory associations and will be excluded.")

    excludedFeatures <- addExcluded(excludedFeatures, targets_without_regulators,
                                    "target", "No prior regulatory associations found")

    targetData <- targetData[rownames(targetData) %in% active_targets_in_priors, , drop = FALSE]
  }

  if (nrow(targetData) == 0) {
    cli::cli_abort("No target features remaining after filtering for prior associations.")
  }

  # Filter regulators without any active prior association in their respective layer
  for (m_name in names(regulatoryData)) {
    active_regs <- unique(priorAssociations[[m_name]][[1]])
    unused_regs <- setdiff(rownames(regulatoryData[[m_name]]), active_regs)

    if (length(unused_regs) > 0) {
      cli::cli_alert_info("Layer {.field {m_name}}: {length(unused_regs)} regulator(s) have no active target associations and will be excluded.")

      excludedFeatures <- addExcluded(excludedFeatures, unused_regs, m_name,
                                      "No active prior target associations found")

      regulatoryData[[m_name]] <- regulatoryData[[m_name]][active_regs, , drop = FALSE]
    }
  }

  cli::cli_alert_success("Data quality filtering and prior matching completed successfully.")

  return(list(
    targetData        = targetData,
    regulatoryData    = regulatoryData,
    priorAssociations = priorAssociations,
    excludedFeatures  = excludedFeatures
  ))
}
