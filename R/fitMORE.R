#' Fit regulatory models for every target feature
#'
#' Main function of the first stage of \pkg{scMore}. For each target feature
#' (for example a gene) it fits a (generalized) linear mixed model at single-cell
#' level in which the expression is explained by its prior regulators (for
#' example ATAC regions and transcription factors), the experimental design
#' (fixed effects such as condition or cell type), the interactions between
#' regulators and design, and random intercepts (for example sample). The fitted
#' models are stored as \code{emmeans} trend lists so that regulatory effects
#' can be contrasted between conditions or cell types in the second stage.
#'
#' @details
#' The pipeline validates and sanitizes the input identifiers
#' (\code{\link{checkInputNames}}), validates the design
#' (\code{\link{validateModelDesign}}), filters low-quality features and priors
#' (\code{\link{filterInputData}}), detects the omic types
#' (\code{\link{detectOmicTypes}}), optionally residualizes the regulators on the
#' random effects (\code{\link{prepareTwoSteps}}) and finally fits one model per
#' target (\code{\link{fitTargetModel}}).
#'
#' For a target \eqn{g} with prior regulators \eqn{R_g}, cell \eqn{c} and sample
#' \eqn{s(c)}, the model with interactions is
#' \deqn{y_{gc} = \beta_0 + \sum_{r \in R_g} \beta_r x_{rc} + \sum_{f} \gamma_f z_{fc} +
#' \sum_{r \in R_g}\sum_{f} \delta_{rf} x_{rc} z_{fc} + u_{s(c)} + \epsilon_{gc}}{
#' y = b0 + sum(b_r x_r) + sum(g_f z_f) + sum(d_rf x_r z_f) + u_s + e}
#' with \eqn{u_s \sim N(0, \sigma_u^2)}{u_s ~ N(0, s_u^2)} and
#' \eqn{\epsilon_{gc} \sim N(0, \sigma^2)}{e ~ N(0, s^2)}, where \eqn{z_{fc}} are
#' the fixed effects (and all their interactions). The number of fixed
#' coefficients grows as \eqn{(k + 1)\prod_f L_f}{(k + 1) * prod(L_f)}, with
#' \eqn{k} regulators and \eqn{L_f} levels per fixed effect, so the number of
#' regulators per target and of design levels should be kept moderate relative to
#' the number of cells. The target is always modelled as Gaussian. The omic type
#' only affects step 1 of the two-step approach.
#'
#' With \code{twoSteps = TRUE} (default) every regulator is first replaced by
#' the residuals of a null model \code{regulator ~ (1 | random effects)}. The
#' regulator effects are then within-sample (cell-level) effects. This step is
#' skipped, with a warning, if no random effects are given.
#'
#' @param targetData Matrix, sparse \code{Matrix} or data.frame with the target
#'   omic (for example log-normalized gene expression). Features (for example
#'   gene IDs) are in rows and cells in columns.
#' @param regulatoryData Named list in which each element is a different
#'   regulatory omic (for example chromatin accessibility or transcription
#'   factor activity). Each element is a matrix, sparse \code{Matrix} or
#'   data.frame with regulators in rows and cells in columns. A single matrix is
#'   also accepted.
#' @param priorAssociations Named list with one data.frame per regulatory omic
#'   (same names as \code{regulatoryData}; unnamed lists are matched by
#'   position). The first column holds the regulators, the second the target
#'   features and any further column is ignored. A single data.frame is also
#'   accepted. \code{NULL} priors (all regulators for all targets) are not
#'   supported in the current version.
#' @param metadata Data.frame with cell metadata. Row names must be the cell IDs
#'   used as column names of the omic matrices.
#' @param fixedEffects Character vector of \code{metadata} columns used as
#'  fixed tearm (for example \code{c("condition", "cell_type")}).
#'   By default, \code{NULL}.
#' @param randomEffects Character vector of \code{metadata} columns used
#'  as random intercepts (for example \code{"sample"}), if null linear models
#'  are fitted. By default, \code{NULL}.
#' @param omicType Numeric vector indicating whether each regulatory omic is
#'   continuous (0) or binary (1). It can be a single value (applied to all
#'   omics), an unnamed vector with one value per omic or a named vector. If
#'   \code{NULL} the types are inferred and printed. If the inferred types are
#'   incorrect, please stop the process and set them manually. By default,
#'   \code{NULL}.
#' @param twoSteps Logical. If \code{TRUE}, regulators are residualized on the
#' random effects before fitting the models. By default, TRUE
#' @param interaction Logical. If \code{TRUE}, regulators are crossed
#'   with all fixed effects and their interactions. By default, \code{TRUE}
#' @param parallel Reserved for future use. Parallelization is not implemented
#'   yet and models are always fitted sequentially.By default, \code{FALSE}
#' @param formula Reserved for future use. Custom formulas are not supported yet
#'   and must be \code{NULL}.
#' @param keepResiduals Logical. If \code{TRUE}, the residualized regulatory
#'   matrices are returned in \code{GlobalSummary$residuals}. By default,
#'   \code{FALSE}, because they can be large.
#' @param seed Numeric. Seed for reproducibility. By default, 1808.
#'
#' @return An object of class \code{scMoreFit}: a list with
#' \itemize{
#'   \item \code{Models}: named list (one element per fitted target) with the
#'     output of \code{\link{fitTargetModel}}.
#'   \item \code{GlobalSummary}: a list with the arguments used (\code{args}),
#'     the \code{omicType}, the \code{nameMap} of sanitized feature names, the
#'     filtered \code{priorAssociations}, the \code{excludedFeatures} table, the
#'     per-target \code{fitSummary} table (number of regulators, convergence,
#'     singularity, warnings and errors), the \code{cells} used and, optionally,
#'     the \code{residuals}.
#' }
#'
#' @examples
#' \dontrun{
#' fit <- fitMORE(
#'   targetData        = targetData,
#'   regulatoryData    = list(ATAC = atacMatrix, TFs = tfActivity),
#'   priorAssociations = list(ATAC = regionGeneLinks, TFs = tfGeneEdges),
#'   metadata          = metadata,
#'   fixedEffects      = c("condition", "cell_type"),
#'   randomEffects     = "sample"
#' )
#' fit$GlobalSummary$fitSummary
#' }
#'
#' @importFrom cli cli_h2 cli_alert_info cli_alert_warning cli_alert_success cli_abort
#' @importFrom pbapply pblapply
#' @export
fitMORE <- function(targetData,
                    regulatoryData,
                    priorAssociations,
                    metadata,
                    fixedEffects  = NULL,
                    randomEffects = NULL,
                    omicType      = NULL,
                    twoSteps      = TRUE,
                    interaction   = TRUE,
                    parallel      = FALSE,
                    formula       = NULL,
                    keepResiduals = FALSE,
                    seed          = 1808) {

  #Set the seed for the reproducibility
  set.seed(seed)

  cli::cli_h1("Running scMORE")
  cli::cli_h2("Fit regulatory models for every target feature")

  # Internal modelling choices
  targetOmicType <- 0
  REML           <- FALSE
  verbose        <- FALSE


  if (!isFALSE(parallel)) {
    cli::cli_alert_warning("{.arg parallel} is not implemented yet: models are fitted sequentially.")
  }

  # ~~~~~~~~~~ 1. Check and prepare data ~~~~~~~~~~ #
  # A single matrix / data.frame is wrapped in a list
  regulatoryData <- asLayerList(regulatoryData)
  if (is.null(priorAssociations)) {
    cli::cli_abort("{.arg priorAssociations} = NULL (all regulators for all targets) is not supported yet.")
  }
  priorAssociations <- asLayerList(priorAssociations)

  if( is.null(priorAssociations) || ! length(priorAssociations) == length(regulatoryData)){
    stop("The number of elements in `regulatoryData` and `priorAssociations` do not match; please ensure that both are lists containing as many elements as there are regulatory omics.")
  }

  # Layer names
  if (is.null(names(regulatoryData)) || any(!nzchar(names(regulatoryData)))) {
    names(regulatoryData) <- paste0("Layer", seq_along(regulatoryData))
    cli::cli_alert_info("{.arg regulatoryData} is unnamed: layers were named {.val {names(regulatoryData)}}.")
  }
  if (anyDuplicated(names(regulatoryData))) {
    cli::cli_abort("Names of {.arg regulatoryData} must be unique.")
  }

  # Match priors to layers
  if (length(priorAssociations) != length(regulatoryData)) {
    cli::cli_abort("The number of elements in {.arg regulatoryData} and {.arg priorAssociations} do not match; both must be lists with as many elements as regulatory omics.")
  }
  if (is.null(names(priorAssociations))) {
    names(priorAssociations) <- names(regulatoryData)
  } else if (!all(names(regulatoryData) %in% names(priorAssociations))) {
    cli::cli_abort("Names of {.arg priorAssociations} must match the names of {.arg regulatoryData}.")
  }
  priorAssociations <- priorAssociations[names(regulatoryData)]

  # # Name the user-provided omic types before layers can be dropped
  # if (!is.null(omicType)) {
  #   omicType <- normalizeOmicType(omicType, names(regulatoryData))
  # }

  # 1.1 Checking feature names and cell alignment
  cli::cli_h2("Validating and Sanitizing Input Names")
  checkOut          <- checkInputNames(targetData, regulatoryData, metadata, priorAssociations)
  targetData        <- checkOut$targetData
  regulatoryData    <- checkOut$regulatoryData
  metadata          <- checkOut$metadata
  priorAssociations <- checkOut$priorAssociations
  nameMap           <- checkOut$nameMap
  #list2env(checkOut, envir = environment())

  # 1.2. Checking metadata and design
  cli::cli_h2("Validating model design and metadata integrity")
  designOut     <- validateModelDesign(metadata, fixedEffects, randomEffects, formula)
  metadata      <- designOut$metadata
  fixedEffects  <- designOut$fixedEffects
  randomEffects <- designOut$randomEffects
  #list2env(metadataOut, envir = environment())

  # 1.3. Filtering features and priors
  cli::cli_h2("Filtering data quality and prior associations")
  filterOut         <- filterInputData(targetData, regulatoryData, priorAssociations)
  targetData        <- filterOut$targetData
  regulatoryData    <- filterOut$regulatoryData
  priorAssociations <- filterOut$priorAssociations
  excludedFeatures  <- filterOut$excludedFeatures
  #list2env(metadataOut, envir = environment())

  # 1.4 Check model variable names
  checkNameCollisions(regulatoryData, metadata)

  # 1.5 Eval OmicType
  cli::cli_h2("Omic Data Type Summary")
  residualized <- FALSE
  omicType <- detectOmicTypes(regulatoryData, omicType)

  # ~~~~~~~~~~ 2. Two-Steps models ~~~~~~~~~~ #
  # Residualize regulatory features across cells
  if (twoSteps) {
    cli::cli_h2("Step 1: Residualizing regulatory features")
    if(length(randomEffects) > 0){
      regulatoryData <- prepareTwoSteps(
      regulatoryData = regulatoryData,
      metadata       = metadata,
      randomEffects  = randomEffects,
      omicType       = omicType,
      REML           = REML,
      verbose        = verbose
    )
    residualized <- TRUE
    } else {
      cli::cli_alert_warning("{.arg randomEffects} not specified: skipping Step 1 (residualization).")
    }
  }

  # ~~~~~~~~~~ 3. Fit models per target ~~~~~~~~~~ #
  cli::cli_h2("Step 2: Fitting models per target")
  #  3.1. Combine all regulatory matrices into a single matrix (total regulators x cells)
  combinedRegData <- do.call(rbind, unname(regulatoryData))

  # 3.2. Cell alignment check
  # Models are fitted by position, so cells must be aligned
  if (!identical(colnames(combinedRegData), colnames(targetData)) ||
      !identical(rownames(metadata), colnames(targetData))) {
    cli::cli_abort("Internal error: cells are not aligned across target data, regulators and metadata.")
  }

  # 3.3 Build lookup list: target -> regulators
  regulatorsByTarget <- preparePriorDict(priorAssociations)

  # 3.4. Define target feature
  target_names <- rownames(targetData)
  cli::cli_alert_info("Fitting models for {length(target_names)} target feature(s).")

  # 3.5. Model fitting
  results <- pbapply::pblapply(target_names, function(target) {

    # Retrieve pre-filtered regulators for the current target
    target_regs <- regulatorsByTarget[[target]]

    # Fit individual target model
    fitTargetModel(
      target         = target,
      targetData     = targetData,
      regData        = combinedRegData,
      metadata       = metadata,
      regulators     = target_regs,
      fixedEffects   = fixedEffects,
      randomEffects  = randomEffects,
      interaction    = interaction,
      targetOmicType = targetOmicType,
      REML           = REML,
      verbose        = verbose,
      asEmmGrid      = TRUE
    )
  })
  names(results) <- target_names

  # ~~~~~~~~~~ 4. Summary ~~~~~~~~~~ #
  fitSummary <- do.call(rbind, lapply(results, function(r) {
    data.frame(
      target       = r$target,
      n_regulators = length(r$regulators),
      converged    = r$converged,
      is_singular  = r$is_singular,
      has_warnings = !is.null(r$warnings),
      error        = if (is.null(r$error)) NA_character_ else r$error,
      stringsAsFactors = FALSE
    )
  }))

  rownames(fitSummary) <- NULL

  n_failed <- sum(!is.na(fitSummary$error))
  n_nonconv <- sum(!fitSummary$converged & is.na(fitSummary$error))
  n_sing <- sum(fitSummary$is_singular %in% TRUE)
  cli::cli_alert_success("Fitted {nrow(fitSummary) - n_failed} of {nrow(fitSummary)} model(s).")
  if (n_failed > 0)  cli::cli_alert_warning("{n_failed} model(s) failed (see {.code GlobalSummary$fitSummary}).")
  if (n_nonconv > 0) cli::cli_alert_warning("{n_nonconv} model(s) did not converge.")
  if (n_sing > 0)    cli::cli_alert_info("{n_sing} model(s) have a singular fit.")

  # Create scMoreFit class
  Out <- list(
    Models = results,
    GlobalSummary = list(
      args = list(
        fixedEffects  = fixedEffects,
        randomEffects = randomEffects,
        twoSteps      = twoSteps,
        residualized  = residualized,
        interaction   = interaction,
        seed          = seed
      ),
      omicType          = omicType,
      nameMap           = nameMap,
      priorAssociations = priorAssociations,
      excludedFeatures  = excludedFeatures,
      fitSummary        = fitSummary,
      cells             = colnames(targetData),
      residuals         = if (keepResiduals && residualized) regulatoryData else NULL
    )
  )
  class(Out) <- c("scMoreFit", "list")

  return(Out)
}
