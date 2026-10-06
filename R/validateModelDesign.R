#' Validate model design and metadata integrity
#'
#' Validates that the fixed and random effect variables exist in the metadata,
#' contain no missing values, have at least two levels, and have syntactically
#' valid names for R formulas. Character and logical variables are converted to
#' factors, random effects are always converted to factors, and a summary of the
#' detected variable types is printed so the user can confirm them.
#'
#' @details
#' Numeric fixed effects are kept numeric (they are modelled as continuous
#' covariates). If a numeric variable has few distinct values a warning is
#' issued, because categorical variables coded as numbers (for example
#' \code{0/1}) must be converted to factors by the user. The function also warns
#' when random effects have very few levels (variance components are poorly
#' estimated) and when some combinations of fixed-effect levels are empty
#' (their contrasts cannot be estimated).
#'
#' @param metadata Data.frame with cell metadata. Row names must be cell IDs.
#' @param fixedEffects Character vector of \code{metadata} columns used as
#'  fixed tearm.
#' @param randomEffects Character vector of \code{metadata} columns used
#'  as random intercepts.
#' @param formula Optional custom formula (reserved for a future release; must
#'   be \code{NULL}).
#'
#' @return A list with:
#'   \item{metadata}{Metadata with sanitized column names and coerced factors.}
#'   \item{fixedEffects}{Sanitized vector of fixed effect variable names.}
#'   \item{randomEffects}{Sanitized vector of random effect variable names.}
#'
#' @importFrom cli cli_abort cli_alert_info cli_alert_warning cli_alert_success cli_ul cli_li cli_end
#' @keywords internal
validateModelDesign <- function(metadata,
                                fixedEffects,
                                randomEffects,
                                formula) {

  # ~~~~~~~~~~ 1. Advanced Custom Formula Check ~~~~~~~~~~ #
  if (!is.null(formula)) {
    cli::cli_abort(c(
      "Custom formula specification via {.arg formula} is currently under active development and will be enabled in a future release.",
      "i" = "Please specify your design using {.arg fixedEffects} and {.arg randomEffects}."
    ))
  }

  if (!is.data.frame(metadata)) {
    cli::cli_abort("{.arg metadata} must be a data.frame.")
  }
  if (!is.null(fixedEffects) && !is.character(fixedEffects)) {
    cli::cli_abort("{.arg fixedEffects} must be a character vector of {.arg metadata} column names.")
  }
  if (!is.null(randomEffects) && !is.character(randomEffects)) {
    cli::cli_abort("{.arg randomEffects} must be a character vector of {.arg metadata} column names.")
  }

  fixedEffects  <- unique(fixedEffects)
  randomEffects <- unique(randomEffects)

  if (length(randomEffects) == 0) {
    cli::cli_alert_info("No random effects specified. Linear models (lm) will be fitted.")
  }

  # Group effects
  all_design_vars <- unique(c(fixedEffects, randomEffects))

  if (length(all_design_vars) == 0) {
    cli::cli_alert_warning("No fixed or random effects specified. Fitting regulation-only models.")
    return(list(metadata = metadata,
                fixedEffects = character(0),
                randomEffects = character(0)))
  }

  overlap <- intersect(fixedEffects, randomEffects)
  if (length(overlap) > 0) {
    cli::cli_abort("Variable(s) {.val {overlap}} cannot be both fixed and random effects.")
  }

  # ~~~~~~~~~~ 2. Verify existence in metadata ~~~~~~~~~~ #
  missing_cols <- setdiff(all_design_vars, colnames(metadata))
  if (length(missing_cols) > 0) {
    cli::cli_abort("The following variable(s) were not found in {.arg metadata}: {.val {missing_cols}}.")
  }

  # ~~~~~~~~~~ 3. Check for missing values (NA) ~~~~~~~~~~ #
  na_cols <- all_design_vars[sapply(all_design_vars, function(col) any(is.na(metadata[[col]])))]
  if (length(na_cols) > 0) {
    cli::cli_abort("Missing values (NA) found in {.arg metadata} column(s): {.val {na_cols}}. Please filter or impute NAs before fitting.")
  }

  # ~~~~~~~~~~ 4. Check for single-level variables ~~~~~~~~~~ #
  single_level_cols <- all_design_vars[sapply(all_design_vars, function(col) length(unique(metadata[[col]])) < 2)]
  if (length(single_level_cols) > 0) {
    cli::cli_abort("The following variable(s) have fewer than 2 unique levels in {.arg metadata}: {.val {single_level_cols}}.")
  }

  # ~~~~~~~~~~ 5. Syntactic name sanitization ~~~~~~~~~~ #
  clean_vars <- make.names(all_design_vars)
  non_syntactic <- all_design_vars[clean_vars != all_design_vars]

  if (length(non_syntactic) > 0) {
    old_colnames <- colnames(metadata)
    new_colnames <- make.names(old_colnames)
    if (anyDuplicated(new_colnames)) {
      cli::cli_abort("Sanitizing {.arg metadata} column names produced duplicates. Please rename the columns manually.")
    }
    name_map <- stats::setNames(new_colnames, old_colnames)

    fixedEffects  <- unname(name_map[fixedEffects])
    randomEffects <- unname(name_map[randomEffects])
    colnames(metadata) <- new_colnames

    cli::cli_alert_warning(
      "Converted design variable names to syntactically valid R names for formula compatibility."
    )
  }

  # ~~~~~~~~~~ 6. Variable types ~~~~~~~~~~ #
  for (fv in fixedEffects) {
    if (is.character(metadata[[fv]]) || is.logical(metadata[[fv]])) {
      metadata[[fv]] <- as.factor(metadata[[fv]])
    }
    if (is.numeric(metadata[[fv]]) && length(unique(metadata[[fv]])) <= 10) {
      cli::cli_alert_warning(
        "Fixed effect {.field {fv}} is numeric with only {length(unique(metadata[[fv]]))} distinct values and will be modelled as a continuous covariate. Convert it to a factor if it is categorical."
      )
    }
  }
  for (rv in randomEffects) {
    if (!is.factor(metadata[[rv]])) {
      metadata[[rv]] <- as.factor(metadata[[rv]])
    }
    n_lev <- nlevels(droplevels(metadata[[rv]]))
    if (n_lev >= nrow(metadata)) {
      cli::cli_abort("Random effect {.field {rv}} has one level per cell and cannot be estimated.")
    }
    if (n_lev < 5) {
      cli::cli_alert_warning(
        "Random effect {.field {rv}} has only {n_lev} levels: its variance will be poorly estimated."
      )
    }
  }

  # Empty combinations of fixed-effect levels
  fixed_factors <- fixedEffects[vapply(fixedEffects, function(v) is.factor(metadata[[v]]), logical(1))]
  if (length(fixed_factors) > 1) {
    tab <- table(metadata[, fixed_factors, drop = FALSE])
    if (any(tab == 0)) {
      cli::cli_alert_warning(
        "Some combinations of {.field {fixed_factors}} have no cells: contrasts involving them cannot be estimated."
      )
    }
  }

  # Report the detected types so the user can confirm them
  cli::cli_text("Detected design variables:")
  cli::cli_ul()
  for (v in fixedEffects) {
    x <- metadata[[v]]
    desc <- if (is.factor(x)) paste0("factor, ", nlevels(droplevels(x)), " levels") else paste0(class(x)[1], " (continuous)")
    cli::cli_li("{.field {v}} [fixed]: {desc}")
  }
  for (v in randomEffects) {
    cli::cli_li("{.field {v}} [random intercept]: factor, {nlevels(droplevels(metadata[[v]]))} levels")
  }
  cli::cli_end()

  cli::cli_alert_success("Model design and metadata successfully validated.")

  list(
    metadata      = metadata,
    fixedEffects  = fixedEffects,
    randomEffects = randomEffects
  )
}}
