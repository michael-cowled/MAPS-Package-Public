#' Identify Redundant Features Based on 3-Tier Hierarchy and Retention Time
#'
#' Flags potentially redundant entries in a dataset based on:
#' 1. CID
#' 2. Identical compound.name
#' 3. A combination of compound.name AND smiles
#' Explicitly sorts by confidence to protect Level 1/2 annotations.
#'
#' @param dataset A data frame containing the data to process.
#' @param rt_column Column name (as string) indicating retention time.
#' @param rt_tolerance Numeric value specifying the maximum allowable RT difference.
#'
#' @return The original dataset with a new logical column \code{redundant}.
#' @export
redundancy_fixer <- function(dataset, rt_column = "rt", rt_tolerance = 1) {

  # --- THE VITAL FIX: Sort by Confidence ---
  # Guarantee Level 1 > Level 2 > Level 3 order so the best hits become Keepers.
  if ("confidence.level" %in% colnames(dataset)) {
    if ("confidence.score" %in% colnames(dataset)) {
      dataset <- dataset %>%
        dplyr::mutate(confidence.score = suppressWarnings(as.numeric(confidence.score))) %>%
        dplyr::arrange(confidence.level, dplyr::desc(confidence.score))
    } else {
      dataset <- dataset %>%
        dplyr::arrange(confidence.level)
    }
  }

  # Ensure RT is numeric
  if (!is.numeric(dataset[[rt_column]])) {
    dataset[[rt_column]] <- suppressWarnings(as.numeric(dataset[[rt_column]]))
  }

  # Initialize all to FALSE
  dataset$redundant <- FALSE

  # --- HELPER FUNCTION ---
  apply_redundancy_pass <- function(data, target_vector) {
    # Get unique, valid targets (ignore NAs, empty strings, and generic placeholders)
    unique_targets <- unique(target_vector)
    unique_targets <- unique_targets[!is.na(unique_targets) &
                                       unique_targets != "" &
                                       toupper(unique_targets) != "N/A" &
                                       unique_targets != "0" &
                                       toupper(unique_targets) != "UNKNOWN"]

    for (target in unique_targets) {
      idx <- which(target_vector == target & !data$redundant)

      if (length(idx) > 1) {
        rt_values <- data[[rt_column]][idx]
        used_indices <- rep(FALSE, length(rt_values))

        for (i in seq_along(rt_values)) {
          if (!used_indices[i] && !is.na(rt_values[i])) {
            used_indices[i] <- TRUE # This is the "Keeper" (now guaranteed to be highest confidence)

            for (j in (i + 1):length(rt_values)) {
              if (!used_indices[j] && !is.na(rt_values[j])) {
                if (abs(rt_values[j] - rt_values[i]) <= rt_tolerance) {
                  used_indices[j] <- TRUE
                  data$redundant[idx[j]] <- TRUE # Flag lower confidence hit as redundant
                }
              }
            }
          }
        }
      }
    }
    return(data)
  }

  # --- TIER 1: Check by CID ---
  dataset <- apply_redundancy_pass(dataset, dataset$CID)

  # --- TIER 2: Check by compound.name ---
  dataset <- apply_redundancy_pass(dataset, dataset$compound.name)

  # --- TIER 3: Check by combined compound.name AND smiles ---
  combined_key <- paste(dataset$compound.name, dataset$smiles, sep = "_")

  invalid_name <- is.na(dataset$compound.name) | dataset$compound.name == "" | toupper(dataset$compound.name) == "N/A"
  invalid_smiles <- is.na(dataset$smiles) | dataset$smiles == "" | toupper(dataset$smiles) == "N/A"

  combined_key[invalid_name | invalid_smiles] <- NA

  dataset <- apply_redundancy_pass(dataset, combined_key)

  return(dataset)
}
