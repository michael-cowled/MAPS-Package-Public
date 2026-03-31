#' @title Identify Redundant Features Based on 3-Tier Hierarchy and Retention Time#' Identify Redundant Features Based on 3-Tier Hierarchy and Retention Time
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

  # --- STEP 0: Immediate Type Standardization ---
  # This ensures this function plays nice with the rest of your pipeline
  if ("confidence.score" %in% colnames(dataset)) {
    dataset$confidence.score <- suppressWarnings(as.numeric(dataset$confidence.score))
  }
  if ("confidence.level" %in% colnames(dataset)) {
    dataset$confidence.level <- as.character(dataset$confidence.level)
  }
  if (rt_column %in% colnames(dataset)) {
    dataset[[rt_column]] <- suppressWarnings(as.numeric(dataset[[rt_column]]))
  }

  # --- STEP 1: Sort by Confidence ---
  # Level 1 > Level 2 > Level 3. Best hits are now at the top of their groups.
  if ("confidence.level" %in% colnames(dataset)) {
    if ("confidence.score" %in% colnames(dataset)) {
      dataset <- dataset %>%
        dplyr::arrange(confidence.level, dplyr::desc(confidence.score))
    } else {
      dataset <- dataset %>%
        dplyr::arrange(confidence.level)
    }
  }

  # Initialize redundant flag
  dataset$redundant <- FALSE

  # --- HELPER FUNCTION (Internal) ---
  apply_redundancy_pass <- function(data, target_vector) {
    # Filter out empty/invalid keys for comparison
    unique_targets <- unique(target_vector)
    unique_targets <- unique_targets[!is.na(unique_targets) &
                                       unique_targets != "" &
                                       toupper(unique_targets) != "N/A" &
                                       unique_targets != "0" &
                                       toupper(unique_targets) != "UNKNOWN"]

    for (target in unique_targets) {
      # Find all rows matching this ID that haven't been marked redundant yet
      idx <- which(target_vector == target & !data$redundant)

      if (length(idx) > 1) {
        rt_values <- data[[rt_column]][idx]
        used_indices <- rep(FALSE, length(rt_values))

        for (i in seq_along(rt_values)) {
          if (!used_indices[i] && !is.na(rt_values[i])) {
            used_indices[i] <- TRUE # The "Keeper" (highest confidence available)

            # Check all subsequent hits against this Keeper's RT
            for (j in (i + 1):length(rt_values)) {
              if (!used_indices[j] && !is.na(rt_values[j])) {
                if (abs(rt_values[j] - rt_values[i]) <= rt_tolerance) {
                  used_indices[j] <- TRUE
                  data$redundant[idx[j]] <- TRUE # Mark as redundant
                }
              }
            }
          }
        }
      }
    }
    return(data)
  }

  # --- TIER 1: CID ---
  if ("CID" %in% colnames(dataset)) {
    dataset <- apply_redundancy_pass(dataset, as.character(dataset$CID))
  }

  # --- TIER 2: compound.name ---
  if ("compound.name" %in% colnames(dataset)) {
    dataset <- apply_redundancy_pass(dataset, as.character(dataset$compound.name))
  }

  # --- TIER 3: Name + SMILES ---
  # We use a separator that is unlikely to be in the strings
  combined_key <- paste(dataset$compound.name, dataset$smiles, sep = "|||")

  # Invalidate keys where either part is missing or generic
  invalid_name <- is.na(dataset$compound.name) | dataset$compound.name == "" | toupper(dataset$compound.name) == "N/A"
  invalid_smiles <- is.na(dataset$smiles) | dataset$smiles == "" | toupper(dataset$smiles) == "N/A"
  combined_key[invalid_name | invalid_smiles] <- NA

  dataset <- apply_redundancy_pass(dataset, combined_key)

  return(dataset)
}
