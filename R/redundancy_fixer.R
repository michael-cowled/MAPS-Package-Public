#' Identify Redundant Features Based on 3-Tier Hierarchy and Retention Time
#'
#' Flags potentially redundant entries in a dataset based on:
#' 1. CID
#' 2. Identical compound.name
#' 3. A combination of compound.name AND smiles
#' It never checks 'smiles' exclusively.
#'
#' @param dataset A data frame containing the data to process.
#' @param rt_column Column name (as string) indicating retention time.
#' @param rt_tolerance Numeric value specifying the maximum allowable RT difference.
#'
#' @return The original dataset with a new logical column \code{redundant}, indicating flagged redundant rows.
#' @export
redundancy_fixer <- function(dataset, rt_column = "rt", rt_tolerance = 1) {

  # --- Error Handling: Check column existence ---
  required_cols <- c(rt_column, "CID", "compound.name", "smiles")
  missing_cols <- setdiff(required_cols, colnames(dataset))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in dataset:", paste(missing_cols, collapse = ", ")))
  }

  # Ensure RT is numeric
  if (!is.numeric(dataset[[rt_column]])) {
    dataset[[rt_column]] <- suppressWarnings(as.numeric(dataset[[rt_column]]))
    if (any(is.na(dataset[[rt_column]]))) {
      warning(paste("Some values in", rt_column, "were coerced to NA during numeric conversion."))
    }
  }

  # Initialize redundancy column (everything starts as FALSE)
  dataset$redundant <- FALSE

  # --- HELPER FUNCTION: Apply Redundancy Logic ---
  # This function takes a vector of target values and flags redundancies
  # for any features within the RT tolerance that share the target value.
  apply_redundancy_pass <- function(data, target_vector) {
    # Get unique, valid targets (ignore NAs and empty strings)
    unique_targets <- unique(target_vector)
    unique_targets <- unique_targets[!is.na(unique_targets) & unique_targets != "" & unique_targets != "N/A"]

    for (target in unique_targets) {
      # Find the row indices for features matching the target that ARE NOT YET REDUNDANT
      idx <- which(target_vector == target & !data$redundant)

      # If there is more than 1 feature matching this target, check RTs
      if (length(idx) > 1) {
        rt_values <- data[[rt_column]][idx]

        # Track which of these specific RTs have been processed/flagged
        used_indices <- rep(FALSE, length(rt_values))

        for (i in seq_along(rt_values)) {
          if (!used_indices[i] && !is.na(rt_values[i])) {
            used_indices[i] <- TRUE # This is the "Keeper" (First occurrence)

            # Check all subsequent matching features against the Keeper's RT
            for (j in (i + 1):length(rt_values)) {
              if (!used_indices[j] && !is.na(rt_values[j])) {

                # If within RT tolerance, flag the 'j' index as redundant!
                if (abs(rt_values[j] - rt_values[i]) <= rt_tolerance) {
                  used_indices[j] <- TRUE
                  data$redundant[idx[j]] <- TRUE # Mark in the main dataframe
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
  # Create a combined vector. If either is NA/Empty, we set the combined key to NA
  # so we don't accidentally match "NA_SMILES1" with "NA_SMILES2".
  combined_key <- paste(dataset$compound.name, dataset$smiles, sep = "_")

  invalid_name <- is.na(dataset$compound.name) | dataset$compound.name == "" | dataset$compound.name == "N/A"
  invalid_smiles <- is.na(dataset$smiles) | dataset$smiles == "" | dataset$smiles == "N/A"

  combined_key[invalid_name | invalid_smiles] <- NA

  dataset <- apply_redundancy_pass(dataset, combined_key)

  return(dataset)
}
