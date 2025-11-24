#' Identify Redundant Features Based on SMILES and Retention Time
#'
#' Flags potentially redundant entries in a dataset based on a shared SMILES value and similar retention time (RT).
#' Groups features with identical SMILES whose RT values are within a specified tolerance.
#' For each RT group, all but the first occurrence are marked as redundant.
#'
#' @param dataset A data frame containing the data to process.
#' @param column_to_check Column name (as string) that contains the identifier for grouping (e.g., SMILES).
#' @param rt_column Column name (as string) indicating retention time.
#' @param rt_tolerance Numeric value specifying the maximum allowable RT difference to consider features redundant.
#'
#' @return The original dataset with a new logical column \code{redundant}, indicating flagged redundant rows.
#' @export
#'
#' @examples
#' \dontrun{
#' cleaned <- redundancy_fixer(mydata, column_to_check = "Best.Annotation.Smiles", rt_column = "rt", rt_tolerance = 1)
#' }
redundancy_fixer <- function(dataset, column_to_check = "Best.Annotation.Smiles", rt_column = "rt", rt_tolerance = 1) {
  # --- Error Handling: Check column existence ---
  if (!column_to_check %in% colnames(dataset)) {
    stop(paste("Column '", column_to_check, "' not found in the dataset.", sep = ""))
  }
  if (!rt_column %in% colnames(dataset)) {
    stop(paste("Column '", rt_column, "' not found in the dataset.", sep = ""))
  }

  # --- Ensure rt_column is numeric ---
  if (!is.numeric(dataset[[rt_column]])) {
    warning(paste("Column '", rt_column, "' is not numeric. Attempting to convert.", sep = ""))
    dataset[[rt_column]] <- as.numeric(as.character(dataset[[rt_column]]))
    if (any(is.na(dataset[[rt_column]]))) {
      stop(paste("Some values in ", rt_column, " could not be converted to numeric. Please inspect data.", sep = ""))
    }
  }

  # Initialize redundancy column
  dataset$redundant <- FALSE

  unique_values <- unique(dataset[[column_to_check]])

  for (value in unique_values) {
    temp_df <- dataset %>% dplyr::filter(!!rlang::sym(column_to_check) == value)

    if (nrow(temp_df) > 1) {
      rt_values <- sort(temp_df[[rt_column]])
      rt_values <- rt_values[!is.na(rt_values)]

      grouped_rts <- list()
      used_indices <- rep(FALSE, length(rt_values))

      for (i in seq_along(rt_values)) {
        if (!used_indices[i]) {
          group <- rt_values[i]
          used_indices[i] <- TRUE
          for (j in (i + 1):length(rt_values)) {
            if (!is.na(rt_values[j]) && abs(rt_values[j] - rt_values[i]) <= rt_tolerance) {
              group <- c(group, rt_values[j])
              used_indices[j] <- TRUE
            }
          }
          grouped_rts <- c(grouped_rts, list(sort(group)))
        }
      }

      for (group in grouped_rts) {
        group_df <- temp_df %>% dplyr::filter(!!rlang::sym(rt_column) %in% group)
        if (nrow(group_df) > 1) {
          # Flag all but the first as redundant
          dataset$redundant[which(dataset[[rt_column]] %in% group_df[[rt_column]][-1] &
                                    dataset[[column_to_check]] == value)] <- TRUE
        }
      }
    }
  }

  return(dataset)
}
