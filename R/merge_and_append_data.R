#' @title Merge and Append Data Frames
#' @description Merges a new data frame of annotations with an existing one, handling column discrepancies, ensuring consistent data types, and specifically standardizing 'mz.diff.ppm' to numeric.
#'
#' @param new_data A data frame containing new annotations to be appended.
#' @param existing_annotations A data frame of existing annotations.
#'
#' @return The combined data frame with new annotations appended.
#' @importFrom dplyr %>% bind_rows mutate
#' @export
merge_and_append_data <- function(new_data, existing_annotations) {

  # 1. Handle columns present in existing_annotations but not in new_data
  missing_cols_in_new <- setdiff(colnames(existing_annotations), colnames(new_data))
  for (col in missing_cols_in_new) {
    # Get the data type from the existing data frame
    col_type <- typeof(existing_annotations[[col]])

    # Initialize the new column in new_data with the correct type
    if (col_type == "character") {
      new_data[[col]] <- as.character(NA)
    } else if (col_type %in% c("integer", "double", "numeric")) {
      new_data[[col]] <- as.numeric(NA)
    } else if (col_type == "logical") {
      new_data[[col]] <- NA
    } else {
      # Fallback for complex types or unknown types
      new_data[[col]] <- NA
    }
  }

  # 2. Handle columns present in new_data but not in existing_annotations
  missing_cols_in_existing <- setdiff(colnames(new_data), colnames(existing_annotations))
  for (col in missing_cols_in_existing) {
    # Get the data type from the new data frame
    col_type <- typeof(new_data[[col]])

    # Initialize the new column in existing_annotations with the correct type
    if (col_type == "character") {
      existing_annotations[[col]] <- as.character(NA)
    } else if (col_type %in% c("integer", "double", "numeric")) {
      existing_annotations[[col]] <- as.numeric(NA)
    } else if (col_type == "logical") {
      existing_annotations[[col]] <- NA
    } else {
      # Fallback for complex types or unknown types
      existing_annotations[[col]] <- NA
    }
  }

  # --- NEW STEP ---
  # 3. Explicitly standardize the 'mz.diff.ppm' column to numeric in both data frames
  # This prevents the <character> and <double> type mismatch error during bind_rows.
  # Note: If 'mz.diff.ppm' contains non-numeric strings (like "N/A"), they will become NA,
  # which is the likely source of the warning you saw previously.
  if ("mz.diff.ppm" %in% colnames(existing_annotations)) {
    existing_annotations <- existing_annotations %>%
      dplyr::mutate(mz.diff.ppm = as.numeric(mz.diff.ppm))
  }
  if ("mz.diff.ppm" %in% colnames(new_data)) {
    new_data <- new_data %>%
      dplyr::mutate(mz.diff.ppm = as.numeric(mz.diff.ppm))
  }

  # 4. Append the new data to the existing annotations
  combined_annotations <- existing_annotations %>%
    dplyr::bind_rows(new_data)

  return(combined_annotations)
}
