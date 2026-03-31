#' @title Merge and Append Data Frames
#' @description Merges a new data frame of annotations with an existing one,
#' handling column discrepancies and ensuring consistent data types (e.g., character for IDs, numeric for scores).
#'
#' @param new_data A data frame containing new annotations to be appended.
#' @param existing_annotations A data frame of existing annotations.
#'
#' @return The combined data frame with new annotations appended.
#' @importFrom dplyr %>% bind_rows
#' @export
merge_and_append_data <- function(new_data, existing_annotations) {

  # 1. Handle columns present in existing_annotations but not in new_data
  missing_cols_in_new <- setdiff(colnames(existing_annotations), colnames(new_data))
  for (col in missing_cols_in_new) {
    col_type <- typeof(existing_annotations[[col]])

    if (col_type == "character") {
      new_data[[col]] <- as.character(NA)
    } else if (col_type %in% c("integer", "double", "numeric")) {
      new_data[[col]] <- as.numeric(NA)
    } else if (col_type == "logical") {
      new_data[[col]] <- NA
    } else {
      new_data[[col]] <- NA
    }
  }

  # 2. Handle columns present in new_data but not in existing_annotations
  missing_cols_in_existing <- setdiff(colnames(new_data), colnames(existing_annotations))
  for (col in missing_cols_in_existing) {
    col_type <- typeof(new_data[[col]])

    if (col_type == "character") {
      existing_annotations[[col]] <- as.character(NA)
    } else if (col_type %in% c("integer", "double", "numeric")) {
      existing_annotations[[col]] <- as.numeric(NA)
    } else if (col_type == "logical") {
      existing_annotations[[col]] <- NA
    } else {
      existing_annotations[[col]] <- NA
    }
  }

  # 3. Explicit Standardization (The Type-Mismatch Fix)

  # A. Force ID and Name columns to Character (prevents Integer vs Character errors)
  char_cols <- c("feature.ID", "compound.name", "smiles", "IUPAC", "Formula")
  for (col in char_cols) {
    if (col %in% colnames(new_data)) {
      new_data[[col]] <- as.character(new_data[[col]])
    }
    if (col %in% colnames(existing_annotations)) {
      existing_annotations[[col]] <- as.character(existing_annotations[[col]])
    }
  }

  # B. Force Score and Mass columns to Numeric (prevents Character vs Double errors)
  num_cols <- c("mz.diff.ppm", "confidence.score", "precursor_mz", "Monoisotopic.Mass", "rt")
  for (col in num_cols) {
    if (col %in% colnames(new_data)) {
      new_data[[col]] <- as.numeric(new_data[[col]])
    }
    if (col %in% colnames(existing_annotations)) {
      existing_annotations[[col]] <- as.numeric(existing_annotations[[col]])
    }
  }

  # 4. Final Safety Check: Match remaining common column types
  # If a column isn't in the lists above but still differs,
  # force new_data to match the existing schema.
  common_cols <- intersect(colnames(new_data), colnames(existing_annotations))
  for (col in common_cols) {
    if (class(new_data[[col]]) != class(existing_annotations[[col]])) {
      suppressWarnings(class(new_data[[col]]) <- class(existing_annotations[[col]]))
    }
  }

  # 5. Append the data frames safely
  combined_annotations <- dplyr::bind_rows(existing_annotations, new_data)

  return(combined_annotations)
}
