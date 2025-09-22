#' @title Merge and Append Data Frames
#' @description Merges a new data frame of annotations with an existing one, handling column discrepancies.
#'
#' @param new_data A data frame containing new annotations to be appended.
#' @param existing_annotations A data frame of existing annotations.
#'
#' @return The combined data frame with new annotations appended.
#' @importFrom dplyr %>% bind_rows
#' @export
merge_and_append_data <- function(new_data, existing_annotations) {
  # Identify columns present in existing_annotations but not in new_data
  missing_cols_in_new <- setdiff(colnames(existing_annotations), colnames(new_data))
  for (col in missing_cols_in_new) {
    new_data[[col]] <- NA
  }

  # Identify columns present in new_data but not in existing_annotations
  missing_cols_in_existing <- setdiff(colnames(new_data), colnames(existing_annotations))
  for (col in missing_cols_in_existing) {
    existing_annotations[[col]] <- NA
  }

  # Append the new data to the existing annotations
  combined_annotations <- existing_annotations %>%
    dplyr::bind_rows(new_data)

  return(combined_annotations)
}
