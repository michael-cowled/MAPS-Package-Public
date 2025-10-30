#' @title Append MS2Query Analogue Annotations
#' @description Processes and appends MS2Query Level 3 analogue data to existing annotations.
#'
#' @param ms2query_data A data frame containing MS2Query Level 3 annotations.
#' @param existing_annotations A data frame of existing annotations (lv1, lv2, lv3).
#'
#' @return The updated annotations data frame with MS2Query data appended.
#' @importFrom dplyr %>% filter select bind_rows
#' @export
append_ms2query_analogues <- function(ms2query_data, existing_annotations) {
  # Convert the mz.diff.ppm column in existing_annotations to numeric
  if ("mz.diff.ppm" %in% names(existing_annotations) && is.character(existing_annotations$mz.diff.ppm)) {
    existing_annotations$mz.diff.ppm <- as.numeric(existing_annotations$mz.diff.ppm)
  }
  print("check3")
  # Modify compound names to indicate they are analogues
  ms2query_data$compound.name <- paste0("Analogue of ", ms2query_data$compound.name)
  print("check4")
  # Filter for unique features not already annotated
  ##TO REMOVE
 # unique_in_ms2query <- setdiff(ms2query_data$feature.ID, existing_annotations$feature.ID)
 # ms2query_data <- ms2query_data %>%
 #   dplyr::filter(feature.ID %in% unique_in_ms2query) %>%
  #  dplyr::select(-mz.diff, -precursor_mz)
  print("check5")
  # Identify and add missing columns to ensure successful binding
  missing_cols <- colnames(existing_annotations)[!colnames(existing_annotations) %in% colnames(ms2query_data)]
  print("existing_annotations")
  print(head(existing_annotations))
  print("ms2query_data")
  print(head(ms2query_data))
  if(length(missing_cols) > 0) {
    for (col in missing_cols) {
      ms2query_data[[col]] <- NA
    }
  }
  print("check6")
  # Append the processed data to the existing annotations
  updated_annotations <- existing_annotations %>%
    dplyr::bind_rows(ms2query_data)
  print("check7")
  return(updated_annotations)
}
