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

  # 1. Standardize column names to match existing annotations
  # This fixes the "Dot vs Underscore" issue (e.g., compound.name vs compound_name)
  common_cols <- intersect(names(ms2query_data), names(existing_annotations))

  # 2. Handle the Compound Name prefixing safely
  # Find which 'name' column exists in MS2Query to avoid creating a second one
  name_col <- if ("compound.name" %in% names(ms2query_data)) "compound.name" else "compound_name"

  if (name_col %in% names(ms2query_data)) {
    ms2query_data[[name_col]] <- paste0("Analogue of ", ms2query_data[[name_col]])
  }

  # 3. Ensure critical joining columns (like feature.ID) are the same type (Character)
  # This prevents bind_rows from failing due to numeric vs character mismatches
  if ("feature.ID" %in% names(ms2query_data)) {
    ms2query_data$feature.ID <- as.numeric(ms2query_data$feature.ID)
  }
  if ("feature.ID" %in% names(existing_annotations)) {
    existing_annotations$feature.ID <- as.numeric(existing_annotations$feature.ID)
  }

  # 4. Handle mz.diff.ppm type safety
  if ("mz.diff.ppm" %in% names(existing_annotations)) {
    existing_annotations$mz.diff.ppm <- as.numeric(as.character(existing_annotations$mz.diff.ppm))
  }

  # 5. Simple Append
  # bind_rows handles the missing columns automatically (much safer than the manual loop)
  updated_annotations <- dplyr::bind_rows(existing_annotations, ms2query_data)

  # Return the dataframe directly (simpler than a list)
  return(updated_annotations)
}
