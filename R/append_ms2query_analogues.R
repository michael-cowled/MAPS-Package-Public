#' @title Format and Append MS2Query Analogue Annotations
#' @description Processes MS2Query Level 3 analogue data, calculates mass differences to identify specific modifications, formats the compound names, and prepares the data for merging.
#'
#' @param ms2query_data A data frame containing MS2Query Level 3 annotations.
#' @param existing_annotations A data frame of existing annotations.
#' @param mod_db The modification dataframe.
#' @param abs_tol The absolute mass tolerance in Daltons (default 0.01).
#'
#' @return The updated MS2Query data frame, formatted and ready for merging.
#' @importFrom dplyr %>% mutate select
#' @importFrom purrr map_chr
#' @export
append_ms2query_analogues <- function(ms2query_data, existing_annotations, mod_db, abs_tol = 0.01) {

  # 1. Standardize column names to match existing annotations
  common_cols <- intersect(names(ms2query_data), names(existing_annotations))

  # Find the correct name column
  name_col <- if ("compound.name" %in% names(ms2query_data)) "compound.name" else "compound_name"

  # 2. Calculate mass deltas and apply the Modification DB
  if (name_col %in% names(ms2query_data) && "mz" %in% names(ms2query_data) && "analogue_mz" %in% names(ms2query_data)) {

    ms2query_data <- ms2query_data %>%
      dplyr::mutate(
        mz_delta = as.double(mz) - as.double(analogue_mz), # Remember to check 'analogue_mz' name!

        mod_name = purrr::map_chr(mz_delta, function(x) {
          if (is.na(x)) return("Probable")

          diffs <- abs(mod_db$Mass.Change - x)
          match_idx <- which(diffs <= abs_tol)

          if (length(match_idx) > 0) {
            best_match <- match_idx[which.min(diffs[match_idx])]
            return(mod_db$Modification[best_match])
          } else {
            return("Probable")
          }
        })
      )

    # Update the compound name using base R to dynamically target 'name_col'
    ms2query_data[[name_col]] <- paste0(ms2query_data$mod_name, " analogue of: ", ms2query_data[[name_col]])

    # Clean up calculation columns
    ms2query_data <- ms2query_data %>% dplyr::select(-mz_delta, -mod_name)
  }

  # 3. Ensure critical joining columns (like feature.ID) are Character/Numeric matched
  if ("feature.ID" %in% names(ms2query_data)) {
    ms2query_data$feature.ID <- as.numeric(ms2query_data$feature.ID)
  }
  if ("feature.ID" %in% names(existing_annotations)) {
    existing_annotations$feature.ID <- as.numeric(existing_annotations$feature.ID)
  }

  # Return the nicely formatted dataframe
  return(ms2query_data)
}
