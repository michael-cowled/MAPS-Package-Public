#' @title Format and Append MS2Query Analogue Annotations
#' @description Processes MS2Query Level 3 analogue data, maps the existing mz.diff against a modification database, unconditionally formats the compound names, and prepares the data for merging.
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

  name_col <- if ("compound.name" %in% names(ms2query_data)) "compound.name" else "compound_name"

  if (name_col %in% names(ms2query_data)) {

    # Set the baseline assumption
    ms2query_data$mod_prefix <- "Probable"

    # Since process_ms2query_data keeps 'mz.diff' for Level 3, we use it directly!
    if ("mz.diff" %in% names(ms2query_data)) {

      ms2query_data <- ms2query_data %>%
        dplyr::mutate(
          mod_prefix = purrr::map_chr(mz.diff, function(x) {
            if (is.na(x)) return("Probable")

            # Compare the absolute magnitude of the shift to the DB
            diffs <- abs(abs(mod_db$Mass.Change) - abs(x))
            match_idx <- which(diffs <= abs_tol)

            if (length(match_idx) > 0) {
              best_match <- match_idx[which.min(diffs[match_idx])]
              return(mod_db$Modification[best_match])
            } else {
              return("Probable")
            }
          })
        )
    } else {
      warning("Column 'mz.diff' not found. Defaulting all to 'Probable'.")
    }

    # UNCONDITIONALLY apply the naming structure
    ms2query_data[[name_col]] <- paste0(ms2query_data$mod_prefix, " analogue of: ", ms2query_data[[name_col]])

    # Clean up the temporary prefix column
    ms2query_data <- ms2query_data %>% dplyr::select(-mod_prefix)

    # (Optional: If you want to drop mz.diff so it matches Level 2, you could add:
    # ms2query_data <- ms2query_data %>% dplyr::select(-mz.diff, -precursor_mz) )
  }

  # Ensure critical joining columns are Character/Numeric matched
  if ("feature.ID" %in% names(ms2query_data)) {
    ms2query_data$feature.ID <- as.numeric(ms2query_data$feature.ID)
  }
  if ("feature.ID" %in% names(existing_annotations)) {
    existing_annotations$feature.ID <- as.numeric(existing_annotations$feature.ID)
  }

  return(ms2query_data)
}
