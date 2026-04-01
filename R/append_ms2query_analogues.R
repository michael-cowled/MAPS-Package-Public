#' @title Format and Append MS2Query Analogue Annotations
#' @description Processes MS2Query Level 3 analogue data. Standardises base annotations offline, maps the existing mz.diff against a modification database, unconditionally formats the compound names, and prepares the data for merging.
#'
#' @param ms2query_data A data frame containing MS2Query Level 3 annotations.
#' @param existing_annotations A data frame of existing annotations.
#' @param mod_db The modification dataframe.
#' @param abs_tol The absolute mass tolerance in Daltons (default 0.01).
#' @param cid_cache_df A data frame used for caching CID lookups.
#' @param lipids.file Path to a lipids file for standardization.
#' @param cid_database_path Path to the PubChem CID database.
#' @param standardise_annotation A function to standardize annotations.
#' @param standardisation Logical; whether to run the standardisation passes.
#' @param cache.location Path to save the updated cache.
#' @param enable_local_db Logical; whether to use the local SQLite DB for property lookups.
#'
#' @return A list containing the updated annotations data frame and the updated CID cache.
#' @importFrom dplyr %>% mutate select
#' @importFrom purrr map_chr
#' @export
append_ms2query_analogues <- function(
    ms2query_data,
    existing_annotations,
    mod_db,
    abs_tol = 0.01,
    cid_cache_df,
    lipids.file,
    cid_database_path,
    standardise_annotation,
    standardisation,
    cache.location,
    enable_local_db = TRUE
) {

  name_col <- if ("compound.name" %in% names(ms2query_data)) "compound.name" else "compound_name"

  if (name_col %in% names(ms2query_data)) {

    # 1. --- Ensure SMILES column exists for standardisation ---
    if (!"smiles" %in% names(ms2query_data)) {
      ms2query_data$smiles <- NA_character_
    }

    # 2. --- STANDARDISATION (Base Compound, Offline Mode) ---
    # We do this FIRST so the DB/Cache lookups match the clean compound name!
    std_result <- standardise_annotation(
      data = ms2query_data,
      name_col = name_col,
      smiles_col = "smiles",
      cid_cache_df = cid_cache_df,
      lipids.file = lipids.file,
      cid_database_path = cid_database_path,
      standardisation = standardisation,
      cache.location = cache.location,
      enable_local_db = enable_local_db,
      enable_api = FALSE  # STRICTLY DISABLED
    )

    ms2query_data <- std_result$data
    updated_cid_cache_df <- std_result$cache

    # 3. --- MODIFICATION MAPPING ---
    # Set the baseline assumption
    ms2query_data$mod_prefix <- "Probable"

    # Since process_ms2query_data keeps 'mz.diff' for Level 3, we use it directly
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

    # 4. --- ANALOGUE FORMATTING ---
    # UNCONDITIONALLY apply the naming structure to the now-standardised name
    ms2query_data[[name_col]] <- paste0(ms2query_data$mod_prefix, " analogue of: ", ms2query_data[[name_col]])

    # Clean up the temporary prefix column
    ms2query_data <- ms2query_data %>% dplyr::select(-mod_prefix)
  } else {
    # If there is no name column to begin with, standardisation skips, pass cache through
    updated_cid_cache_df <- cid_cache_df
  }

  # Ensure critical joining columns are Character/Numeric matched
  if ("feature.ID" %in% names(ms2query_data)) {
    ms2query_data$feature.ID <- as.numeric(ms2query_data$feature.ID)
  }
  if ("feature.ID" %in% names(existing_annotations)) {
    existing_annotations$feature.ID <- as.numeric(existing_annotations$feature.ID)
  }

  # Return both the formatted data and the updated cache
  return(list(annotations = ms2query_data, cache = updated_cid_cache_df))
}
