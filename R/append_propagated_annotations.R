#' @title Append Propagated Annotations to a Data Frame
#' @description Merges propagated annotations into the main annotation data frame, calculates mass differences to identify specific modifications, and updates compound names.
#'
#' @param full.annotation.data The original full annotation data frame. Must contain 'feature.ID' and 'mz'.
#' @param propagated_df The data frame of newly propagated annotations. Must contain 'feature.ID' (the target) and 'Propagated.Feature.ID' (the source/parent feature ID).
#' @param mod_db The modification dataframe (defaults to the modification_db defined above).
#' @param ppm_tol The mass tolerance in ppm (default 10).
#' @param abs_tol The absolute mass tolerance in Daltons (default 0.01).
#'
#' @return A single data frame with the propagated annotations appended and specific modifications named.
#' @importFrom dplyr %>% left_join mutate select filter rename
#' @importFrom purrr map_chr
#' @export
append_propagated_annotations <- function(full.annotation.data,
                                          propagated_df,
                                          mod_db = modification_db,
                                          ppm_tol = 10,
                                          abs_tol = 0.01) {

  # --- DEFENSIVE CHECK: Early Exit with Schema Consistency ---
  if (is.null(propagated_df) || nrow(propagated_df) == 0 || ncol(propagated_df) == 0) {
    message("Notice: No propagation data found. Returning original data with empty placeholders.")

    # Ensure columns exist so downstream functions don't crash
    cols_to_add <- c("Propagated.Feature.ID", "Propagated.Annotation.Type", "Propagated.Annotation.Class", "Propagated.Annotation.Smiles")
    for (col in cols_to_add) {
      if (!(col %in% colnames(full.annotation.data))) full.annotation.data[[col]] <- NA
    }
    return(full.annotation.data)
  }

  # 1. Retrieve the Source/Parent MZ
  mz_lookup <- full.annotation.data %>%
    dplyr::select(feature.ID, mz) %>%
    dplyr::rename(parent_mz = mz, Propagated.Feature.ID = feature.ID) %>%
    dplyr::mutate(Propagated.Feature.ID = as.character(Propagated.Feature.ID))

  # Join parent M/Z into propagated_df
  propagated_df_w_mass <- propagated_df %>%
    dplyr::mutate(Propagated.Feature.ID = as.character(Propagated.Feature.ID)) %>%
    dplyr::left_join(mz_lookup, by = "Propagated.Feature.ID")

  # 2. Main Join and Calculation
  propagated_data <- full.annotation.data %>%
    dplyr::left_join(propagated_df_w_mass, by = "feature.ID") %>%
    dplyr::mutate(
      mz = as.double(mz),
      parent_mz = as.double(parent_mz),
      confidence.level = as.character(confidence.level),

      # Determine if a row is a candidate for propagation name update
      propagation_mask = is.na(compound.name) & !is.na(Probable.Analogue.Of),

      # --- CALCULATION BLOCK ---
      mz_delta = ifelse(propagation_mask, mz - parent_mz, NA),

      mod_name = purrr::map_chr(mz_delta, function(x) {
        if (is.na(x)) return(NA_character_)
        diffs <- abs(mod_db$Mass.Change - x)
        match_idx <- which(diffs <= abs_tol)
        if (length(match_idx) > 0) {
          best_match <- match_idx[which.min(diffs[match_idx])]
          return(mod_db$Modification[best_match])
        } else {
          return("Probable")
        }
      }),

      # --- TEXT REPLACEMENT BLOCK ---
      compound.name = ifelse(
        propagation_mask,
        paste0(mod_name, " analogue of: ", Probable.Analogue.Of),
        compound.name
      ),
      confidence.level = ifelse(propagation_mask, "3", confidence.level),
      NPC.pathway = as.character(NPC.pathway),
      NPC.superclass = ifelse(
        propagation_mask & !is.na(Propagated.Annotation.Class),
        Propagated.Annotation.Class,
        NPC.superclass
      ),
      annotation.type = ifelse(propagation_mask, "GNPS Propagation", annotation.type)
    ) %>%
    dplyr::select(-propagation_mask, -mz_delta, -mod_name, -parent_mz)

  return(propagated_data)
}
