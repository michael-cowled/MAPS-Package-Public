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

  # 1. Retrieve the Source/Parent MZ

  # Create a temporary lookup for M/Z from the full data
  mz_lookup <- full.annotation.data %>%
    dplyr::select(feature.ID, mz) %>%
    dplyr::rename(parent_mz = mz, Propagated.Feature.ID = feature.ID) %>%
    # FIX 1: Convert feature ID to character for consistent joining
    dplyr::mutate(Propagated.Feature.ID = as.character(Propagated.Feature.ID))

  # Join parent M/Z into propagated_df
  propagated_df_w_mass <- propagated_df %>%
    dplyr::left_join(mz_lookup, by = "Propagated.Feature.ID")


  # 2. Main Join and Calculation
  propagated_data <- full.annotation.data %>%
    dplyr::left_join(propagated_df_w_mass, by = "feature.ID") %>%
    dplyr::mutate(
      # FIX 2: Ensure mz and parent_mz are numeric before calculation
      mz = as.double(mz),
      parent_mz = as.double(parent_mz),

      confidence.level = as.character(confidence.level),
      propagation_mask = is.na(compound.name) & !is.na(Probable.Analogue.Of),

      # --- CALCULATION BLOCK ---
      # Calculate Delta: Target (full.mz) - Source (parent_mz)
      mz_delta = ifelse(propagation_mask, mz - parent_mz, NA),

      # Identify Modification
      mod_name = purrr::map_chr(mz_delta, function(x) {
        if (is.na(x)) return(NA_character_)

        # Calculate absolute mass differences against the modification database
        diffs <- abs(mod_db$Mass.Change - x)

        # Check tolerance (Absolute)
        match_idx <- which(diffs <= abs_tol)

        if (length(match_idx) > 0) {
          # Return the modification name with the smallest error
          best_match <- match_idx[which.min(diffs[match_idx])]
          return(mod_db$Modification[best_match])
        } else {
          # If no specific modification is matched within tolerance
          return("Probable")
        }
      }),

      # --- TEXT REPLACEMENT BLOCK ---
      compound.name = ifelse(
        propagation_mask,
        paste0(mod_name, " analogue of: ", Probable.Analogue.Of),
        compound.name
      ),

      confidence.level = ifelse(
        propagation_mask,
        "3", # Level 3 for propagated
        confidence.level
      ),

      NPC.pathway = as.character(NPC.pathway), # Coerce to character for consistent updating

      NPC.superclass = ifelse(
        propagation_mask & !is.na(Propagated.Annotation.Class),
        Propagated.Annotation.Class,
        NPC.superclass
      ),

      annotation.type = ifelse(
        propagation_mask,
        "GNPS Propagation",
        annotation.type
      )
    ) %>%
    # Clean up helper columns. Note: Propagated.Feature.ID is no longer needed.
    dplyr::select(-propagation_mask, -mz_delta, -mod_name, -parent_mz)

  return(propagated_data)
}
