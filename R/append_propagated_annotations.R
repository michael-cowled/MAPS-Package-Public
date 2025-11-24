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

  message("full.annotations columns: ", paste(names(full.annotation.data), collapse=", "))
  message("propagated columns: ", paste(names(propagated_df), collapse=", "))

  # 1. Retrieve the Source/Parent MZ
  # We look up the 'mz' of the 'Propagated.Feature.ID' (the parent/source)

  # Create a temporary lookup for M/Z from the full data
  mz_lookup <- full.annotation.data %>%
    dplyr::select(feature.ID, mz) %>%
    # Rename the feature.ID from the full data to match the PARENT ID column
    # in propagated_df, which is 'Propagated.Feature.ID'
    dplyr::rename(parent_mz = mz, Propagated.Feature.ID = feature.ID)

  # Join parent M/Z into propagated_df
  # Join on 'Propagated.Feature.ID' to link the parent's mass
  propagated_df_w_mass <- propagated_df %>%
    dplyr::left_join(mz_lookup, by = "Propagated.Feature.ID")

  # Ensure the necessary parent_mz column was successfully added
  if (!"parent_mz" %in% names(propagated_df_w_mass)) {
    stop("Failed to join 'parent_mz' to propagated_df. Check column names and data.")
  }

  # 2. Main Join and Calculation
  propagated_data <- full.annotation.data %>%
    # Main join on 'feature.ID' (the target feature)
    dplyr::left_join(propagated_df_w_mass, by = "feature.ID") %>%
    dplyr::mutate(
      confidence.level = as.character(confidence.level),
      # Mask for features that need propagation: unannotated (NA compound.name)
      # but have a probable analogue (not NA Probable.Analogue.Of)
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
    dplyr::select(-propagation_mask, -mz_delta, -mod_name, -parent_mz, -Propagated.Feature.ID)

  return(propagated_data)
}
