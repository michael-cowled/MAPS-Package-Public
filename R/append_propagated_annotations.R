#' @title Append Propagated Annotations to a Data Frame
#' @description Merges propagated annotations into the main annotation data frame,
#' calculates mass differences to identify specific modifications, and updates compound names.
#'
#' @export
append_propagated_annotations <- function(full.annotation.data,
                                          propagated_df,
                                          mod_db = modification_db,
                                          ppm_tol = 10,
                                          abs_tol = 0.01) {

  # --- FIX 1: Ensure feature.ID and confidence.score are consistent types before starting ---
  full.annotation.data <- full.annotation.data %>%
    dplyr::mutate(
      feature.ID = as.character(feature.ID),
      # Force confidence.score to numeric (turning "N/A" into real NA)
      confidence.score = suppressWarnings(as.numeric(confidence.score))
    )

  propagated_df <- propagated_df %>%
    dplyr::mutate(
      feature.ID = as.character(feature.ID),
      Propagated.Feature.ID = as.character(Propagated.Feature.ID),
      # Force Propagated score to numeric
      Propagated.Confidence.Score = suppressWarnings(as.numeric(Propagated.Confidence.Score))
    )

  # 1. Retrieve the Source/Parent MZ
  mz_lookup <- full.annotation.data %>%
    dplyr::select(feature.ID, mz) %>%
    dplyr::rename(parent_mz = mz, Propagated.Feature.ID = feature.ID)

  # Join parent M/Z into propagated_df
  propagated_df_w_mass <- propagated_df %>%
    dplyr::left_join(mz_lookup, by = "Propagated.Feature.ID")

  # 2. Main Join and Calculation
  propagated_data <- full.annotation.data %>%
    dplyr::left_join(propagated_df_w_mass, by = "feature.ID") %>%
    dplyr::mutate(
      # Clean up any lingering "N/A" strings in character columns
      dplyr::across(tidyselect::where(is.character), ~ dplyr::na_if(.x, "N/A")),

      mz = as.double(mz),
      parent_mz = as.double(parent_mz),
      confidence.level = as.character(confidence.level),

      propagation_mask = (is.na(compound.name) | confidence.level == "3") & !is.na(Probable.Analogue.Of),
      mz_delta = ifelse(propagation_mask, mz - parent_mz, NA),

      # Identify Modification
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

      smiles = ifelse(
        propagation_mask & !is.na(Propagated.Smiles),
        as.character(Propagated.Smiles),
        as.character(smiles)
      ),

      confidence.level = ifelse(
        propagation_mask,
        "3",
        as.character(confidence.level)
      ),

      # --- FIX 2: Ensure ifelse outputs consistent numeric types ---
      confidence.score = ifelse(
        propagation_mask & !is.na(Propagated.Confidence.Score),
        as.numeric(Propagated.Confidence.Score),
        as.numeric(confidence.score)
      ),

      NPC.pathway = as.character(NPC.pathway),

      annotation.type = ifelse(
        propagation_mask,
        "GNPS Propagation",
        as.character(annotation.type)
      )
    ) %>%
    # Clean up helper columns
    dplyr::select(-propagation_mask, -mz_delta, -mod_name, -parent_mz)

  return(propagated_data)
}
