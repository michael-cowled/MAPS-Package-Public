#' @title Append Propagated Annotations to a Data Frame
#' @description Merges propagated annotations into the main annotation data frame, calculates mass differences to identify specific modifications, and updates compound names.
#'
#' @param full.annotation.data The original full annotation data frame. Must contain 'feature.ID' and 'mz'.
#' @param propagated_df The data frame of newly propagated annotations. Must contain 'feature.ID' (the target) and 'connected_to_ID' (the source/parent feature ID).
#' @param mod_db The modification dataframe (defaults to the modification_db defined above).
#' @param ppm_tol The mass tolerance in ppm (default 10).
#' @param abs_tol The absolute mass tolerance in Daltons (default 0.01).
#'
#' @return A single data frame with the propagated annotations appended and specific modifications named.
#' @importFrom dplyr %>% left_join mutate select filter rename
#' @export
append_propagated_annotations <- function(full.annotation.data,
                                          propagated_df,
                                          mod_db = modification_db,
                                          ppm_tol = 10,
                                          abs_tol = 0.01) {

  message("full.annotations columns: ", paste(names(full.annotation.data), collapse=", "))
  message("propagated columns: ", paste(names(propagated_df), collapse=", "))

  # 1. Retrieve the Source/Parent MZ
  # We look up the 'mz' of the 'connected_to_ID' from the full dataset
  # Assumes propagated_df has a column 'connected_to_ID' representing the parent node
  if(!"connected_to_ID" %in% names(propagated_df)){
    stop("propagated_df must contain a column 'connected_to_ID' to retrieve the parent M/Z.")
  }

  # Create a temporary lookup for M/Z
  mz_lookup <- full.annotation.data %>%
    dplyr::select(feature.ID, mz) %>%
    dplyr::rename(parent_mz = mz, connected_to_ID = feature.ID)

  # Join parent M/Z into propagated_df
  propagated_df_w_mass <- propagated_df %>%
    dplyr::rename(connected_to_ID = feature.ID)
    dplyr::left_join(mz_lookup, by = "connected_to_ID")


  # 2. Main Join and Calculation
  propagated_data <- full.annotation.data %>%
    dplyr::left_join(propagated_df_w_mass, by = "feature.ID") %>%
    dplyr::mutate(
      confidence.level = as.character(confidence.level),
      propagation_mask = is.na(compound.name) & !is.na(Probable.Analogue.Of),

      # --- CALCULATION BLOCK ---
      # Calculate Delta: Target (full) - Source (propagated)
      mz_delta = ifelse(propagation_mask, mz - parent_mz, NA),

      # Identify Modification
      # We use a custom logic to scan the mod_db for the closest match
      # We map over the rows because fuzzy matching is hard to vectorize simply
      mod_name = purrr::map_chr(mz_delta, function(x) {
        if (is.na(x)) return(NA_character_)

        # Calculate errors
        diffs <- abs(mod_db$Mass.Change - x)

        # Check tolerance (Absolute or PPM)
        # We use absolute here as it is safer for small mass differences (2Da),
        # but you can switch to ppm if preferred: (diffs / abs(x) * 1e6) <= ppm_tol
        match_idx <- which(diffs <= abs_tol)

        if (length(match_idx) > 0) {
          # Return the modification name with the smallest error
          best_match <- match_idx[which.min(diffs[match_idx])]
          return(mod_db$Modification[best_match])
        } else {
          return("Probable") # Default if no specific mass match found
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

      NPC.pathway = as.character(NPC.pathway),

      NPC.superclass = ifelse(
        propagation_mask,
        Propagated.Annotation.Class,
        NPC.superclass
      ),

      annotation.type = ifelse(
        propagation_mask,
        "GNPS Propagation",
        annotation.type
      )
    ) %>%
    # Clean up helper columns
    dplyr::select(-propagation_mask, -mz_delta, -mod_name, -parent_mz)

  return(propagated_data)
}
