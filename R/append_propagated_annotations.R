#' @title Append Propagated Annotations to a Data Frame
#' @description Merges propagated annotations into the main annotation data frame and updates relevant columns.
#'
#' @param full.annotation.data The original full annotation data frame.
#' @param propagated_df The data frame of newly propagated annotations.
#'
#' @return A single data frame with the propagated annotations appended.
#' @importFrom dplyr %>% left_join mutate select
#' @export
append_propagated_annotations <- function(full.annotation.data, propagated_df) {
  propagated_data <- full.annotation.data %>%
    dplyr::left_join(propagated_df, by = "feature.ID") %>%
    dplyr::mutate(
      confidence.level = as.character(confidence.level),  # ensure editable
      propagation_mask = is.na(compound.name) & !is.na(Probable.Analogue.Of), # condition mask
      compound.name = ifelse(
        propagation_mask,
        paste0("Probable analogue of: ", Probable.Analogue.Of),
        compound.name
      ),
      confidence.level = ifelse(
        propagation_mask,
        "3",
        confidence.level
      ),
      NPC.pathway = as.character(NPC.pathway), # Ensure the type is consistent
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
    dplyr::select(-propagation_mask)

  return(propagated_data)
}
