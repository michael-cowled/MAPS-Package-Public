#' @title Append Level 4 and 5 Annotations
#' @description Appends annotations based on a provided mask and updates relevant columns.
#'
#' @param data The data frame to be annotated.
#' @param mask_condition A logical vector (or an expression that evaluates to one) indicating which rows to update.
#' @param compound_name The new compound name to be assigned.
#' @param confidence_level The confidence level to be assigned.
#' @param npc_pathway The NPC pathway to be assigned.
#' @param npc_superclass The NPC superclass to be assigned.
#' @param annotation_type The annotation type to be assigned.
#'
#' @return The data frame with the new annotations appended.
#' @importFrom dplyr %>% mutate select
#' @export
append_annotations <- function(data, mask_condition, compound_name, confidence_level, npc_pathway, npc_superclass, annotation_type) {
  data %>%
    dplyr::mutate(
      confidence.level = as.character(confidence.level),
      propagation_mask = {{ mask_condition }},
      compound.name = ifelse(
        propagation_mask,
        compound_name,
        compound.name
      ),
      confidence.level = ifelse(
        propagation_mask,
        confidence_level,
        confidence.level
      ),
      NPC.pathway = ifelse(
        propagation_mask,
        npc_pathway,
        NPC.pathway
      ),
      NPC.superclass = ifelse(
        propagation_mask,
        npc_superclass,
        NPC.superclass
      ),
      annotation.type = ifelse(
        propagation_mask,
        annotation_type,
        annotation.type
      )
    ) %>%
    dplyr::select(-propagation_mask)
}
