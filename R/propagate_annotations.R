#' @title Propagate Annotations Based on GNPS Cluster Pairs
#' @description Finds and propagates annotations for unknown features based on GNPS cluster relationships.
#'
#' @param full.annotation.data A data frame containing all existing annotations.
#' @param gnps.cluster.pairs A data frame of feature ID pairs from GNPS clusters.
#' @param paired_feature_finder A function to find paired features (e.g., from MAPS.Package).
#' @param get_result A function to retrieve annotation details for a feature ID.
#'
#' @return A data frame with propagated annotations and structural metadata.
#' @export
propagate_annotations <- function(full.annotation.data, gnps.cluster.pairs, paired_feature_finder, get_result) {

  # 1. Identify "Targets": Features that need a better ID.
  # We include Level 3 here so they can be "upgraded" by a Level 1/2 neighbor.
  na.rows <- dplyr::filter(full.annotation.data,
                           is.na(compound.name) |
                             confidence.level == "3")

  na.feature.ids <- na.rows$feature.ID

  if (length(na.feature.ids) == 0) {
    return(tibble::tibble())
  }

  pb <- progress::progress_bar$new(
    format = "Propagating annotations [:bar] :percent eta: :eta",
    total = length(na.feature.ids),
    width = 60
  )

  # 2. Main Propagation Loop
  propagated_df <- purrr::map_dfr(na.feature.ids, function(i) {
    pb$tick()

    # Find features clustered with this ID
    paired_values <- paired_feature_finder(i, gnps.cluster.pairs)

    selected_paired_value <- NA
    final_result_data <- list(value = NA, column = NA, superclass = NA)
    parent_smiles <- NA_character_

    # Get the dataframe of paired features and their cosine scores
    paired_df <- paired_feature_finder(i, gnps.cluster.pairs)

    selected_cosine_score <- NA_real_ # Initialize empty score

    # Iterate over the rows of the new dataframe
    for (row_idx in seq_len(nrow(paired_df))) {
      value <- paired_df$paired_value[row_idx]
      current_cosine <- paired_df$Cosine[row_idx]

      # Extract metadata for the potential parent
      parent_meta <- full.annotation.data %>%
        dplyr::filter(feature.ID == value) %>%
        dplyr::select(annotation.type, confidence.level, smiles) %>%
        dplyr::slice(1)

      is_valid_source <- !is.na(parent_meta$confidence.level) &&
        parent_meta$confidence.level != "3"

      if (is_valid_source) {
        result_data <- get_result(value, full.annotation.data)

        if (!is.na(result_data$value)) {
          selected_paired_value <- value
          final_result_data <- result_data
          parent_smiles <- parent_meta$smiles
          final_result_data$column <- parent_meta$annotation.type
          selected_cosine_score <- current_cosine # Save the winning cosine score
          break
        }
      }
    }

    # 3. Construct the output row (add confidence.score here)
    if (!is.na(selected_paired_value)) {
      tibble::tibble(
        feature.ID = i,
        Probable.Analogue.Of = final_result_data$value,
        Propagated.Feature.ID = selected_paired_value,
        Propagated.Smiles = parent_smiles,
        Propagated.Annotation.Type = final_result_data$column, # Changed from .Class to .Type
        Propagated.Confidence.Score = selected_cosine_score
      )
    } else {
      tibble::tibble()
    }
  })

  return(propagated_df)
}
