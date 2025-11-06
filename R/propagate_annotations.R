#' @title Propagate Annotations Based on GNPS Cluster Pairs
#' @description Finds and propagates annotations for unknown features based on GNPS cluster relationships.
#'
#' @param full.annotation.data A data frame containing all existing annotations.
#' @param gnps.cluster.pairs A data frame of feature ID pairs from GNPS clusters.
#' @param paired_feature_finder A function to find paired features (e.g., from MAPS.Package).
#' @param get_result A function to retrieve annotation details for a feature ID.
#'
#' @return A data frame with propagated annotations.
#' @importFrom dplyr %>% filter pull first
#' @importFrom progress progress_bar
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @export
propagate_annotations <- function(full.annotation.data, gnps.cluster.pairs, paired_feature_finder, get_result) {
  # Identify the unknown features that need annotation
  na.rows <- dplyr::filter(full.annotation.data, is.na(compound.name))
  na.feature.ids <- na.rows$feature.ID

  if (length(na.feature.ids) == 0) {
    return(tibble())
  }

  # Initialize progress bar
  pb <- progress::progress_bar$new(
    format = "Propagating annotations [:bar] :percent eta: :eta",
    total = length(na.feature.ids),
    width = 60
  )

  # Propagate annotations and collect results
  propagated_df <- purrr::map_dfr(na.feature.ids, function(i) {
    pb$tick()
    paired_values <- paired_feature_finder(i, gnps.cluster.pairs)

    selected_paired_value <- NA
    final_result_data <- list(value = NA, column = NA, superclass = NA)

    for (value in paired_values) {
      result_data <- get_result(value, full.annotation.data)
      if (!is.na(result_data$value)) {
        selected_paired_value <- value
        final_result_data <- result_data
        annotation_type_value <- full.annotation.data %>%
          dplyr::filter(feature.ID == selected_paired_value) %>%
          dplyr::pull(annotation.type) %>%
          first()
        final_result_data$column <- annotation_type_value
        break
      }
    }

    if (!is.na(selected_paired_value)) {
      tibble(
        feature.ID = i,
        Probable.Analogue.Of = final_result_data$value,
        Propagated.Feature.ID = selected_paired_value,
        Propagated.Annotation.Type = final_result_data$column,
        Propagated.Annotation.Class = final_result_data$superclass
      )
    } else {
      tibble() # Return an empty tibble if no annotation is found
    }
  })

  return(propagated_df)
}
