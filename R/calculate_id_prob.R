#' Compute Identification Probability and Deduplicate
#'
#' Calculates identification probability for annotations based on a score
#' threshold and then deduplicates the data frame, keeping only the best
#' scoring annotation per feature.
#'
#' @param df A data frame with annotation data.
#' @param score_col A string specifying the name of the score column.
#' @param threshold A numeric value for the confidence score threshold.
#'
#' @return A data frame with a new 'id.prob' column and one row per 'feature.ID'.
#' @export
#' @importFrom rlang .data sym
#' @importFrom dplyr group_by mutate ungroup select arrange slice desc
#' @importFrom magrittr %>%
#'
#' @examples
#' my_data <- data.frame(feature.ID = c(1, 1, 2, 2, 2),
#'                      confidence.score = c(0.8, 0.9, 0.7, 0.6, 0.9),
#'                      name = c("A", "B", "C", "D", "E"))
#' result <- compute_id_prob(my_data, "confidence.score", 0.7)
#' print(result)
compute_id_prob <- function(df, score_col, threshold) {
  score_col_sym <- rlang::sym(score_col)
  df %>%
    dplyr::group_by(.data$feature.ID) %>%
    dplyr::mutate(
      n_above_thresh = sum(.data[[score_col]] >= threshold, na.rm = TRUE),
      id.prob = ifelse(.data[[score_col]] >= threshold & n_above_thresh > 0, 1 / n_above_thresh, 0)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$n_above_thresh) %>%
    dplyr::group_by(.data$feature.ID) %>%
    dplyr::arrange(dplyr::desc(.data[[score_col]])) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
}
