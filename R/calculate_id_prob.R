#' Calculate Identification Probability
#'
#' Calculates a new 'id.prob' column. For each ID, it checks how many
#' annotations are above a specified probability threshold and assigns an
#' identification probability as the inverse of that count. Rows below the
#' threshold get an ID probability of 0.
#'
#' @param df A data frame containing annotation data.
#' @param id_col A bare column name or string specifying the ID column to group by.
#' @param score_col A bare column name or string specifying the numeric score column.
#' @param prob_threshold A numeric value representing the probability cutoff.
#'
#' @return The input data frame with a new 'id.prob' column.
#' @importFrom dplyr group_by mutate ungroup select
#' @importFrom rlang `[[`
#' @export
#'
#' @examples
#' my_df <- data.frame(
#'   feature_id = c(1, 1, 2, 2),
#'   score = c(0.8, 0.6, 0.7, 0.9)
#' )
#' df_with_prob <- calculate_id_prob(my_df, id_col = feature_id, score_col = score, prob_threshold = 0.7)
#' print(df_with_prob)
calculate_id_prob <- function(df, id_col, score_col, prob_threshold) {
  df %>%
    group_by(across({{ id_col }})) %>%
    mutate(
      n_above_thresh = sum(get(as.character(substitute(score_col))) >= prob_threshold),
      id.prob = ifelse(get(as.character(substitute(score_col))) >= prob_threshold & n_above_thresh > 0, 1 / n_above_thresh, 0)
    ) %>%
    ungroup() %>%
    select(-n_above_thresh)
}
