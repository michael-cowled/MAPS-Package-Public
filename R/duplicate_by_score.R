#' Deduplicate a Data Frame by Score
#'
#' This function groups a data frame by a specified ID column and
#' keeps only the row with the maximum value in a score column.
#' In case of ties, it keeps the first row encountered.
#'
#' @param df A data frame to be processed.
#' @param id_col A bare column name or string specifying the ID column to group by.
#' @param score_col A bare column name or string specifying the numeric score column to filter by.
#'
#' @return A data frame with deduplicated rows based on the highest score.
#' @export
#'
#' @examples
#' my_df <- data.frame(
#'   feature_id = c(1, 1, 2, 2, 2),
#'   score = c(0.8, 0.9, 0.7, 0.7, 0.8),
#'   name = c("A", "B", "C", "D", "E")
#' )
#' deduplicated_df <- deduplicate_by_score(my_df, id_col = feature_id, score_col = score)
#' print(deduplicated_df)
deduplicate_by_score <- function(df, id_col, score_col) {
  df %>%
    group_by(across({{ id_col }})) %>%
    filter(get(as.character(substitute(score_col))) == max(get(as.character(substitute(score_col))), na.rm = TRUE)) %>%
    slice(1) %>%
    ungroup()
}
