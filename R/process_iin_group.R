#' Process Ion Identity Network (IIN) Group
#'
#' For a given group of rows sharing the same ion identity network (IIN), this function determines
#' which features (rows) to retain based on the presence of non-NA values in prioritized columns.
#'
#' It retains:
#' \itemize{
#'   \item All rows except the one with the highest \code{feature.ID}, if multiple non-NA entries exist in the column.
#'   \item The single row with a non-NA value if only one exists.
#'   \item As a fallback, the row with the smallest \code{feature.ID}.
#' }
#'
#' @param group_data A data frame representing one IIN group.
#' @param feature_columns A character vector of column names to check for non-NA values (in priority order).
#'
#' @return A subset of \code{group_data} containing only the retained rows.
#' @export
#'
#' @examples
#' \dontrun{
#' filtered <- process_iin_group(group_df, c("annotation", "smiles"))
#' }
process_iin_group <- function(group_data, feature_columns) {
  # Step 1: try by lowest confidence.level
  if ("confidence.level" %in% names(group_data) && any(!is.na(group_data$confidence.level))) {
    min_level <- min(group_data$confidence.level, na.rm = TRUE)
    level_candidates <- group_data[group_data$confidence.level == min_level, ]
  } else {
    level_candidates <- group_data
  }

  # Step 2: within candidates, try by highest confidence.score
  if ("confidence.score" %in% names(level_candidates) && any(!is.na(level_candidates$confidence.score))) {
    max_score <- max(level_candidates$confidence.score, na.rm = TRUE)
    score_candidates <- level_candidates[level_candidates$confidence.score == max_score, ]
  } else {
    score_candidates <- level_candidates
  }

  # Step 3: if still multiple rows, keep lowest feature.ID
  chosen <- score_candidates[which.min(score_candidates$feature.ID), , drop = FALSE]

  return(chosen)
}
