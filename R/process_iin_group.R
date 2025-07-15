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
  rows_to_retain <- logical(nrow(group_data))

  for (col in feature_columns) {
    if (any(!is.na(group_data[[col]]))) {
      values <- group_data$feature.ID[!is.na(group_data[[col]])]
      if (length(values) > 1) {
        # Retain all but the highest feature.ID
        highest_value <- max(values, na.rm = TRUE)
        rows_to_retain <- group_data$feature.ID %in% values & group_data$feature.ID != highest_value
      } else {
        # Only one match, retain it
        rows_to_retain <- group_data$feature.ID == values
      }
      break  # Stop at the first matching column
    }
  }

  if (all(!rows_to_retain)) {
    # Fallback: retain row with smallest feature.ID
    rows_to_retain <- group_data$feature.ID == min(group_data$feature.ID, na.rm = TRUE)
  }

  return(group_data[rows_to_retain, ])
}
