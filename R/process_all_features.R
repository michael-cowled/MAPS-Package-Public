#' Process All Features by Ion Identity Network (IIN)
#'
#' Applies \code{process_iin_group()} to all non-NA groups in a dataset based on a specified ID column
#' (typically an ion identity network identifier). Returns a list of data frames with redundant features removed.
#'
#' @param data A data frame containing the full dataset.
#' @param id_column A string giving the name of the column that identifies ion identity groups (e.g., \code{"ion.identity.ID"}).
#' @param feature_columns A character vector of column names to check for annotation presence (used for pruning).
#'
#' @return A named list of data frames, one per unique non-NA ID, with filtered feature rows.
#' @export
#'
#' @seealso \code{\link{process_iin_group}}
#'
#' @examples
#' \dontrun{
#' cleaned_groups <- process_all_features(mydata, "ion.identity.ID", c("annotation", "smiles"))
#' }
process_all_features <- function(data, id_column, feature_columns) {
  unique_ids <- unique(data[[id_column]][!is.na(data[[id_column]])])

  results <- lapply(unique_ids, function(id) {
    group_data <- dplyr::filter(data, !!rlang::sym(id_column) == id)
    process_iin_group(group_data, feature_columns)
  })

  names(results) <- as.character(unique_ids)
  return(results)
}
