#' Update Data Frame Based on Results
#'
#' Updates a data frame by replacing rows with non-NA IDs with those in a list of result data frames,
#' while retaining rows with NA in the specified ID column.
#'
#' @param df A data frame to update.
#' @param results A list of data frames containing updates (must include `feature.ID` column).
#' @param id_column A string specifying the column name used to distinguish rows (e.g., "ion.identity.ID").
#'
#' @return A data frame combining updated rows from `results` with the original rows that had `NA` in `id_column`.
#' @export
#'
#' @examples
#' \dontrun{
#' updated <- update_data_frame(main_df, list(results1, results2), "ion.identity.ID")
#' }
update_data_frame <- function(df, results, id_column) {
  df_non_na <- df %>% filter(!is.na(!!sym(id_column)))
  results_df <- bind_rows(results)

  updated_df <- df_non_na %>% filter(feature.ID %in% results_df$feature.ID)
  df_na <- df %>% filter(is.na(!!sym(id_column)))

  final_df <- bind_rows(updated_df, df_na)
  return(final_df)
}
