## 9. update_data_frame
# Define a function to update the main data frame
update_data_frame <- function(df, results, id_column) {
  df_non_na <- df %>% filter(!is.na(!!sym(id_column)))
  results_df <- bind_rows(results)
  
  updated_df <- df_non_na %>% filter(feature.ID %in% results_df$feature.ID)
  df_na <- df %>% filter(is.na(!!sym(id_column)))
  
  final_df <- bind_rows(updated_df, df_na)
  return(final_df)
}