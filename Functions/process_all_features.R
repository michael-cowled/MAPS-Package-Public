## 10. process_all_features
# Main processing function for collapsing by ion identity networking
process_all_features <- function(data, id_column, feature_columns) {
  unique_ids <- unique(data[[id_column]][!is.na(data[[id_column]])])
  results <- lapply(unique_ids, function(id) {
    group_data <- filter(data, !!sym(id_column) == id)
    process_iin_group(group_data, feature_columns)
  })
  names(results) <- as.character(unique_ids)
  return(results)
}