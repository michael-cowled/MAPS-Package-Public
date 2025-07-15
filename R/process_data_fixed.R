## 11. process_data_fixed
redundancy_fixer <- function(dataset, column_to_check = "Best.Annotation.Smiles", rt_column = "rt", rt_tolerance = 1) {
  # --- Error Handling: Check column existence ---
  if (!column_to_check %in% colnames(dataset)) {
    stop(paste("Column '", column_to_check, "' not found in the dataset.", sep = ""))
  }
  if (!rt_column %in% colnames(dataset)) {
    stop(paste("Column '", rt_column, "' not found in the dataset.", sep = ""))
  }
  
  # --- Ensure rt_column is numeric ---
  if (!is.numeric(dataset[[rt_column]])) {
    warning(paste("Column '", rt_column, "' is not numeric. Attempting to convert.", sep = ""))
    dataset[[rt_column]] <- as.numeric(as.character(dataset[[rt_column]]))
    if (any(is.na(dataset[[rt_column]]))) {
      stop(paste("Some values in ", rt_column, " could not be converted to numeric. Please inspect data", sep = ""))
    }
  }
  
  # Add a new column to mark redundancy (initialized as FALSE)
  dataset$redundant <- FALSE
  
  unique_values <- unique(dataset[[column_to_check]]) # corrected line
  
  for (value in unique_values) {
    temp_df <- dataset %>% filter(!!sym(column_to_check) == value)
    
    if (nrow(temp_df) > 1) {
      rt_values <- sort(temp_df[[rt_column]])
      rt_values <- rt_values[!is.na(rt_values)]
      
      grouped_rts <- list()
      used_indices <- rep(FALSE, length(rt_values))
      
      for (i in seq_along(rt_values)) {
        if (!used_indices[i]) {
          group <- rt_values[i]
          used_indices[i] <- TRUE
          for (j in (i + 1):length(rt_values)) {
            if (!is.na(rt_values[i]) && !is.na(rt_values[j])) {
              if (abs(rt_values[j] - rt_values[i]) <= rt_tolerance) {
                group <- c(group, rt_values[j])
                used_indices[j] <- TRUE
              }
            }
          }
          grouped_rts <- c(grouped_rts, list(sort(group)))
        }
      }
      
      for (group in grouped_rts) {
        group_df <- temp_df %>% filter(!!sym(rt_column) %in% group)
        if (nrow(group_df) > 1) {
          dataset$redundant[which(dataset[[rt_column]] %in% group_df[[rt_column]][-1])] <- TRUE
        }
      }
    }
  }
  
  return(dataset)
}