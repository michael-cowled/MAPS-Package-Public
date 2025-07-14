## 6. get_result
# Function to get the first non-NA value from specified columns, the column name, and the superclass
get_result <- function(paired_value, summary_data) {
  # Define the compound name column
  compound_col <- "compound.name"
  
  # Define the superclass columns to check
  superclass_columns <- c("NPC.superclass", "canopus.NPC.superclass")
  
  # Filter and ensure only a single row is returned
  data_subset <- summary_data %>%
    filter(feature.ID == paired_value) %>%
    select(all_of(compound_col), all_of(superclass_columns), smiles) %>%
    slice(1)  # Take the first row if duplicates exist
  
  compound_value <- data_subset[[compound_col]]
  
  if (!is.na(compound_value) && compound_value != "null") {
    # Find the first valid superclass
    superclass_value <- NA
    for (col in superclass_columns) {
      sc_val <- data_subset[[col]]
      if (!is.na(sc_val) && sc_val != "null") {
        superclass_value <- sc_val
        break
      }
    }
    return(list(value = compound_value, column = compound_col, superclass = superclass_value))
  }
  
  # Fallback to smiles
  csi_value <- data_subset$smiles
  if (!is.na(csi_value) && csi_value != "null") {
    return(list(value = csi_value, column = "smiles", superclass = NA))
  }
  
  return(list(value = NA, column = NA, superclass = NA))
}