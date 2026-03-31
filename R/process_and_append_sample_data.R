#' @title Process and Append Sample Data
#' @description Cleans and processes sample data to create a presence/absence list,
#' then appends it to a main annotation data frame.
#'
#' @param sample.data A data frame containing sample information, including peak areas.
#' @param propagated.annotation.data A data frame of existing annotations.
#'
#' @export
process_and_append_sample_data <- function(sample.data, propagated.annotation.data) {

  # 1. Standardize the Join Key (feature.ID) in the annotation data
  propagated.annotation.data <- propagated.annotation.data %>%
    dplyr::mutate(feature.ID = as.character(feature.ID))

  # 2. Select and Rename columns in sample data
  sample.data <- sample.data %>%
    dplyr::select(id, contains(".area"))

  # Rename 'id' to 'feature.ID'
  colnames(sample.data)[1] <- "feature.ID"

  # --- THE FIX: Force sample IDs to character immediately ---
  sample.data$feature.ID <- as.character(sample.data$feature.ID)

  # 3. Clean up column names (Remove prefixes/suffixes)
  colnames_sample <- colnames(sample.data)
  colnames_sample <- sub("^datafile\\.", "", colnames_sample)
  colnames_sample <- sub("\\.mzML\\.area$", "", colnames_sample)
  colnames(sample.data) <- colnames_sample

  # Make a copy for intermediate output (this will now have character IDs too)
  sample.data2 <- sample.data

  # 4. Convert peak areas to binary (1/0)
  # We start from index 2 to skip the 'feature.ID' column
  sample.data[, 2:ncol(sample.data)] <- lapply(sample.data[, 2:ncol(sample.data)], function(x) {
    x <- as.numeric(x) # Ensure areas are numeric
    x[is.na(x)] <- 0
    x[x > 0] <- 1
    return(x)
  })

  # 5. Create 'Samples' summary column
  # Identify which samples had a '1' and collapse them into a string
  sample_cols <- colnames(sample.data)[2:ncol(sample.data)]

  sample.data$Samples <- apply(sample.data[, 2:ncol(sample.data)], 1, function(row) {
    present_samples <- sample_cols[which(row == 1)]
    paste(present_samples, collapse = "; ")
  })

  # 6. Keep only the Join Key and the new summary column
  sample.data_final <- sample.data[, c("feature.ID", "Samples")]

  # 7. Join the data
  # Since both are now <character>, the join will succeed
  combined_data <- propagated.annotation.data %>%
    dplyr::full_join(sample.data_final, by = "feature.ID")

  return(list(
    combined_data = combined_data,
    sample.data2 = sample.data2
  ))
}
