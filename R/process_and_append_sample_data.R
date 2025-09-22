#' @title Process and Append Sample Data
#' @description Cleans and processes sample data to create a presence/absence list, then appends it to a main annotation data frame. Also returns the intermediate `sample.data2` data frame.
#'
#' @param sample.data A data frame containing sample information, including peak areas.
#' @param propagated.annotation.data A data frame of existing annotations to which the sample data will be appended.
#'
#' @return A list containing two data frames: `combined_data` (the merged annotation and sample data) and `sample.data2` (the processed sample data before presence/absence conversion).
#' @importFrom dplyr %>% select full_join contains
#' @export
process_and_append_sample_data <- function(sample.data, propagated.annotation.data) {

  # Select columns containing ".area" and rename the first column
  sample.data <- sample.data %>%
    dplyr::select(id, contains(".area"))

  # Rename the first column to "feature.ID"
  colnames(sample.data)[1] <- "feature.ID"

  # Remove prefixes and suffixes from column names
  colnames_sample <- colnames(sample.data)
  colnames_sample <- sub("^datafile\\.", "", colnames_sample)
  colnames_sample <- sub("\\.mzML\\.area$", "", colnames_sample)
  colnames(sample.data) <- colnames_sample

  # Make a copy for potential other uses
  sample.data2 <- sample.data

  # Convert peak areas to a binary presence/absence format (1 for positive area, 0 otherwise)
  sample.data[, 2:ncol(sample.data)] <- lapply(sample.data[, 2:ncol(sample.data)], function(x) {
    x[x > 0] <- 1
    return(x)
  })

  # Create a new column 'Samples' with a semicolon-separated list of present samples
  sample.data$Samples <- apply(sample.data[, 2:ncol(sample.data)], 1, function(row) {
    paste(colnames(sample.data[, 2:ncol(sample.data)])[which(row == 1)], collapse = "; ")
  })

  # Select only the 'feature.ID' and 'Samples' columns
  sample.data <- sample.data[, c(1, ncol(sample.data))]

  # Append the processed sample data to the annotations data frame
  combined_data <- propagated.annotation.data %>%
    dplyr::full_join(sample.data, by = "feature.ID")

  # Return a list containing both the combined data and sample.data2
  return(list(
    combined_data = combined_data,
    sample.data2 = sample.data2
  ))
}
