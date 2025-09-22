#' @title Process ZODIAC Data
#' @description Loads, cleans, and processes ZODIAC annotations from a TSV file.
#' @param zodiac.data_path Path to the ZODIAC annotations TSV file.
#' @return A data frame containing the processed ZODIAC annotations.
#' @importFrom dplyr %>% group_by filter slice ungroup
#' @export
process_zodiac_data <- function(zodiac.data_path) {
  # Load the data from the provided file path
  zodiac.data <- read_checked_tsv(zodiac.data_path)

  # Select and rename columns
  zodiac.data <- zodiac.data[, c(2, 5, 21)]
  names(zodiac.data) <- c("zodiac.formula", "zodiac.confidence.score", 'feature.ID')

  # Filter data to keep the highest confidence score entry for each feature.ID
  zodiac.data <- zodiac.data %>%
    group_by(feature.ID) %>%
    filter(!all(zodiac.confidence.score == 0)) %>%
    filter(zodiac.confidence.score == max(zodiac.confidence.score, na.rm = TRUE)) %>%
    slice(1) %>%
    ungroup()

  return(zodiac.data)
}
