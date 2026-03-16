#' @title Process ZODIAC Data
#' @description Loads, cleans, and processes ZODIAC annotations from a TSV file.
#' @param zodiac.data_path Path to the ZODIAC annotations TSV file.
#' @return A data frame containing the processed ZODIAC annotations.
#' @importFrom dplyr %>% group_by filter slice ungroup select
#' @export
process_zodiac_data <- function(zodiac.data_path) {
  # 1. Load the data
  zodiac.data <- read_checked_tsv(zodiac.data_path)

  # 2. Select and rename using column names for stability
  zodiac.data <- zodiac.data %>%
    dplyr::select(
      zodiac.formula = molecularFormula,
      zodiac.confidence.score = ZodiacScore,
      feature.ID = mappingFeatureId
    )

  # 3. Ensure numeric types for logical operations
  zodiac.data <- zodiac.data %>%
    dplyr::mutate(
      zodiac.confidence.score = as.numeric(zodiac.confidence.score),
      feature.ID = as.numeric(feature.ID)
    )

  # 4. Filter data to keep the highest confidence score entry for each feature.ID
  zodiac.data <- zodiac.data %>%
    dplyr::group_by(feature.ID) %>%
    dplyr::filter(!all(zodiac.confidence.score == 0 | is.na(zodiac.confidence.score))) %>%
    dplyr::filter(zodiac.confidence.score == max(zodiac.confidence.score, na.rm = TRUE)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  return(zodiac.data)
}
