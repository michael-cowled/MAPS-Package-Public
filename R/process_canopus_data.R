#' @title Process CANOPUS Data
#' @description Loads, cleans, and processes CANOPUS annotations from a TSV file.
#' @param canopus.data_path Path to the CANOPUS annotations TSV file.
#' @return A data frame containing the processed CANOPUS annotations.
#' @importFrom dplyr %>% group_by filter slice ungroup
#' @export
process_canopus_data <- function(canopus.data_path) {
  # Load the data from the provided file path
  canopus.data <- read_checked_tsv(canopus.data_path)

  # Select and rename columns
  canopus.data <- canopus.data[, c(5:8, 27)]
  names(canopus.data) <- c(
    "canopus.NPC.pathway", "canopus.NPC.pathway.probability",
    "canopus.NPC.superclass", "canopus.NPC.superclass.probability",
    'feature.ID'
  )

  # Filter data to keep the highest probability entry for each feature.ID
  canopus.data <- canopus.data %>%
    group_by(feature.ID) %>%
    filter(!all(canopus.NPC.pathway.probability == 0)) %>%
    filter(canopus.NPC.pathway.probability == max(canopus.NPC.pathway.probability, na.rm = TRUE)) %>%
    slice(1) %>%
    ungroup()

  return(canopus.data)
}
