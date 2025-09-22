#' @title Prepare and Append Cytoscape File
#' @description Reads the GNPS Cytoscape file, appends the best annotations, and writes the final file to disk.
#'
#' @param folder The base folder path for file input/output.
#' @param propagated.annotation.data.with.samples The data frame containing propagated annotations and samples.
#' @param read_checked_csv A function to read checked CSV files.
#' @param write_large_csv A function to write large CSV files.
#' @return The processed Cytoscape data frame.
#' @importFrom dplyr %>% select full_join
#' @export
prepare_cytoscape_file <- function(
    folder,
    propagated.annotation.data.with.samples,
    read_checked_csv,
    write_large_csv
) {
  # Read the original Cytoscape file
  cytoscape <- read_checked_csv(paste0(folder, "/gnps/cytoscape.csv"))

  # Prepare a data frame with feature.ID and compound.name for merging
  cytoscape.annotations <- propagated.annotation.data.with.samples %>%
    dplyr::select(feature.ID, compound.name)
  names(cytoscape.annotations) <- c("shared.name", "library_compound_name2")

  # Tidy up the compound names
  cytoscape.annotations$library_compound_name2 <- ifelse(
    is.na(cytoscape.annotations$library_compound_name2) | cytoscape.annotations$library_compound_name2 == "NA",
    "",
    cytoscape.annotations$library_compound_name2
  )
  cytoscape.annotations$shared.name <- as.numeric(cytoscape.annotations$shared.name)

  # Join the annotations with the Cytoscape data
  cytoscape <- cytoscape %>%
    dplyr::full_join(cytoscape.annotations, by = "shared.name")
  colnames(cytoscape)[colnames(cytoscape) == 'shared.name'] <- 'shared name'

  # Write the final Cytoscape file
  write_large_csv(cytoscape, paste0(folder, "/cytoscape-v2.csv"))

  return(cytoscape)
}
