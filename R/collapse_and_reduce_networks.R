#' @title Collapse Ion Identity Networks and Reduce Redundancy
#' @description This function collapses ion identity networks to retain the best annotation and then reduces redundancy based on SMILES. It processes the full annotation data, handles samples, and tidies the final output.
#'
#' @param propagated.annotation.data.with.samples The main data frame containing all annotations.
#' @param sample.data2 The processed sample data with peak areas, before presence/absence conversion.
#' @param process_all_features A function to process all features (e.g., from an external package).
#' @param update_data_frame A function to update the data frame by removing rows (e.g., from an external package).
#' @param fix_compound_names A function to fix compound names (e.g., from an external package).
#' @param redundancy_fixer A function to reduce redundancy (e.g., from an external package).
#' @return A list containing the final, processed annotation data frame and a samples data frame for visualization.
#' @importFrom dplyr %>% filter select mutate full_join
#' @importFrom tidyr pivot_longer
#' @export
collapse_and_reduce_networks <- function(
    propagated.annotation.data.with.samples,
    sample.data2,
    process_all_features,
    update_data_frame,
    fix_compound_names,
    redundancy_fixer
) {

  # Step 1: Collapsing Ion Identity Networks
  # This part of the code assumes that process_all_features and update_data_frame are functions you have defined elsewhere.
  # They are passed as arguments to make the function more robust and flexible.

  # Filter for features with Ion Identity IDs
  iin.features <- dplyr::filter(propagated.annotation.data.with.samples, !is.na(ion.identity.ID) & ion.identity.ID != "")

  # Process all features to get the best annotations
  results <- process_all_features(iin.features, "ion.identity.ID", c("feature1", "feature2"))

  # Update the final data frame by removing redundant rows
  final_annotation_df <- update_data_frame(propagated.annotation.data.with.samples, results, "ion.identity.ID") %>%
    dplyr::select(-ion.identity.ID, -Probable.Analogue.Of)

  # Step 2: Create Samples Data Frame
  # Remove the 'Samples' column if it exists and pivot the data to long format
  if ("Samples" %in% names(sample.data2)) {
    sample.data2 <- sample.data2 %>% dplyr::select(-Samples)
  }

  long_df <- sample.data2 %>%
    tidyr::pivot_longer(
      cols = -feature.ID,
      names_to = "samples",
      values_to = "area"
    )

  # Join with the final annotation data
  samples_df <- final_annotation_df %>%
    dplyr::select(feature.ID, feature.usi, compound.name, smiles, Formula, IUPAC, Monoisotopic.Mass)
  samples_df$feature.ID <- as.numeric(samples_df$feature.ID)

  samples_df <- long_df %>%
    dplyr::full_join(samples_df, by = "feature.ID")

  # Step 3: Fix Compound Names
  final_annotation_df <- fix_compound_names(final_annotation_df, "compound.name")

  # Step 4: Redundancy Reduction
  # Separate data into those with and without SMILES strings
  dataset <- final_annotation_df %>%
    dplyr::filter(!is.na(smiles) & smiles != "N/A")

  other_data <- final_annotation_df %>%
    dplyr::filter(is.na(smiles) | smiles == "N/A")

  # Add this block to explicitly convert the 'rt' column
  if("rt" %in% names(dataset)) {
    dataset$rt <- as.numeric(dataset$rt)
  }
  if("rt" %in% names(other_data)) {
    other_data$rt <- as.numeric(other_data$rt)
  }

  if(nrow(dataset) > 0) {
    dataset$rt <- as.numeric(dataset$rt)
    dataset$redundant <- FALSE
    dataset <- redundancy_fixer(dataset, column_to_check = "smiles")
  }

  # Combine the processed data frames
  final_annotation_df <- dplyr::bind_rows(dataset, other_data)

  return(list(
    final_annotation_df = final_annotation_df,
    samples_df = samples_df
  ))
}
