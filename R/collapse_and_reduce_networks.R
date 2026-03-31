#' @title Collapse Ion Identity Networks and Reduce Redundancy
#' @export
collapse_and_reduce_networks <- function(
    propagated.annotation.data.with.samples,
    sample.data2,
    process_all_features,
    update_data_frame,
    fix_compound_names,
    redundancy_fixer
) {

  # --- PRE-STEP: Global ID Standardization ---
  # Ensure all joining keys are characters immediately
  propagated.annotation.data.with.samples$feature.ID <- as.character(propagated.annotation.data.with.samples$feature.ID)
  sample.data2$feature.ID <- as.character(sample.data2$feature.ID)

  # --- Step 1: Collapsing Ion Identity Networks ---
  iin.features <- dplyr::filter(propagated.annotation.data.with.samples,
                                !is.na(ion.identity.ID) & ion.identity.ID != "")

  results <- process_all_features(iin.features, "ion.identity.ID", c("feature1", "feature2"))

  final_annotation_df <- update_data_frame(propagated.annotation.data.with.samples, results, "ion.identity.ID") %>%
    dplyr::select(-ion.identity.ID, -Probable.Analogue.Of)

  # --- Step 2: Create Samples Data Frame ---
  if ("Samples" %in% names(sample.data2)) {
    sample.data2 <- sample.data2 %>% dplyr::select(-Samples)
  }

  long_df <- sample.data2 %>%
    tidyr::pivot_longer(
      cols = -feature.ID,
      names_to = "samples",
      values_to = "area"
    )

  samples_df <- final_annotation_df %>%
    dplyr::select(feature.ID, feature.usi, compound.name, smiles, Formula, IUPAC, Monoisotopic.Mass)

  # FIX: Keep feature.ID as character to match long_df
  samples_df$feature.ID <- as.character(samples_df$feature.ID)

  samples_df <- long_df %>%
    dplyr::full_join(samples_df, by = "feature.ID")

  # --- Step 3: Fix Compound Names ---
  final_annotation_df <- fix_compound_names(final_annotation_df, "compound.name")

  # --- Step 4: Redundancy Reduction ---

  # 1. Protected Data
  protected_data <- final_annotation_df %>%
    dplyr::filter(annotation.type == "ms2query" & confidence.level == "3")

  # 2. Dataset for Redundancy Check
  dataset <- final_annotation_df %>%
    dplyr::filter(!(feature.ID %in% protected_data$feature.ID))

  # Ensure RT and IDs are consistent types for processing
  if("rt" %in% names(dataset)) dataset$rt <- as.numeric(dataset$rt)
  if("rt" %in% names(protected_data)) protected_data$rt <- as.numeric(protected_data$rt)

  dataset$feature.ID <- as.character(dataset$feature.ID)
  protected_data$feature.ID <- as.character(protected_data$feature.ID)

  # Run redundancy fixer
  if(nrow(dataset) > 0) {
    dataset <- redundancy_fixer(dataset, rt_column = "rt", rt_tolerance = 1)
  }

  # Ensure protected data gets the 'redundant' column
  if(nrow(protected_data) > 0) {
    protected_data$redundant <- FALSE
  }

  # --- SAFETY: Final Type Match before bind_rows ---
  # If redundancy_fixer changed anything, force dataset to match protected_data
  common_cols <- intersect(colnames(dataset), colnames(protected_data))
  for (col in common_cols) {
    if (class(dataset[[col]]) != class(protected_data[[col]])) {
      class(dataset[[col]]) <- class(protected_data[[col]])
    }
  }

  # Combine streams back together
  final_annotation_df <- dplyr::bind_rows(dataset, protected_data)

  return(list(
    final_annotation_df = final_annotation_df,
    samples_df = samples_df
  ))
}
