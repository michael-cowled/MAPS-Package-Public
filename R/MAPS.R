#' @title Metabolomics ANnotation Processing and Standardisation (MAPS)
#' @description This is the main function for processing and integrating
#'   metabolomics annotation data from various sources (GNPS, SIRIUS, MS2Query).
#' @name MAPS
#' @param dataset.id Character string identifying the dataset (for file naming).
#' @param folder Character string. **Required**: Path to the specified Processed Data Folder containing results from MZmine, GNPS, SIRIUS, etc.
#' @param gnps.task.id Character string. **Required**: External GNPS ID for the task (can be extracted from metadata).
#' @param gnps.prob Numeric. Minimum probability/score for high-confidence GNPS annotations (default: 0.7).
#' @param canopus.prob Numeric. Minimum probability for CANOPUS (NPC Pathway) annotations (default: 0.7).
#' @param csi.prob Numeric. Minimum probability/score for CSI:FingerID annotations (default: 0.64).
#' @param ms2query.prob Numeric. Minimum probability/score for MS2Query annotations (default: 0.7).
#' @param rt.tol Numeric. Retention time tolerance, e.g., for C18 column (default: 0.1).
#' @param cid_database_path Character string. Path to the indexed PubChem SQLite database (default: "Y:/MA_BPA_Microbiome/Databases/PubChem/PubChem_Indexed.sqlite").
#' @return An invisible data frame containing the final, processed, and propagated annotations.
#' @importFrom dplyr filter select mutate
#' @importFrom readr write_csv
#' @export
MAPS <- function(
    dataset.id,
    folder, # Required: Specified Processed Data Folder
    gnps.task.id, # Required: External GNPS ID (can be extracted from metadata)
    gnps.prob = 0.7,
    canopus.prob = 0.7,
    csi.prob = 0.64,
    ms2query.prob = 0.7,
    rt.tol = 0.1, # Min for C18
    cid_database_path = NULL
) {

  # NOTE: All required libraries (dplyr, MAPS.Package, etc.) are assumed to be loaded.

  #-----------------------------------------------------------------------------------------------------------------------#
  ## 1. Data Check (MAPS.Package::validate_and_get_paths must be available)
  if (folder == "") {
    stop("Processed Data Folder 'folder' cannot be empty. Please specify a path or ensure metadata extraction is successful.")
  }

  paths <- MAPS.Package::validate_and_get_paths(folder)

  mzmine.data <- paths$mzmine_data
  mzmine.annotations <- paths$mzmine_annotations
  canopus.data <- paths$canopus_data
  csi.data <- paths$csi_data
  zodiac.data <- paths$zodiac_data
  ms2query.data <- paths$ms2query_data
  cytoscape <- paths$cytoscape

  #-----------------------------------------------------------------------------------------------------------------------#
  ## 3. Process MZMine Data
  message("Processing MZMINE data")
  processed_data <- MAPS.Package::process_mzmine_data(mzmine.annotations, mzmine.data, gnps.prob)

  mzmine.data <- processed_data$processed.data
  mzmine.annotations.final <- processed_data$annotations.data
  sample.data <- processed_data$initial.mzmine.data
  cid_cache_df <- processed_data$cid.cache
  lipids.file <- processed_data$lipids.file


  #-----------------------------------------------------------------------------------------------------------------------#
  ## 4. Process GNPS Data
  message("Processing GNPS data:")
  if (gnps.task.id == "") {
    warning("GNPS Task ID is empty. Skipping GNPS processing.")
  } else {
    gnps.processed.data <- MAPS.Package::process_gnps_data(gnps.task.id, gnps.prob, mzmine.annotations.final)

    gnps.data.lv2.high.conf <- gnps.processed.data$lv2.high.conf.data
    gnps.data.lv2.low.conf <- gnps.processed.data$lv2.low.conf.data
    gnps.data.lv3 <- gnps.processed.data$lv3.data
    gnps.cluster.data <- gnps.processed.data$gnps.cluster.data
    gnps.cluster.pairs <- gnps.processed.data$gnps.cluster.pairs

    gnps_lv2_results <- MAPS.Package::standardise_and_compute_gnps(
      gnps.data.lv2.high.conf, cid_cache_df, lipids.file,
      cid_database_path = cid_database_path,
      ms2query.prob
    )

    gnps.data.lv2.high.conf <- gnps_lv2_results$data
    gnps.data.lv2.high.conf$feature.ID <- as.numeric(gnps.data.lv2.high.conf$feature.ID)
    gnps.data.lv2.high.conf$confidence.score <- as.numeric(gnps.data.lv2.high.conf$confidence.score)
    cid_cache_df <- gnps_lv2_results$cache

    lv1.and.lv2.annotations <- MAPS.Package::merge_and_append_data(
      new_data = gnps.data.lv2.high.conf,
      existing_annotations = mzmine.annotations.final
    )
  }

  # Initialize if skipped (to avoid errors in later steps)
  if (!exists("lv1.and.lv2.annotations")) {
    lv1.and.lv2.annotations <- mzmine.annotations.final
  }


  #-----------------------------------------------------------------------------------------------------------------------#
  ## 5. Process ms2query Data
  message("Processing ms2query data:")
  ms2query_results <- MAPS.Package::process_ms2query_data(
    ms2query.data = paths$ms2query_data,
    lv1.and.lv2.annotations = lv1.and.lv2.annotations,
    ms2query.prob = ms2query.prob
  )

  ms2query.data.lv2 <- ms2query_results$lv2.data
  ms2query.data.lv3 <- ms2query_results$lv3.data

  lv2_processed_results <- MAPS.Package::standardize_and_compute_all_lv2(
    gnps.data.lv2.low.conf = if (exists("gnps.data.lv2.low.conf")) gnps.data.lv2.low.conf else NULL,
    ms2query.data.lv2 = ms2query.data.lv2,
    cid_cache_df = cid_cache_df,
    lipids.file = lipids.file,
    cid_database.path = cid_database_path,
    ms2query.prob = ms2query.prob
  )

  lv2.annotations <- lv2_processed_results$data
  lv2.annotations$feature.ID <- as.numeric(lv2.annotations$feature.ID)
  lv2.annotations$confidence.score <- as.numeric(lv2.annotations$confidence.score)
  cid_cache_df <- lv2_processed_results$cache

  lv1.and.lv2.annotations <- MAPS.Package::merge_and_append_data(
    new_data = lv2.annotations,
    existing_annotations = lv1.and.lv2.annotations
  )

  #-----------------------------------------------------------------------------------------------------------------------#
  ## 6. Process lv3 in silico matches from GNPS
  if (exists("gnps.data.lv3") && nrow(gnps.data.lv3) > 0) {
    unique_in_lv3 <- setdiff(gnps.data.lv3$feature.ID, lv1.and.lv2.annotations$feature.ID)
    gnps.data.lv3 <- dplyr::filter(gnps.data.lv3, feature.ID %in% unique_in_lv3)

    gnps.data.lv3 <- MAPS.Package::compute_id_prob(gnps.data.lv3, "confidence.score", gnps.prob)

    missing_cols <- colnames(lv1.and.lv2.annotations)[!colnames(lv1.and.lv2.annotations) %in% colnames(gnps.data.lv3)]
    for (col in missing_cols) {
      gnps.data.lv3[[col]] <- NA
    }

    gnps.data.lv3 <- gnps.data.lv3 %>% dplyr::select(colnames(lv1.and.lv2.annotations))

    lv1.lv2.lv3.annotations <- lv1.and.lv2.annotations %>%
      rbind(gnps.data.lv3)
  } else {
    lv1.lv2.lv3.annotations <- lv1.and.lv2.annotations
  }

  #-----------------------------------------------------------------------------------------------------------------------#
  ## 7. Process SIRIUS data
  message("Processing SIRIUS data:")
  canopus.data <- MAPS.Package::process_canopus_data(canopus.data)
  zodiac.data <- MAPS.Package::process_zodiac_data(zodiac.data)

  csi_results <- MAPS.Package::process_and_append_csi(
    csi.data, existing_annotations = lv1.lv2.lv3.annotations,
    csi.prob, cid_cache_df, lipids.file,
    cid_database_path = cid_database_path,
    compute_id_prob = MAPS.Package::compute_id_prob,
    deduplicate_data = MAPS.Package::deduplicate_data,
    standardize_annotation = MAPS.Package::standardise_annotation
  )

  lv1.lv2.lv3.annotations <- csi_results$annotations
  lv1.lv2.lv3.annotations$mz.diff.ppm <- as.numeric(lv1.lv2.lv3.annotations$mz.diff.ppm)
  cid_cache_df <- csi_results$cache

  #-----------------------------------------------------------------------------------------------------------------------#
  ## 8. Appending all other features and annotations
  lv1.lv2.lv3.annotations <- MAPS.Package::append_ms2query_analogues(
    ms2query_data = ms2query.data.lv3,
    existing_annotations = lv1.lv2.lv3.annotations
  )

  full.annotation.data <- MAPS.Package::create_full_annotation_table(
    mzmine_data = mzmine.data,
    lv_annotations = lv1.lv2.lv3.annotations,
    canopus_data = canopus.data,
    zodiac_data = zodiac.data,
    gnps_data = if (exists("gnps.cluster.data")) gnps.cluster.data else NULL,
    gnps_task_id = gnps.task.id
  )

  #-----------------------------------------------------------------------------------------------------------------------#
  ## 9. Propagation of annotations
  propagated_df <- MAPS.Package::propagate_annotations(
    full.annotation.data,
    gnps.cluster.pairs = if (exists("gnps.cluster.pairs")) gnps.cluster.pairs else data.frame(),
    paired_feature_finder = MAPS.Package::paired_feature_finder,
    get_result = MAPS.Package::get_result
  )

  propagated.annotation.data <- MAPS.Package::append_propagated_annotations(
    full.annotation.data, propagated_df
  )

  #-----------------------------------------------------------------------------------------------------------------------#
  ## 10. Append level 4 and 5 annotations
  propagated.annotation.data <- MAPS.Package::append_annotations(
    data = propagated.annotation.data,
    mask_condition = is.na(compound.name) & !is.na(canopus.NPC.pathway) & canopus.NPC.pathway.probability >= canopus.prob,
    compound_name = paste0("Predicted NPC Pathway: ", propagated.annotation.data$canopus.NPC.pathway),
    confidence_level = "4",
    npc_pathway = propagated.annotation.data$canopus.NPC.pathway,
    npc_superclass = propagated.annotation.data$canopus.NPC.superclass,
    annotation_type = "canopus"
  )  %>%
    dplyr::mutate(confidence.level = ifelse(is.na(confidence.level), "5", confidence.level))

  #-----------------------------------------------------------------------------------------------------------------------#
  ## 11. Append Sample List Data and Collapse Ion Identity Networking
  processed_results <- MAPS.Package::process_and_append_sample_data(sample.data, propagated.annotation.data)
  propagated.annotation.data.with.samples <- processed_results$combined_data
  sample.data2 <- processed_results$sample.data2

  processed_results <- MAPS.Package::collapse_and_reduce_networks(
    propagated.annotation.data.with.samples, sample.data2,
    process_all_features = MAPS.Package::process_all_features,
    update_data_frame = MAPS.Package::update_data_frame,
    fix_compound_names = MAPS.Package::fix_compound_names,
    redundancy_fixer = MAPS.Package::redundancy_fixer
  )
  final.annotation.df <- processed_results$final_annotation_df
  samples.df <- processed_results$samples_df

  final.annotation.df <- get_hmdb_from_cid(final.annotation.df, cid_database_path)

  #-----------------------------------------------------------------------------------------------------------------------#
  ## 12. Writing Final Files
  message("Writing files to disk: PLEASE WAIT...")
  MAPS.Package::write_final_files(
    final.annotation.df, samples.df, folder, dataset.id, mzmine.data, cid_cache_df,
    write_large_csv = MAPS.Package::write_large_csv
  )

  cytoscape_df <- MAPS.Package::prepare_cytoscape_file(
    folder, propagated.annotation.data.with.samples,
    read_checked_csv = MAPS.Package::read_checked_csv,
    write_large_csv = MAPS.Package::write_large_csv
  )

  tryCatch({
    readr::write_csv(cid_cache_df, "~/MAPS/cid_cache.csv")
    message("Cache updated and saved successfully.")
  }, error = function(e) {
    warning("Failed to save cache: ", e$message)
  })

  closeAllConnections()

  # Return the final annotation data frame
  invisible(final.annotation.df)
}
