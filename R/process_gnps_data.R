#' Process GNPS annotation and cluster data
#'
#' Loads, tidies, and standardizes GNPS2 annotation data, and splits into
#' Level 2 high-confidence, Level 2 low-confidence, and Level 3 annotation tables.
#'
#' @param gnps.task.id GNPS task ID (character)
#' @param gnps.prob Minimum confidence score (numeric)
#' @param mzmine.annotations.final Data frame of MZmine annotations
#'
#' @return A list with five data frames including levels 2 and 3 annotations and networking data.
#' @export
process_gnps_data <- function(gnps.task.id, gnps.prob, mzmine.annotations.final) {

  # 1. Load and Standardize GNPS2 annotation data
  gnps.url <- paste0("https://gnps2.org/resultfile?task=", gnps.task.id, "&file=nf_output/library/merged_results_with_gnps.tsv")
  gnps.annotation.data <- read_checked_tsv(gnps.url)

  # Select, Rename, and Initial Clean in one pipe
  gnps.annotation.data <- gnps.annotation.data %>%
    dplyr::select(
      feature.ID = `#Scan#`,
      compound.name = Compound_Name,
      smiles = Smiles,
      confidence.score = MQScore,
      mz.diff.ppm = MZErrorPPM,
      gnps.shared.peaks = SharedPeaks,
      library.name = LibraryName,
      library.quality = LibraryQualityString,
      NPC.superclass = npclassifier_superclass, # matched to your previous request
      NPC.pathway = npclassifier_pathway,
      gnps.library.usi = library_usi
    ) %>%
    dplyr::mutate(
      # Standardize USI format
      gnps.library.usi = stringr::str_replace(gnps.library.usi, "mzspec:GNPS:GNPS-LIBRARY:(.*)", "mzspec:GNPS:GNPS-LIBRARY:accession:\\1"),
      # Extract bile acid info
      gnps.in.silico.bile.acid.info = ifelse(grepl("Candidate ", compound.name), compound.name, NA_character_),
      annotation.type = "gnps",
      confidence.score = as.numeric(confidence.score),
      mz.diff.ppm = as.numeric(mz.diff.ppm)
    ) %>%
    # Initial Filter
    dplyr::filter(mz.diff.ppm <= 5, confidence.score >= gnps.prob)

  # 2. Load and process GNPS cluster data and pairs
  gnps.cluster.url <- paste0("https://gnps2.org/resultfile?task=", gnps.task.id, "&file=nf_output/networking/clustersummary_with_network.tsv")
  gnps.cluster.data <- read_checked_tsv(gnps.cluster.url) %>%
    dplyr::select(feature.ID = `cluster index`, gnps.cluster.ID = component)

  gnps.pairs.url <- paste0("https://gnps2.org/resultfile?task=", gnps.task.id, "&file=nf_output/networking/filtered_pairs.tsv")
  gnps.cluster.pairs <- read_checked_tsv(gnps.pairs.url)

  # 3. Confidence Level Splitting

  # Remove features already identified in Level 1 (MZmine)
  gnps.annotation.data <- gnps.annotation.data %>%
    dplyr::filter(!(feature.ID %in% mzmine.annotations.final$feature.ID))

  # --- Level 2 (Library Matches) ---
  gnps.data.lv2 <- gnps.annotation.data %>%
    dplyr::filter(library.quality != "Insilico") %>%
    dplyr::mutate(confidence.level = "2") %>%
    fix_compound_names("compound.name")

  # Split Lv2 into High (Gold/Silver) and Low (Bronze)
  gnps.data.lv2.high.conf <- gnps.data.lv2 %>% dplyr::filter(library.quality != "Bronze")
  gnps.data.lv2.low.conf  <- gnps.data.lv2 %>%
    dplyr::filter(library.quality == "Bronze") %>%
    dplyr::filter(!(feature.ID %in% gnps.data.lv2.high.conf$feature.ID))

  # --- Level 3 (In-silico / Putative) ---
  gnps.data.lv3 <- gnps.annotation.data %>%
    dplyr::filter(library.quality == "Insilico") %>%
    dplyr::mutate(
      confidence.level = "3",
      # Clean Candidate names
      compound.name = ifelse(grepl("Candidate ", compound.name), sub("\\s*\\(.*$", "", compound.name), compound.name),
      # Remove quotes and add placeholders
      compound.name = gsub('"', "", compound.name),
      gnps.in.silico.bile.acid.info = gsub('"', "", gnps.in.silico.bile.acid.info),
      CID = NA,
      id.prob = NA
    )

  return(list(
    lv2.high.conf.data = gnps.data.lv2.high.conf,
    lv2.low.conf.data = gnps.data.lv2.low.conf,
    lv3.data = gnps.data.lv3,
    gnps.cluster.data = gnps.cluster.data,
    gnps.cluster.pairs = gnps.cluster.pairs
  ))
}
