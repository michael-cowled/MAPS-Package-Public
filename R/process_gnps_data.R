#' Process GNPS annotation and cluster data
#'
#' Loads, tidies, and standardizes GNPS2 annotation data, and splits into
#' Level 2 high-confidence, Level 2 low-confidence, and Level 3 annotation tables.
#'
#' @param gnps.task.id GNPS task ID (character)
#' @param gnps.prob Minimum confidence score (numeric)
#' @param mzmine.annotations.final Data frame of MZmine annotations
#'
#' @return A list with three data frames:
#' \describe{
#'   \item{lv2.high.conf.data}{Level 2 high-confidence GNPS annotations}
#'   \item{lv2.low.conf.data}{Level 2 low-confidence GNPS annotations}
#'   \item{lv3.data}{Level 3 GNPS annotations}
#' }
#' @export
process_gnps_data <- function(gnps.task.id, gnps.prob, mzmine.annotations.final) {

  # 1. Load, tidy, and standardize GNPS2 annotation data
  gnps.annotation.data <- read_checked_tsv(paste0("https://gnps2.org/resultfile?task=", gnps.task.id, "&file=nf_output/library/merged_results_with_gnps.tsv"))

  # Select and rename columns for clarity
  gnps.annotation.data <- gnps.annotation.data[, c(2, 4, 5, 8, 9, 15, 27, 35, 43, 45, 46)]
  names(gnps.annotation.data) <- c("feature.ID", "library.name", "confidence.score",
                                   "mz.diff.ppm", "gnps.shared.peaks", "compound.name",
                                   "smiles", "library.quality", "NPC.superclass", "NPC.pathway", "gnps.library.usi")

  # Standardize GNPS library USI format
  gnps.annotation.data$gnps.library.usi <- str_replace(
    gnps.annotation.data$gnps.library.usi,
    "mzspec:GNPS:GNPS-LIBRARY:(.*)",
    "mzspec:GNPS:GNPS-LIBRARY:accession:\\1"
  )

  # Extract bile acid info into a new column
  gnps.annotation.data$gnps.in.silico.bile.acid.info <- NA
  gnps.annotation.data$gnps.in.silico.bile.acid.info[grepl("Candidate ", gnps.annotation.data$compound.name)] <-
    gnps.annotation.data$compound.name[grepl("Candidate ", gnps.annotation.data$compound.name)]

  # Reorder and filter annotations
  gnps.annotation.data <- gnps.annotation.data[, c(1, 6, 7, 3:5, 2, 8, 10, 9, 11, 12)]
  gnps.annotation.data <- gnps.annotation.data %>%
    filter(mz.diff.ppm <= 5) %>%
    filter(confidence.score >= gnps.prob)

  # 2. Load and process GNPS cluster data and pairs
  gnps.cluster.data <- read_checked_tsv(paste0("https://gnps2.org/resultfile?task=", gnps.task.id, "&file=nf_output/networking/clustersummary_with_network.tsv"))
  gnps.cluster.data <- select(gnps.cluster.data, 'cluster index', 'component')
  names(gnps.cluster.data) <- c('feature.ID', "gnps.cluster.ID")

  gnps.cluster.pairs <- read_checked_tsv(paste0("https://gnps2.org/resultfile?task=", gnps.task.id, "&file=nf_output/networking/filtered_pairs.tsv"))

  # 3. Create different confidence levels (Levels 2 and 3)
  gnps.annotation.data$annotation.type <- "gnps"

  # Filter out features that are already in MZMINE Level 1
  unique_in_gnps <- setdiff(gnps.annotation.data$feature.ID, mzmine.annotations.final$feature.ID)
  gnps.annotation.data <- filter(gnps.annotation.data, feature.ID %in% unique_in_gnps)

  # Create level 2 data frames
  gnps.data.lv2 <- filter(gnps.annotation.data, library.quality != "Insilico")
  gnps.data.lv2$confidence.level <- "2"
  # This relies on a user-defined function, assumed to be available
  gnps.data.lv2 <- fix_compound_names(gnps.data.lv2, "compound.name")

  # Separate high and low confidence Level 2 data
  gnps.data.lv2.high.conf <- filter(gnps.data.lv2, library.quality != "Bronze")
  gnps.data.lv2.low.conf <- filter(gnps.data.lv2, library.quality == "Bronze")

  # Filter low confidence data to remove duplicates with high confidence
  unique_in_gnps <- setdiff(gnps.data.lv2.low.conf$feature.ID, gnps.data.lv2.high.conf$feature.ID)
  gnps.data.lv2.low.conf <- filter(gnps.data.lv2.low.conf, feature.ID %in% unique_in_gnps)

  # Create level 3 data frame
  gnps.data.lv3 <- filter(gnps.annotation.data, library.quality == "Insilico")
  gnps.data.lv3$confidence.level <- "3"

  # Clean up compound names for Level 3
  gnps.data.lv3 <- gnps.data.lv3 %>%
    mutate(
      compound.name = ifelse(
        grepl("Candidate ", compound.name),
        sub("\\s*\\(.*$", "", compound.name),
        compound.name
      )
    )

  # Tidy up additional columns
  gnps.data.lv3$gnps.in.silico.bile.acid.info <- gsub('"', "", gnps.data.lv3$gnps.in.silico.bile.acid.info)
  gnps.data.lv3$compound.name <- gsub('"', "", gnps.data.lv3$compound.name)
  gnps.data.lv3$CID <- NA
  gnps.data.lv3$id.prob <- NA

  # Return the processed data frames
  return(list(
    lv2.high.conf.data = gnps.data.lv2.high.conf,
    lv2.low.conf.data = gnps.data.lv2.low.conf,
    lv3.data = gnps.data.lv3,
    gnps.cluster.data = gnps.cluster.data,
    gnps.cluster.pairs = gnps.cluster.pairs
  ))
}
