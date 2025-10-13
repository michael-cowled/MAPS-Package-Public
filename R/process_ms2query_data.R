#' Process MS2Query Data and Split by Confidence Level
#'
#' Reads raw MS2Query output, cleans and filters the annotations,
#' and splits them into Level 2 (high-confidence) and Level 3 (lower-confidence)
#' annotation tables based on confidence scores and m/z differences.
#'
#' @param ms2query.data.path Character. Path to the raw MS2Query CSV file.
#' @param lv1.and.lv2.annotations Data frame of existing Level 1 and 2 annotations,
#'   used to filter out duplicate features.
#' @param ms2query.prob Numeric. The probability threshold for keeping annotations.
#'
#' @return A list with two data frames:
#' \describe{
#'   \item{lv2.data}{Level 2 high-confidence MS2Query annotations}
#'   \item{lv3.data}{Level 3 lower-confidence MS2Query annotations}
#' }
#'
#' @importFrom readr read_csv
#' @importFrom dplyr %>% mutate filter select
#'
#' @export
process_ms2query_data <- function(ms2query.data.path, lv1.and.lv2.annotations, ms2query.prob) {
  # Read raw MS2Query results
  ms2query.data <- read_checked_csv(ms2query.data.path)

  # Select and rename relevant columns
  ms2query.data <- ms2query.data[, c(2, 3, 4, 7, 8, 10, 17, 18)]
  names(ms2query.data) <- c(
    "confidence.score", "mz.diff", "precursor_mz",
    "compound.name", "smiles", "feature.ID",
    "NPC.superclass", "NPC.pathway"
  )

  # Standardise compound names
  ms2query.data <- fix_compound_names(ms2query.data, "compound.name")
  ms2query.data$mz.diff.ppm <- as.numeric(ms2query.data$mz.diff.ppm)

  # Calculate ppm difference and filter by tolerance
  ms2query.data <- ms2query.data %>%
    mutate(mz.diff.ppm = mz.diff / precursor_mz * 1e6) %>%
    filter(mz.diff.ppm <= 5)

  # Mark annotation source
  ms2query.data$annotation.type <- "ms2query"

  # Exclude features already in existing Level 1/2 annotations
  unique_in_ms2query <- setdiff(ms2query.data$feature.ID, lv1.and.lv2.annotations$feature.ID)
  ms2query.data <- filter(ms2query.data, feature.ID %in% unique_in_ms2query)

  # -------------------------
  # Create Level 2 annotations
  # -------------------------
  ms2query.data.lv2 <- ms2query.data %>%
    filter(confidence.score > ms2query.prob) %>%
    filter(mz.diff <= 0.001) %>%
    select(-mz.diff, -precursor_mz)

  ms2query.data.lv2$confidence.level <- "2"
  ms2query.data.lv2$gnps.shared.peaks <- NA
  ms2query.data.lv2$library.name <- "ms2query"
  ms2query.data.lv2$library.quality <- "ms2query"
  ms2query.data.lv2$gnps.library.usi <- NA
  ms2query.data.lv2$gnps.in.silico.bile.acid.info <- NA

  # -------------------------
  # Create Level 3 annotations
  # -------------------------
  ms2query.data.lv3 <- ms2query.data %>%
    filter(confidence.score > ms2query.prob) %>%
    filter(mz.diff >= 0.001) %>%
    filter(!feature.ID %in% lv1.and.lv2.annotations$feature.ID)

  ms2query.data.lv3$confidence.level <- "3"

  # Return both levels
  return(list(
    lv2.data = ms2query.data.lv2,
    lv3.data = ms2query.data.lv3
  ))
}
