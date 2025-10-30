#' @title Create Full Annotation Table
#' @description Joins all annotation data frames into a single, comprehensive table.
#'
#' @param mzmine_data The initial MzMine data frame.
#' @param lv_annotations The consolidated lv1, lv2, and lv3 annotations.
#' @param canopus_data The CANOPUS annotation data.
#' @param zodiac_data The ZODIAC annotation data.
#' @param gnps_data The GNPS cluster data.
#' @param gnps_task_id The GNPS task ID for creating USIs.
#'
#' @return A single data frame with all annotation sources joined.
#' @export
create_full_annotation_table <- function(mzmine_data, lv_annotations, canopus_data, zodiac_data, gnps_data, gnps_task_id) {
  print("check3")
   # Convert feature.ID to numeric for all data frames before joining
  lv_annotations$feature.ID <- as.numeric(lv_annotations$feature.ID)
  canopus_data$feature.ID <- as.numeric(canopus_data$feature.ID)
  zodiac_data$feature.ID <- as.numeric(zodiac_data$feature.ID)
  gnps_data$feature.ID <- as.numeric(gnps_data$feature.ID)
  print("check4")
  print("mzmine")
  print(head(mzmine_data))
  print("lv annotations")
  print(head(lv_annotations))
  print("canopus")
  print(head(canopus_data))
  print("zodiac")
  print(head(zodiac_data))
  print("gnps")
  print(head(gnps_data))
  # Perform full joins to combine all data
  full.annotation.data <- mzmine_data[, 1:4] %>%
    dplyr::full_join(lv_annotations, by = "feature.ID") %>%
    dplyr::full_join(canopus_data, by = "feature.ID") %>%
    dplyr::full_join(zodiac_data, by = "feature.ID") %>%
    dplyr::full_join(gnps_data, by = "feature.ID")
  print("check5")
  # Add ion.identity.ID column and generate USIs
  full.annotation.data$ion.identity.ID <- NA_character_
  full.annotation.data$feature.usi <- paste0("mzspec:GNPS2:TASK-", gnps_task_id, "-nf_output/clustering/spectra_reformatted.mgf:scan:", full.annotation.data$feature.ID)
  print("check6")
  # Tidy up compound names
  full.annotation.data$compound.name[grepl("CID", full.annotation.data$compound.name)] <-
    full.annotation.data$IUPAC[grepl("CID", full.annotation.data$compound.name)]
  print("check7")
  return(full.annotation.data)
}
