#' @title Extract GNPS Information from Metadata
#' @description This function reads an Excel metadata file, converts its sheets to
#'     CSV files, and extracts a GNPS task ID and retention time tolerance
#'     from the resulting data.
#'
#' @param excel.metadata.file A character string. The path to the Excel metadata file.
#' @param dataset.id A character string. The ID used to filter the dataset in the metadata.
#' @param output.dir A character string. The directory where CSV files will be saved. Defaults to "HGM".
#'
#' @return A list containing two elements: `gnps.task.id` (a character string) and `rt.tol` (a numeric value).
#'
#' @importFrom readxl excel_sheets read_excel
#' @importFrom readr read_csv write_csv
#' @importFrom dplyr filter
#'
#' @export
#'
extract_gnps_info <- function(excel.metadata.file, dataset.id, output.dir = "HGM") {

  # Create the output directory if it doesn't exist
  if (!dir.exists(output.dir)) {
    dir.create(output.dir)
  }

  # Read all sheets from the Excel file and save them as separate CSVs
  sheet.names <- readxl::excel_sheets(excel.metadata.file)
  for (sheet in sheet.names) {
    data <- readxl::read_excel(excel.metadata.file, sheet = sheet)
    readr::write_csv(data, file = file.path(output.dir, paste0(sheet, ".csv")))
  }

  # Read the "D - Dataset.csv" and filter to find the GNPS task ID
  dataset <- readr::read_csv(file.path(output.dir, "D - Dataset.csv")) %>%
    dplyr::filter(HGMD.ID == dataset.id)

  gnps.task.id <- dataset$gnps.task.ID[1]

  # Determine the retention time tolerance based on column type
  if (!is.na(dataset$column.type[1]) && dataset$column.type[1] == "HILIC") {
    rt.tol <- 0.2
  } else {
    rt.tol <- 0.1
  }

  # Stop if the GNPS task ID is missing
  if (is.na(gnps.task.id)) {
    stop("gnps.task.ID is missing. Add to 'D - Dataset.csv' before re-running.")
  }

  # Return the extracted information
  return(list(gnps.task.id = gnps.task.id, rt.tol = rt.tol))
}
