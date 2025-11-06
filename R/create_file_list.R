#' @title Create a list of files from a specified folder
#' @description Reads dataset information from a CSV, identifies a target folder,
#'     and generates a CSV file listing all files within that folder.
#'
#' @param dataset.id The ID used to filter the dataset.
#' @param output.directory The path to the directory containing the dataset metadata.
#'     Defaults to "HGM".
#'
#' @return A data frame containing the filtered dataset information.
#'
#' @importFrom readr read_csv write_csv
#' @importFrom dplyr filter
#'
#' @export
#'
create_file_list <- function(dataset.id, output.directory = "HGM") {

  # Construct the path to the dataset CSV file
  dataset.path <- file.path(output.directory, "D - Dataset.csv")

  # Read the dataset and filter based on the provided ID
  dataset.info <- readr::read_csv(dataset.path) %>%
    dplyr::filter(HGMD.ID == dataset.id)

  # Get the output directory and construct the mzML folder path
  output.folder <- dataset.info$Processed.Data.Folder[1]
  mzml.folder <- paste0(output.folder, "\\mzml")

  print(mzml.folder)

  # Check if the folder exists and process the files
  if (dir.exists(mzml.folder)) {
    files.to.process <- list.files(mzml.folder)
    full.paths <- file.path(mzml.folder, files.to.process)
    windows.paths <- gsub("/", "\\\\", full.paths)

    file.summary <- data.frame(
      folder = mzml.folder,
      filename = files.to.process,
      full_path = windows.paths
    )

    print(file.summary)

    # Define the output file name and path
    output.filename <- paste0(dataset.id, "_file-list.csv")
    output.path <- file.path(output.folder, output.filename)

    # Write the file list to the specified path
    readr::write_csv(file.summary, file = output.path)
    print(paste("File list saved to:", output.path))

  } else {
    print("Folder does not exist. Check the path.")
  }

  # Return the dataset object
  return(invisible(dataset.info))
}
