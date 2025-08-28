#' Read a Checked CSV File
#'
#' This function reads a CSV file, first checking for file existence and ensuring
#' the file is not empty. If checks pass, it reads the data into a data frame.
#'
#' @param path A string specifying the path to the CSV file.
#' @param ... Additional arguments to pass to `read.csv`.
#'
#' @return A data frame containing the data from the CSV file.
#' @export
#' @importFrom utils read.csv
#'
#' @examples
#' \dontrun{
#'   # Assuming 'data.csv' exists and is not empty
#'   my_data <- read_checked_csv("data.csv")
#' }
read_checked_csv <- function(path, ...) {
  if (!file.exists(path)) stop("File does not exist: ", path)
  if (file.info(path)$size == 0) stop("File exists but is empty: ", path)
  read.csv(path, stringsAsFactors = FALSE, ...)
}
