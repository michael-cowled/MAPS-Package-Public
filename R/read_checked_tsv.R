#' Read a Checked TSV File
#'
#' This function reads a TSV file from either a local path or a remote URL.
#' For local files, it checks for existence first. It uses `readr` for efficient
#' parsing.
#'
#' @param path A string specifying the path or URL to the TSV file.
#' @param ... Additional arguments to pass to `readr::read_tsv`.
#'
#' @return A tibble containing the data from the TSV file.
#' @export
#' @importFrom readr read_tsv cols
#' @importFrom dplyr select
#'
#' @examples
#' \dontrun{
#'   # Assuming 'data.tsv' exists
#'   my_data <- read_checked_tsv("data.tsv")
#'
#'   # Reading from a remote URL
#'   remote_data <- read_checked_tsv("http://example.com/data.tsv")
#' }
read_checked_tsv <- function(path, ...) {
  if (!grepl("^https?://", path) && !file.exists(path)) stop("File does not exist: ", path)
  readr::read_tsv(path, col_types = readr::cols(.default = "c"), ...)
}
