#' Clean and Standardize Compound Names
#'
#' Applies regex-based cleanup and casing corrections to specified column(s) of a data frame.
#' Useful for cleaning up messy annotations from GNPS or similar metabolomics tools.
#'
#' @param df A data frame containing compound names to clean.
#' @param column A character vector of column names in \code{df} to clean.
#'
#' @return A cleaned version of \code{df} with standardized compound name columns.
#' @export
#'
#' @examples
#' \dontrun{
#' cleaned_df <- fix_compound_names(my_annotations, column = "compound.name")
#' }
fix_compound_names <- function(df, column) {
  # Convert to data.table for efficient updates
  df <- data.table::as.data.table(df)

  # Define patterns and replacements
  patterns <- c(
    "\\.(alpha|beta|gamma|omega)\\.-", "\\.(DELTA)\\.",
    "[[:space:]]*\\[IIN-based(?:[[:space:]]+on)?[[:space:]]*:[[:space:]]*[^)]*\\]",
    "NCG[[:alnum:]-]+_[[:alnum:]]+_", "NCG[^\\s!]+[\\s!]?", "Spectral Match to ",
    " from NIST14", "40.0 eV", "70.0 eV", "50.0 eV", ".00 eV", "\"", "CollisionEnergy",
    "\\|.*", "_", "\\[", "\\]", "\\{", "\\}", " :102040", " :205060", " - ",
    "30.0 eV", "Massbank:[[:alnum:]]+\\s? ", "MassbankEU:[[:alnum:]]+\\s? ",
    "MoNA:[[:alnum:]]+\\s? ", " \\(Chimeric precursor selection\\)", "; \\(M\\+.*"
  )
  replacements <- c(
    "\\1-", "delta-", "", "", "", "", "", "", "", "", "", "", "", "", " ", "(", ")", "(", ")",
    "", "", "", "", "", "", "", "", "", ""
  )

  # Apply regex replacements and formatting to each specified column
  for (col in column) {
    for (i in seq_along(patterns)) {
      df[grepl(patterns[i], get(col), ignore.case = TRUE),
         (col) := gsub(patterns[i], replacements[i], get(col), ignore.case = TRUE)]
    }

    # Capitalize if ALL CAPS or ALL LOWERCASE
    df[, (col) := ifelse(
      grepl("^[[:upper:]\\s\\d\\W]*$", get(col)) & grepl("[[:alpha:]]", get(col)),
      paste0(toupper(substring(get(col), 1, 1)), tolower(substring(get(col), 2))),
      get(col)
    )]

    df[, (col) := ifelse(
      grepl("^[[:lower:]\\s\\d\\W]*$", get(col)) & grepl("[[:alpha:]]", get(col)),
      paste0(toupper(substring(get(col), 1, 1)), substring(get(col), 2)),
      get(col)
    )]
  }

  return(as.data.frame(df))
}
