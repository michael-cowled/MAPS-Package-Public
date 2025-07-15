#' Standardise Compound Annotations Using PubChem
#'
#' For each row in a data frame, attempts to resolve compound annotations using PubChem.
#' Updates the compound name, SMILES, and adds columns for CID, molecular formula, and monoisotopic mass.
#' Rows containing "candidate" (case-insensitive) in the name column are excluded.
#'
#' @param data A data frame containing compound annotations.
#' @param name_col Name of the column in `data` containing compound names (as a string).
#' @param smiles_col Name of the column in `data` containing SMILES strings (as a string).
#'
#' @return A data frame with standardised annotations and additional PubChem info columns:
#'   \itemize{
#'     \item \code{CID}
#'     \item \code{MolecularFormula}
#'     \item \code{MonoisotopicMass}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' df <- standardise_annotation(data, "compound.name", "smiles")
#' }
standardise_annotation <- function(data, name_col, smiles_col) {
  # Filter out rows with "candidate"
  data <- data[!grepl("candidate", data[[name_col]], ignore.case = TRUE), ]

  n <- nrow(data)
  pb <- txtProgressBar(min = 0, max = n, style = 3)

  # Ensure output columns exist
  data[["CID"]] <- NA
  data[["MolecularFormula"]] <- NA
  data[["MonoisotopicMass"]] <- NA

  for (i in seq_len(n)) {
    current_name <- data[[name_col]][i]
    current_smiles <- data[[smiles_col]][i]

    if (!is.na(current_name) && nzchar(current_name)) {
      message(paste("Processing:", current_name))

      info <- get_cid_with_fallbacks(current_name, current_smiles)

      if (!is.null(info)) {
        data[["CID"]][i] <- info$CID
        if (!is.na(info$ResolvedName) && info$ResolvedName != current_name) {
          info$ResolvedName <- sub(" \\|.*", "", info$ResolvedName)
          data[[name_col]][i] <- info$ResolvedName
          message(paste("  Updated", name_col, ":", current_name, "→", info$ResolvedName))
        }
        data[[smiles_col]][i] <- info$SMILES
        data[["MolecularFormula"]][i] <- info$MolecularFormula
        data[["MonoisotopicMass"]][i] <- info$MonoisotopicMass
      } else {
        message(paste("  No CID found for", current_name))
      }
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)
  return(data)
}
