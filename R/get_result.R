#' Get First Valid Compound Annotation for a Feature
#'
#' Extracts the first non-NA and non-"null" value from the `compound.name` column
#' of a feature row, along with the first valid superclass (from multiple columns),
#' or falls back to the `smiles` column if no compound name is found.
#'
#' @param paired_value A single feature ID (value from `feature.ID`) to look up.
#' @param summary_data A data frame containing compound annotation columns,
#'   including `feature.ID`, `compound.name`, `smiles`, `NPC.superclass`, and `canopus.NPC.superclass`.
#'
#' @return A named list with three elements:
#'   \item{value}{The resolved compound name or fallback value}
#'   \item{column}{The column from which the value was taken (`compound.name` or `smiles`)}
#'   \item{superclass}{The first non-NA superclass value found, if any}
#' @export
#'
#' @examples
#' \dontrun{
#' get_result(123, summary_data)
#' }
get_result <- function(paired_value, summary_data) {
  compound_col <- "compound.name"
  superclass_columns <- c("NPC.superclass", "canopus.NPC.superclass")

  # Filter to feature row and select only relevant columns
  data_subset <- summary_data %>%
    dplyr::filter(feature.ID == paired_value) %>%
    dplyr::select(dplyr::all_of(compound_col), dplyr::all_of(superclass_columns), smiles) %>%
    dplyr::slice(1)  # In case of duplicates

  compound_value <- data_subset[[compound_col]]

  if (!is.na(compound_value) && compound_value != "null") {
    # Try to extract a valid superclass
    superclass_value <- NA
    for (col in superclass_columns) {
      sc_val <- data_subset[[col]]
      if (!is.na(sc_val) && sc_val != "null") {
        superclass_value <- sc_val
        break
      }
    }
    return(list(value = compound_value, column = compound_col, superclass = superclass_value))
  }

  # Fallback: use SMILES
  csi_value <- data_subset$smiles
  if (!is.na(csi_value) && csi_value != "null") {
    return(list(value = csi_value, column = "smiles", superclass = NA))
  }

  return(list(value = NA, column = NA, superclass = NA))
}
