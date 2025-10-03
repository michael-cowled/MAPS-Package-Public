#' @title Update Compound Names from PubChem
#' @description Queries the PubChem PUG-REST API using CIDs from a specified column
#'   and overwrites the values in a target column with the official PubChem compound title.
#'
#' @param data A data frame or tibble containing compound data.
#' @param name_col A character string specifying the name of the column to be
#'   overwritten with the fetched PubChem title. If the column does not exist, it will be created.
#' @param cid_col A character string specifying the name of the column containing
#'   the PubChem Compound IDs (CIDs). Defaults to "CID".
#'
#' @return The modified data frame with the \code{name_col} updated.
#' @export
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' \dontrun{
#' # Dummy data frame for demonstration
#' compounds_df <- data.frame(
#'   ID = 1:4,
#'   My_CID_Col = c(2244, 440955, 0, NA), # Aspirin, Caffeine, Invalid, NA
#'   Custom_Name = c("Aspirin_Old", "Caffeine_Old", "Zero_CID", "NA_Comp"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Run the function
#' updated_df <- update_compound_names(
#'   data = compounds_df,
#'   name_col = "Custom_Name",
#'   cid_col = "My_CID_Col"
#' )
#' print(updated_df)
#' # Note: The row with CID = 0 and CID = NA will keep their original Custom_Name values.
#' }
update_compound_names <- function(data, name_col, cid_col = "CID") {

  # 1. Input Checks
  if (!cid_col %in% names(data)) {
    stop(paste("CID column '", cid_col, "' not found in data."))
  }

  # If name_col is missing, create it
  if (!name_col %in% names(data)) {
    message(paste("Creating new column '", name_col, "'."))
    data[[name_col]] <- NA_character_
  }

  # Ensure CIDs are numeric for filtering (results in NA for non-numeric CIDs)
  cids <- suppressWarnings(as.numeric(data[[cid_col]]))

  # 2. Identify the rows where the CID is valid (NOT NA AND >= 1)
  valid_indices <- which(!is.na(cids) & cids >= 1)

  # 3. Apply the helper function only to valid CIDs
  if (length(valid_indices) > 0) {

    cids_to_query <- cids[valid_indices]
    num_queries <- length(cids_to_query)

    # Initialize progress bar using utils::txtProgressBar
    message("Starting PubChem lookups...")
    pb <- utils::txtProgressBar(min = 0, max = num_queries, style = 3)

    fetched_titles <- character(num_queries)

    # Iterate through CIDs to query and update progress
    for (i in seq_along(cids_to_query)) {

      # Fetch the title using the internal helper function
      fetched_titles[i] <- fetch_pubchem_title(cids_to_query[i])

      # Update progress bar after each query
      utils::setTxtProgressBar(pb, i)
    }

    # Close progress bar
    close(pb)
    message("PubChem lookups complete.")

    # 4. Overwrite the name_col in the original data frame
    data[[name_col]][valid_indices] <- fetched_titles
  } else {
    message("No valid CIDs (NA or < 1) found to query. Returning original data.")
  }

  return(data)
}
