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
      fetched_titles[i] <- .fetch_pubchem_title(cids_to_query[i])

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
