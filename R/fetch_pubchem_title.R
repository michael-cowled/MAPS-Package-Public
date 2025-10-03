fetch_pubchem_title <- function(cid) {

  # 1. Input Validation: Skip if NA, NULL, or not a positive integer
  if (is.na(cid) || is.null(cid) || suppressWarnings(as.numeric(cid)) < 1) {
    return(NA_character_)
  }

  # Ensure CID is converted to character for URL construction
  cid <- as.character(cid)

  # 2. Construct the PUG-REST URL
  url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                cid, "/description/JSON")

  # 3. Fetch and Parse with Error Handling
  tryCatch({

    # **MESSAGE FOR LOOKUP**
    message(paste("Querying PubChem for CID:", cid))

    # Use jsonlite::fromJSON (imported via roxygen tag)
    data <- jsonlite::fromJSON(url)

    # Navigate the nested JSON structure to get the Title
    title <- data$InformationList$Information$Title

    # 4. Implement a 0.2 second delay to respect PubChem API usage policy
    Sys.sleep(0.2)

    # Return the title if valid
    if (length(title) > 0 && !is.null(title) && !is.na(title)) {
      return(title[1])
    } else {
      # Handle case where CID is queried but no title is returned
      return(NA_character_)
    }
  }, error = function(e) {
    # Implement a delay even on error to avoid rapid retries
    Sys.sleep(0.2)
    # Handle network or parsing errors
    warning(paste("Error fetching title for CID", cid, ". Skipping. Message:", e$message),
            call. = FALSE)
    return(NA_character_)
  })
}
