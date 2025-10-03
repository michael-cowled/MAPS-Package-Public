#' @title Fetch PubChem Title
#' @description Internal helper function to query the PubChem PUG-REST API
#'   for a compound's official title based on its CID.
#'
#' @param cid An integer representing the PubChem Compound ID (CID).
#' @return A character string containing the compound title, or \code{NA_character_}
#'   if the lookup fails, the CID is invalid (\code{< 1}), or the title is not found.
#' @keywords internal
#' @importFrom jsonlite fromJSON
#' @export
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
    # Use jsonlite::fromJSON (imported via roxygen tag)
    data <- jsonlite::fromJSON(url)

    # Navigate the nested JSON structure to get the Title
    title <- data$InformationList$Information$Title

    # Return the title if valid
    if (length(title) > 0 && !is.null(title) && !is.na(title)) {
      return(title[1])
    } else {
      # Handle case where CID is queried but no title is returned
      return(NA_character_)
    }
  }, error = function(e) {
    # Handle network or parsing errors
    warning(paste("Error fetching title for CID", cid, ". Skipping. Message:", e$message),
            call. = FALSE)
    return(NA_character_)
  })
}
