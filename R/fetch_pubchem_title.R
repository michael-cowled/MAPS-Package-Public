#' @title Fetch PubChem Title
#' @description Internal helper function to query the PubChem PUG-REST API
#'   for a compound's official title based on its CID, using XML parsing.
#'
#' @param cid An integer representing the PubChem Compound ID (CID).
#' @return A character string containing the compound title, or \code{NA_character_}
#'   if the lookup fails, the CID is invalid (\code{< 1}), or the title is not found.
#' @keywords internal
#' @importFrom xml2 read_xml xml_find_first xml_text
#' @importFrom httr GET status_code content
#' @export
fetch_pubchem_title <- function(cid) {

  # 1. Input Validation: Skip if NA, NULL, or not a positive integer
  if (is.na(cid) || is.null(cid) || suppressWarnings(as.numeric(cid)) < 1) {
    return(NA_character_)
  }

  # Ensure CID is converted to character for URL construction
  cid <- as.character(cid)

  # 2. Construct the PUG-REST URL, requesting XML format
  url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                cid, "/description/XML")

  # 3. Fetch and Parse with Error Handling
  tryCatch({

    message(paste("Querying PubChem for CID:", cid))

    # Use httr::GET to fetch the response reliably
    response <- httr::GET(url)

    # 3a. Check HTTP status code explicitly
    if (httr::status_code(response) != 200) {
      warning(paste("HTTP Error:", httr::status_code(response),
                    "received for CID", cid, ". Skipping."))
      Sys.sleep(0.2)
      return(NA_character_)
    }

    # 3b. Extract content and parse the XML
    xml_content <- httr::content(response, as = "text", encoding = "UTF-8")
    xml_doc <- xml2::read_xml(xml_content)
    print(xml_doc)

    # Use XPath to find the first occurrence of the <Title> tag
    # FIX: Using the local-name() XPath function to ignore the XML namespace,
    # which caused the previous parsing failure.
    title_node <- xml2::xml_find_first(xml_doc, "//*[local-name()='Title']")
    print(title_node)

    # Extract the text content from the node
    retrieved_title <- xml2::xml_text(title_node)
    print(retrieved_title)

    # 4. Implement a 0.2 second delay to respect PubChem API usage policy
    Sys.sleep(0.2)

    # Check if a valid title was retrieved (not empty or NA)
    if (length(retrieved_title) == 1 && retrieved_title != "") {
      message(paste("Successfully retrieved title for CID", cid, ":", retrieved_title))
      return(retrieved_title)
    } else {
      # Handle case where CID is queried but no title is returned
      message(paste("No title found for CID:", cid))
      return(NA_character_)
    }

  }, error = function(e) {
    # Implement a delay even on error to avoid rapid retries
    Sys.sleep(0.2)
    # Handle network or parsing errors
    warning(paste("Network/Parsing Error fetching title for CID", cid, ". Skipping. Message:", e$message),
            call. = FALSE)
    return(NA_character_)
  })
}
