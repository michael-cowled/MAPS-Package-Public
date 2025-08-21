#' Query PubChem for CIDs
#'
#' Queries the PubChem REST API using a compound name, SMILES, or synonym
#' to retrieve a single compound CID. This is a simplified version of the
#' original function, designed to only perform CID lookups.
#'
#' @param query The input string to search for (compound name, SMILES, or synonym).
#' @param type A character string: one of \code{"name"}, \code{"smiles"}, or \code{"synonym"}.
#' @param property Always "cids" in this simplified version.
#'
#' @return A numeric CID or \code{NA_real_} on failure or no result.
#' @export
#' @import httr
#' @import jsonlite
get_pubchem_lite <- function(query, type, property = "cids") {

  # Base URL for the PubChem PUG REST API
  base_url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

  # Warn if the property parameter is not "cids" as this function is simplified
  if (property != "cids") {
    warning("This function is designed to only retrieve CIDs. 'property' parameter will be ignored.")
  }

  # --- 1. Validate the query type ---
  if (!type %in% c("name", "smiles", "synonym")) {
    message("Invalid query type. Must be 'name', 'smiles', or 'synonym'.")
    return(NA_real_)
  }

  # --- 2. Construct the request body and URL based on type ---
  # Using POST for all types to handle special characters and standardize API calls
  body_data <- list(
    "query" = query
  )

  # Determine the specific API endpoint based on the query type
  if (type %in% c("name", "smiles")) {
    url <- paste0(base_url, "/compound/", type, "/cids/JSON")
  } else if (type == "synonym") {
    # Synonym search uses the /name/ endpoint but with a different property
    url <- paste0(base_url, "/compound/name/", "/synonyms/JSON")
    body_data <- list("name" = query) # Synonym search is a name lookup
  }

  # --- 3. Execute the API request with error handling ---
  response <- tryCatch({
    httr::POST(url, body = body_data, encode = "form")
  }, error = function(e) {
    message(paste("  [get_pubchem ERROR]", type, "CID lookup failed for '", query, "':", e$message))
    return(NULL)
  })

  # Check for connection errors or non-200 HTTP status codes
  if (is.null(response) || httr::http_error(response)) {
    return(NA_real_)
  }

  # --- 4. Parse the JSON response and extract the CID ---
  response_json <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  Sys.sleep(0.2) # Be polite to the API

  # Extract CID based on the query type
  if (type %in% c("name", "smiles")) {
    if (!is.null(response_json$IdentifierList$CID)) {
      return(response_json$IdentifierList$CID[1])
    }
  } else if (type == "synonym") {
    if (!is.null(response_json$InformationList$Information[[1]]$CID)) {
      return(response_json$InformationList$Information[[1]]$CID[1])
    }
  }

  # Return NA if no CID was found
  return(NA_real_)
}
