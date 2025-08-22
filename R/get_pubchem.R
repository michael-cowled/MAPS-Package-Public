#' Query PubChem for CIDs
#'
#' Queries the PubChem REST API using a compound name, SMILES, or synonym
#' to retrieve a single compound CID.
#'
#' @param query The input string to search for (compound name, SMILES, or synonym).
#' @param type A character string: one of \code{"name"}, \code{"smiles"}, or \code{"synonym"}.
#'
#' @return A numeric CID or \code{NA_real_} on failure or no result.
#' @export
#' @import httr
#' @import jsonlite
get_pubchem_lite <- function(query, type) {

  # Base URL for the PubChem PUG REST API
  base_url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

  # --- 1. Validate the query type ---
  if (!type %in% c("name", "smiles", "synonym")) {
    message("Invalid query type. Must be 'name', 'smiles', or 'synonym'.")
    return(NA_real_)
  }

  # --- 2. Construct and execute the API request ---
  if (type %in% c("name", "smiles")) {
    # Use POST for name and smiles queries to handle special characters
    url <- paste0(base_url, "/compound/", type, "/cids/JSON")
    body_data <- list("query" = query)

    response <- tryCatch({
      httr::POST(url, body = body_data, encode = "form")
    }, error = function(e) {
      message(paste("  [get_pubchem ERROR] POST failed for '", query, "':", e$message))
      return(NULL)
    })

    # Check for errors and parse the response
    if (is.null(response) || httr::http_error(response)) {
      return(NA_real_)
    }
    response_json <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))

    if (!is.null(response_json$IdentifierList$CID)) {
      Sys.sleep(0.2)
      return(response_json$IdentifierList$CID[1])
    }
  } else if (type == "synonym") {
    # Use a GET request for synonym searches as per API design
    url <- paste0(base_url, "/compound/name/", URLencode(query, reserved = TRUE), "/synonyms/JSON")

    # Check for errors and parse the response
    if (is.null(response) || httr::http_error(response)) {
      return(NA_real_)
    }

    response_json <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))

    # The fix: Check if the response is a list and contains the expected structure
    if (is.list(response_json) &&
        !is.null(response_json$InformationList) &&
        is.list(response_json$InformationList) &&
        length(response_json$InformationList$Information) > 0 &&
        is.list(response_json$InformationList$Information[[1]]) &&
        !is.null(response_json$InformationList$Information[[1]]$CID)) {

      Sys.sleep(0.2)
      return(response_json$InformationList$Information[[1]]$CID[1])
    }
  }

  # Return NA if no CID was found
  return(NA_real_)
}
