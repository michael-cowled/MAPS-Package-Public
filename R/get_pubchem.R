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
get_pubchem_lite <- function(query, type, property = "cids") {

  base_url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

  if (property != "cids") {
    warning("This function is designed to only retrieve CIDs. 'property' parameter will be ignored.")
  }

  if (type == "name") {
    url <- paste0(base_url, "/compound/name/", URLencode(query), "/cids/JSON")
  } else if (type == "smiles") {
    url <- paste0(base_url, "/compound/smiles/", URLencode(query), "/cids/JSON")
  } else if (type == "synonym") {
    # The synonym API is different; it returns a list of information with the CID
    url <- paste0(base_url, "/compound/name/", URLencode(query), "/synonyms/JSON")
    response <- tryCatch(jsonlite::fromJSON(url), error = function(e) {
      message(paste("  [get_pubchem ERROR] Synonym CID lookup failed for '", query, "':", e$message))
      return(NULL)
    })

    if (!is.null(response$InformationList$Information[[1]]$CID)) {
      Sys.sleep(0.2)
      # Extract the first CID from the list
      return(response$InformationList$Information[[1]]$CID[1])
    } else {
      return(NA_real_)
    }
  } else {
    message("Invalid query type.")
    return(NA_real_)
  }

  # For 'name' and 'smiles' lookups
  if (type %in% c("name", "smiles")) {
    response <- tryCatch(jsonlite::fromJSON(url), error = function(e) {
      message(paste("  [get_pubchem ERROR]", type, "CID lookup failed for '", query, "':", e$message))
      return(NULL)
    })

    if (!is.null(response$IdentifierList$CID)) {
      Sys.sleep(0.2)
      return(response$IdentifierList$CID[1])
    } else {
      return(NA_real_)
    }
  }

  return(NA_real_)
}
