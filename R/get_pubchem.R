library(httr)
library(jsonlite)

#' Query PubChem for CIDs (Lite version)
#'
#' Queries the PubChem REST API using a compound name, SMILES, or synonym
#' to retrieve a single compound CID.
#'
#' @param query The input string to search for (compound name, SMILES, or synonym).
#' @param type A character string: one of "name", "smiles", or "synonym".
#' @param property Always "cids" in this simplified version.
#'
#' @return A numeric CID or NA_real_ on failure or no result.
#' @export
get_pubchem_lite <- function(query, type, property = "cids") {
  base_url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

  if (property != "cids") {
    warning("This function is designed to only retrieve CIDs. 'property' parameter will be ignored.")
  }

  # Name and synonym lookups still use GET
  if (type == "name") {
    url <- paste0(base_url, "/compound/name/", URLencode(query), "/cids/JSON")
    response <- tryCatch(fromJSON(url), error = function(e) {
      message("[get_pubchem ERROR] Name CID lookup failed for '", query, "': ", e$message)
      return(NULL)
    })
    if (!is.null(response$IdentifierList$CID)) return(response$IdentifierList$CID[1])
    return(NA_real_)
  }

  if (type == "synonym") {
    url <- paste0(base_url, "/compound/name/", URLencode(query), "/synonyms/JSON")
    response <- tryCatch(fromJSON(url), error = function(e) {
      message("[get_pubchem ERROR] Synonym CID lookup failed for '", query, "': ", e$message)
      return(NULL)
    })
    if (!is.null(response$InformationList$Information[[1]]$CID)) {
      return(response$InformationList$Information[[1]]$CID[1])
    } else {
      return(NA_real_)
    }
  }

  # SMILES lookups use POST
  if (type == "smiles") {
    url <- paste0(base_url, "/compound/smiles/cids/JSON")
    response <- tryCatch({
      res <- POST(url, body = list(smiles = query), encode = "form")
      if (status_code(res) != 200) stop("HTTP error ", status_code(res))
      fromJSON(content(res, as = "text", encoding = "UTF-8"))
    }, error = function(e) {
      message("[get_pubchem ERROR] SMILES CID lookup failed for '", query, "': ", e$message)
      return(NULL)
    })

    if (!is.null(response$IdentifierList$CID)) return(response$IdentifierList$CID[1])
    return(NA_real_)
  }

  message("Invalid query type.")
  return(NA_real_)
}
