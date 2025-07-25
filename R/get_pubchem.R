#' Query PubChem for CIDs, Titles, or Properties
#'
#' Queries the PubChem REST API using a compound name, SMILES, or CID
#' to retrieve CIDs, compound titles, or molecular properties.
#'
#' @param query The input string to search for (compound name, SMILES, or CID).
#' @param type A character string: one of \code{"name"}, \code{"smiles"}, \code{"synonym"}, or \code{"cid"}.
#' @param property A character string specifying what to retrieve. Options:
#'   \itemize{
#'     \item \code{"cids"}: to retrieve the compound ID(s)
#'     \item \code{"title"}: to get the PubChem title from a CID
#'     \item \code{"properties"}: to fetch SMILES (MolecularFormula and MonoisotopicMass are no longer retrieved)
#'   }
#'
#' @return Depending on the request type:
#'   \itemize{
#'     \item A numeric CID
#'     \item A character string (compound title)
#'     \item A list with \code{SMILES}
#'     \item Or \code{NA} on failure or no result
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' get_pubchem("glucose", type = "name", property = "cids")
#' get_pubchem("50-99-7", type = "name", property = "cids")
#' get_pubchem("CCO", type = "smiles", property = "cids")
#' get_pubchem("5793", type = "cid", property = "title")
#' get_pubchem("5793", type = "cid", property = "properties")
#' }
get_pubchem <- function(query, type, property = NULL) {
  base_url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

  # Helper function to safely extract text, returning NA if node not found or empty
  safe_xml_text <- function(xml_obj, xpath) {
    # Check if xml_obj is valid before trying to find nodes
    if (is.null(xml_obj) || !inherits(xml_obj, "xml_document")) {
      return(NA_character_)
    }

    node <- xml2::xml_find_first(xml_obj, xpath)
    # Check if node exists (length > 0)
    if (length(node) == 0) {
      return(NA_character_)
    }

    text <- xml2::xml_text(node)
    # Check for empty string result, which can happen if node exists but has no text
    if (length(text) == 0 || nchar(text) == 0) {
      return(NA_character_)
    }
    return(text)
  }

  if (type == "name" && property == "cids") {
    url <- paste0(base_url, "/compound/name/", URLencode(query), "/cids/JSON")
    response <- tryCatch(jsonlite::fromJSON(url), error = function(e) {
      message(paste("  [get_pubchem ERROR] Name CID lookup failed for '", query, "':", e$message))
      return(NULL)
    })
    if (!is.null(response$IdentifierList$CID)) {
      Sys.sleep(0.2)
      return(response$IdentifierList$CID[1])
    } else return(NA_real_) # Return numeric NA for CID

  } else if (type == "smiles" && property == "cids") {
    url <- paste0(base_url, "/compound/smiles/", URLencode(query), "/cids/JSON")
    response <- tryCatch(jsonlite::fromJSON(url), error = function(e) {
      message(paste("  [get_pubchem ERROR] SMILES CID lookup failed for '", query, "':", e$message))
      return(NULL)
    })
    if (!is.null(response$IdentifierList$CID)) {
      Sys.sleep(0.2)
      return(response$IdentifierList$CID[1])
    } else return(NA_real_) # Return numeric NA for CID

  } else if (type == "synonym" && property == "cids") {
    url <- paste0(base_url, "/compound/name/", URLencode(query), "/synonyms/JSON")
    response <- tryCatch(jsonlite::fromJSON(url), error = function(e) {
      message(paste("  [get_pubchem ERROR] Synonym CID lookup failed for '", query, "':", e$message))
      return(NULL)
    })
    if (!is.null(response$InformationList$Information[[1]]$CID)) {
      Sys.sleep(0.2)
      return(response$InformationList$Information[[1]]$CID)
    } else {
      message("    Synonym search returned no CID.")
      return(NA_real_) # Return numeric NA for CID
    }

  } else if (type == "cid" && property == "title") {
    page_url <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", query)
    page <- tryCatch(xml2::read_html(page_url), error = function(e) {
      message(paste("  [get_pubchem ERROR] Failed to retrieve title page for CID", query, ":", e$message))
      return(NULL)
    })
    if (!is.null(page)) {
      title_node <- rvest::html_node(page, "title")
      if (length(title_node) > 0) {
        title <- rvest::html_text(title_node)
        return(gsub(" - PubChem", "", title, fixed = TRUE))
      } else {
        return(NA_character_) # Return character NA if title node not found
      }
    } else return(NA_character_) # Return character NA on page fetch failure

  } else if (type == "cid" && property == "properties") {
    # Only retrieve SMILES, MolecularFormula and MonoisotopicMass are no longer needed
    url <- paste0(base_url, "/compound/cid/", query, "/property/SMILES/XML")
    response <- tryCatch({
      xml <- xml2::read_xml(url)

      # Explicitly check if xml object is NULL or an error after read_xml
      # If so, return NULL from this tryCatch block, which will then be handled by the caller
      if (is.null(xml) || inherits(xml, "try-error")) {
        stop("XML object is NULL or an error after read_xml, forcing NULL return.")
      }

      result <- list(
        SMILES = safe_xml_text(xml, ".//SMILES")
      )
      Sys.sleep(0.2)
      return(result)
    }, error = function(e) {
      message(paste("  [get_pubchem ERROR] Failed to retrieve properties for CID", query, ":", e$message))
      # This return NULL handles cases where read_xml fails completely, or safe_xml_text fails
      return(NULL)
    })
    return(response) # This 'response' can be NULL if the tryCatch block returns NULL

  } else {
    message("Invalid query type or property.")
    return(NA)
  }
}
