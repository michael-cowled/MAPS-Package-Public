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
#'     \item \code{"properties"}: to fetch SMILES, molecular formula, and mass
#'   }
#'
#' @return Depending on the request type:
#'   \itemize{
#'     \item A numeric CID
#'     \item A character string (compound title)
#'     \item A list with \code{SMILES}, \code{MolecularFormula}, and \code{MonoisotopicMass}
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

  if (type == "name" && property == "cids") {
    url <- paste0(base_url, "/compound/name/", URLencode(query), "/cids/JSON")
    response <- tryCatch(jsonlite::fromJSON(url), error = function(e) return(NULL))
    if (!is.null(response$IdentifierList$CID)) {
      Sys.sleep(0.2)
      return(response$IdentifierList$CID[1])
    } else return(NA)

  } else if (type == "smiles" && property == "cids") {
    url <- paste0(base_url, "/compound/smiles/", URLencode(query), "/cids/JSON")
    response <- tryCatch(jsonlite::fromJSON(url), error = function(e) return(NULL))
    if (!is.null(response$IdentifierList$CID)) {
      Sys.sleep(0.2)
      return(response$IdentifierList$CID[1])
    } else return(NA)

  } else if (type == "synonym" && property == "cids") {
    url <- paste0(base_url, "/compound/name/", URLencode(query), "/synonyms/JSON")
    response <- tryCatch(jsonlite::fromJSON(url), error = function(e) return(NULL))
    if (!is.null(response$InformationList$Information[[1]]$CID)) {
      Sys.sleep(0.2)
      return(response$InformationList$Information[[1]]$CID)
    } else {
      message("    Synonym search returned no CID.")
      return(NA)
    }

  } else if (type == "cid" && property == "title") {
    page_url <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", query)
    page <- tryCatch(xml2::read_html(page_url), error = function(e) return(NULL))
    if (!is.null(page)) {
      title <- rvest::html_text(rvest::html_node(page, "title"))
      return(gsub(" - PubChem", "", title, fixed = TRUE))
    } else return(NA)

  } else if (type == "cid" && property == "properties") {
    url <- paste0(base_url, "/compound/cid/", query, "/property/SMILES,MolecularFormula,MonoisotopicMass/XML")
    response <- tryCatch({
      xml <- xml2::read_xml(url)
      result <- list(
        SMILES = xml2::xml_text(xml2::xml_find_first(xml, ".//SMILES")),
        MolecularFormula = xml2::xml_text(xml2::xml_find_first(xml, ".//MolecularFormula")),
        MonoisotopicMass = as.numeric(xml2::xml_text(xml2::xml_find_first(xml, ".//MonoisotopicMass")))
      )
      Sys.sleep(0.2)
      return(result)
    }, error = function(e) return(NULL))
    return(response)

  } else {
    message("Invalid query type or property.")
    return(NA)
  }
}
