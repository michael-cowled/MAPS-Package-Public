# get_pubchem.R
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
#'     \item \code{"properties"}: to fetch SMILES
#'   }
#' @param db_con An active DBI connection object to your SQLite database, or NULL if not available.
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
#' # Assuming a db_con is established for properties lookup
#' # get_pubchem("glucose", type = "name", property = "cids")
#' # get_pubchem("50-99-7", type = "name", property = "cids")
#' # get_pubchem("CCO", type = "smiles", property = "cids")
#' # get_pubchem("5793", type = "cid", property = "title")
#' # get_pubchem("5793", type = "cid", property = "properties", db_con = my_db_connection)
#' }
get_pubchem <- function(query, type, property = NULL, db_con = NULL) {
  base_url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

  # Early return for invalid CID input
  if (type == "cid" && (!suppressWarnings(!is.na(as.numeric(query))) || as.numeric(query) <= 0)) {
    message(paste("[get_pubchem SKIP] Non-positive or invalid CID:", query, "- skipping lookup."))
    return(NULL)
  }

  # Helper function to safely extract text, returning NA if node not found or empty
  safe_xml_text <- function(xml_obj, xpath) {
    if (is.null(xml_obj) || !inherits(xml_obj, "xml_document")) {
      return(NA_character_)
    }

    node <- xml2::xml_find_first(xml_obj, xpath)
    if (length(node) == 0) {
      return(NA_character_)
    }

    text <- xml2::xml_text(node)
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
    } else return(NA_real_)

  } else if (type == "smiles" && property == "cids") {
    url <- paste0(base_url, "/compound/smiles/", URLencode(query), "/cids/JSON")
    response <- tryCatch(jsonlite::fromJSON(url), error = function(e) {
      message(paste("  [get_pubchem ERROR] SMILES CID lookup failed for '", query, "':", e$message))
      return(NULL)
    })
    if (!is.null(response$IdentifierList$CID)) {
      Sys.sleep(0.2)
      return(response$IdentifierList$CID[1])
    } else return(NA_real_)

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
      return(NA_real_)
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
        return(NA_character_)
      }
    } else return(NA_character_)

  } else if (type == "cid" && property == "properties") {
    # Check if a database connection is provided and valid
    if (!is.null(db_con) && DBI::dbIsValid(db_con)) {
      tryCatch({
        query_db <- sprintf("SELECT SMILES FROM pubchem_data WHERE CID = %s", sQuote(query))
        message(paste("  [get_pubchem DEBUG] Executing DB query:", query_db))
        db_result <- DBI::dbGetQuery(db_con, query_db)

        if (nrow(db_result) > 0 && !is.na(db_result$SMILES[1])) {
          message(paste("  [get_pubchem INFO] Retrieved SMILES from local DB for CID", query))
          return(list(SMILES = db_result$SMILES[1]))
        } else {
          message(paste("  [get_pubchem INFO] SMILES not found in local DB for CID", query, ". Falling back to PubChem API."))
        }
      }, error = function(e) {
        message(paste("  [get_pubchem ERROR] Database lookup failed for CID", query, ":", e$message))
      })
    } else {
      message(paste("  [get_pubchem INFO] No valid database connection for properties lookup. Querying PubChem API for CID", query))
    }

    # Fallback to PubChem API if no database, DB query failed, or DB didn't contain the info
    # Ensure query is not an invalid CID like 0 or -1 before proceeding
    if (!is.na(as.numeric(query)) && as.numeric(query) > 0) {
      url <- paste0(base_url, "/compound/cid/", query, "/property/IsomericSMILES/XML")
      xml_response <- tryCatch({
        xml2::read_xml(url)
      }, error = function(e) {
        message(paste("  [get_pubchem ERROR] Failed to retrieve properties from PubChem API for CID", query, ":", e$message))
        return(NULL)
      })

      if (!is.null(xml_response) && !inherits(xml_response, "try-error")) {
        smiles_from_pubchem <- safe_xml_text(xml_response, ".//IsomericSMILES")
        Sys.sleep(0.2)
        return(list(SMILES = smiles_from_pubchem))
      }
    } else {
      message("  [get_pubchem INFO] Invalid CID for properties lookup. Skipping API call.")
    }

    return(NULL)
  } else {
    message("Invalid query type or property.")
    return(NA)
  }
}
