## 3. get_pubchem
get_pubchem <- function(query, type, property = NULL) {
  base_url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
  
  if (type == "name" && property == "cids") {
    url <- paste0(base_url, "/compound/name/", URLencode(query), "/cids/JSON")
    response <- tryCatch(fromJSON(url), error = function(e) return(NULL))
    if (!is.null(response$IdentifierList$CID)) {
      Sys.sleep(0.2)
      return(response$IdentifierList$CID[1])
    } else return(NA)
    
  } else if (type == "smiles" && property == "cids") {
    url <- paste0(base_url, "/compound/smiles/", URLencode(query), "/cids/JSON")
    response <- tryCatch(fromJSON(url), error = function(e) return(NULL))
    if (!is.null(response$IdentifierList$CID)) {
      Sys.sleep(0.2)
      return(response$IdentifierList$CID[1])
    } else return(NA)
    
  } else if (type == "synonym" && property == "cids") {
    url <- paste0(base_url, "/compound/name/", URLencode(query), "/synonyms/JSON")
    response <- tryCatch(fromJSON(url), error = function(e) return(NULL))
    if (!is.null(response$InformationList$Information[[1]]$CID)) {
      Sys.sleep(0.2)
      return(response$InformationList$Information[[1]]$CID)
    } else {
      message("    Synonym search returned no CID.")
      return(NA)
    }
    
  } else if (type == "cid" && property == "title") {
    page_url <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", query)
    page <- tryCatch(read_html(page_url), error = function(e) return(NULL))
    if (!is.null(page)) {
      title <- html_text(html_node(page, "title"))
      return(gsub(" - PubChem", "", title, fixed = TRUE))
    } else return(NA)
    
  } else if (type == "cid" && property == "properties") {
    url <- paste0(base_url, "/compound/cid/", query, "/property/SMILES,MolecularFormula,MonoisotopicMass/XML")
    response <- tryCatch({
      xml <- read_xml(url)
      result <- list(
        SMILES = xml_text(xml_find_first(xml, ".//SMILES")),
        MolecularFormula = xml_text(xml_find_first(xml, ".//MolecularFormula")),
        MonoisotopicMass = as.numeric(xml_text(xml_find_first(xml, ".//MonoisotopicMass")))
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