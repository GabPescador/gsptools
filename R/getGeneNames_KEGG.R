#' Gets gene names from a KEGG pathway ID
#'
#' This function takes a KEGG pathway ID and outputs the gene names associated with them. Requires KEGGREST library.
#'
#' @param pathwayID pathwayID should be either number or with organism specified. e.g., dre03010 where "dre" refers to Danio rerio
#' @return Generates a vector with all your gene name associated to a KEGG pathway ID.
#' @export
# Create a function that gives me gene names based on the KEGG IDs
getGeneNames_KEGG <- function(pathwayID){
  
  # pathwayID should be either number or with organism specified
  # e.g., dre03010 where "dre" refers to Danio rerio
  pathway_genes <- keggGet(pathwayID)
  
  # Extract gene information for the pathway
  genes_info <- character()
  for(i in 1:length(pathway_genes)){
    temp <- pathway_genes[[i]]$GENE # Get gene annotation from pathway
    temp <- temp[seq(2, length(temp), by = 2)] # get only gene names
    temp <- sapply(strsplit(temp, " "), `[`, 1) # Keep only the first word, which is the gene name
    temp <- sub("[^a-zA-Z0-9:]", "", temp) # Get rid of ; from the name
    genes_info <- append(genes_info, temp)
  }
  print(genes_info)
}