#' Nuclear protein IDs from UNIPROT from humans
#'
#' Dataframe with nuclear protein IDs downloaded from uniprot and GO term GO:00058 on 2025/09/23
#'
#' @format A data frame:
#' \describe{
#'   \item{ProteinID}{Accession ID from uniprot}
#'   \item{Reviewed}{If protein ID is from reviewed list or not}
#'   \item{Entry Name}{Uniprot entry name}
#'   \item{Protein names}{Long format of protein name}
#'   \item{Gene Names}{Short protein names assosicated with ID, might have multiple}
#'   \item{Organism}{Which organism it is from}
#'   \item{Length}{Aminoacid length of protein}
#'   \item{ProteinName}{Gene names taken from the first name of Gene Names column, should be unique}
#'   \item{Ensembl}{Ensembl transcript ID mapped to the protein ID, might have multiple}
#'   ...
#' }
#' @source Downloaded from Uniprot on 2025/09/23.
"nucleusProteinsHs"
