#' E3 ligase proteins and their substrates
#'
#' Dataframe with E3 ligase proteins and their substrates downloaded from UbiBrowser on 2025/09/26.
#'
#' @format A data frame:
#' \describe{
#'   \item{E3_protein_id}{Protein ID from E3 ligase}
#'   \item{E3_protein_name}{Protein name from E3 ligase}
#'   \item{substrate_protein_id}{Protein ID from substrate}
#'   \item{substrate_protein_name}{Protein name from substrate}
#'   \item{SwissProt ID (E3)}{SwissProt ID}
#'   \item{SwissProt ID (Substrate)}{SwissProt ID}
#'   \item{SOURCE}{Source of substrate link}
#'   \item{SOURCEID}{Source ID of substrate link}
#'   \item{E3TYPE}{Type of E3 ligase, e.g., RING domain}
#'   ...
#' }
#' @source Downloaded from UbiBrowser on 2025/09/26.
"E3ligaseSubstrate"
