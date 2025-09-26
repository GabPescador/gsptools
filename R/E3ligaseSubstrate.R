#' E3 ligase proteins and their substrates
#'
#' Dataframe with E3 ligase proteins and their substrates downloaded from UbiBrowser on 2025/09/26.
#'
#' @format A data frame:
#' \describe{
#'   \item{ProteinName}{Protein name}
#'   \item{SwissProt ID (E3)}{SwissProt ID}
#'   \item{SwissProt ID (Substrate)}{SwissProt ID}
#'   \item{E3_ProteinID}{Protein ID from E3 ligase}
#'   \item{Substrate_ProteinID}{Protein ID from subtrate}
#'   \item{Substrate_ProteinName}{Protein name from substrate}
#'   \item{SOURCE}{Source of substrate link}
#'   \item{SOURCEID}{Source ID of substrate link}
#'   \item{E3TYPE}{Type of E3 ligase, e.g., RING domain}
#'   ...
#' }
#' @source Downloaded from UbiBrowser on 2025/09/26.
"E3ligaseSubstrate"
