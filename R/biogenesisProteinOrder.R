#' Ribosome biogenesis factors and groupings in step-wise manner
#'
#' Dataframe with ribosome biogenesis factors from https://doi.org/10.15252/embj.2022112699 in a step-wise manner
#'
#' @format A data frame:
#' \describe{
#'   \item{protein_id}{Protein ID from Uniprot}
#'   \item{protein_name}{Protein name}
#'   \item{Yeat}{Yeast protein name}
#'   \item{Category}{Which step-wise category it is participating}
#'   \item{Type}{Sub-category, if exists}
#'   \item{Function in rRNA transcription}{Functions if annotated}
#'   \item{Citation}{Citation of function description}
#'   ...
#' }
#' @source Inferred from https://doi.org/10.15252/embj.2022112699.
"biogenesisProteinOrder"
