#' Ribosome biogenesis factors and groupings in step-wise manner
#'
#' Dataframe with ribosome biogenesis factors from https://doi.org/10.15252/embj.2022112699 in a step-wise manner
#'
#' @format A data frame:
#' \describe{
#'   \item{protein_id}{Protein ID from Uniprot}
#'   \item{protein_name}{Protein name}
#'   \item{Dr_protein_id}{Protein ID from Uniprot for zebrafish matched by OrthoFinder}
#'   \item{Dr_protein_name}{Protein name from Uniprot when available, if it was NA then uses the human protein name in lower case. If Dr_protein_id is NA, it means it is using the human protein name in its place}
#'   \item{Yeat}{Yeast protein name}
#'   \item{Category}{Which step-wise category it is participating}
#'   \item{Type}{Sub-category, if exists}
#'   \item{Function in rRNA transcription}{Functions if annotated}
#'   \item{Citation}{Citation of function description}
#'   ...
#' }
#' @source Inferred from https://doi.org/10.15252/embj.2022112699 and matched with OrthoFinder by Eric Ross (Stowers CompBio).
"biogenesisProteinOrder"
