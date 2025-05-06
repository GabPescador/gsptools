#' ERCC Annotation
#'
#' Dataframe with ERCC annotation used by the sequencing core at SIMR. Useful for edgeR library 
#' normalization.
#'
#' @format A data frame with 92 rows and 5 columns:
#' \describe{
#'   \item{ERCC_ID}{ERCC ID}
#'   \item{GenBank}{GenBank ID}
#'   \item{5prime_assay}{5prime}
#'   \item{3prime_assay}{3prime}
#'   \item{Sequence}{Sequence}
#'   ...
#' }
#' @source Downloaded from supplier website: https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095047.txt
"ERCCannotation"