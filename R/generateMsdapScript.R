#' Generates msdap R script
#'
#' This function generates an R script to run msdap.
#' It is a helper function for the processDiannMSdap.R function.
#' 
#' Settings for MSdap are the defaults for DIA data, if necessary change inside the first
#' glue::glue block.
#'
#' NOTE: sample_metadata.xlsx can be generated manually or with the function
#' write_template_for_sample_metadata() from library(msdap).
#'
#' @param inputPath Path where input files report.tsv (diann), sample_metadata.xlsx, FASTA used in the search and contrast.csv are located.
#' @param outputPath Path where msdap results will be saved.
#' @return Creates an R script to run msdap.
#' @export

generateMsdapScript <- function(inputPath, outputPath) {
  
  # Resolve file paths
  metadata   <- here::here(inputPath, "sample_metadata.xlsx")
  report     <- here::here(inputPath, "report.tsv")
  fasta      <- Sys.glob(here::here(inputPath, "*.fasta"))
  if (length(fasta) == 0) fasta <- Sys.glob(here::here(inputPath, "*.fas"))
  fasta_str  <- paste(fasta, collapse = '", "')
  msdap_path <- file.path(outputPath, "msdap_results")
  
  # Build contrast_pairs string for the generated script
  contrast_pairs     <- data.table::fread(here::here(inputPath, "contrasts.csv"))
  contrast_pairs     <- Map(c, contrast_pairs$contrast1, contrast_pairs$contrast2)
  contrast_pairs_str <- paste0(
    "list(",
    paste(
      sapply(contrast_pairs, function(p) paste0('c("', p[1], '", "', p[2], '")')),
      collapse = ", "
    ),
    ")"
  )
  
  glue::glue('
library(msdap)

dataset = import_dataset_diann(filename = "{report}")
dataset = import_fasta(dataset, files = c("{fasta_str}"))
dataset = import_sample_metadata(dataset, filename = "{metadata}")
dataset = setup_contrasts(dataset, contrast_list = {contrast_pairs_str})

dataset = analysis_quickstart(
  dataset,
  filter_min_detect = 3,
  filter_min_quant = 3,
  filter_fraction_detect = 0.75,
  filter_fraction_quant = 0.75,
  filter_min_peptide_per_prot = 1,
  filter_by_contrast = TRUE,
  norm_algorithm = c("vsn", "modebetween_protein"),
  dea_algorithm = c("deqms", "msempire", "msqrob"),
  dea_qvalue_threshold = 0.01,
  dea_log2foldchange_threshold = NA,
  output_qc_report = TRUE,
  output_abundance_tables = TRUE,
  output_dir = "{msdap_path}",
  output_within_timestamped_subdirectory = FALSE,
  dump_all_data = TRUE
)
  ')
}