#' Import diaNN output into MSdap for processing
#'
#' This function takes report.tsv from diaNN, sample_metadata.xlsx from MSdap, a contrasts.csv file
#' and the FASTA file used in the proteomics search and creates an .R script to run MSdap and
#' a SLURM script to run in HPC cluster.
#'
#' Settings for MSdap are the defaults for DIA data, if necessary change inside the first
#' glue::glue block.
#'
#' NOTE: sample_metadata.xlsx can be generated manually or with the function
#' write_template_for_sample_metadata() from MSdap.
#'
#' @param inputPath Path where the 4 files necessary are located
#' @param outputPath Path to where you want to save the MSdap output
#' @return Creates 2 scripts to run MSdap on HPC.
#' @export

importDiannMSdap <- function(inputPath, outputPath) {

  # Check that tables are in the input folder ----
  if (file.exists(here(inputPath, "sample_metadata.xlsx")) &
      file.exists(here(inputPath, "report.tsv")) &
      file.exists(here(inputPath, "contrasts.csv")) &
      any(file.exists(c(Sys.glob(here(inputPath, "*.fas")),
                        Sys.glob(here(inputPath, "*.fasta")))))
  ) {

    print("Input tables exists, proceeding with MSdap...")

  } else {

    stop("No sample_metadata.xlsx, contrasts.csv, report.tsv and/or .fasta/.fas found, stopping execution.")

  }

  # Read paths to files ----

  metadata <- here::here(inputPath, "sample_metadata.xlsx")
  report <- here::here(inputPath, "report.tsv")
  contrast_pairs <- data.table::fread(here(inputPath, "contrasts.csv"))
  contrast_pairs <- Map(c, contrast_pairs$contrast1, contrast_pairs$contrast2)
  contrast_pairs_str <- paste0(
    "list(",
    paste(
      sapply(contrast_pairs, function(p) paste0('c("', p[1], '", "', p[2], '")')),
      collapse = ", "
    ),
    ")"
  )
  fasta <- Sys.glob(here::here(inputPath, "*.fasta"))
  if (length(fasta) == 0){
    fasta <- Sys.glob(here::here(inputPath, "*.fas"))
  }
  msdap_path <- paste0(outputPath, '/msdap_results')
  if(!dir.exists(msdap_path)){
    dir.create(msdap_path, recursive = TRUE)
  }

  # Create the msdap script as .R ----

  script <- glue::glue('
    # Auto-generated script
    library(here)
    library(tidyverse)
    library(data.table)
    library(readxl)
    library(rstatix)
    library(SummarizedExperiment)
    library(msdap)
    library(plotrix)
    library(stringr)
    library(dplyr)
    library(gsptools)

    # MSdap code

    dataset = import_dataset_diann(filename = "{report}")
    dataset = import_fasta(dataset, files = "{fasta}")
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
      diffdetect_min_peptides_observed = 2,
      diffdetect_min_samples_observed = 3,
      diffdetect_min_fraction_observed = 0.5,
      output_qc_report = TRUE,
      output_abundance_tables = TRUE,
      output_dir = "{msdap_path}",
      output_within_timestamped_subdirectory = TRUE
    )

  ')

  script_path <- file.path(here::here(dirname(inputPath), "0_scripts/"), paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_msdap_script.R"))
  writeLines(script, script_path)

  # Create the slurm script as .sh ----

  script_slurm <- glue::glue(r'(
    #!/bin/bash
    #SBATCH --cpus-per-task=15
    #SBATCH --mem=100G
    #SBATCH --time=01-00:00:00
    #SBATCH --output="/home/gd2417/logs/%A_%a.out"
    #SBATCH --error="/home/gd2417/logs/%A_%a.err"
    #SBATCH --mail-type=ALL

    cd <<dirname(inputPath)>>
    mkdir -p <<dirname(inputPath)>>/logs/

    # This will copy logs to where our files are being saved
    trap "cp /home/gd2417/logs/$SLURM_JOB_ID* <<dirname(inputPath)>>/logs/" EXIT

    ml R/4.4.1

    Rscript <<script_path>>

    echo "========================================="
    echo "Done!"
    echo "End time: $(date)"
    echo "To see resource usage, run job_stats $SLURM_JOB_ID"
    echo "========================================="
    echo " "
    echo "========================================="
    echo "Resource usage:"
    sacct -j $SLURM_JOB_ID \
        --format=JobID,Elapsed,CPUTime,MaxRSS,State \
        --units=G

    # also run job_stats if available
    command -v job_stats && job_stats $SLURM_JOB_ID | sed "s/\x1B\[[0-9;]*m//g"
    echo "========================================="
  )', .open = "<<", .close = ">>")

  script_path_slurm <- file.path(here::here(dirname(inputPath), "0_scripts/"), paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_msdap_slurm.sh"))
  writeLines(script_slurm, script_path_slurm)
  message("Scripts written to: ", here::here(dirname(inputPath), "0_scripts/"))
}
