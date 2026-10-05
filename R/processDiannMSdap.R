#' Generates diann processing pipeline
#'
#' This function generates the full diann processing pipeline and puts all files on its places.
#'
#' @param inputPath Path where input files report.tsv (diann), sample_metadata.xlsx, FASTA used in the search and contrast.csv are located.
#' @param outputPath Path where results will be saved.
#' @param rawFiles Path where your raw files are located.
#' @return Creates all scripts to import DIA from diann to msdap and post-process.
#' @export

processDiannMSdap <- function(inputPath, outputPath, rawFiles) {

  # Validate inputs
  if (!file.exists(here::here(inputPath, "sample_metadata.xlsx")) |
      !file.exists(here::here(inputPath, "report.tsv")) |
      !file.exists(here::here(inputPath, "contrasts.csv")) |
      !any(file.exists(c(Sys.glob(here::here(inputPath, "*.fas")),
                         Sys.glob(here::here(inputPath, "*.fasta")))))) {
    stop("Missing input files: sample_metadata.xlsx, report.tsv, contrasts.csv, and/or .fasta/.fas")
  }

  # Create directory structure
  scripts_dir <- file.path(outputPath, "scripts")
  dir.create(scripts_dir,                              recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outputPath, "msdap_results"),   recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outputPath, "plots"),           recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dirname(inputPath), "logs"),    recursive = TRUE, showWarnings = FALSE)

  # Step 1: MSdap
  msdap_r_path    <- file.path(scripts_dir, "msdap_analysis.R")
  msdap_slurm_path <- file.path(scripts_dir, "msdap_analysis.sh")
  writeLines(generateMsdapScript(inputPath, outputPath), msdap_r_path)
  writeLines(generateSlurmScript(msdap_r_path, inputPath,
                                   jobName = "msdap",
                                   cpus = 15, mem = "100G"), msdap_slurm_path)

  # Step 2: mzML conversion
  mzML_slurm_path <- file.path(scripts_dir, "msconvert.sh")
  writeLines(generateMsconvertScript(rawFiles = rawFiles,
                                     inputPath = inputPath,
                                     mzML_slurm_path))

  # Step 3: Post-processing
  post_r_path     <- file.path(scripts_dir, "postprocessing.R")
  post_slurm_path <- file.path(scripts_dir, "postprocessing.sh")
  writeLines(generatePostprocessingScript(outputPath), post_r_path)
  writeLines(generateSlurmScript(post_r_path, outputPath,
                                   jobName = "postprocess",
                                   cpus = 4, mem = "32G"), post_slurm_path)

  # Master submission script
  pipeline_sh_path <- file.path(scripts_dir, "run_pipeline.sh")
  writeLines(
    generatePipelineSh(list(
      msdap       = msdap_slurm_path,
      mzML = mzML_slurm_path,
      postprocess = post_slurm_path
      ),
    inputPath = inputPath),
    pipeline_sh_path
  )
  system2("chmod", c("+x", pipeline_sh_path))

  message("=== All scripts generated ===")
  message("Output directory structure:")
  message("  ", scripts_dir, "/")
  message("    msdap_analysis.R")
  message("    msdap_analysis.sh")
  message("    msconvert.sh")
  message("    postprocessing.R")
  message("    postprocessing.sh")
  message("")
  message("To submit, run in terminal:")
  message("  sbatch ", pipeline_sh_path)
}
