#' Build the full DIA-NN -> MSDAP -> Shiny pipeline
#'
#' Generates every script needed to take a DIA-NN search + raw files through
#' MSDAP and into a ready-to-serve Shiny app, and chains them together with
#' slurm dependencies so a single `sbatch run_pipeline.sh` runs the whole
#' thing in the correct order:
#'   1. msconvert   - raw files -> .mzML, staged next to report.tsv
#'   2. msdap       - import_dataset_diann() + analysis_quickstart()
#'   3. postprocess - MSDAP output -> app.R's CSVs/PNGs/xlsx/zip
#'   4. assemble    - lays those files out at finalPath/{data,www}/ and
#'                    writes a job-specific app.R
#'
#' @param inputPath Path where report.tsv (diann), sample_metadata.xlsx, contrasts.csv and the FASTA live.
#' @param outputPath Path where msdap_results/, plots/ and scripts/ are written.
#' @param rawFiles Path where the raw vendor files are located.
#' @param jobname Job identifier, e.g. "PROT-1256". Used for output filenames and the final deploy path.
#' @param species Species code used by app.R's gene-group logic. Defaults to "Hs".
#' @param templateAppPath Path to the master app.R template that gets customized per job.
#' @param finalPath Where the finished Shiny app is deployed. Defaults to /home/gd2417/<jobname>.
#' @return Creates all scripts to import DIA-NN into msdap, post-process, and deploy the Shiny app.
#' @export

processDiannMSdap <- function(inputPath, outputPath, rawFiles, jobname, species = "Hs",
                               templateAppPath = system.file("templates", "app.r", package = "gsptools"),
                               finalPath = file.path("/home/gd2417/ShinyApps", jobname)) {

  # --- Validate inputs ---
  if (!file.exists(here::here(inputPath, "sample_metadata.xlsx")) |
      !file.exists(here::here(inputPath, "report.tsv")) |
      !file.exists(here::here(inputPath, "contrasts.csv")) |
      !any(file.exists(c(Sys.glob(here::here(inputPath, "*.fas")),
                          Sys.glob(here::here(inputPath, "*.fasta")))))) {
    stop("Missing input files: sample_metadata.xlsx, report.tsv, contrasts.csv, and/or .fasta/.fas")
  }
  if (!nzchar(templateAppPath) || !file.exists(templateAppPath)) {
    stop("Shiny template not found: '", templateAppPath, "'. ",
         "Is gsptools installed with inst/templates/app_template.R?", call. = FALSE)
  }

  # --- Create directory structure ---
  scripts_dir <- file.path(outputPath, "scripts")
  dir.create(scripts_dir,                            recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outputPath, "msdap_results"),  recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outputPath, "plots"),          recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dirname(inputPath), "logs"),   recursive = TRUE, showWarnings = FALSE)

  # === Step 1: mzML conversion (MUST run before msdap) ===
  mzML_slurm_path <- file.path(scripts_dir, "msconvert.sh")
  writeLines(generateMsconvertScript(rawFiles = rawFiles,
                                      inputPath = inputPath,
                                      jobName = "msconvert",
                                      cpus = 15, mem = "100G"),
             mzML_slurm_path)

  # === Step 2: MSdap ===
  msdap_r_path     <- file.path(scripts_dir, "msdap_analysis.R")
  msdap_slurm_path <- file.path(scripts_dir, "msdap_analysis.sh")
  writeLines(generateMsdapScript(inputPath, outputPath), msdap_r_path)
  writeLines(generateSlurmScript(msdap_r_path, inputPath,
                                  jobName = "msdap",
                                  cpus = 15, mem = "100G"), msdap_slurm_path)

  # === Step 3: Post-processing (builds Shiny-ready CSVs/PNGs/xlsx/zip) ===
  post_r_path     <- file.path(scripts_dir, "postprocessing.R")
  post_slurm_path <- file.path(scripts_dir, "postprocessing.sh")
  writeLines(generatePostprocessingScript(inputPath, outputPath, jobname), post_r_path)
  writeLines(generateSlurmScript(post_r_path, inputPath,
                                  jobName = "postprocess",
                                  cpus = 4, mem = "32G"), post_slurm_path)

  # === Step 4: Assemble the deployable Shiny app at finalPath ===
  assemble_r_path     <- file.path(scripts_dir, "assemble_shinyapp.R")
  assemble_slurm_path <- file.path(scripts_dir, "assemble_shinyapp.sh")
  writeLines(glue::glue('
library(gsptools)

gsptools::generateShinyAppScript(
  plotsPath       = "{file.path(outputPath, "plots")}",
  templateAppPath = "{templateAppPath}",
  jobname         = "{jobname}",
  species         = "{species}",
  finalPath       = "{finalPath}"
)
  '), assemble_r_path)
  writeLines(generateSlurmScript(assemble_r_path, outputPath,
                                  jobName = "assemble",
                                  cpus = 2, mem = "8G"), assemble_slurm_path)

  # === Master submission script, steps chained in execution order ===
  pipeline_sh_path <- file.path(scripts_dir, "run_pipeline.sh")
  writeLines(
    generatePipelineSh(list(
      msconvert   = mzML_slurm_path,
      msdap       = msdap_slurm_path,
      postprocess = post_slurm_path,
      assemble    = assemble_slurm_path
      ),
    inputPath = inputPath),
    pipeline_sh_path
  )
  system2("chmod", c("+x", pipeline_sh_path))

  message("=== All scripts generated ===")
  message("Output directory structure:")
  message("  ", scripts_dir, "/")
  message("    msconvert.sh")
  message("    msdap_analysis.R / msdap_analysis.sh")
  message("    postprocessing.R / postprocessing.sh")
  message("    assemble_shinyapp.R / assemble_shinyapp.sh")
  message("")
  message("Deployed app (after the pipeline finishes) will be at: ", finalPath)
  message("  ", finalPath, "/data  - CSVs, results.xlsx, results.zip")
  message("  ", finalPath, "/www   - chromatogram PNGs")
  message("  ", finalPath, "/app.R - customized for jobname = ", jobname)
  message("")
  message("To submit, run in terminal:")
  message("  sbatch ", pipeline_sh_path)
}
