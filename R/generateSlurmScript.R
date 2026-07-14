#' Generates slurm script for submitting R scripts to HPC
#'
#' This function generates a slurm script to run R scripts on HPC cluster.
#' It is a helper function for the processDiannMSdap.R function.
#'
#' @param scriptPath Path to your R script.
#' @param inputPath Path where input files are stored for processDiannMSdap.R function. Uses parent directory to this to save logs on a log folder.
#' @param jobName Name for the job submission. Defaults to "job".
#' @param cpus Integer for how many cpus should be used for the submission. Defaults to 15.
#' @param mem How much memory should be allocated for the submission. This should be an integer followed by G, like 50G. Defaults to 100G.
#' @param time How much time should be requested for the submission. This should be integers with the following pattern: Days-Hours:Minutes:Seconds. Defaults to 1 day (01-00:00:00).
#' @return Creates slurm scripts to run MSdap and post processing on HPC.
#' @keywords internal
#' @noRd

generateSlurmScript <- function(scriptPath, inputPath,
                                  jobName = "job",
                                  cpus = 15, mem = "100G",
                                  time = "01-00:00:00") {

  log_dir <- file.path(dirname(inputPath), "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  glue::glue(r'(
#!/bin/bash
#SBATCH --job-name=<<jobName>>
#SBATCH --cpus-per-task=<<cpus>>
#SBATCH --mem=<<mem>>
#SBATCH --time=<<time>>
#SBATCH --output="<<log_dir>>/%j_%A_%a.out.out"
#SBATCH --error="<<log_dir>>/%j_%A_%a.out.err"
#SBATCH --mail-type=ALL

cd <<dirname(inputPath)>>
mkdir -p <<log_dir>>

trap "cp <<log_dir>>/$SLURM_JOB_ID* <<log_dir>>/" EXIT

ml R/4.4.1

Rscript <<scriptPath>>

echo "========================================="
echo "Done!"
echo "End time: $(date)"
echo "========================================="
echo "Resource usage:"
sacct -j $SLURM_JOB_ID \
--format=JobID,Elapsed,CPUTime,MaxRSS,State \
--units=G

command -v job_stats && job_stats $SLURM_JOB_ID | sed "s/\x1B\[[0-9;]*m//g"
echo "========================================="
  )', .open = "<<", .close = ">>")
}
