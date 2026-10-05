#' Generates slurm script to submit all jobs to HPC
#'
#' Chains a list of slurm scripts together using --dependency=afterok, so
#' each step only starts once the previous one succeeds. Step order is taken
#' directly from the order of `scriptPaths`, so callers must pass steps in the
#' order they should run (msconvert -> msdap -> postprocessing -> shinyapp).
#' It is a helper function for the processDiannMSdap.R function.
#'
#' @param scriptPaths Named list of paths to the slurm scripts to submit, in execution order.
#' @param inputPath Path used to resolve the shared log directory (dirname(inputPath)/logs).
#' @return Creates slurm script to run whole pipeline on HPC.
#' @keywords internal
#' @noRd

generatePipelineSh <- function(scriptPaths, inputPath) {

  # --- Build one sbatch submission line per step, chained via dependencies ---
  steps <- lapply(seq_along(scriptPaths), function(i) {
    name     <- toupper(names(scriptPaths)[i])
    var_name <- paste0("JOB", i)
    dep      <- if (i == 1) "" else paste0("--dependency=afterok:$JOB", i - 1, " ")

    paste0(
      "# Step ", i, ": ", name, "\n",
      var_name, "=$(sbatch ", dep, scriptPaths[[i]], " | awk '{print $4}')\n",
      'echo "', name, ' job submitted: $', var_name, '"'
    )
  })

  log_dir  <- file.path(dirname(inputPath), "logs")
  job_vars <- paste(paste0("$JOB", seq_along(scriptPaths)), collapse = ",")

  glue::glue(r'(
#!/bin/bash
#SBATCH --job-name=pipeline_sh
#SBATCH --cpus-per-task=15
#SBATCH --mem=100
#SBATCH --time=01-00:00:00
#SBATCH --output="<<log_dir>>/%j_%A_%a.out"
#SBATCH --error="<<log_dir>>/%j_%A_%a.err"
#SBATCH --mail-type=ALL

echo "Submitting pipeline..."
echo "Start time: $(date)"
echo "========================================="

<<paste(steps, collapse = "\n\n")>>

echo "========================================="
echo "Pipeline submitted successfully"
echo "Track with: squeue -j <<job_vars>>"
echo "========================================="
  )', .open = "<<", .close = ">>")
}
