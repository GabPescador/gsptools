#' Generates slurm script for submitting msconvert to HPC
#'
#' This function generates a slurm script to run msconvert on HPC cluster.
#' It is a helper function for the processDiannMSdap.R function.
#'
#' @param jobname Name of your PROT-XXXX folder.
#' @param baseDir Path to where you want things created. Defaults to /home/gd2417/ShinyApps
#' @return Creates shiny app folder structure.
#' @keywords internal
#' @noRd
#' 
createShinyStructure <- function(jobname, baseDir = "/home/gd2417/ShinyApps") {
  job_path <- file.path(baseDir, jobname)
  
  subdirs <- c("data", "www")
  paths <- file.path(job_path, subdirs)
  
  for (p in paths) {
    if (!dir.exists(p)) {
      dir.create(p, recursive = TRUE)
      message("Created: ", p)
    } else {
      message("Already exists: ", p)
    }
  }
  
  invisible(job_path)
}