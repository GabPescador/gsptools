#' Creates folder structure for PROT-XXXX projects
#'
#' This function creates the folder structure inside the path you input that is used for all other
#' mass spec related functions from gsptools package.
#'
#' @param path Path to your PROT-XXXX folder
#' @return Creates directories inside PROT-XXXX folder
#' @export
folderStructure <- function (path) {

  outputPath <- here(path, "output/")

  # Create output directories in case they don't exist
  paths <- c(
    file.path(outputPath),
    file.path(outputPath, "plots"),
    file.path(outputPath, "tables"),
    file.path(outputPath, "normalyzerDE"),
    file.path(path, "0_scripts"),
    file.path(path, "1_raw"),
    file.path(path, "2_fasta"),
    file.path(path, "3_fragpipe"),
    file.path(path, "4_input")
  )
  if (matrixStats::count(dir.exists(paths)) < 9) {
    invisible(sapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))
  }
}
