#' Assembles the final Shiny app folder structure
#'
#' Takes the tables/images produced by the post-processing step and lays them
#' out the way app.R expects: CSVs + xlsx/zip in `<finalPath>/data`, chromatogram
#' PNGs in `<finalPath>/www`, and a copy of app.R with `jobname`/`species`
#' swapped in for this job. This is the last step in the pipeline, run after
#' generatePostprocessingScript.R.
#'
#' @param outputPath Path where postprocessing.R wrote its CSVs/PNGs/xlsx/zip (outputPath).
#' @param templateAppPath Path to the master app.R (the reusable Shiny app template).
#' @param jobname Job identifier, e.g. "PROT-1256".
#' @param species Species code used by app.R's gene-group logic, e.g. "Hs" or "Dr".
#' @param finalPath Destination root for the deployed app. Defaults to /home/gd2417/<jobname>.
#' @return Invisibly returns `finalPath`. Called for its side effect of writing files to disk.
#' @keywords internal
#' @noRd

generateShinyAppScript <- function(outputPath,
                                   templateAppPath = system.file("templates", "app.r", package = "gsptools"),
                                   jobname,
                                   species = "Hs",
                                   finalPath = file.path("/home/gd2417/ShinyApps", jobname)) {

  # --- 1) Create destination folder structure ---
  data_dir <- file.path(finalPath, "data")
  www_dir  <- file.path(finalPath, "www")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(www_dir,  recursive = TRUE, showWarnings = FALSE)

  # --- 2) Copy tables + combined results into data/ ---
  table_files <- list.files(
    paste0(outputPath, "/tables/"),
    pattern    = paste0("^", jobname, "_(BasicStats|NormalizedAbundances|DEA_long|UniqueProteins|results)\\.(csv|xlsx|zip)$"),
    full.names = TRUE
  )
  if (length(table_files) == 0) {
    warning("No matching data files found in ", outputPath, " - check generatePostprocessingScript.R ran successfully.")
  }
  file.copy(table_files, data_dir, overwrite = TRUE)

  # --- 3) Copy chromatogram PNGs into www/ ---
  chromatogram_files <- list.files(
    paste0(outputPath, "/plots/"),
    pattern = paste0("^", jobname, "_chromatograms_[[:alnum:]]+_page[0-9]+\\.png$"),
    full.names = TRUE
  )
  file.copy(chromatogram_files, www_dir, overwrite = TRUE)

  # --- 4) Write a job-specific app.R from the shared template ---
  app_lines <- readLines(templateAppPath)

  # Guard: if the template doesn't have a line starting with "jobname <- "
  # (e.g. someone reformatted app_template.R), sub() below would silently
  # match nothing and this app would deploy with whatever jobname it already
  # had - fail loudly instead.
  if (!any(grepl('^jobname <- ', app_lines))) {
    stop("templateAppPath has no line starting with 'jobname <- ' - cannot customize it. Check ", templateAppPath)
  }
  if (!any(grepl('^species <- ', app_lines))) {
    stop("templateAppPath has no line starting with 'species <- ' - cannot customize it. Check ", templateAppPath)
  }

  app_lines <- sub('^jobname <- .*',
                   sprintf('jobname <- "%s"', jobname), app_lines)
  app_lines <- sub('^species <- .*',
                   sprintf('species <- "%s" # set per job by generateShinyAppScript', species), app_lines)

  app_out_path <- file.path(finalPath, "app.R")
  writeLines(app_lines, app_out_path)

  # Verify what actually landed on disk matches what we intended to write -
  # catches silent mismatches from stale/regenerated scripts, wrong
  # templateAppPath, etc. rather than deploying a wrong jobname unnoticed.
  written <- readLines(app_out_path)
  jobname_line <- written[grepl('^jobname <- ', written)][1]
  if (!identical(jobname_line, sprintf('jobname <- "%s"', jobname))) {
    stop("app.R was written but its jobname line reads '", jobname_line,
         "' instead of the expected 'jobname <- \"", jobname, "\"'. ",
         "Re-run processDiannMSdap() for this jobname before deploying, ",
         "or check for a stale assemble_shinyapp.R in outputPath/scripts/.")
  }

chmod_recursive <- function(path, mode = "0777") {
  targets <- c(
    path,
    list.files(path, recursive = TRUE, full.names = TRUE,
               all.files = TRUE, include.dirs = TRUE, no.. = TRUE)
  )
  Sys.chmod(targets, mode = mode, use_umask = FALSE)
}

chmod_recursive(file.path("/home/gd2417/ShinyApps", jobname))

  message("=== Shiny app assembled at ", finalPath, " ===")
  message("  data/ : ", length(table_files), " file(s)")
  message("  www/  : ", length(chromatogram_files), " chromatogram(s)")
  message("  app.R written for jobname = ", jobname, ", species = ", species)

  invisible(finalPath)
}
