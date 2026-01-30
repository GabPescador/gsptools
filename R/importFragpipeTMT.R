#' Imports TMT data from fragpipe and normalizes with NormalyzerDE
#'
#' This function takes an inputPath with "metadata.csv", "contrast.csv" and
#' "abundance_protein_MD.tsv", performs normalization with NormalyzerDE and saves
#' all outputs in the ouputPath. Examples of input files can be found in inst/extdata/importFragpipeTMT_examples.
#'
#' @param inputPath Input path where "metadata.csv", "contrast.csv" and "abundance_protein_MD.tsv" are located.
#' @param jobname Specifies the folder name for NormalyzerDE, for convenience I always use the PROT-XXXX.
#' @param outputpath Output path to save all the outputs.
#' @param force Default to FALSE. Forces normalyzer DE to run.
#' @return Generates the NormalyzerDE report from fragpipe input.
#' @export

importFragpipeTMT <- function(inputPath, jobname, outputPath, force = FALSE){

  # Check that tables are in the input folder
  if (file.exists(paste0(inputPath, "metadata.csv")) &
      file.exists(paste0(inputPath, "contrasts.csv")) &
      file.exists(paste0(inputPath, "abundance_protein_MD.tsv"))) {

    print("Input tables exists, proceeding with pipeline...")

  } else {

    stop("No metadata.csv, contrasts.csv and/or abundance_protein_MD.tsv found, stopping execution.")

  }

  # Create output directories in case they don't exist
      paths <- c(
        file.path(outputPath, "plots"),
        file.path(outputPath, "tables"),
        file.path(outputPath, "normalyzerDE")
            )
      invisible(sapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))

  # Read input file and remove contaminants and rev proteins
      tmt <- fread(here(inputPath, "abundance_protein_MD.tsv")) %>%
             rename("ProteinID" = "Index",
                    "ProteinName" = "Gene",
                    # "EntryID" = "Protein ID"
                    ) %>%
        select(-`Protein ID`) %>%
        mutate(ProteinName = tolower(ProteinName)) %>%
        filter(!str_detect(ProteinID, "Cont")) %>%
        filter(!str_detect(ProteinID, "rev")) %>%
        setNames(snakecase::to_snake_case(names(.)))

  # Normalization
  ## Creating the summarizedExperiment object
  ### Matrix table
      matrix <- tmt %>%
        select(protein_id, reference_intensity:last_col(), -reference_intensity) %>%
        as.data.frame() %>%
        tibble::column_to_rownames("protein_id") %>%
        as.matrix()

  ### Design table
      design <- fread(here(inputPath, "metadata.csv")) %>%
        mutate(across(where(is.character), snakecase::to_snake_case))

  ### colData
      colData <- design

  ### rowData
      rowData <- tmt %>%
        select("protein_id","protein_name", "reference_intensity")

  ### SumarizedObject
      se <- SummarizedExperiment(assay=list(raw=matrix),
                                 colData = colData,
                                 rowData = rowData)

  ### Normalization with NormalizerDE
      if (force == TRUE) {

        print("Running NormalizerDE...")

        normalyzer(
          jobName = jobname,
          #designPath = NULL,
          #dataPath = NULL,
          experimentObj = se,
          outputDir = here(outputPath, "normalyzerDE"),
          forceAllMethods = FALSE,
          omitLowAbundSamples = FALSE,
          #sampleAbundThres = 0,
          #tinyRunThres = 50,
          requireReplicates = TRUE,
          #normalizeRetentionTime = TRUE,
          #plotRows = 3,
          #plotCols = 4,
          zeroToNA = TRUE,
          sampleColName = "sample",
          groupColName = "group",
          #inputFormat = "default",
          skipAnalysis = FALSE,
          quiet = FALSE,
          noLogTransform = TRUE,
          #writeReportAsPngs = FALSE,
          #rtStepSizeMinutes = 1,
          #rtWindowMinCount = 100,
          #rtWindowShifts = 1,
          #rtWindowMergeMethod = "mean"
        )

      } else if (dir.exists(here(outputPath, "normalyzerDE", jobname))) { # simpler way to check directory already exists

        print(paste0("Normalization already performed and can be found in ",
                     here(outputPath, "normalyzerDE", jobname)))

      } else {
        print("Running NormalizerDE...")

        normalyzer(
          jobName = jobname,
          #designPath = NULL,
          #dataPath = NULL,
          experimentObj = se,
          outputDir = here(outputPath, "normalyzerDE"),
          forceAllMethods = FALSE,
          omitLowAbundSamples = FALSE,
          #sampleAbundThres = 0,
          #tinyRunThres = 50,
          requireReplicates = TRUE,
          #normalizeRetentionTime = TRUE,
          #plotRows = 3,
          #plotCols = 4,
          zeroToNA = TRUE,
          sampleColName = "sample",
          groupColName = "group",
          #inputFormat = "default",
          skipAnalysis = FALSE,
          quiet = FALSE,
          noLogTransform = TRUE,
          #writeReportAsPngs = FALSE,
          #rtStepSizeMinutes = 1,
          #rtWindowMinCount = 100,
          #rtWindowShifts = 1,
          #rtWindowMergeMethod = "mean"
        )
      }

}
