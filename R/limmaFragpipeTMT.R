#' Imports normalized TMT data from NormalyzerDE and performs limma contrasts.
#'
#' This function takes an inputPath with "metadata.csv", "contrast.csv" and
#' "abundance_protein_MD.tsv", and the chosen normalization method from the
#' NormalizerDE output folder, then perform limma differential gene expression
#' and saves all outputs in the ouputPath. jobname should be the same used from
#' importFragpipeTMT.
#'
#' @param inputPath Input path where "metadata.csv", "contrast.csv" and "abundance_protein_MD.tsv" are located.
#' @param jobname Specifies the folder name for NormalyzerDE, for convenience I always use the PROT-XXXX.
#' @param method Specifies the chosen normalization method, defaults to "quantile". Options are: "cycloess", "gi", "log2", "mean", "median", "quantile", "rlr"
#' @param outputpath Output path to save all the outputs.
#' @return Generates the limma contrasts and saves them in the outputPath.
#' @export

limmaFragpipeTMT <- function(inputPath, jobname, method = "quantile", outputPath){

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
    file.path(outputPath, "output", "plots"),
    file.path(outputPath, "output", "tables"),
    file.path(outputPath, "output", "normalyzerDE")
  )
  invisible(sapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))

  if (method == "cycloess") {

    norm <- data.table::fread(here::here(outputPath, "output", "normalyzerDE", jobname, "CycLoess-normalized.txt"))

  } else if (method == "gi") {

    norm <- data.table::fread(here::here(outputPath, "output", "normalyzerDE", jobname, "GI-normalized.txt"))

  } else if (method == "log2") {

    norm <- data.table::fread(here::here(outputPath, "output", "normalyzerDE", jobname, "log2-normalized.txt"))

  } else if (method == "mean") {

    norm <- data.table::fread(here::here(outputPath, "output", "normalyzerDE", jobname, "mean-normalized.txt"))

  } else if (method == "median") {

    norm <- data.table::fread(here::here(outputPath, "output", "normalyzerDE", jobname, "median-normalized.txt"))

  } else if (method == "quantile") {

    norm <- data.table::fread(here::here(outputPath, "output", "normalyzerDE", jobname, "Quantile-normalized.txt"))

  } else if (method == "rlr") {

    norm <- data.table::fread(here::here(outputPath, "output", "normalyzerDE", jobname, "RLR-normalized.txt"))

  }

  # QC for before and after normalization
  # Before normalization
  tmt <- data.table::fread(here::here(inputPath, "abundance_protein_MD.tsv")) %>%
    dplyr::rename("ProteinID" = "Index",
           "ProteinName" = "Gene") %>%
    mutate(ProteinName = tolower(ProteinName)) %>%
    filter(!stringr::str_detect(.data$ProteinID, "Cont")) %>%
    filter(!stringr::str_detect(.data$ProteinID, "rev"))

  p1 <- tmt %>%
    reshape2::melt(id.vars = colnames(tmt)[1:11]) %>%
    ggplot(aes(x=.data$variable, y=.data$value)) +
    geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("") +
    ylab("Log2(Abundance)") +
    ggtitle("Before Normalization")

  # After normalization
  p2 <- norm %>%
    reshape2::melt(id.vars = colnames(norm)[1:7]) %>%
    ggplot(aes(x=.data$variable, y=.data$value)) +
    geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("") +
    ylab("Log2(Abundance)") +
    ggtitle("After Normalization")

  # Plots together
  grDevices::pdf(here::here(outputPath, "output", "plots", paste0(jobname, "_normalization_boxplot.pdf")))
  print(cowplot::plot_grid(p1, p2, ncol = 2))
  grDevices::dev.off()

  # Matrix for limma
  ### Matrix table
  matrix_norm <- norm %>%
    select(ProteinID, `Indistinguishable Proteins`:last_col(), -`Indistinguishable Proteins`) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("ProteinID") %>%
    as.matrix()

  ### Design table
  design <- fread(here::here(inputPath, "metadata.csv"))
  design_limma <- stats::model.matrix(~0 + factor(design$group))
  colnames(design_limma) <- levels(factor(design$group))
  rownames(design_limma) <- design$sample

  ### Contrast table
  contrast_table <- fread(here::here(inputPath, "contrasts.csv"))

  contrast_formulas <- vector()
  for(i in 1:nrow(contrast_table)){
    temp <- c(paste0(contrast_table[i,1], "-", contrast_table[i,2]))
    contrast_formulas <- c(contrast_formulas, temp)
  }

  contrast.matrix <- limma::makeContrasts(
    contrasts = contrast_formulas,
    levels = design_limma
  )

  # Fit the data and create contrasts
  fit <- limma::lmFit(matrix_norm, design_limma)
  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit2)

  grDevices::pdf(here::here(outputPath, "output", "plots", paste0(jobname, "limma_mds-plot.pdf")))
  limma::plotMDS(fit)
  grDevices::dev.off()

  # Saving contrasts
  df <- data.frame()
  for(i in colnames(contrast.matrix)){
    temp <- data.frame()
    temp <- limma::topTable(fit2, coef=i, number=Inf)
    temp$contrast <- i
    temp$ProteinID <- rownames(temp)
    temp <- dplyr::relocate(temp, "ProteinID", .before = "logFC")

    df <- rbind(df, temp)
  }

  # Putting protein names back
  df <- merge(df, norm[,c(1,2)], by = "ProteinID") %>%
    relocate(ProteinName, .after = ProteinID)

  if (file.exists(paste0(outputPath, jobname, "_limmaContrasts.csv"))){

    print(paste0("Contrasts already performed and can be found in ",
                 here::here(outputPath, "output", "tables", paste0(jobname, "_limmaContrasts.csv"))))

  } else {

    write_csv(df, here::here(outputPath, "output", "tables", paste0(jobname, "_limmaContrasts.csv")))

  }

}
