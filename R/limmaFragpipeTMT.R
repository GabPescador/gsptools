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
#' @param replicateFilter Defaults to TRUE. Filters proteins that are not present in at least 2 replicates.
#' @param groups Defaults to FALSE. If TRUE, will use the type column in metadata.csv to search unique proteins based on those groupings.
#' @param force Defaults to FALSE. Forces saving the results table in case it already exists.
#' @return Generates the limma contrasts and saves them in the outputPath.
#' @export

limmaFragpipeTMT <- function(inputPath, jobname, method = "quantile", outputPath, replicateFilter = TRUE, groups = FALSE, force = FALSE){

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

    norm <- data.table::fread(here::here(outputPath, "output", "normalyzerDE", jobname, "CycLoess-normalized.txt")) %>%
      setNames(snakecase::to_snake_case(names(.)))

  } else if (method == "gi") {

    norm <- data.table::fread(here::here(outputPath, "output", "normalyzerDE", jobname, "GI-normalized.txt")) %>%
      setNames(snakecase::to_snake_case(names(.)))

  } else if (method == "log2") {

    norm <- data.table::fread(here::here(outputPath, "output", "normalyzerDE", jobname, "log2-normalized.txt")) %>%
      setNames(snakecase::to_snake_case(names(.)))

  } else if (method == "mean") {

    norm <- data.table::fread(here::here(outputPath, "output", "normalyzerDE", jobname, "mean-normalized.txt")) %>%
      setNames(snakecase::to_snake_case(names(.)))

  } else if (method == "median") {

    norm <- data.table::fread(here::here(outputPath, "output", "normalyzerDE", jobname, "median-normalized.txt")) %>%
      setNames(snakecase::to_snake_case(names(.)))

  } else if (method == "quantile") {

    norm <- data.table::fread(here::here(outputPath, "output", "normalyzerDE", jobname, "Quantile-normalized.txt")) %>%
      setNames(snakecase::to_snake_case(names(.)))

  } else if (method == "rlr") {

    norm <- data.table::fread(here::here(outputPath, "output", "normalyzerDE", jobname, "RLR-normalized.txt")) %>%
      setNames(snakecase::to_snake_case(names(.)))

  }

  metadata <- fread(here::here(inputPath, "metadata.csv")) %>%
    mutate(across(where(is.character), snakecase::to_snake_case))

  # QC for before and after normalization
  # Before normalization
  tmt <- data.table::fread(here::here(inputPath, "abundance_protein_MD.tsv")) %>%
    dplyr::rename("ProteinID" = "Index",
           "ProteinName" = "Gene",
           "EntryID" = "Protein ID") %>%
    mutate(ProteinName = tolower(ProteinName)) %>%
    filter(!stringr::str_detect(.data$ProteinID, "Cont")) %>%
    filter(!stringr::str_detect(.data$ProteinID, "rev")) %>%
    setNames(snakecase::to_snake_case(names(.)))

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

  # Taking out proteins that are not present in at least 2 replicates if replicateFilter == TRUE
  if (replicateFilter == TRUE){

    norm_long <- norm %>%
      reshape2::melt(id.vars = colnames(.)[1:7]) %>%
      # select(protein_id, variable, value) %>%
      merge(., metadata[,-3], by.x = "variable", by.y = "sample")

    norm_filtered <- norm_long %>%
      filter(!is.na(value)) %>%
      group_by(protein_id, group) %>%
      summarize(counts =n())

    norm_long <- merge(norm_long, norm_filtered, by = c("protein_id", "group"))

    norm_long <- norm_long %>%
      mutate(value = if_else(counts < 2, NA_real_, value))
      # group_by(protein_id, protein_name, group, protein, entry_id, entry_name, protein_description, indistinguishable_proteins) %>%
      # summarize(mean = mean(value, na.rm=TRUE))

    norm <- norm_long %>%
      pivot_wider(
        id_cols = c(protein_id, protein_name, protein, entry_id, entry_name, protein_description, indistinguishable_proteins),
        names_from = variable,
        values_from = c(value),
        names_sep = "_"
        )

  } else {

    norm <- norm

  }

  if (groups == FALSE) {
  # Creating a list of unique proteins based on groups
  uniquePerGroup <- norm %>%
    reshape2::melt(id.vars = colnames(.)[1:7]) %>%
    select(protein_id, variable, value) %>%
    merge(., metadata[,-3], by.x = "variable", by.y = "sample") %>%
    filter(!is.na(value)) %>%
    group_by(protein_id) %>%
    filter(n_distinct(group) == 1) %>%
    ungroup() %>%
    distinct(protein_id, group)

  uniquePerGroup <- merge(uniquePerGroup, norm[,c("protein_id", "protein_name")]) %>%
    relocate(protein_name, .before = group)

  } else {

    uniquePerGroup <- norm %>%
      reshape2::melt(id.vars = colnames(.)[1:7]) %>%
      select(protein_id, variable, value) %>%
      merge(., metadata[,-3], by.x = "variable", by.y = "sample") %>%
      filter(!is.na(value)) %>%
      group_by(protein_id) %>%
      filter(n_distinct(group) == 1) %>%
      ungroup() %>%
      distinct(protein_id, group)

    # Creating a list of unique proteins based on groups
    uniquePerType <- norm %>%
      reshape2::melt(id.vars = colnames(.)[1:7]) %>%
      select(protein_id, variable, value) %>%
      merge(., metadata[,-3], by.x = "variable", by.y = "sample") %>%
      filter(!is.na(value)) %>%
      group_by(protein_id) %>%
      filter(n_distinct(type) == 1) %>%
      ungroup() %>%
      distinct(protein_id, type) %>%
      rename("group" = "type")

    uniquePerGroup <- rbind(uniquePerGroup, uniquePerType)
    uniquePerGroup <- merge(uniquePerGroup, norm[,c("protein_id", "protein_name")]) %>%
      relocate(protein_name, .before = group)

  }

  # Matrix for limma
  ### Matrix table
  matrix_norm <- norm %>%
    select(protein_id, indistinguishable_proteins:last_col(), -indistinguishable_proteins) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("protein_id") %>%
    as.matrix()

  ### Design table
  design <- metadata
  design_limma <- stats::model.matrix(~0 + factor(design$group))
  colnames(design_limma) <- levels(factor(design$group))
  rownames(design_limma) <- design$sample

  ### Contrast table
  contrast_table <- fread(here::here(inputPath, "contrasts.csv")) %>%
    mutate(across(where(is.character), snakecase::to_snake_case))

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
    temp$protein_id <- rownames(temp)
    temp <- dplyr::relocate(temp, "protein_id", .before = "logFC")

    df <- rbind(df, temp)
  }

  # Putting protein names back
  df <- merge(df, norm[,c(1,2)], by = "protein_id") %>%
    relocate(protein_name, .after = protein_id)

  df2 <- df %>%
    pivot_wider(
      id_cols = c(protein_id, protein_name),
      names_from = contrast,
      values_from = c(logFC, P.Value, adj.P.Val),
      names_sep = "_"
    )

  sheets <- list(
    "Metadata" = metadata,
    "NormalizedAbundances" = norm,
    "Contrasts" = as.data.frame(contrast_formulas),
    "limmaDEA" = df2,
    "UniqueProteins" = uniquePerGroup
  )

  if (force == TRUE) {

    writexl::write_xlsx(x = sheets,
                        path = here::here(outputPath, "output", "tables", paste0(jobname, "_results.xlsx")),
                        format_headers = FALSE)

    for(i in names(sheets)) {
    readr::write_csv(x = sheets[[i]],
                     path = here::here(outputPath, "output", "tables", paste0(jobname, "_", i, ".csv")))
    }

  } else if (file.exists(here::here(outputPath, "output", "tables", paste0(jobname, "_results.xlsx")))) {

    print(paste0("Contrasts already performed and can be found in ",
                 here::here(outputPath, "output", "tables", paste0(jobname, "_results.xlsx"))))

  } else {

    writexl::write_xlsx(x = sheets,
                        path = here::here(outputPath, "output", "tables", paste0(jobname, "_results.xlsx")),
                        format_headers = FALSE)

    for(i in names(sheets)) {
      readr::write_csv(x = sheets[[i]],
                       path = here::here(outputPath, "output", "tables", paste0(jobname, "_", i, ".csv")))

    }
  }
}
