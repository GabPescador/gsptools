#' Imports normalized TMT diGly data from NormalyzerDE and performs limma contrasts.
#'
#' This function takes an inputPath with "metadata.csv", "contrast.csv",
#' "abundance_single-site_MD.tsv" and normalized protein abundances, and the chosen normalization
#' method from the NormalizerDE output folder, then normalizes each site with its paired input
#' sample, lastly performs limma differential gene expression
#' and saves all outputs in the ouputPath. jobname should be the same used from
#' importFragpipeTMTdiGly.
#'
#' @param inputPath Input path where "metadata.csv", "contrast.csv", "abundance_single-site_MD.tsv" and "abundance_protein_MD.tsv" are located.
#' @param jobname Specifies the folder name for NormalyzerDE, for convenience I always use the PROT-XXXX.
#' @param method Specifies the chosen normalization method, defaults to "quantile". Options are: "cycloess", "gi", "log2", "mean", "median", "quantile", "rlr"
#' @param outputpath Output path to save all the outputs.
#' @param replicateFilter Defaults to TRUE. Filters proteins that are not present in at least 2 replicates.
#' @param proteinNormalization Defaults to FALSE. Normalizes diGly values for each site based on the input values of respective proteins.
#' @param groups Defaults to FALSE. If TRUE, will use the type column in metadata.csv to search unique proteins based on those groupings.
#' @param force Defaults to FALSE. Forces saving the results table in case it already exists.
#' @param exclude Defaults to NULL. Defines a character vector of columns to be excluded from the pipeline.
#' @return Generates the limma contrasts and saves them in the outputPath.
#' @export

limmaFragpipeTMTdiGly <- function(inputPath,
                             jobname,
                             method = "quantile",
                             outputPath,
                             replicateFilter = TRUE,
                             proteinNormalization = FALSE,
                             groups = FALSE,
                             force = FALSE,
                             exclude = NULL){

  if (proteinNormalization == TRUE){
    # Check that tables are in the input folder
    if (file.exists(paste0(inputPath, "metadata.csv")) &
        file.exists(paste0(inputPath, "contrasts.csv")) &
        file.exists(paste0(inputPath, "abundance_single-site_MD.tsv")) &
        file.exists(paste0(inputPath, "abundance_protein_MD.tsv")) &
        file.exists(Sys.glob(paste0(inputPath, "input_*normalized.txt")))) {

      print("Input tables exists, proceeding with pipeline...")

    } else {

      stop("No metadata.csv, contrasts.csv, abundance_single-site_MD.tsv and/or normalized.txt found, stopping execution.")

    }
  } else {

  # Check that tables are in the input folder
  if (file.exists(paste0(inputPath, "metadata.csv")) &
      file.exists(paste0(inputPath, "contrasts.csv")) &
      file.exists(paste0(inputPath, "abundance_single-site_MD.tsv"))) {

    print("Input tables exists, proceeding with pipeline...")

  } else {

    stop("No metadata.csv, contrasts.csv, abundance_single-site_MD.tsv and/or normalized.txt found, stopping execution.")

    }
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

  if (!is.null(exclude)) {

    norm <- norm %>%
      select(-exclude)

    metadata <- fread(here::here(inputPath, "metadata.csv")) %>%
      mutate(across(where(is.character), snakecase::to_snake_case)) %>%
      filter(!sample %in% exclude)

  } else {

    norm <- norm

    metadata <- fread(here::here(inputPath, "metadata.csv")) %>%
      mutate(across(where(is.character), snakecase::to_snake_case))

  }

  if (proteinNormalization == TRUE){
  # Normalization to the input proteins
  # Import normalized protein input
  input <- data.table::fread(here::here(Sys.glob(paste0(inputPath, "input_*normalized.txt")))) %>%
    setNames(snakecase::to_snake_case(names(.)))

  long_norm <- norm %>%
    reshape2::melt(id.vars = colnames(norm)[1:8])

  long_input <- input %>%
    reshape2::melt(id.vars = colnames(input)[1:3]) %>%
    select(-protein_name, -reference_intensity) %>%
    filter(!is.na(value)) %>%
    rename("input" = "value")

  long_norm_2 <- merge(long_norm, long_input,
                       by = c("protein_id", "variable")) %>%
    mutate(norm_value = value-input)

  prot_norm <- long_norm_2 %>%
    select(-value, -input) %>%
    pivot_wider(
      id_cols = c(modified_site, protein_id, protein_name, peptide, sequence_window, start, end, max_pep_prob),
      names_from = variable,
      values_from = c(norm_value),
      names_sep = "_")
  }
  # QC for before and after normalization
  # Before normalization
  tmt <- data.table::fread(here::here(inputPath, "abundance_single-site_MD.tsv")) %>%
    dplyr::rename("ModifiedSite" = "Index",
           "ProteinName" = "Gene",
           "ProteinID" = "ProteinID") %>%
    mutate(ProteinName = tolower(ProteinName)) %>%
    filter(!stringr::str_detect(.data$ProteinID, "Cont")) %>%
    filter(!stringr::str_detect(.data$ProteinID, "rev")) %>%
    setNames(snakecase::to_snake_case(names(.)))

  p1 <- tmt %>%
    reshape2::melt(id.vars = colnames(tmt)[1:9]) %>%
    ggplot(aes(x=.data$variable, y=.data$value)) +
    geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
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

  if (proteinNormalization == TRUE){
  p3 <- long_norm_2 %>%
    ggplot(aes(x=.data$variable, y=.data$norm_value)) +
    geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("") +
    ylab("Log2(diGly) - Log2(Input)") +
    ggtitle("Protein-level Normalization")

  # Plots together
  grDevices::pdf(here::here(outputPath, "output", "plots", paste0(jobname, "_normalization_boxplot.pdf")))
  print(cowplot::plot_grid(p1, p2, p3, ncol = 3))
  grDevices::dev.off()

  } else {
    # Plots together
    grDevices::pdf(here::here(outputPath, "output", "plots", paste0(jobname, "_normalization_boxplot.pdf")))
    print(cowplot::plot_grid(p1, p2, ncol = 2))
    grDevices::dev.off()
  }

  if (proteinNormalization == TRUE){

    norm <- prot_norm

  } else {

    norm <- norm

  }

  # Taking out sites that are not present in at least 2 replicates if replicateFilter == TRUE
  if (replicateFilter == TRUE){

    norm_long <- norm %>%
      reshape2::melt(id.vars = colnames(.)[1:8]) %>%
      # select(protein_id, variable, value) %>%
      merge(., metadata[,-3], by.x = "variable", by.y = "sample")

    norm_filtered <- norm_long %>%
      filter(!is.na(value)) %>%
      group_by(modified_site, group) %>%
      summarize(counts =n())

    norm_long <- merge(norm_long, norm_filtered, by = c("modified_site", "group"))

    norm_long <- norm_long %>%
      mutate(value = if_else(counts < 2, NA_real_, value))

    norm <- norm_long %>%
      pivot_wider(
        id_cols = c(modified_site, protein_id, protein_name, peptide, sequence_window, start, end, max_pep_prob),
        names_from = variable,
        values_from = c(value),
        names_sep = "_"
      )

  } else {

    norm <- norm

  }

norm <- norm %>%
  mutate(protein_name = toupper(protein_name))

  if (groups == FALSE) {
  # Creating a list of unique proteins based on groups
  uniquePerGroup <- norm %>%
    reshape2::melt(id.vars = colnames(.)[1:8]) %>%
    select(modified_site, protein_id, variable, value) %>%
    merge(., metadata[,-3], by.x = "variable", by.y = "sample") %>%
    filter(!is.na(value)) %>%
    group_by(modified_site) %>%
    filter(n_distinct(group) == 1) %>%
    ungroup() %>%
    distinct(modified_site, group)

  uniquePerGroup <- merge(uniquePerGroup, norm[,c("modified_site", "protein_id", "protein_name")]) %>%
    relocate(protein_name, .before = group)

  } else {

    uniquePerGroup <- norm %>%
      reshape2::melt(id.vars = colnames(.)[1:8]) %>%
      select(modified_site, protein_id, variable, value) %>%
      merge(., metadata[,-3], by.x = "variable", by.y = "sample") %>%
      filter(!is.na(value)) %>%
      group_by(modified_site) %>%
      filter(n_distinct(group) == 1) %>%
      ungroup() %>%
      distinct(modified_site, group)

    # Creating a list of unique proteins based on groups
    uniquePerType <- norm %>%
      reshape2::melt(id.vars = colnames(.)[1:8]) %>%
      select(modified_site, protein_id, variable, value) %>%
      merge(., metadata[,-3], by.x = "variable", by.y = "sample") %>%
      filter(!is.na(value)) %>%
      group_by(modified_site) %>%
      filter(n_distinct(type) == 1) %>%
      ungroup() %>%
      distinct(modified_site, type) %>%
      rename("group" = "type")

    uniquePerGroup <- rbind(uniquePerGroup, uniquePerType)
    uniquePerGroup <- merge(uniquePerGroup, norm[,c("modified_site", "protein_id", "protein_name")]) %>%
      relocate(protein_name, .before = group)

  }

  # Matrix for limma
  ### Matrix table
  matrix_norm <- norm %>%
    select(modified_site, max_pep_prob:last_col(), -max_pep_prob) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("modified_site") %>%
    as.matrix()

  ### Design table
  #### This comes from Claude, for some reason limma was not being able to keep the same order
  #### of columns form the design and the abundance matrix. Although it worked fine outisde
  #### the function, it was not working inside the function. This should now fix it.
  design <- metadata %>%
    mutate(across(where(is.character), snakecase::to_snake_case))

  # Match matrix_norm columns to design samples to get their groups
  sample_order <- colnames(matrix_norm)
  group_order_from_matrix <- as.character(design$group[match(sample_order, design$sample)])

  # Get unique groups in the order they first appear in matrix_norm
  group_order <- character()
  for(sample in colnames(matrix_norm)) {
    group <- design$group[design$sample == sample]
    if(length(group) > 0 && !(group %in% group_order)) {
      group_order <- c(group_order, group)
    }
  }

  # Set factor levels based on this order
  design$group <- factor(design$group, levels = group_order)

  design_limma <- stats::model.matrix(~0 + design$group)
  colnames(design_limma) <- levels(design$group)
  rownames(design_limma) <- design$sample

  # Reorder matrix to match design
  matrix_norm <- matrix_norm[, rownames(design_limma)]

  ### Contrast table
  contrast_table <- fread(here::here(inputPath, "contrasts.csv")) %>%
    mutate(across(where(is.character), snakecase::to_snake_case))

  contrast_formulas <- vector()
  for(i in 1:nrow(contrast_table)){
    temp <- c(paste0(contrast_table[i,1], "-", contrast_table[i,2]))
    contrast_formulas <- c(contrast_formulas, temp)
  }

  contrast.matrix <- do.call(
    limma::makeContrasts,
    c(as.list(setNames(contrast_formulas, contrast_formulas)),
      list(levels = design_limma))
  )

  # Fit the data and create contrasts
  fit <- limma::lmFit(matrix_norm, design_limma)
  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit2)

  if (dim(fit$design)[2] <= 2) {

    print("Less than 2 samples, use a scatter plot to compare both samples.")

  } else {

    grDevices::pdf(here::here(outputPath, "output", "plots", paste0(jobname, "limma_mds-plot.pdf")))
    limma::plotMDS(fit)
    grDevices::dev.off()

  }

  # Saving contrasts
  df <- data.frame()
  for(i in colnames(contrast.matrix)){
    temp <- data.frame()
    temp <- limma::topTable(fit2, coef=i, number=Inf)
    temp$contrast <- i
    temp$modified_site <- rownames(temp)
    temp <- dplyr::relocate(temp, "modified_site", .before = "logFC")

    df <- rbind(df, temp)
  }

  # Putting protein names back
  df <- merge(df, norm[,c(1,2,3)], by = "modified_site") %>%
    relocate(protein_name, .after = modified_site) %>%
    relocate(protein_id, .before = protein_name) %>%
    mutate(protein_name = toupper(protein_name)) %>%
    mutate(site_name = paste0(protein_name,
                              "_",
                              str_split_fixed(modified_site, "_", n=Inf)[,2])) %>%
    relocate(site_name, .after = protein_id)

  df2 <- df %>%
    pivot_wider(
      id_cols = c(modified_site, protein_id, site_name, protein_name),
      names_from = contrast,
      values_from = c(logFC, P.Value, adj.P.Val),
      names_sep = "_"
    )

  if (proteinNormalization == TRUE){
  sheets <- list(
    "Metadata" = metadata,
    "NormalizedAbundances" = norm,
    "ProteinNormalizedAbundances" = prot_norm,
    "Contrasts" = as.data.frame(contrast_formulas),
    "limmaDEA_long" = df,
    "limmaDEA_wide" = df2,
    "UniqueSites" = uniquePerGroup
  )
  } else {
    sheets <- list(
      "Metadata" = metadata,
      "NormalizedAbundances" = norm,
      "Contrasts" = as.data.frame(contrast_formulas),
      "limmaDEA_long" = df,
      "limmaDEA_wide" = df2,
      "UniqueSites" = uniquePerGroup
    )
  }
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
