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
#' @param outputpath Output path to save all the outputs.
#' @param replicateFilter Defaults to TRUE. Filters proteins that are not present in at least 2 replicates.
#' @param proteinInput Defaults to TRUE. Includes protein input for limma contrasts.
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
                             proteinInput = FALSE,
                             groups = FALSE,
                             force = FALSE,
                             exclude = NULL){

  if (proteinInput == TRUE){
    # Check that tables are in the input folder
    if (file.exists(paste0(inputPath, "metadata.csv")) &
        file.exists(paste0(inputPath, "contrasts.csv")) &
        file.exists(paste0(inputPath, "abundance_single-site_MD.tsv")) &
        file.exists(Sys.glob(paste0(inputPath, "diGly_*normalized.txt"))) &
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
      file.exists(Sys.glob(paste0(inputPath, "diGly_*normalized.txt"))) &
      file.exists(paste0(inputPath, "abundance_single-site_MD.tsv"))) {

    print("Input tables exists, proceeding with pipeline...")

  } else {

    stop("No metadata.csv, contrasts.csv, abundance_single-site_MD.tsv and/or normalized.txt found, stopping execution.")

    }
  }
  # Create output directories in case they don't exist
  paths <- c(
    file.path(outputPath, "output", "plots"),
    file.path(outputPath, "output", "tables")
  )
  invisible(sapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))

  # Import normalized values
  diGly <- data.table::fread(here::here(Sys.glob(paste0(inputPath, "diGly_*normalized.txt")))) %>%
        setNames(snakecase::to_snake_case(names(.)))

  if (proteinInput == TRUE){
  input <- data.table::fread(here::here(Sys.glob(paste0(inputPath, "input_*normalized.txt")))) %>%
    setNames(snakecase::to_snake_case(names(.)))
  }

  # Exclude replicates from exclude parameter
    if (!is.null(exclude) & proteinInput == TRUE) {

      diGly <- diGly %>%
        select(-exclude)

      input <- input %>%
        select(-exclude)

      metadata <- data.table::fread(here::here(inputPath, "metadata.csv")) %>%
        mutate(across(where(is.character), snakecase::to_snake_case)) %>%
        filter(!sample %in% exclude)

    } else if (!is.null(exclude) & proteinInput == FALSE) {

      diGly <- diGly %>%
        select(-exclude)

      metadata <- data.table::fread(here::here(inputPath, "metadata.csv")) %>%
        mutate(across(where(is.character), snakecase::to_snake_case)) %>%
        filter(!sample %in% exclude)

    } else {

      metadata <- data.table::fread(here::here(inputPath, "metadata.csv")) %>%
        mutate(across(where(is.character), snakecase::to_snake_case))

    }

  # QC for before and after normalization
  # Before normalization diGly
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
    ggplot2::ggplot(aes(x=.data$variable, y=.data$value)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggplot2::xlab("") +
    ggplot2::ylab("Log2(Abundance)") +
    ggplot2::ggtitle("Before Normalization")

  # After normalization
  p2 <- diGly %>%
    reshape2::melt(id.vars = colnames(diGly)[1:7]) %>%
    ggplot2::ggplot(aes(x=.data$variable, y=.data$value)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    ggplot2::xlab("") +
    ggplot2::ylab("Log2(Abundance)") +
    ggplot2::ggtitle("After Normalization")

    # Plots together
    grDevices::pdf(here::here(outputPath, "output", "plots", paste0(jobname, "_diGly_normalization_boxplot.pdf")))
    print(cowplot::plot_grid(p1, p2, ncol = 2))
    grDevices::dev.off()

    # Before normalization diGly
    tmt <- data.table::fread(here::here(inputPath, "abundance_protein_MD.tsv")) %>%
      dplyr::rename("ProteinID" = "Index",
                    "ProteinName" = "Gene") %>%
      mutate(ProteinName = tolower(ProteinName)) %>%
      filter(!stringr::str_detect(.data$ProteinID, "Cont")) %>%
      filter(!stringr::str_detect(.data$ProteinID, "rev")) %>%
      setNames(snakecase::to_snake_case(names(.)))

    p1 <- tmt %>%
      reshape2::melt(id.vars = colnames(tmt)[1:5]) %>%
      ggplot2::ggplot(aes(x=.data$variable, y=.data$value)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggplot2::xlab("") +
      ggplot2::ylab("Log2(Abundance)") +
      ggplot2::ggtitle("Before Normalization")

    # After normalization
    p2 <- input %>%
      reshape2::melt(id.vars = colnames(input)[1:3]) %>%
      ggplot2::ggplot(aes(x=.data$variable, y=.data$value)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      ggplot2::xlab("") +
      ggplot2::ylab("Log2(Abundance)") +
      ggplot2::ggtitle("After Normalization")

    # Plots together
    grDevices::pdf(here::here(outputPath, "output", "plots", paste0(jobname, "_input_normalization_boxplot.pdf")))
    print(cowplot::plot_grid(p1, p2, ncol = 2))
    grDevices::dev.off()

  # Taking out sites that are not present in at least 2 replicates if replicateFilter == TRUE
  if (replicateFilter == TRUE){

    # diGly
    diGly_long <- diGly %>%
      reshape2::melt(id.vars = colnames(.)[1:8]) %>%
      # select(protein_id, variable, value) %>%
      merge(., metadata[,-3], by.x = "variable", by.y = "sample")

    diGly_filtered <- diGly_long %>%
      filter(!is.na(value)) %>%
      group_by(modified_site, group) %>%
      summarize(counts =n())

    diGly_long <- merge(diGly_long, diGly_filtered, by = c("modified_site", "group"))

    diGly_long <- diGly_long %>%
      mutate(value = if_else(counts < 2, NA_real_, value))

    diGly <- diGly_long %>%
      tidyr::pivot_wider(
        id_cols = c(modified_site, protein_id, protein_name, peptide, sequence_window, start, end, max_pep_prob),
        names_from = variable,
        values_from = c(value),
        names_sep = "_"
      )

    # Input
    input_long <- input %>%
      reshape2::melt(id.vars = colnames(.)[1:3]) %>%
      # select(protein_id, variable, value) %>%
      merge(., metadata[,-3], by.x = "variable", by.y = "sample")

    input_filtered <- input_long %>%
      filter(!is.na(value)) %>%
      group_by(protein_id, group) %>%
      summarize(counts =n())

    input_long <- merge(input_long, input_filtered, by = c("protein_id", "group"))

    input_long <- input_long %>%
      mutate(value = if_else(counts < 2, NA_real_, value))

    input <- input_long %>%
      tidyr::pivot_wider(
        id_cols = c(protein_id, protein_name, reference_intensity),
        names_from = variable,
        values_from = c(value),
        names_sep = "_"
      )

  }

  if (groups == FALSE) {
  # Creating a list of unique proteins based on groups
    # diGly
  uniquePerGroup_diGly <- diGly %>%
    reshape2::melt(id.vars = colnames(.)[1:8]) %>%
    select(modified_site, protein_id, variable, value) %>%
    merge(., metadata[,-3], by.x = "variable", by.y = "sample") %>%
    filter(!is.na(value)) %>%
    group_by(modified_site) %>%
    filter(n_distinct(group) == 1) %>%
    ungroup() %>%
    distinct(modified_site, group)

  uniquePerGroup_diGly <- merge(uniquePerGroup_diGly, diGly[,c("modified_site", "protein_id", "protein_name")]) %>%
    relocate(protein_name, .before = group)

  # input
  uniquePerGroup_input <- input %>%
    reshape2::melt(id.vars = colnames(.)[1:3]) %>%
    select(protein_id, variable, value) %>%
    merge(., metadata[,-3], by.x = "variable", by.y = "sample") %>%
    filter(!is.na(value)) %>%
    group_by(protein_id) %>%
    filter(n_distinct(group) == 1) %>%
    ungroup() %>%
    distinct(protein_id, group)

  uniquePerGroup_input <- merge(uniquePerGroup_input, input[,c("protein_id", "protein_name")]) %>%
    relocate(protein_name, .before = group)

  } else {

    # diGly
    uniquePerGroup_diGly <- diGly %>%
      reshape2::melt(id.vars = colnames(.)[1:8]) %>%
      select(modified_site, protein_id, variable, value) %>%
      merge(., metadata[,-3], by.x = "variable", by.y = "sample") %>%
      filter(!is.na(value)) %>%
      group_by(modified_site) %>%
      filter(n_distinct(group) == 1) %>%
      ungroup() %>%
      distinct(modified_site, group)

    # Creating a list of unique proteins based on groups
    uniquePerType_diGly <- diGly %>%
      reshape2::melt(id.vars = colnames(.)[1:8]) %>%
      select(modified_site, protein_id, variable, value) %>%
      merge(., metadata[,-3], by.x = "variable", by.y = "sample") %>%
      filter(!is.na(value)) %>%
      group_by(modified_site) %>%
      filter(n_distinct(type) == 1) %>%
      ungroup() %>%
      distinct(modified_site, type) %>%
      rename("group" = "type")

    uniquePerGroup_diGly <- rbind(uniquePerGroup_diGly, uniquePerType_diGly)
    uniquePerGroup_diGly <- merge(uniquePerGroup_diGly, diGly[,c("modified_site", "protein_id", "protein_name")]) %>%
      relocate(protein_name, .before = group)

    # input
    uniquePerGroup_input <- input %>%
      reshape2::melt(id.vars = colnames(.)[1:3]) %>%
      select(protein_id, variable, value) %>%
      merge(., metadata[,-3], by.x = "variable", by.y = "sample") %>%
      filter(!is.na(value)) %>%
      group_by(protein_id) %>%
      filter(n_distinct(group) == 1) %>%
      ungroup() %>%
      distinct(protein_id, group)

    # Creating a list of unique proteins based on groups
    uniquePerType_input <- input %>%
      reshape2::melt(id.vars = colnames(.)[1:3]) %>%
      select(protein_id, variable, value) %>%
      merge(., metadata[,-3], by.x = "variable", by.y = "sample") %>%
      filter(!is.na(value)) %>%
      group_by(protein_id) %>%
      filter(n_distinct(type) == 1) %>%
      ungroup() %>%
      distinct(protein_id, type) %>%
      rename("group" = "type")

    uniquePerGroup_input <- rbind(uniquePerGroup_input, uniquePerType_input)
    uniquePerGroup_input <- merge(uniquePerGroup_input, input[,c("protein_id", "protein_name")]) %>%
      relocate(protein_name, .before = group)
  }

  # Matrix for limma
  ### Matrix table
    # diGly
  matrix_diGly <- diGly %>%
    select(modified_site, max_pep_prob:last_col(), -max_pep_prob) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("modified_site") %>%
    as.matrix()

  ### Design table
  #### This comes from Claude, for some reason limma was not being able to keep the same order
  #### of columns from the design and the abundance matrix. Although it worked fine outisde
  #### the function, it was not working inside the function. This should now fix it.
  design <- metadata %>%
    mutate(across(where(is.character), snakecase::to_snake_case))

  # Match matrix_norm columns to design samples to get their groups
  sample_order <- colnames(matrix_diGly)
  group_order_from_matrix <- as.character(design$group[match(sample_order, design$sample)])

  # Get unique groups in the order they first appear in matrix_norm
  group_order <- character()
  for(sample in colnames(matrix_diGly)) {
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
  matrix_diGly <- matrix_diGly[, rownames(design_limma)]

  ### Contrast table
  contrast_table <- data.table::fread(here::here(inputPath, "contrasts.csv")) %>%
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
  fit <- limma::lmFit(matrix_diGly, design_limma)
  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit2)

  grDevices::pdf(here::here(outputPath, "output", "plots", paste0(jobname, "_diGly_limma_mds-plot.pdf")))
  limma::plotMDS(matrix_diGly)
  grDevices::dev.off()

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
  df_diGly <- merge(df, diGly[,c(1,2,3)], by = "modified_site") %>%
    relocate(protein_name, .after = modified_site) %>%
    relocate(protein_id, .before = protein_name) %>%
    mutate(protein_name = toupper(protein_name)) %>%
    mutate(site_name = paste0(protein_name,
                              "_",
                              stringr::str_split_fixed(modified_site, "_", n=Inf)[,2])) %>%
    relocate(site_name, .after = protein_id)

  df2_diGly <- df_diGly %>%
    split(.$contrast)

  if (proteinInput == TRUE) {
    # Matrix for limma
    ### Matrix table
    # input
    matrix_input <- input %>%
      select(protein_id, reference_intensity:last_col(), -reference_intensity) %>%
      as.data.frame() %>%
      tibble::column_to_rownames("protein_id") %>%
      as.matrix()

    # Reorder matrix to match design
    matrix_input <- matrix_input[, rownames(design_limma)]

    # Fit the data and create contrasts
    fit <- limma::lmFit(matrix_input, design_limma)
    fit2 <- limma::contrasts.fit(fit, contrast.matrix)
    fit2 <- limma::eBayes(fit2)


    grDevices::pdf(here::here(outputPath, "output", "plots", paste0(jobname, "_input_limma_mds-plot.pdf")))
    limma::plotMDS(matrix_input)
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
    df_input <- merge(df, input[,c(1,2)], by = "protein_id") %>%
      relocate(protein_name, .after = protein_id) %>%
      mutate(protein_name = toupper(protein_name))

    df2_input <- df_input %>%
      split(.$contrast)
  }

  if (proteinInput == TRUE){
  sheets_input <- c(
    list(data.table::fread(system.file("inst", "extdata", "importFragpipeTMT_examples", "outputDescription.csv",
                                       package = "gsptools"))),
    list(
      "Metadata" = metadata,
      "NormalizedAbundances" = input,
      "Contrasts" = as.data.frame(contrast_formulas),
      "limmaDEA_long" = df_input,
      "UniqueSites" = uniquePerGroup_input
    ),
  df2_input)

  sheets_diGly <- c(
    list(data.table::fread(system.file("inst", "extdata", "importFragpipeTMT_examples", "outputDescription.csv",
                                       package = "gsptools"))),
    list(
      "Metadata" = metadata,
      "NormalizedAbundances" = diGly,
      "Contrasts" = as.data.frame(contrast_formulas),
      "limmaDEA_long" = df_diGly,
      "UniqueSites" = uniquePerGroup_diGly
    ),
  df2_diGly)

  } else {
    sheets_diGly <- c(
      list(data.table::fread(system.file("inst", "extdata", "importFragpipeTMT_examples", "outputDescription.csv",
                                         package = "gsptools"))),
      list(
        "Metadata" = metadata,
        "NormalizedAbundances" = diGly,
        "Contrasts" = as.data.frame(contrast_formulas),
        "limmaDEA_long" = df_diGly,
        "UniqueSites" = uniquePerGroup_diGly
      ),
    df2_diGly)
  }

  if (force == TRUE) {
    if(proteinInput == TRUE){
      #input
      writexl::write_xlsx(x = sheets_input,
                          path = here::here(outputPath, "output", "tables", paste0(jobname, "_input_results.xlsx")),
                          format_headers = FALSE)

      for(i in names(sheets_input)) {
      readr::write_csv(x = sheets_input[[i]],
                       file = here::here(outputPath, "output", "tables", paste0(jobname, "_input_", i, ".csv")))
      }

      #diGly
      writexl::write_xlsx(x = sheets_diGly,
                          path = here::here(outputPath, "output", "tables", paste0(jobname, "_diGly_results.xlsx")),
                          format_headers = FALSE)

      for(i in names(sheets_diGly)) {
        readr::write_csv(x = sheets_diGly[[i]],
                         file = here::here(outputPath, "output", "tables", paste0(jobname, "_diGly_", i, ".csv")))
      }

    } else {
      writexl::write_xlsx(x = sheets_diGly,
                          path = here::here(outputPath, "output", "tables", paste0(jobname, "_diGly_results.xlsx")),
                          format_headers = FALSE)

      for(i in names(sheets_diGly)) {
        readr::write_csv(x = sheets_diGly[[i]],
                         file = here::here(outputPath, "output", "tables", paste0(jobname, "_diGly_", i, ".csv")))
      }
    }
  } else if (file.exists(here::here(outputPath, "output", "tables", paste0(jobname, "_diGly_results.xlsx")))) {

    print(paste0("Contrasts already performed and can be found in ",
                 here::here(outputPath, "output", "tables", paste0(jobname, "_diGly_results.xlsx"))))

    }
}
