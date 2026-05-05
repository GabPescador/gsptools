#' Imports SILAC diGly data from Fragpipe and performs analysis.
#'
#' This function takes an inputPath with "metadata.csv" and "combined_modified_peptide_label_quant.tsv"
#' generated, imports search results, combines site localization, perform t.test comparisons based on replicates
#' and outputs result files as an excel table for later import in ShinyApp.
#'
#' @param inputPath Input path where "metadata.csv", "contrast.csv", "digly_combined_modified_peptide_label_quant.tsv" and "input_combined_modified_peptide_label_quant.tsv" are located.
#' @param psmPath Path to where all psm.tsv files are for the digly search.
#' @param jobname Specifies the folder name for convenience I always use the PROT-XXXX.
#' @param outputpath Output path to save all the outputs.
#' @param replicateFilter Defaults to FALSE. Filters proteins that are not present in at least 2 replicates. Since majority of SILAC experiments only have 2 replicates (swaps) this is FALSE by default.
#' @param proteinInput Defaults to TRUE. Includes protein input for comparisons.
#' @param groups Defaults to FALSE. If TRUE, will use the type column in metadata.csv to search unique proteins based on those groupings.
#' @param force Defaults to FALSE. Forces saving the results table in case it already exists.
#' @param exclude Defaults to NULL. Defines a character vector of columns to be excluded from the pipeline.
#' @return Generates the results files and saves them in the outputPath.
#' @export

importFragpipeSILACdiGly <- function(inputPath,
                                     psmPath,
                                     jobname,
                                     outputPath,
                                     replicateFilter = FALSE,
                                     proteinInput = TRUE,
                                     groups = FALSE,
                                     force = FALSE,
                                     exclude = NULL){

  if (proteinInput == TRUE){
    # Check that tables are in the input folder
    if (file.exists(paste0(inputPath, "metadata.csv")) &
        # file.exists(paste0(inputPath, "contrasts.csv")) &
        file.exists(paste0(inputPath, "combined_modified_peptide_label_quant.tsv")) &
        # file.exists(Sys.glob(paste0(inputPath, "diGly_*normalized.txt"))) &
        file.exists(paste0(inputPath, "combined_protein_label_quant.tsv"))
        # file.exists(Sys.glob(paste0(inputPath, "input_*normalized.txt"))))
        ) {

      print("Input tables exists, proceeding with pipeline...")

    } else {

      stop("No metadata.csv, contrasts.csv, combined_protein_label_quant.tsv and/or combined_modified_peptide_label_quant.tsv found, stopping execution.")

    }
  } else {

  # Check that tables are in the input folder
  if (file.exists(paste0(inputPath, "metadata.csv")) &
      # file.exists(paste0(inputPath, "contrasts.csv")) &
      file.exists(Sys.glob(paste0(inputPath, "combined_modified_peptide_label_quant.tsv")))) {

    print("Input tables exists, proceeding with pipeline...")

  } else {

    stop("No metadata.csv, contrasts.csv and/or combined_modified_peptide_label_quant.tsv found, stopping execution.")

    }
  }
  # Create output directories in case they don't exist
  paths <- c(
    file.path(outputPath, "plots"),
    file.path(outputPath, "tables")
  )
  invisible(sapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))

  # Import digly
  diGly <- data.table::fread(here::here(inputPath, "combined_modified_peptide_label_quant.tsv")) %>%
    setNames(snakecase::to_snake_case(names(.))) %>%
    dplyr::filter(!str_detect(protein, "^contam"))

  # Find all psm.tsv files recursively in a directory
  files <- list.files(
    path = here(psmPath),
    pattern = "^psm\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  # Read and bind all files into one dataframe and filter only rows that have K(114)
  df <- files %>%
    lapply(readr::read_tsv) %>%
    dplyr::bind_rows() %>%
    setNames(snakecase::to_snake_case(names(.))) %>%
    #filter(str_detect(k_114_0429, "\\(")) %>%
    dplyr::select(peptide, extended_peptide, k_114_0429) %>%
    dplyr::distinct()

  diGly <- merge(diGly, df, by.x="peptide_sequence", by.y="peptide") %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      # 1. Find the starting index of the modification string "[114.0429]"
      # We use fixed() because "[" is a special character in regex
      mod_index = str_locate(k_114_0429, fixed("(1.0000"))[, 1],
      # 2. The residue index in the sequence is the character right before the "["
      # This is mod_index - 1
      residue_idx_in_peptide = mod_index - 1,
      # 3. Calculate absolute position: start + (index - 1)
      # This simplifies to: start + mod_index - 2
      absolute_mod_pos = start + (residue_idx_in_peptide - 1)
    ) %>%
    dplyr::mutate(modified_site = paste0(gene, "_K", absolute_mod_pos)) %>%
    dplyr::relocate(modified_site, .before = peptide_sequence) %>%
    dplyr::filter(!is.na(k_114_0429)) %>% # filter out any psm that does not have K(114)
    select(modified_site, protein_id, gene, peptide_sequence, modified_peptide, contains("max_lfq"), contains("log_2_ratio"))

  if (proteinInput == TRUE){
  input <- data.table::fread(here::here(Sys.glob(paste0(inputPath, "combined_protein_label_quant.tsv")))) %>%
    setNames(snakecase::to_snake_case(names(.))) %>%
    filter(!str_detect(protein, "contam"))
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

  # Taking out sites that are not present in at least 2 replicates if replicateFilter == TRUE
  if (replicateFilter == TRUE){

    # diGly
    diGly_long <- diGly %>%
      select(modified_site:modified_peptide, contains("log_2")) %>%
      reshape2::melt(id.vars = colnames(.)[1:5]) %>%
      mutate(variable = sub("_log_2.*", "", variable)) %>%
      base::merge(., metadata[,-4], by.x = "variable", by.y = "sample")

    diGly_filtered <- diGly_long %>%
      filter(!is.na(value)) %>%
      group_by(modified_site, group) %>%
      summarize(counts =n())

    diGly_long <- base::merge(diGly_long, diGly_filtered, by = c("modified_site", "group"))

    diGly_long <- diGly_long %>%
      mutate(value = if_else(counts < 2, NA_real_, value))

    diGly <- diGly_long %>%
      tidyr::pivot_wider(
        id_cols = c(modified_site, protein_id, gene, peptide_sequence, modified_peptide, group, experiment),
        names_from = variable,
        values_from = c(value),
        names_sep = "_"
      )

    # Input
    input_long <- input %>%
      select(protein_id, gene, contains("median_log_2")) %>%
      reshape2::melt(id.vars = colnames(.)[1:2]) %>%
      mutate(variable = sub("_median_log_2.*", "", variable)) %>%
      base::merge(., metadata[,-4], by.x = "variable", by.y = "sample")

    input_filtered <- input_long %>%
      filter(!is.na(value)) %>%
      group_by(protein_id, group) %>%
      summarize(counts = n())

    input_long <- base::merge(input_long, input_filtered, by = c("protein_id", "group"))

    input_long <- input_long %>%
      mutate(value = if_else(counts < 2, NA_real_, value))

    input <- input_long %>%
      tidyr::pivot_wider(
        id_cols = c(protein_id, gene),
        names_from = variable,
        values_from = c(value),
        names_sep = "_"
      )

  }

  if (groups == FALSE) {
  # Creating a list of unique proteins based on groups
    # diGly
  uniquePerGroup_diGly <- diGly %>%
    select(modified_site:modified_peptide, contains("log_2")) %>%
    reshape2::melt(id.vars = colnames(.)[1:5]) %>%
    mutate(variable = sub("_log_2.*", "", variable)) %>%
    base::merge(., metadata[,-4], by.x = "variable", by.y = "sample") %>%
    filter(!is.na(value)) %>%
    group_by(modified_site) %>%
    filter(n_distinct(group) == 1) %>%
    ungroup() %>%
    distinct(modified_site, group)

  uniquePerGroup_diGly <- base::merge(uniquePerGroup_diGly, diGly[,c("modified_site", "protein_id", "gene")]) %>%
    relocate(gene, .before = group)

  # input
  uniquePerGroup_input <- input %>%
    select(protein_id, gene, contains("median_log_2")) %>%
    reshape2::melt(id.vars = colnames(.)[1:2]) %>%
    mutate(variable = sub("_median_log_2.*", "", variable)) %>%
    base::merge(., metadata[,-4], by.x = "variable", by.y = "sample") %>%
    filter(!is.na(value)) %>%
    group_by(protein_id) %>%
    filter(n_distinct(group) == 1) %>%
    ungroup() %>%
    distinct(protein_id, group)

  uniquePerGroup_input <- base::merge(uniquePerGroup_input, input[,c("protein_id", "gene")]) %>%
    relocate(gene, .before = group)

  } else {

    # diGly
    uniquePerGroup_diGly <- diGly %>%
      select(modified_site:modified_peptide, contains("log_2")) %>%
      reshape2::melt(id.vars = colnames(.)[1:5]) %>%
      mutate(variable = sub("_log_2.*", "", variable)) %>%
      base::merge(., metadata[,-4], by.x = "variable", by.y = "sample") %>%
      filter(!is.na(value)) %>%
      group_by(modified_site) %>%
      filter(n_distinct(group) == 1) %>%
      ungroup() %>%
      distinct(modified_site, group)

    # Creating a list of unique proteins based on type
    uniquePerType_diGly <- diGly %>%
      select(modified_site:modified_peptide, contains("log_2")) %>%
      reshape2::melt(id.vars = colnames(.)[1:5]) %>%
      mutate(variable = sub("_log_2.*", "", variable)) %>%
      base::merge(., metadata[,-4], by.x = "variable", by.y = "sample") %>%
      filter(!is.na(value)) %>%
      group_by(modified_site) %>%
      filter(n_distinct(type) == 1) %>%
      ungroup() %>%
      distinct(modified_site, type) %>%
      rename("group" = "type")

    uniquePerGroup_diGly <- rbind(uniquePerGroup_diGly, uniquePerType_diGly)
    uniquePerGroup_diGly <- base::merge(uniquePerGroup_diGly, diGly[,c("modified_site", "protein_id", "gene")]) %>%
      relocate(gene, .before = group)

    # input
    uniquePerGroup_input <- input %>%
      select(protein_id, gene, contains("median_log_2")) %>%
      reshape2::melt(id.vars = colnames(.)[1:2]) %>%
      mutate(variable = sub("_median_log_2.*", "", variable)) %>%
      base::merge(., metadata[,-4], by.x = "variable", by.y = "sample") %>%
      filter(!is.na(value)) %>%
      group_by(protein_id) %>%
      filter(n_distinct(group) == 1) %>%
      ungroup() %>%
      distinct(protein_id, group)

    # Creating a list of unique proteins based on groups
    uniquePerType_input <- input %>%
      select(protein_id, gene, contains("median_log_2")) %>%
      reshape2::melt(id.vars = colnames(.)[1:2]) %>%
      mutate(variable = sub("_median_log_2.*", "", variable)) %>%
      base::merge(., metadata[,-4], by.x = "variable", by.y = "sample") %>%
      filter(!is.na(value)) %>%
      group_by(protein_id) %>%
      filter(n_distinct(type) == 1) %>%
      ungroup() %>%
      distinct(protein_id, type) %>%
      rename("group" = "type")

    uniquePerGroup_input <- rbind(uniquePerGroup_input, uniquePerType_input)
    uniquePerGroup_input <- base::merge(uniquePerGroup_input, input[,c("protein_id", "gene")]) %>%
      relocate(gene, .before = group)
  }

  # t.tests based on number of replicates
  # diGly
  diGly_long <- diGly %>%
    select(modified_site:modified_peptide, contains("log_2")) %>%
    reshape2::melt(id.vars = colnames(.)[1:5]) %>%
    mutate(variable = sub("_log_2.*", "", variable)) %>%
    base::merge(., metadata, by.x = "variable", by.y = "sample") %>%
    mutate(value = case_when(
      negate == TRUE ~ value*-1,
      negate == FALSE ~ value
    )) %>%
    group_by(group, experiment, modified_site) %>%
    summarize(
      p.val = tryCatch(
        t.test(value, mu = 0)$p.value,
        error = function(e) NA_real_
      ),
      mean_ratio = mean(value, na.rm = TRUE),
      n = sum(!is.na(value)),
      .groups = "drop"
    ) %>%
    dplyr::filter(n >= 1)

# diGly_df <- diGly_long %>%
#     pivot_wider(
#       id_cols = c(modified_site),
#       names_from = c(experiment, group),
#       values_from = c(mean_ratio, p.val),
#       names_sep = "_"
#     )

df2_digly <- diGly_long %>%
  filter(n >=1) %>%
  split(.$group)

  if (proteinInput == TRUE) {

    # diGly
    input_long <- input %>%
      select(protein_id, gene, contains("log_2")) %>%
      reshape2::melt(id.vars = colnames(.)[1:2]) %>%
      mutate(variable = sub("_median_log_2.*", "", variable)) %>%
      base::merge(., metadata, by.x = "variable", by.y = "sample") %>%
      mutate(value = case_when(
        negate == TRUE ~ value*-1,
        negate == FALSE ~ value
      )) %>%
      group_by(group, experiment, protein_id, gene) %>%
      summarize(
        p.val = tryCatch(
          t.test(value, mu = 0)$p.value,
          error = function(e) NA_real_
        ),
        mean_ratio = mean(value, na.rm = TRUE),
        n = sum(!is.na(value)),
        .groups = "drop"
      ) %>%
      dplyr::filter(n >= 1)

    # input_df <- input_long %>%
    #   pivot_wider(
    #     id_cols = c(protein_id, gene),
    #     names_from = c(experiment, group),
    #     values_from = c(mean_ratio, p.val),
    #     names_sep = "_"
    #   )

    df2_input <- input_long %>%
      dplyr::filter(n >= 1) %>%
      split(.$group)

  }

  # Description of result tables
  # description <- data.table::fread(system.file("extdata", "importFragpipeSILACdigly_examples", "outputDescription.csv",
  #                               package = "gsptools"))
  description <- fread("/home/gd2417/useful_things/Scripts/gsptools/inst/extdata/importFragpipeSILACdigly_examples/outputDescription.csv")

  if (proteinInput == TRUE){
  sheets_input <- c(
    list(
      "Description" = description,
      "Metadata" = metadata,
      "Abundances" = input,
      "OneSample-ttest" = input_long,
      "UniqueSites" = uniquePerGroup_input
    ),
  df2_input)

  sheets_diGly <- c(
    list(
      "Description" = description,
      "Metadata" = metadata,
      "Abundances" = diGly,
      "OneSample-ttest" = diGly_long,
      "UniqueSites" = uniquePerGroup_diGly
    ),
    df2_digly)

  } else {
    sheets_diGly <- c(
      list(
        "Description" = description,
        "Metadata" = metadata,
        "Abundances" = diGly,
        "OneSample-ttest" = diGly_long,
        "UniqueSites" = uniquePerGroup_diGly
      ),
      df2_digly)
  }

  if (force == TRUE) {
    if(proteinInput == TRUE){
      #input
      writexl::write_xlsx(x = sheets_input,
                          path = here::here(outputPath, "tables", paste0(jobname, "_input_results.xlsx")),
                          format_headers = FALSE)

      for(i in names(sheets_input)) {
      readr::write_csv(x = sheets_input[[i]],
                       file = here::here(outputPath, "tables", paste0(jobname, "_input_", i, ".csv")))
      }

      #diGly
      writexl::write_xlsx(x = sheets_diGly,
                          path = here::here(outputPath, "tables", paste0(jobname, "_diGly_results.xlsx")),
                          format_headers = FALSE)

      for(i in names(sheets_diGly)) {
        readr::write_csv(x = sheets_diGly[[i]],
                         file = here::here(outputPath, "tables", paste0(jobname, "_diGly_", i, ".csv")))
      }

    } else {
      writexl::write_xlsx(x = sheets_diGly,
                          path = here::here(outputPath, "tables", paste0(jobname, "_diGly_results.xlsx")),
                          format_headers = FALSE)

      for(i in names(sheets_diGly)) {
        readr::write_csv(x = sheets_diGly[[i]],
                         file = here::here(outputPath, "tables", paste0(jobname, "_diGly_", i, ".csv")))
      }
    }
  } else if (force == FALSE & file.exists(here::here(outputPath, "tables", paste0(jobname, "_diGly_results.xlsx")))) {

    print(paste0("Contrasts already performed and can be found in ",
                 here::here(outputPath, "tables", paste0(jobname, "_diGly_results.xlsx"))))

    } else if (force == FALSE & !file.exists(here::here(outputPath, "tables", paste0(jobname, "_diGly_results.xlsx")))) {
      if(proteinInput == TRUE){
      #input
      writexl::write_xlsx(x = sheets_input,
                          path = here::here(outputPath, "tables", paste0(jobname, "_input_results.xlsx")),
                          format_headers = FALSE)

      for(i in names(sheets_input)) {
      readr::write_csv(x = sheets_input[[i]],
                       file = here::here(outputPath, "tables", paste0(jobname, "_input_", i, ".csv")))
      }

      #diGly
      writexl::write_xlsx(x = sheets_diGly,
                          path = here::here(outputPath, "tables", paste0(jobname, "_diGly_results.xlsx")),
                          format_headers = FALSE)

      for(i in names(sheets_diGly)) {
        readr::write_csv(x = sheets_diGly[[i]],
                         file = here::here(outputPath, "tables", paste0(jobname, "_diGly_", i, ".csv")))
      }

      } else {
        writexl::write_xlsx(x = sheets_diGly,
                            path = here::here(outputPath, "tables", paste0(jobname, "_diGly_results.xlsx")),
                            format_headers = FALSE)

        for(i in names(sheets_diGly)) {
          readr::write_csv(x = sheets_diGly[[i]],
                          file = here::here(outputPath, "tables", paste0(jobname, "_diGly_", i, ".csv")))
        }
      }

    }

}

