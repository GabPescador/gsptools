#' Imports SILAC PSM tables from Fragpipe-DDA mode and calculates labeling efficiency
#'
#' This function takes an inputPath as the 3_fragpipe/results folder and reads all "psm.tsv" to
#' calculate SILAC labeling efficiency. It outputs a dataframe with all stats and saves a histogram
#' and a scatter plot for each sample found.
#'
#' @param inputPath Input path directing to 3_fragpipe/results.
#' @param outputpath Output path to save all the outputs.
#' @return Generates the labeling efficiency stats and saves QC plots in outputPath.
#' @export

silacDDALabelingEfficiency <- function(inputPath,
                                       outputPath){

  # Peptide modification information
  file_paths <- list.files(inputPath,
                           pattern = "^psm\\.tsv$",
                           recursive = TRUE,
                           full.names = TRUE)

  df_list <- map(file_paths, fread)
  names(df_list) <- basename(dirname(file_paths))

  df_list2 <- list()
  for(i in names(df_list)){
    temp <- df_list[[i]] %>%
      filter(Qvalue <= 0.01) # FDR of 1% (0.01)

    df_list2[[i]] <- temp
  }

  results <- list()
  ratios <- list()
  for(i in names(df_list2)){
    # Identify heavy-labeled peptides (typically K+8 and R+10 for SILAC)
    data <- df_list2[[i]] %>%
      mutate(
        has_heavy_K = str_detect(`Assigned Modifications`, fixed("K(8.0")),
        has_heavy_R = str_detect(`Assigned Modifications`, fixed("R(10.0")),
        has_heavy_P = str_detect(`Assigned Modifications`, fixed("P(10.0")),
        is_heavy = has_heavy_K | has_heavy_R,
        is_heavy_P = has_heavy_P,
        is_light = !is_heavy
      )

    # Count heavy vs light peptides
    labeling_summary <- data %>%
      summarise(
        Total_PSMs = n(),
        Heavy_PSMs = sum(is_heavy),
        Light_PSMs = sum(is_light),
        Heavy_P_PSMs = sum(is_heavy_P),
        Heavy_Percentage = (Heavy_PSMs / Total_PSMs) * 100,
        Light_Percentage = (Light_PSMs / Total_PSMs) * 100,
        R_to_P_Percentage = (Heavy_P_PSMs / Total_PSMs) *100
      )

    # For peptides containing K
    k_peptides <- data %>%
      filter(str_detect(Peptide, "K")) %>%
      summarise(
        Total_with_K = n(),
        Heavy_K = sum(has_heavy_K),
        K_labeling_efficiency = (Heavy_K / Total_with_K) * 100
      )

    # For peptides containing R
    r_peptides <- data %>%
      filter(str_detect(Peptide, "R")) %>%
      summarise(
        Total_with_R = n(),
        Heavy_R = sum(has_heavy_R),
        R_labeling_efficiency = (Heavy_R / Total_with_R) * 100
      )

    # For peptides containing P
    p_peptides <- data %>%
      filter(str_detect(Peptide, "P")) %>%
      summarise(
        Total_with_P = n(),
        Heavy_P = sum(has_heavy_P),
        P_labeling_efficiency = (Heavy_P / Total_with_P) * 100
      )

    labeling_efficiency <- bind_rows(
      tibble(Amino_Acid = "Lysine (K)", Efficiency = k_peptides$K_labeling_efficiency),
      tibble(Amino_Acid = "Arginine (R)", Efficiency = r_peptides$R_labeling_efficiency),
      tibble(Amino_Acid = "Proline (P)", Efficiency = p_peptides$P_labeling_efficiency)
    )

    # Using intensity data
    ratio_analysis <- data %>%
      group_by(Protein = `Protein ID`) %>%  # Adjust column name as needed
      summarise(
        sample = i,
        Heavy_Intensity = sum(Intensity[is_heavy], na.rm = TRUE),
        Light_Intensity = sum(Intensity[is_light], na.rm = TRUE),
        H_L_Ratio = Heavy_Intensity / Light_Intensity,
        Log2_Ratio = log2(H_L_Ratio)
      ) %>%
      filter(Heavy_Intensity > 0 & Light_Intensity > 0)

    # Summary statistics
    ratio_stats <- ratio_analysis %>%
      summarise(
        Median_Ratio = median(H_L_Ratio, na.rm = TRUE),
        Mean_Ratio = mean(H_L_Ratio, na.rm = TRUE),
        Median_Log2_Ratio = median(Log2_Ratio, na.rm = TRUE)
      )

    # Plot 2: Heavy/Light ratio distribution
    p1 <- ggplot(ratio_analysis, aes(x = Log2_Ratio)) +
      geom_histogram(bins = 50, fill = "black", alpha = 0.7) +
      geom_vline(xintercept = 0, color = "#AA4499", size = 1) +
      geom_vline(xintercept = c(-1.22,1.22), linetype = "dashed", color = "#AA4499", size = 1) +
      labs(title = paste0(i, " Heavy/Light Ratio Distribution"),
           subtitle = paste0("Observed median = ",
                             round(ratio_stats$Median_Log2_Ratio, 2)),
           x = "Log2(Heavy/Light)",
           y = "Count") +
      theme_minimal()

    # Plot 3: Heavy vs Light intensity scatter
    p2 <- ggplot(ratio_analysis, aes(x = log2(Light_Intensity), y = log2(Heavy_Intensity))) +
      geom_point() +
      geom_abline(slope = 1, intercept = 0, color = "#AA4499") +
      geom_abline(slope = 1, intercept = c(-1.22,1.22), color = "#AA4499", linetype = "dashed") +
      labs(title = paste0(i, " Heavy vs Light Intensity"),
           x = "Log10(Light Intensity)",
           y = "Log10(Heavy Intensity)") +
      theme_minimal()

    # Saving LE plots
    pdf(here(outputPath, "plots", paste0(i, "_LE.pdf")), width = 8, height = 3)
    print(p1 + p2)
    dev.off()

    rm(p1,p2)

    # Find proteins with mixed heavy/light peptides
    mixed_labeling <- data %>%
      group_by(Protein = `Protein ID`) %>%
      summarise(
        Has_Heavy = any(is_heavy),
        Has_Light = any(is_light),
        Total_Peptides = n(),
        Heavy_Peptides = sum(is_heavy),
        Light_Peptides = sum(is_light)
      ) %>%
      filter(Has_Heavy & Has_Light) %>%
      mutate(Labeling_Completeness = (Heavy_Peptides / Total_Peptides) * 100)

    # Proteins with poor labeling
    poor_labeling <- mixed_labeling %>%
      filter(Labeling_Completeness < 90) %>%
      arrange(Labeling_Completeness)

    # Create comprehensive summary report
    qc_report <- tibble(
      Metric = c(
        "Total PSMs",
        "Heavy PSMs (%)",
        "Light PSMs (%)",
        "K Labeling Efficiency (%)",
        "R Labeling Efficiency (%)",
        "R>P Conversion Rate (%)",
        "Median H/L Ratio"
      ),
      !!i := c(
        labeling_summary$Total_PSMs,
        round(labeling_summary$Heavy_Percentage, 2),
        round(labeling_summary$Light_Percentage, 2),
        round(k_peptides$K_labeling_efficiency, 2),
        round(r_peptides$R_labeling_efficiency, 2),
        round(p_peptides$P_labeling_efficiency, 2),
        round(log2(ratio_stats$Median_Ratio), 2)
      )
    )

    results[[i]] <- qc_report
    ratios[[i]] <- ratio_analysis

  }

  ratio_bind <- ratios %>%
    map(select, sample, Log2_Ratio) %>%
    bind_rows()

  # Merge all dataframes by the Metric column
  combined_df <- purrr::reduce(results, left_join, by = "Metric")

  # Create a lookup for expected ratios
  expected_ratios <- tibble(
    Metric = "Expected H/L Ratio"
  ) %>%
    bind_cols(
      map_dfc(names(combined_df)[-1], function(col_name) {
        ratio <- case_when(
          str_detect(tolower(col_name), "30h.*70l|30.*70") ~ log2(30/70),
          str_detect(tolower(col_name), "50.*50") ~ log2(50/50),
          str_detect(tolower(col_name), "70h.*30l|70.*30") ~ log2(70/30),
          TRUE ~ NA_real_
        )
        tibble(!!col_name := ratio)
      })
    )

  combined_df <- bind_rows(combined_df, expected_ratios)

  # Extract the two rows as numeric vectors
  observed_values <- combined_df %>%
    filter(Metric == "Median H/L Ratio") %>%
    select(-Metric) %>%
    as.numeric()

  expected_values <- combined_df %>%
    filter(Metric == "Expected H/L Ratio") %>%
    select(-Metric) %>%
    as.numeric()

  # Calculate deviations
  deviations <- abs(observed_values - expected_values)

  # Create new row
  deviation_row <- tibble(
    Metric = "Ratio Deviation",
    !!!setNames(as.list(deviations), names(combined_df)[-1])
  )

  combined_df <- bind_rows(combined_df, deviation_row)

  write_csv(combined_df, here(outputPath, "tables", "SILAC_labeling_efficiency_stats.csv"))

  p3 <- combined_df %>%
    reshape2::melt() %>%
    filter(Metric %in% c("Heavy PSMs (%)", "Light PSMs (%)")) %>%
    ggplot(aes(x = variable, y = value, fill = Metric)) +
    geom_col(position = "dodge") +
    geom_hline(yintercept = c(30,50,70), linetype = "dashed") +
    scale_fill_manual(values = c("#AA4499", "grey")) +
    labs(title = paste0("SILAC Labeling Efficiency"),
         y = "Labeling Efficiency (%)",
         x = "") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, 100)

  p4 <- combined_df %>%
    reshape2::melt() %>%
    filter(Metric %in% c("K Labeling Efficiency (%)", "R Labeling Efficiency (%)", "R>P Conversion Rate (%)")) %>%
    ggplot(aes(x = variable, y = value, fill = Metric)) +
    geom_col(position = "dodge") +
    geom_hline(yintercept = c(30,50,70), linetype = "dashed") +
    scale_fill_manual(values = c("black", "grey", "#AA4499")) +
    labs(title = paste0("SILAC Labeling Efficiency"),
         y = "Labeling Efficiency (%)",
         x = "") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, 100)

  p5 <- combined_df %>%
    reshape2::melt() %>%
    filter(Metric %in% c("Median H/L Ratio")) %>%
    ggplot(aes(x = variable, y = value)) +
    geom_col(fill="black") +
    geom_hline(yintercept = c(-1.22,0,1.22), linetype = "dashed") +
    labs(title = "SILAC Ratio",
         y = "Heavy/Light Ratio",
         x = "") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p6 <- ratio_bind %>%
    ggplot(aes(x = sample, y = Log2_Ratio)) +
    geom_boxplot() +
    geom_hline(yintercept = c(-1.22,0,1.22), linetype = "dashed") +
    labs(title = "SILAC Ratio",
         y = "Heavy/Light Log2 Ratio",
         x = "") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  pdf(here(outputPath, "plots", "SILAC_stats.pdf"), width = 5, height = 4)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  dev.off()

  return(list(
    "DF" = combined_df,
    "LE1" = p3,
    "LE2" = p4,
    "Ratios" = p5,
    "ratio_boxplot" = p6)
    )

}
