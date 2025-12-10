#' Generates an html report based on outputs from limmaFragpipeTMTdiGly function
#'
#' This function takes the output from limmaFragpipeTMTdiGly function and generates a report with
#' all results and tables.
#'
#' @param path Main path of the project, usually as PROT-XXXX and at the top level.
#' @param jobname Specifies the project name, for convenience I always use the PROT-XXXX.
#' @return Generates a report with all standard analyses for TMT diGly data.
#' @export

reportTMTdiGly <- function(path,
                           jobname){

  ########################################
  #        Quality Control Plots         #
  ########################################

  print("Generating quality control plots...")

  #       Normalization Boxplots         #
  ########################################

  # Define theme to use in boxplots
  norm_boxplot <- list(
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    )

  # load input data
  input_raw <- data.table::fread(Sys.glob(paste0(path, "/4_input/abundance_protein_MD.tsv"))) %>%
    select(Index, ReferenceIntensity:last_col(), -ReferenceIntensity) %>%
    reshape2::melt()

  input_norm <- data.table::fread(Sys.glob(paste0(path, "4_input/input*normalized.txt"))) %>%
    select(protein_id, reference_intensity:last_col(), -reference_intensity) %>%
    reshape2::melt()

  p1 <- input_raw %>%
    ggplot(aes(x=variable, y=value)) +
    geom_boxplot() +
    ylab("Log2(Abundance)") +
    xlab("") +
    ggtitle("Input Before Normalization") +
    norm_boxplot

  p2 <- input_norm %>%
    ggplot(aes(x=variable, y=value)) +
    geom_boxplot() +
    ylab("Log2(Abundance)") +
    xlab("") +
    ggtitle("Input After Normalization") +
    norm_boxplot

  ########################################
  QCp1 <- cowplot::plot_grid(p1, p2, ncol = 2)
  ########################################

  # Load diGly data
  diGly_raw <- data.table::fread(Sys.glob(paste0(path, "/4_input/abundance_single-site_MD.tsv"))) %>%
    select(Index, ReferenceIntensity:last_col(), -ReferenceIntensity) %>%
    reshape2::melt()

  diGly_norm <- data.table::fread(Sys.glob(paste0(path, "/4_input/diGly*normalized.txt"))) %>%
    select(modified_site, max_pep_prob:last_col(), -max_pep_prob) %>%
    reshape2::melt()

  p1 <- diGly_raw %>%
    ggplot(aes(x=variable, y=value)) +
    geom_boxplot() +
    ylab("Log2(Abundance)") +
    xlab("") +
    ggtitle("diGly Before Normalization") +
    norm_boxplot

  p2 <- diGly_norm %>%
    ggplot(aes(x=variable, y=value)) +
    geom_boxplot() +
    ylab("Log2(Abundance)") +
    xlab("") +
    ggtitle("diGly After Normalization") +
    norm_boxplot

  ########################################
  QCp2 <- cowplot::plot_grid(p1, p2, ncol = 2)
  ########################################

  #           Correlation Plots          #
  ########################################

  input_omit <- data.table::fread(Sys.glob(paste0(path, "/4_input/input*normalized.txt"))) %>%
    select(protein_id, reference_intensity:last_col(), -reference_intensity) %>%
    tibble::column_to_rownames("protein_id") %>%
    na.omit()

  correlation <- cor(input_omit, use = "complete.obs")

  p1 <- ComplexHeatmap::Heatmap(correlation,
                                name = "Correlation",
                                column_title = "Input Correlation",
                                col = circlize::colorRamp2(c(min(correlation), 1), c("white", "#882255")),
                                show_row_names = TRUE,
                                show_column_names = TRUE,
                                border = TRUE)

  diGly_omit <- data.table::fread(Sys.glob(paste0(path, "/4_input/diGly*normalized.txt"))) %>%
    select(modified_site, max_pep_prob:last_col(), -max_pep_prob) %>%
    tibble::column_to_rownames("modified_site") %>%
    na.omit()

  correlation <- cor(diGly_omit, use = "complete.obs")

  p2 <- ComplexHeatmap::Heatmap(correlation,
                                name = "Correlation",
                                column_title = "diGly Correlation",
                                col = circlize::colorRamp2(c(min(correlation), 1), c("white", "#882255")),
                                show_row_names = TRUE,
                                show_column_names = TRUE,
                                border = TRUE)

  ########################################
  QCp3 <- p1
  QCp4 <- p2
  ########################################

  #               PCA Plots              #
  ########################################

  #compute the principal components
  data_filtered_pca <- FactoMineR::PCA(t(input_omit), graph = FALSE)
  p1 <- factoextra::fviz_eig(data_filtered_pca, addlabels=TRUE)

  #extract the PCA scores for dim 1 and dim 2 so we can make our own plots.
  scores <- as.data.frame(data_filtered_pca$ind$coord) |>
    tibble::rownames_to_column(var = "Replicates") %>%
    mutate(Group = str_remove(Replicates, "_[^_]*$"))
  dim1_score <- round(p1$data$eig[1], digits = 1)
  dim2_score <- round(p1$data$eig[2], digits = 1)

  #plots for PCA
  p2 <- ggplot(scores, aes(x=Dim.1, y = Dim.2, color=Group, label = Replicates)) +
    geom_point(size = 4) +
    ggrepel::geom_text_repel(show.legend = F) +
    scale_color_manual(values = ColorPalette$Hex) +
    ylab(paste0("PCA Dim.2 (", dim2_score, "%)")) +
    xlab(paste0("PCA Dim.1 (", dim1_score, "%)")) +
    theme_classic() +
    theme(legend.position = "bottom") +
    ggtitle("Input PCA")

  #compute the principal components
  data_filtered_pca <- FactoMineR::PCA(t(diGly_omit), graph = FALSE)
  p1 <- factoextra::fviz_eig(data_filtered_pca, addlabels=TRUE)

  #extract the PCA scores for dim 1 and dim 2 so we can make our own plots.
  scores <- as.data.frame(data_filtered_pca$ind$coord) |>
    tibble::rownames_to_column(var = "Replicates") %>%
    mutate(Group = str_remove(Replicates, "_[^_]*$"))
  dim1_score <- round(p1$data$eig[1], digits = 1)
  dim2_score <- round(p1$data$eig[2], digits = 1)

  #plots for PCA
  p3 <- ggplot(scores, aes(x=Dim.1, y = Dim.2, color=Group, label = Replicates)) +
    geom_point(size = 4) +
    ggrepel::geom_text_repel(show.legend = F) +
    scale_color_manual(values = ColorPalette$Hex) +
    ylab(paste0("PCA Dim.2 (", dim2_score, "%)")) +
    xlab(paste0("PCA Dim.1 (", dim1_score, "%)")) +
    theme_classic() +
    theme(legend.position = "bottom") +
    ggtitle("diGly PCA")

  ########################################
  QCp5 <- cowplot::plot_grid(p2, p3, ncol = 2)
  ########################################

  ########################################
  #         Data Visualization           #
  ########################################

  print("Generating visualization plots...")

  contrasts <- data.table::fread(
    Sys.glob(
      paste0(
        path,
        "/output/tables/*diGly_Contrasts.csv"
        )
      )[which.max(file.info(Sys.glob(paste0(path, "/output/tables/*diGly_Contrasts.csv")))$mtime)]) # To get the most recent file in case there's more than one

  #        Unique Sites/Proteins         #
  ########################################

  diGlyUnique <- data.table::fread(
    Sys.glob(
      paste0(
        path,
        "/output/tables/*diGly_UniqueSites.csv"
        )
      )[which.max(file.info(Sys.glob(paste0(path, "/output/tables/*diGly_UniqueSites.csv")))$mtime)])

  inputUnique <- data.table::fread(
    Sys.glob(
      paste0(
        path,
        "/output/tables/*input_UniqueSites.csv"
      )
    )[which.max(file.info(Sys.glob(paste0(path, "/output/tables/*input_UniqueSites.csv")))$mtime)])

  #            Volcano Plots             #
  ########################################

  inputFC <- data.table::fread(
    Sys.glob(
      paste0(
        path,
        "/output/tables/*input_limmaDEA_long.csv"
      )
    )[which.max(file.info(Sys.glob(paste0(path, "/output/tables/*input_limmaDEA_long.csv")))$mtime)]) %>%
    mutate(Sig = ifelse(abs(logFC) >= 1 & adj.P.Val <= 0.05,
                        "Significant",
                        "Non-Significant")) %>%
    mutate(Color = ifelse(Sig == "Significant",
                          "#AA4499",
                          "#BBBBBB")) %>%
    mutate(Biogenesis = ifelse(protein_id %in% c(ribosomeProteinsHs$protein_id,biogenesisProteinOrder$protein_id),
                               "Biogenesis/Ribosome",
                               ""))

  plots <- list()
  for(i in unique(inputFC$contrast)){
    p1 <- inputFC %>%
      filter(contrast == i) %>%
      ggplot(aes(x=logFC, y=-log10(adj.P.Val), color = Color)) +
      geom_point() +
      geom_label_repel(data = filter(inputFC,
                                     contrast == i &
                                       Color == "#AA4499" &
                                       protein_id %in% c(ribosomeProteinsHs$protein_id,
                                                         biogenesisProteinOrder$protein_id)
                                     ),
                       aes(label = protein_name),
                       size = 2,
                       label.padding = 0.1
                      ) +
      scale_color_identity() +
      theme_classic() +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = c(-1,1), linetype = "dashed") +
      ggtitle(paste0("Input ", i))

    plots[[i]] <- p1
  }

  ########################################
  Vp1 <- patchwork::wrap_plots(plots, ncol = 2, widths = 3, heights = 3)
  ########################################

  diGlyFC <- data.table::fread(
    Sys.glob(
      paste0(
        path,
        "/output/tables/*diGly_limmaDEA_long.csv"
      )
    )[which.max(file.info(Sys.glob(paste0(path, "/output/tables/*diGly_limmaDEA_long.csv")))$mtime)]) %>%
    mutate(Sig = ifelse(abs(logFC) >= 1 & adj.P.Val <= 0.05,
                        "Significant",
                        "Non-Significant")) %>%
    mutate(Color = ifelse(Sig == "Significant",
                          "#AA4499",
                          "#BBBBBB")) %>%
    mutate(Biogenesis = ifelse(protein_id %in% c(ribosomeProteinsHs$protein_id,biogenesisProteinOrder$protein_id),
                               "Biogenesis/Ribosome",
                               ""))

  plots <- list()
  for(i in unique(diGlyFC$contrast)){
    p1 <- diGlyFC %>%
      filter(contrast == i) %>%
      ggplot(aes(x=logFC, y=-log10(adj.P.Val), color = Color)) +
      geom_point() +
      geom_label_repel(data = filter(diGlyFC,
                                     contrast == i &
                                       Color == "#AA4499" &
                                       protein_id %in% c(ribosomeProteinsHs$protein_id,
                                                         biogenesisProteinOrder$protein_id)
      ),
      aes(label = protein_name),
      size = 2,
      label.padding = 0.1
      ) +
      scale_color_identity() +
      theme_classic() +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = c(-1,1), linetype = "dashed") +
      ggtitle(paste0("diGly ", i))

    plots[[i]] <- p1
  }

  ########################################
  Vp2 <- patchwork::wrap_plots(plots, ncol = 2, widths = 3, heights = 3)
  ########################################

  #          Significant Genes           #
  ########################################

  Vdf1 <- inputFC %>%
    filter(Sig == "Significant")

  Vdf2 <- diGlyFC %>%
    filter(Sig == "Significant")

  #            Scatter Plots             #
  ########################################

  df1 <- inputFC %>%
    select(-AveExpr, -t, -P.Value, -B, -Sig, -Color, -Biogenesis) %>%
    rename("input_logFC" = "logFC",
           "input_adj.P.Val" = "adj.P.Val")

  df2 <- diGlyFC %>%
    select(-AveExpr, -t, -P.Value, -B) %>%
    rename("diGly_logFC" = "logFC",
           "diGly_adj.P.Val" = "adj.P.Val")

  df3 <- merge(df2, df1,
               by = c("protein_id", "protein_name", "contrast"))

  plots <- list()
  for(i in unique(df3$contrast)){
    p1 <- df3 %>%
      filter(contrast == i) %>%
      ggplot(aes(x=input_logFC, y=diGly_logFC, color = Color)) +
      geom_point() +
      geom_label_repel(data = filter(df3,
                                     contrast == i &
                                       Color == "#AA4499" &
                                       protein_id %in% c(ribosomeProteinsHs$protein_id,
                                                         biogenesisProteinOrder$protein_id)
      ),
      aes(label = protein_name),
      size = 2,
      label.padding = 0.1
      ) +
      scale_color_identity() +
      theme_classic() +
      geom_hline(yintercept = c(-1,1), linetype = "dashed") +
      geom_vline(xintercept = c(-1,1), linetype = "dashed") +
      ggtitle(paste0(i))

    plots[[i]] <- p1
  }

  ########################################
  Vp3 <- patchwork::wrap_plots(plots, ncol = 2, widths = 3, heights = 3)
  ########################################

  #          Significant Genes           #
  ########################################

  Vdf3 <- df3 %>%
    filter(Sig == "Significant" &
             abs(input_logFC) <= 1)

  ########################################
  #     E3 Ligase Enrichment Analysis    #
  ########################################

  ube3APA <- Sys.glob(
    paste0(
      path,
      "/output/ube3APA/E3enrichment_protein_level_normalized_by_protein_grouped*.csv"
    )
  ) %>%
    purrr::set_names(basename(.)) %>%
    purrr::map(data.table::fread) %>%
    bind_rows(.id = "filename") %>%
    mutate(contrast = str_extract(filename, "(?<=grouped_)[A-Za-z0-9]+")) %>%
    mutate(Sig = ifelse(abs(site_avg) >= 0.5 & p_value <= 0.05,
                        "Significant",
                        "Non-Significant")) %>%
    mutate(Color = ifelse(Sig == "Significant",
                          "#AA4499",
                          "#BBBBBB"))

  plots <- list()
  for(i in unique(ube3APA$contrast)){
    p1 <- ube3APA %>%
      filter(contrast == i) %>%
      ggplot(aes(x=site_avg, y=-log10(p_value), color = Color)) +
      geom_point() +
      scale_color_identity() +
      theme_classic() +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed") +
      geom_text_repel(data = filter(ube3APA, contrast == i & Sig == "Significant"),
                      aes(label = leading_e3ligase),
                      size = 2) +
      ggtitle(paste0(i))

    plots[[i]] <- p1
  }

  ########################################
  E3p1 <- patchwork::wrap_plots(plots, ncol = 2, widths = 3, heights = 3)
  ########################################

  #          Significant Genes           #
  ########################################

  E3df1 <- ube3APA %>%
    filter(Sig == "Significant")

  ########################################
  #           List of Results            #
  ########################################

  results <- list(
    "Quality Control" = list(
      plots = list(
        "InputNormalization" = QCp1,
        "diGlyNormalization" = QCp2,
        "InputCor" = QCp3,
        "diGlyCor" = QCp4,
        "PCA" = QCp5
      ),
      tables = list(
        "Contrasts" = contrasts,
        "InputUnique" = inputUnique,
        "diGlyUnique" = diGlyUnique
      )
    ),
    "Visualization" = list(
      plots = list(
        "InputVolcanoPlots" = Vp1,
        "diGlyVolcanoPlots" = Vp2,
        "ScatterPlots" = Vp3
      ),
      tables = list(
        "InputSig" = Vdf1,
        "diGlySig" = Vdf2,
        "ScatterSig" = Vdf3
      )
    ),
    "E3 Ligase Analysis" = list(
      plots = list(
        "VolcanoPlots" = E3p1
      ),
      tables = list(
        "E3Sig" = E3df1
      )
    )
  )

  ########################################
  #          Render HTML Report          #
  ########################################

  if (!dir.exists(paste0(path, "/output/report"))) {
    dir.create(paste0(path, "/output/report"))
  }

  template <- system.file("rmd", "diGly_html_report_template.Rmd", package = "gsptools")
  outDir <- paste0(path, "/output/report")

  suppressMessages({
  rmarkdown::render(
    input = template,
    output_file = paste0(jobname, "_report.html"),
    output_dir = normalizePath(outDir),
    params = list(results = results,
                  report_title = jobname),
    envir = new.env(parent = globalenv())
    )
  }
  )

  print("########################################")
  print("Done! HTML report saved to: ", normalizePath(outDir))
  print("########################################")

}
