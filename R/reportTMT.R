#' Generates an html report based on outputs from limmaFragpipeTMT function
#'
#' This function takes the output from limmaFragpipeTMT function and generates a report with
#' all results and tables.
#'
#' @param path Main path of the project, usually as PROT-XXXX and at the top level.
#' @param jobname Specifies the project name, for convenience I always use the PROT-XXXX.
#' @return Generates a report with all standard analyses for TMT data.
#' @export

reportTMT <- function(path,
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
  input_raw <- data.table::fread(Sys.glob(paste0(path, "4_input/abundance_protein_MD.tsv"))) %>%
    select(Index, ReferenceIntensity:last_col(), -ReferenceIntensity) %>%
    reshape2::melt()

  input_norm <- data.table::fread(Sys.glob(paste0(path, "4_input/*normalized.txt"))) %>%
    select(protein_id, reference_intensity:last_col(), -reference_intensity) %>%
    reshape2::melt()

  p1 <- input_raw %>%
    ggplot(aes(x=variable, y=value)) +
    geom_boxplot() +
    ylab("Log2(Abundance)") +
    xlab("") +
    ggtitle("Before Normalization") +
    norm_boxplot

  p2 <- input_norm %>%
    ggplot(aes(x=variable, y=value)) +
    geom_boxplot() +
    ylab("Log2(Abundance)") +
    xlab("") +
    ggtitle("After Normalization") +
    norm_boxplot

  ########################################
  QCp1 <- cowplot::plot_grid(p1, p2, ncol = 2)
  ########################################

  #           Correlation Plots          #
  ########################################

  input_omit <- data.table::fread(Sys.glob(paste0(path, "4_input/*normalized.txt"))) %>%
    select(protein_id, reference_intensity:last_col(), -reference_intensity) %>%
    tibble::column_to_rownames("protein_id") %>%
    na.omit()

  correlation <- cor(input_omit, use = "complete.obs")

  p1 <- ComplexHeatmap::Heatmap(correlation,
                                name = "Correlation",
                                column_title = "Correlation",
                                col = circlize::colorRamp2(c(min(correlation), 1), c("white", "#882255")),
                                show_row_names = TRUE,
                                show_column_names = TRUE,
                                border = TRUE)

  ########################################
  QCp3 <- p1
  ########################################

  #               PCA Plots              #
  ########################################

  #compute the principal components
  data_filtered_pca <- FactoMineR::PCA(t(input_omit), graph = FALSE)
  p1 <- factoextra::fviz_eig(data_filtered_pca, addlabels=TRUE)

  #extract the PCA scores for dim 1 and dim 2 so we can make our own plots.
  scores <- as.data.frame(data_filtered_pca$ind$coord) |>
    tibble::rownames_to_column(var = "Replicates") %>%
    mutate(Group = stringr::str_remove(Replicates, "_[^_]*$"))
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
    ggtitle("PCA")

  ########################################
  QCp5 <- p2
  ########################################

  ########################################
  #         Data Visualization           #
  ########################################

  print("Generating visualization plots...")

  contrasts <- data.table::fread(
    Sys.glob(
      paste0(
        path,
        "output/tables/*_Contrasts.csv"
        )
      )[which.max(file.info(Sys.glob(paste0(path, "output/tables/*_Contrasts.csv")))$mtime)]) # To get the most recent file in case there's more than one

  #        Unique Sites/Proteins         #
  ########################################

  inputUnique <- data.table::fread(
    Sys.glob(
      paste0(
        path,
        "output/tables/*_UniqueSites.csv"
      )
    )[which.max(file.info(Sys.glob(paste0(path, "output/tables/*_UniqueProteins.csv")))$mtime)])

  #            Volcano Plots             #
  ########################################

  inputFC <- data.table::fread(
    Sys.glob(
      paste0(
        path,
        "output/tables/*_limmaDEA_long.csv"
      )
    )[which.max(file.info(Sys.glob(paste0(path, "output/tables/*_limmaDEA_long.csv")))$mtime)]) %>%
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
      ggrepel::geom_label_repel(data = filter(inputFC,
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
      ggtitle(paste0(i))

    plots[[i]] <- p1
  }

  ########################################
  Vp1 <- patchwork::wrap_plots(plots, ncol = 2, widths = 3, heights = 3)
  ########################################

  #          Significant Genes           #
  ########################################

  Vdf1 <- inputFC %>%
    filter(Sig == "Significant")

  ########################################
  #           List of Results            #
  ########################################

  results <- list(
    "Quality Control" = list(
      plots = list(
        "Normalization" = QCp1,
        "Cor" = QCp3,
        "PCA" = QCp5
      ),
      tables = list(
        "Contrasts" = contrasts,
        "Unique" = inputUnique
      )
    ),
    "Visualization" = list(
      plots = list(
        "VolcanoPlots" = Vp1
      ),
      tables = list(
        "Sig" = Vdf1
      )
    )
  )

  ########################################
  #        Path to Results files         #
  ########################################

  excel <- Sys.glob(
    paste0(
      path,
      "output/tables/*_results.xlsx"
    )
  )
  excel <- excel[!grepl("~\\$", excel)] # Get rid of temporary files that might be present

  zip_file <- here::here(path, "output", "report", paste0(jobname, "_results.zip"))
  zip::zipr(zip_file,
           files = excel,
           mode = "cherry-pick")

  ########################################
  #          Render HTML Report          #
  ########################################

  if (!dir.exists(paste0(path, "output/report"))) {
    dir.create(paste0(path, "output/report"))
  }

  template <- system.file("rmd", "TMT_html_report_template.Rmd", package = "gsptools")
  outDir <- paste0(path, "output/report")

  rmarkdown::render(
    input = template,
    output_file = paste0(jobname, "_report.html"),
    output_dir = normalizePath(outDir),
    params = list(results = results,
                  report_title = jobname,
                  excel_files = zip_file),
    envir = new.env(parent = globalenv())
    )

  message("########################################")
  message("Done! HTML report saved to: ", normalizePath(outDir))
  message("########################################")

}
