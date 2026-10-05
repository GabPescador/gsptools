# Load libraries ----
library(shiny)
library(bslib)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(data.table)
library(DT)
library(gsptools)
library(patchwork)
library(plotly)

# Configuration ----
# NOTE: these two lines are rewritten per job by generateShinyAppScript.R -
# do not remove or reformat them, they are matched by a regex on "^jobname <-"
# and "^species <-".
jobname <- "PROT-1293_test"
species <- "Hs" # set per job by generateShinyAppScript

# Load data ----
# generateShinyAppScript.R stages these under <finalPath>/data before the app
# is ever launched, so these paths are always relative to the app root.
stats <- fread(paste0("data/", jobname, "_BasicStats.csv"))
norm <- fread(paste0("data/", jobname, "_NormalizedAbundances.csv"))
df <- fread(paste0("data/", jobname, "_DEA_long.csv"))
uniqueProt <- fread(paste0("data/", jobname, "_UniqueProteins.csv"))


# Get unique contrasts
contrast_list <- unique(df$contrast)

# Detect chromatogram PNGs per group, sorted numerically by page number
# within each group (so page10 doesn't sort before page2). Filenames follow
# the postprocessing script's <jobname>_chromatograms_<group>_page<N>.png
# convention, where group is QC (hela) / Blank / Samples - varies by dataset
# (number of raw/mzML fractions per group), so this must not be hardcoded.
chromatogram_groups <- c("QC", "Blank", "Samples")

chromatogram_files_by_group <- lapply(chromatogram_groups, function(grp) {
  files <- list.files("www", pattern = paste0("^", jobname, "_chromatograms_", grp, "_page[0-9]+\\.png$"))
  files[order(as.integer(str_extract(files, "(?<=page)[0-9]+(?=\\.png$)")))]
})
names(chromatogram_files_by_group) <- chromatogram_groups

# List of protein ----

# Define your gene groups
if(species == "Dr"){
gene_groups <- list(
  group1 = unique(tolower(ribosomeProteinsZf$protein_name)),
  group2 = unique(tolower(biogenesisProteinOrder$Dr_protein_name)),
  group3 = unique(tolower(filter(translationProteinsZf, !str_detect(protein_name, "^mrp|^rpl|^rps"), !protein_name == "")$protein_name)),
  group4 = unique(tolower(E3LigaseKK$protein_name))
)
} else {
  gene_groups <- list(
    group1 = unique(tolower(ribosomeProteinsHs$protein_name)),
    group2 = unique(tolower(biogenesisProteinOrder$protein_name)),
    group3 = unique(tolower(filter(translationProteinsHs, !str_detect(protein_name, "^MRP|^RPL|^RPS"), !protein_name == "")$protein_name)),
    group4 = unique(tolower(E3LigaseKK$protein_name))
  )
}

# Helper functions ----

basicStats <- function(df){
  
  # Define theme to use in boxplots
  theming <- list(
    theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust=0.5))
  )
  
  
  p1 <- df %>%
    ggplot(aes(x = spectrum_file, y = psms)) +
    geom_bar(stat = "identity", fill = "black") +
    xlab("") +
    ylab("PSMs") +
    ggtitle(paste0(jobname, " PSMs per fraction")) +
    theming
  
  p2 <- df %>%
    ggplot(aes(x = spectrum_file, y = protein_groups)) +
    geom_bar(stat = "identity", fill = "black") +
    xlab("") +
    ylab("Protein Groups") +
    ggtitle(paste0(jobname, " Protein Groups per fraction")) +
    theming
  
  p <- cowplot::plot_grid(plotlist = list(p1, p2), rows = 1, cols = 2)
  
  return(p)
  
}

boxplotNormalization <- function(df){
  
  # Define theme to use in boxplots
  norm_boxplot <- list(
    theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  
  p1 <- df %>%
    reshape2::melt() %>%
    ggplot(aes(x=variable, y=value)) +
    geom_boxplot() +
    ylab("Log2(Abundance)") +
    xlab("") +
    ggtitle("After Normalization") +
    norm_boxplot
  
  return(p1)
}

correlationPlot <- function(df){
  
  input_omit <- df %>%
    tibble::column_to_rownames("protein_id") %>%
    select(-protein_name) %>%
    na.omit()
  
  correlation <- cor(input_omit, use = "complete.obs")
  
  p1 <- ComplexHeatmap::Heatmap(correlation,
                                name = "Correlation",
                                column_title = "Correlation",
                                col = circlize::colorRamp2(c(min(correlation), 1), c("white", "#882255")),
                                show_row_names = TRUE,
                                show_column_names = TRUE,
                                border = TRUE)
  
  return(p1)
}

pcaPlot <- function(df){
  
  input_omit <- df %>%
    tibble::column_to_rownames("protein_id") %>%
    select(-protein_name) %>%
    na.omit()
  
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
  
  return(p2)
}

# Boxplot of normalized abundance for one or more proteins, faceted by
# protein, colored/grouped by sample group (derived from column names, same
# convention as pcaPlot()'s Group derivation).
proteinBoxplot <- function(df, genes, ncol = 3) {
 
  # Keep only the selected proteins (case-insensitive match)
  plot_data <- df %>%
    filter(tolower(protein_name) %in% tolower(genes))
 
  if (nrow(plot_data) == 0) return(NULL)
 
  # Reshape wide (one column per sample) -> long (one row per sample/protein)
  plot_data <- plot_data %>%
    select(-protein_id) %>%
    pivot_longer(cols = -protein_name, names_to = "sample", values_to = "abundance") %>%
    # Derive group by stripping the trailing replicate number, with or
    # without an underscore before it - e.g. "C1"/"C2" -> "C", and
    # "KO_1"/"KO_2" -> "KO". (pcaPlot()'s "_[^_]*$" regex only strips an
    # underscore-delimited suffix, so it leaves "C1"/"C2" as two separate
    # groups instead of pooling them - this version handles both schemes.)
    mutate(Group = stringr::str_remove(sample, "_?[0-9]+$"))
 
  # Force alphabetical ordering of groups on the x-axis (by sample name,
  # not by first-appearance order in the data) - ggplot's default factor
  # conversion for a character x is already alphabetical, but this makes
  # it explicit and immune to any future upstream reordering of `df`.
  plot_data <- plot_data %>%
    mutate(Group = factor(Group, levels = sort(unique(Group))))
 
  p <- ggplot(plot_data, aes(x = Group, y = abundance, color = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
    scale_color_manual(values = ColorPalette$Hex) +
    facet_wrap(~ protein_name, ncol = ncol, scales = "free_y") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      strip.text = element_text(face = "bold")
    ) +
    ylab("Log2(Normalized Abundance)") +
    xlab("")
 
  return(p)
}

volcanoPlot <- function(df, genes_to_highlight, genes_to_label, contrast_name, n_top = 10) {
  
  # Only plot msqrob DEA + MSDAP's differential_detection results - other
  # dea_algorithm rows (e.g. deqms, msempire) present in the long table are
  # excluded here regardless of what else is in the data.
  plot_data <- df %>%
    filter(contrast == contrast_name,
           dea_algorithm %in% c("msqrob", "differential_detection"),
           !is.na(logFC), 
           !is.na(adj.P.Val))
  
  # ALWAYS calculate top proteins
  top_down <- plot_data %>%
    filter(logFC < 0) %>%
    arrange(adj.P.Val) %>%
    slice_min(logFC, n = n_top) %>%
    pull(protein_name)
  
  top_up <- plot_data %>%
    filter(logFC > 0) %>%
    arrange(adj.P.Val) %>%
    slice_max(logFC, n = n_top) %>%
    pull(protein_name)
  
  # Color ALL genes from the list (even non-significant)
  if(length(genes_to_highlight) == 0) {
    plot_data <- plot_data %>%
      mutate(
        highlight = case_when(
          protein_name %in% top_down ~ "Top",
          protein_name %in% top_up ~ "Top",
          TRUE ~ "Other"
        )
      )
  } else {
    # User has selected proteins - color all of them
    # Case-insensitive match (see validated_genes() for why: gene_groups
    # entries are lowercase, protein_name keeps its original case).
    plot_data <- plot_data %>%
      mutate(
        highlight = case_when(
          tolower(protein_name) %in% tolower(genes_to_highlight) ~ "User Selected",
          TRUE ~ "Other"
        )
      )
  }
  
  # Separate data by category
  other_data <- plot_data %>% filter(highlight == "Other")
  top_data <- plot_data %>% filter(highlight == "Top")
  user_data <- plot_data %>% filter(highlight == "User Selected")
  
  # Base plot with grey points first
  p <- ggplot(plot_data, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(data = other_data, color = "#BBBBBB", alpha = 0.6, size = 1.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    theme_classic() +
    theme(legend.position = "right",
          plot.title = element_text(size = 14, hjust = 0.5),
          axis.text = element_text(size = 12, hjust = 0.5)) +
    ggtitle(contrast_name) +
    ylab("-Log10(Adj.P.Val)") +
    xlab("Log2(FC)")
  
  # Add top proteins if present
  if(nrow(top_data) > 0) {
    p <- p + geom_point(data = top_data, color = "#117733", alpha = 0.8, size = 2.5)
  }
  
  # Add user-selected proteins if present (on top!) - ALL of them colored
  if(nrow(user_data) > 0) {
    p <- p + geom_point(data = user_data, color = "#AA4499", alpha = 0.9, size = 3)
  }
  
  # ONLY LABEL the significant ones from genes_to_label
  labeled_data <- plot_data %>%
    filter(protein_name %in% genes_to_label | highlight == "Top")
  
  if(nrow(labeled_data) > 0) {
    
    if(species == "Dr"){
    p <- p +
      geom_label_repel(
        data = labeled_data,
        aes(label = tools::toTitleCase(protein_name)),
        color = ifelse(labeled_data$highlight == "Top", "#117733", "#AA4499"),
        min.segment.length = 0,
        max.overlaps = 20,
        size = 3,
        show.legend = FALSE
      )
    } else {
      p <- p +
        geom_label_repel(
          data = labeled_data,
          aes(label = protein_name),
          color = ifelse(labeled_data$highlight == "Top", "#117733", "#AA4499"),
          min.segment.length = 0,
          max.overlaps = 20,
          size = 3,
          show.legend = FALSE
        )
    }
  }
  
  return(p)
}

# Interactive counterpart to volcanoPlot() - no top-N labeling logic here,
# since this version is only driven by the gene-group radio buttons (no text
# box, no top-N slider). Hover tooltip shows protein name, logFC, adj.P.Val.
volcanoPlotly <- function(df, genes_to_highlight, contrast_name) {
  
  plot_data <- df %>%
    filter(contrast == contrast_name,
           dea_algorithm %in% c("msqrob", "differential_detection"),
           !is.na(logFC),
           !is.na(adj.P.Val)) %>%
    mutate(
      highlight = if (length(genes_to_highlight) == 0) {
        "Other"
      } else {
        case_when(
          tolower(protein_name) %in% tolower(genes_to_highlight) ~ "Selected",
          TRUE ~ "Other"
        )
      },
      tooltip_text = paste0(
        "Protein: ", protein_name,
        "<br>Log2FC: ", round(logFC, 3),
        "<br>Adj. P-Val: ", formatC(adj.P.Val, format = "e", digits = 2)
      )
    )
  
  color_map <- c("Other" = "#BBBBBB", "Selected" = "#AA4499")
  
  plotly::plot_ly(
    data = plot_data,
    x = ~logFC,
    y = ~-log10(adj.P.Val),
    text = ~tooltip_text,
    hoverinfo = "text",
    type = "scatter",
    mode = "markers",
    color = ~highlight,
    colors = color_map,
    marker = list(size = 8, opacity = 0.7)
  ) %>%
    plotly::layout(
      title = list(text = contrast_name, x = 0.5),
      xaxis = list(title = "Log2(FC)", zeroline = FALSE),
      yaxis = list(title = "-Log10(Adj.P.Val)", zeroline = FALSE),
      showlegend = FALSE,
      shapes = list(
        list(type = "line", x0 = -1, x1 = -1, y0 = 0, y1 = 1, yref = "paper",
             line = list(dash = "dash", color = "grey")),
        list(type = "line", x0 = 1, x1 = 1, y0 = 0, y1 = 1, yref = "paper",
             line = list(dash = "dash", color = "grey")),
        list(type = "line", x0 = 0, x1 = 1, xref = "paper",
             y0 = -log10(0.05), y1 = -log10(0.05),
             line = list(dash = "dash", color = "grey"))
      )
    )
}

# UI ----
ui <- page_navbar(
  title = jobname,
  bg = "#44AA99",
  inverse = TRUE,
  
  # Add custom CSS
  tags$style(HTML("
    .nav-pills .nav-link {
      background-color: white;
      color: #333333;
      margin-right: 5px;
    }
    .nav-pills .nav-link.active {
      background-color: #44AA99;
      color: white;
    }
    .nav-pills .nav-link:hover {
      background-color: #66CCBB;
      color: white;
    }
  ")),
  
  ##########################
  # First page with QC plots
  ##########################
  nav_panel(
    title = "QC",
    navset_pill(
        nav_panel(title = "Basic MS Stats",
                  p(),
                  p("For TMT experiments we expect between 8k-12k PSMs and between 1.5k-2.5k protein groups per fraction."),
                  p(),
                  plotOutput("basic_stats", height = "600px", width = "1200px"),
                  div(
                    style = "text-align: left; margin-top: 15px;",
                    downloadButton("stats_png", "PNG"),
                    downloadButton("stats_pdf", "PDF")
                  )
                ),
        nav_panel(title = "Chromatograms",
                  navset_pill(
                    !!!lapply(chromatogram_groups, function(grp) {
                      files <- chromatogram_files_by_group[[grp]]

                      # Fallback message instead of a blank page if a group
                      # has no chromatograms (e.g. a run with no blanks).
                      content <- if (length(files) == 0) {
                        list(p(paste0("No chromatograms found for ", grp, ".")))
                      } else {
                        lapply(files, function(f) img(src = f, width = "50%"))
                      }

                      nav_panel(title = grp, !!!content)
                    })
                  )
                ),
        nav_panel(title = "Normalized Abundance",
                  p(),
                  p("Ideally, all boxplots have a similar distribution between replicates after normalization."),
                  p(),
                  plotOutput("norm_boxplot", height = "600px", width = "600px"),
                  div(
                    style = "text-align: left; margin-top: 15px;",
                    downloadButton("normalization_png", "PNG"),
                    downloadButton("normalization_pdf", "PDF")
                    )
                  ),
        nav_panel(title = "Correlation Plot",
                  p(),
                  p("Ideally, all replicates from the same group should cluster together and have > 0.9 correlation coefficient. For proteomics data it is not unusual to consider > 0.8 or even > 0.7 as good correlation between replicates."),
                  p(),
                  plotOutput("cor_plot", height = "600px", width = "600px"),
                  div(
                    style = "text-align: left; margin-top: 15px;",
                    downloadButton("correlation_png", "PNG"),
                    downloadButton("correlation_pdf", "PDF")
                    )
                  ),
        nav_panel(title = "PCA",
                  p(),
                  p("Ideally, groups should cluster together. If there is technical variance you might see replicates clustering together instead, which represents a problem for downstream analyses.
                    If that’s the case, one might decide to get rid of a weird replicate. Alternatively, one might want to correct the bacth effect by applying something like proBatch."),
                  p(),
                  plotOutput("pca_plot", height = "600px", width = "600px"),
                  div(
                    style = "text-align: left; margin-top: 15px;",
                    downloadButton("PCA_png", "PNG"),
                    downloadButton("PCA_pdf", "PDF")
                    )
                  )
    ),
  ),
  
  ##########################
  # Second page with Visualization plots
  ##########################
  nav_panel(
    title = "Visualization",
    
    # Sidebar settings
    layout_sidebar(
      sidebar = sidebar(
        width = 350,
        
        # Text box to highlight speific proteins
        textAreaInput(
          "protein",
          value = "",
          label = "Enter proteins (one per line)",
          placeholder = "RPS19\nRPS9\nRPL10",
          rows = 10
        ),

        # Radio buttons for predefined gene groups
        radioButtons(
          "gene_group",
          "Highlight gene group:",
          choices = c(
            "None" = "none",
            "Ribosomal Proteins" = "group1",
            "Biogenesis Factors" = "group2",
            "Translation Factors" = "group3",
            "E3 Ligases" = "group4"
          ),
          selected = "none"
        ),
        
        # Slider to select how many Top proteins to show in volcano plot
        sliderInput(
          "n_proteins",
          "Top proteins to show:",
          min = 5,
          max = 50,
          value = 10,
          step = 5
        ),
        
        uiOutput("warning_message")
        
      ),
      
      # Main part of page settings
      navset_pill(
        
        # Recursively checks all contrast names and creates a separate panel for each
        !!!lapply(contrast_list, function(contrast_name) {
          safe_name <- gsub("[^[:alnum:]]", "_", contrast_name)
          nav_panel(
            title = contrast_name,
            # Volcano plots
            plotOutput(paste0("volcano_", safe_name), height = "600px", width = "600px"),
            # Plot download buttons
            div(
              style = "text-align: left; margin-top: 15px;",
              downloadButton(paste0("download_png_", safe_name), "PNG"),
              downloadButton(paste0("download_pdf_", safe_name), "PDF")
            ),
            
            p("This table updates based on the highlighted proteins."),
            p("To download the full table with all data points go to the 'Download results' tab above."),
            
            # Table of Top proteins
            h4("Highlighted Proteins", style = "margin-top: 10px;"),
            div(
              style = "margin-bottom: 3px;",
              selectInput(paste0("table_format_", safe_name), 
                          "Format:", 
                          choices = c("Excel" = "xlsx", "CSV" = "csv"),
                          width = "100px"),
              downloadButton(paste0("download_table_", safe_name), "Download Table")
            ),
            # Table
            tableOutput(paste0("table_", safe_name))
          )
        })
      )
      )
    ),
  
  ##########################
  # Interactive volcano plot page
  # Same gene-group lists as the Visualization tab, but no text box and no
  # top-N slider - only the radio buttons drive highlighting here.
  ##########################
  nav_panel(
    title = "Interactive Volcano",
    
    layout_sidebar(
      sidebar = sidebar(
        width = 350,
        
        # Radio buttons for predefined gene groups (same groups/order as
        # the Visualization tab, but with an "_ly" suffixed input ID so it
        # doesn't share state with input$gene_group)
        radioButtons(
          "gene_group_ly",
          "Highlight gene group:",
          choices = c(
            "None" = "none",
            "Ribosomal Proteins" = "group1",
            "Biogenesis Factors" = "group2",
            "Translation Factors" = "group3",
            "E3 Ligases" = "group4"
          ),
          selected = "none"
        ),
        
        uiOutput("warning_message_ly")
      ),
      
      navset_pill(
        !!!lapply(contrast_list, function(contrast_name) {
          safe_name <- gsub("[^[:alnum:]]", "_", contrast_name)
          nav_panel(
            title = contrast_name,
            plotly::plotlyOutput(paste0("volcano_ly_", safe_name), height = "600px", width = "600px")
          )
        })
      )
    )
  ),
  
  ##########################
  # Protein Boxplots page
  # Same gene-group/manual-list priority pattern as the Visualization tab
  # (group beats manual text list), but plotted against normalized
  # abundance per sample/replicate rather than DEA fold-change stats.
  ##########################
  nav_panel(
    title = "Protein Boxplots",
    
    layout_sidebar(
      sidebar = sidebar(
        width = 350,
        
        # Manual protein list (used only when no gene group is selected)
        textAreaInput(
          "protein_box",
          value = "",
          label = "Enter proteins (one per line)",
          placeholder = "RPS19\nRPS9\nRPL10",
          rows = 10
        ),
        
        # Pre-defined gene group selector (reuses the same gene_groups list
        # already built at app startup for the Visualization tab)
        radioButtons(
          "gene_group_box",
          "Or select a gene group:",
          choices = c(
            "None (use text list)" = "none",
            "Ribosomal Proteins"   = "group1",
            "Biogenesis Factors"   = "group2",
            "Translation Factors"  = "group3",
            "E3 Ligases"           = "group4"
          ),
          selected = "none"
        ),
        
        # Layout control for the facet grid
        sliderInput(
          "ncol_box",
          "Plots per row:",
          min = 1, max = 6, value = 3, step = 1
        ),
        
        uiOutput("warning_message_box")
      ),
      
      plotOutput("protein_boxplot", height = "700px", width = "900px"),
      div(
        style = "text-align: left; margin-top: 15px;",
        downloadButton("boxplot_png", "PNG"),
        downloadButton("boxplot_pdf", "PDF")
      )
    )
  ),
  
  ##########################
  # Third page with unique proteins
  ##########################
  
  nav_panel(
    title = "Unique Proteins",
    p("All comparative plots exclude unique modified sites/proteins since they cannot perform statistics on things that are present in one sample but not in another. These are still interesting candidates and are summarized in the list below."),
    p("Be mindful that the pipeline defines sites/proteins that are not found in at least 2 replicates as not expressed in the samples. Therefore, unique sites/proteins might come from things that were present in just one replicate rather than completely not present in a certain sample."),
    div(
      style = "text-align: left; margin-top: 15px;",
      downloadButton("unique_table", "Download Table")
      ),
    DT::dataTableOutput("unique_prot")
  ),
  
  ##########################
  # Fourth page with full results download
  ##########################
  
  nav_panel(
    title = "Download results",
    p(paste0("This ShinyApp uses as a base the ", jobname, "_results.xlsx file, where all data is compiled after analysis and can be downloaded below.")),
    p("Descriptions of each sheet included in the document can be found in the first page “Description” of the excel file."),
    div(
      style = "text-align: left; margin-top: 15px;",
      downloadButton("download_zip", "Download all results")
      )
  ),
  
  nav_spacer()
)

# Server ----
server <- function(input, output, session) {
  
  ##########################
  # Reactives for QC plots
  ##########################
  
  # Render normalization boxplot
  output[["basic_stats"]] <- renderPlot({
    basicStats(stats)
  })
  
  # Download handler for PNG
  output[["stats_png"]] <- downloadHandler(
    filename = function() {
      paste0(jobname, "_stats_", Sys.Date(), ".png")
    },
    content = function(file) {
      p <- basicStats(stats)
      ggsave(file, plot = p, width = 6, height = 4, dpi = 300)
    }
  )
  
  # Download handler for PDF
  output[["stats_pdf"]] <- downloadHandler(
    filename = function() {
      paste0(jobname, "_stats_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      p <- basicStats(stats)
      ggsave(file, plot = p, width = 6, height = 4, dpi = 300)
    }
  )
  
  # Render normalization boxplot
  output[["norm_boxplot"]] <- renderPlot({
    boxplotNormalization(norm)
  })
  
  # Download handler for PNG
  output[["normalization_png"]] <- downloadHandler(
    filename = function() {
      paste0(jobname, "_normalized_", Sys.Date(), ".png")
    },
    content = function(file) {
      p <- boxplotNormalization(norm)
      ggsave(file, plot = p, width = 6, height = 4, dpi = 300)
    }
  )
  
  # Download handler for PDF
  output[["normalization_pdf"]] <- downloadHandler(
    filename = function() {
      paste0(jobname, "_normalized_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      p <- boxplotNormalization(norm)
      ggsave(file, plot = p, width = 6, height = 4, dpi = 300)
    }
  )
  
  # Render correlation plot
  output[["cor_plot"]] <- renderPlot({
    correlationPlot(norm)
  })
  
  # Download handler for PNG
  output[["correlation_png"]] <- downloadHandler(
    filename = function() {
      paste0(jobname, "_correlation_", Sys.Date(), ".png")
    },
    content = function(file) {
      png(filename = file, width = 6, height = 4, units = "in", res=300)
      print(correlationPlot(norm))
      dev.off()
    },
    contentType = "image/png" 
  )
  
  # Download handler for PDF
  output[["correlation_pdf"]] <- downloadHandler(
    filename = function() {
      paste0(jobname, "_correlation_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      pdf(file = file, width = 6, height = 4)
      print(correlationPlot(norm))  # Add print() in case it's a ggplot
      dev.off()
    },
    contentType = "application/pdf"  # Explicitly set content type
  )
  
  # Render PCA plot
  output[["pca_plot"]] <- renderPlot({
    pcaPlot(norm)
  })
  
  # Download handler for PNG
  output[["PCA_png"]] <- downloadHandler(
    filename = function() {
      paste0(jobname, "_PCA_", Sys.Date(), ".png")
    },
    content = function(file) {
      p <- pcaPlot(norm)
      ggsave(file, plot = p, width = 6, height = 4, dpi = 300)
    }
  )
  
  # Download handler for PDF
  output[["PCA_pdf"]] <- downloadHandler(
    filename = function() {
      paste0(jobname, "_PCA_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      p <- pcaPlot(norm)
      ggsave(file, plot = p, width = 6, height = 4, dpi = 300)
    }
  )
  
  ##########################
  # Reactives for volcano plots and tables
  ##########################
  
  # Create selected_genes reactive
  selected_genes <- reactive({
    # Priority 1: Check if a gene group is selected
    if (input$gene_group != "none") {
      return(gene_groups[[input$gene_group]])
    }
    
    # Priority 2: Check manual text input
    if (!is.null(input$protein) && input$protein != "") {
      gene_list <- trimws(unlist(strsplit(input$protein, "\n")))
      gene_list <- gene_list[gene_list != ""]
      return(gene_list)
    }
    
    # Priority 3: Return empty if using top N
    return(character(0))
  })
  
  # Validate against available proteins
  validated_genes <- reactive({
    genes <- selected_genes()
    
    if(length(genes) == 0) {
      return(character(0))
    }
    
    available_genes <- unique(df$protein_name)
    # Case-insensitive match: gene_groups are stored lowercase, but
    # df$protein_name keeps its original case - a plain %in% here would
    # silently fail every gene-group lookup (they'd never validate, and
    # volcanoPlot() would fall back to its "no selection" top-N branch).
    valid <- genes[tolower(genes) %in% tolower(available_genes)]
    invalid <- genes[!tolower(genes) %in% tolower(available_genes)]
    
    # Create reactive warning message
    output$warning_message <- renderUI({
      # Your code that generates 'invalid' vector
      
      if(length(invalid) > 0) {
        div(
          style = "color: black; font-weight: regular; margin-top: 10px; font-size: 10px",
          paste("⚠ Proteins not found:", paste(invalid, collapse = ", "))
        )
      } else {
        NULL  # Show nothing if no invalid proteins
      }
    })
    
    return(valid)
  })
  
  # Generate individual contrast plots
  lapply(contrast_list, function(contrast_name) {
    safe_name <- gsub("[^[:alnum:]]", "_", contrast_name)
    
    # Render plot with n_top parameter
    output[[paste0("volcano_", safe_name)]] <- renderPlot({
      genes_to_color <- validated_genes()  # ALL genes from list (colored)
      genes_to_label <- highlighted_proteins() %>%  # Only significant (labeled)
        filter(Regulation != "Non-significant") %>%
        pull(protein_name)
      
      n <- input$n_proteins
      volcanoPlot(df, genes_to_color, genes_to_label, contrast_name, n_top = n)
    })
    
    # Download handler for PNG
    output[[paste0("download_png_", safe_name)]] <- downloadHandler(
      filename = function() {
        paste0(jobname, "_", contrast_name, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        genes_to_color <- validated_genes()
        genes_to_label <- highlighted_proteins() %>% 
          filter(Regulation != "Non-significant") %>%
          pull(protein_name)
        n <- input$n_proteins
        p <- volcanoPlot(df, genes_to_color, genes_to_label, contrast_name, n_top = n)
        ggsave(file, plot = p, width = 10, height = 8, dpi = 300)
      }
    )
    
    # Download handler for PDF
    output[[paste0("download_pdf_", safe_name)]] <- downloadHandler(
      filename = function() {
        paste0(jobname, "_", contrast_name, "_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        genes_to_color <- validated_genes()
        genes_to_label <- highlighted_proteins() %>% 
          filter(Regulation != "Non-significant") %>%
          pull(protein_name)
        n <- input$n_proteins
        p <- volcanoPlot(df, genes_to_color, genes_to_label, contrast_name, n_top = n)
        ggsave(file, plot = p, width = 10, height = 8)
      }
    )
    
    # Modified reactive expression
    highlighted_proteins <- reactive({
      # Same restriction as volcanoPlot(): keep only msqrob + differential_detection
      contrast_data <- df %>%
        filter(contrast == contrast_name,
               dea_algorithm %in% c("msqrob", "differential_detection")) %>%
        arrange(adj.P.Val)
      
      # Priority 1: Check if a gene group is selected
      if (input$gene_group != "none") {
        gene_list <- gene_groups[[input$gene_group]]
        
        contrast_data %>%
          filter(tolower(protein_name) %in% tolower(gene_list)) %>%
          select(protein_name, logFC, adj.P.Val) %>%
          mutate(Regulation = case_when(
            abs(logFC) >= 1 & adj.P.Val <= 0.05 & logFC > 0 ~ "Up",
            abs(logFC) >= 1 & adj.P.Val <= 0.05 & logFC < 0 ~ "Down",
            TRUE ~ "Non-significant"
          ))
        
        # Priority 2: Check manual text input
      } else if (!is.null(input$protein) && input$protein != "") {
        gene_list <- trimws(unlist(strsplit(input$protein, "\n")))
        gene_list <- gene_list[gene_list != ""]
        
        contrast_data %>%
          filter(tolower(protein_name) %in% tolower(gene_list)) %>%
          select(protein_name, logFC, adj.P.Val) %>%
          mutate(Regulation = ifelse(logFC > 0, "Up", "Down"))
        
        # Priority 3: Default to top N from slider
      } else {
        n <- input$n_proteins
        
        top_down <- contrast_data %>%
          filter(logFC < 0,
                 abs(logFC) >= 1,
                 adj.P.Val <= 0.05) %>%
          slice_min(logFC, n = n) %>%
          select(protein_name, logFC, adj.P.Val) %>%
          mutate(Regulation = "Down")
        
        top_up <- contrast_data %>%
          filter(logFC > 0,
                 abs(logFC) >= 1,
                 adj.P.Val <= 0.05) %>%
          slice_max(logFC, n = n) %>%
          select(protein_name, logFC, adj.P.Val) %>%
          mutate(Regulation = "Up")
        
        bind_rows(top_down, top_up)
      }
    })
    
    # Rendering Table
    output[[paste0("table_", safe_name)]] <- renderTable({
      highlighted_proteins() %>%
        arrange(desc(Regulation), desc(abs(logFC))) %>%
        mutate(
          logFC = round(logFC, 3),
          adj.P.Val = formatC(adj.P.Val, format = "e", digits = 2)
        ) %>%
        rename(
          `Protein` = protein_name,
          `Log2 FC` = logFC,
          `Adj. P-Value` = adj.P.Val
        )
    }, striped = TRUE, hover = TRUE)
    
    # Download handler for table
    output[[paste0("download_table_", safe_name)]] <- downloadHandler(
      filename = function() {
        paste0(jobname, "_", contrast_name, "_highlighted_proteins_", Sys.Date(), ".csv")
      },
      content = function(file) {
        table_data <- highlighted_proteins() %>%
          arrange(desc(Regulation), desc(abs(logFC))) %>%
          mutate(
            logFC = round(logFC, 3),
            adj.P.Val = formatC(adj.P.Val, format = "e", digits = 2)
          ) %>%
          rename(
            `Protein` = protein_name,
            `Log2 FC` = logFC,
            `Adj. P-Value` = adj.P.Val
          )
        
        write.csv(table_data, file, row.names = FALSE)
      }
    )
    
  })
  
  ##########################
  # Reactives for interactive volcano plots
  ##########################
  
  # Validate the selected gene group against available proteins, same
  # case-insensitive logic as validated_genes() above. Only one input to
  # check here (input$gene_group_ly) - no manual text priority chain needed.
  validated_genes_ly <- reactive({
    if (input$gene_group_ly == "none") {
      return(character(0))
    }
    
    genes <- gene_groups[[input$gene_group_ly]]
    available_genes <- unique(df$protein_name)
    valid <- genes[tolower(genes) %in% tolower(available_genes)]
    invalid <- genes[!tolower(genes) %in% tolower(available_genes)]
    
    output$warning_message_ly <- renderUI({
      if (length(invalid) > 0) {
        div(
          style = "color: black; font-weight: regular; margin-top: 10px; font-size: 10px",
          paste("⚠ Proteins not found:", paste(invalid, collapse = ", "))
        )
      } else {
        NULL
      }
    })
    
    return(valid)
  })
  
  # Generate one interactive volcano plot per contrast
  lapply(contrast_list, function(contrast_name) {
    safe_name <- gsub("[^[:alnum:]]", "_", contrast_name)
    
    output[[paste0("volcano_ly_", safe_name)]] <- plotly::renderPlotly({
      volcanoPlotly(df, validated_genes_ly(), contrast_name)
    })
  })
  
  ##########################
  # Reactives for protein boxplots
  ##########################
  
  # Priority 1: gene group selection, Priority 2: manual text list
  selected_genes_box <- reactive({
    if (input$gene_group_box != "none") {
      return(gene_groups[[input$gene_group_box]])
    }
    if (!is.null(input$protein_box) && input$protein_box != "") {
      gene_list <- trimws(unlist(strsplit(input$protein_box, "\n")))
      return(gene_list[gene_list != ""])
    }
    return(character(0))
  })
  
  # Validate typed/selected genes against what's actually in the normalized
  # abundance table, and surface any typos/unmatched names as a warning
  validated_genes_box <- reactive({
    genes <- selected_genes_box()
    if (length(genes) == 0) return(character(0))
    
    available_genes <- unique(norm$protein_name)
    valid   <- genes[tolower(genes) %in% tolower(available_genes)]
    invalid <- genes[!tolower(genes) %in% tolower(available_genes)]
    
    output$warning_message_box <- renderUI({
      if (length(invalid) > 0) {
        div(
          style = "color: black; font-weight: regular; margin-top: 10px; font-size: 10px",
          paste("⚠ Proteins not found:", paste(invalid, collapse = ", "))
        )
      } else {
        NULL
      }
    })
    
    return(valid)
  })
  
  # Render the boxplot(s)
  output$protein_boxplot <- renderPlot({
    genes <- validated_genes_box()
    req(length(genes) > 0)  # avoid rendering an empty/invalid plot
    proteinBoxplot(norm, genes, ncol = input$ncol_box)
  })
  
  # Download handler for PNG
  output$boxplot_png <- downloadHandler(
    filename = function() paste0(jobname, "_protein_boxplots_", Sys.Date(), ".png"),
    content = function(file) {
      p <- proteinBoxplot(norm, validated_genes_box(), ncol = input$ncol_box)
      ggsave(file, plot = p, width = 10, height = 8, dpi = 300)
    }
  )
  
  # Download handler for PDF
  output$boxplot_pdf <- downloadHandler(
    filename = function() paste0(jobname, "_protein_boxplots_", Sys.Date(), ".pdf"),
    content = function(file) {
      p <- proteinBoxplot(norm, validated_genes_box(), ncol = input$ncol_box)
      ggsave(file, plot = p, width = 10, height = 8)
    }
  )

  ##########################
  # Reactives for unique proteins
  ##########################
  
  # Unique protein table
  output$unique_prot <- DT::renderDataTable({
    # Your reactive data here
    uniqueProt
  })
  
  output$unique_table <- downloadHandler(
    filename = paste0(jobname, "_unique_proteins.csv"),
    content = function(file) {
      file.copy(paste0("data/", jobname, "_UniqueProteins.csv"), file)
    },
    contentType = "text/csv"
  )
  
  ##########################
  # Reactives for downloading full results
  ##########################
  
  output$download_zip <- downloadHandler(
    filename = paste0(jobname, "_results.zip"),
    content = function(file) {
      file.copy(paste0("data/", jobname, "_results.zip"), file)
    },
    contentType = "application/zip"
  )
  
  
# ending bracket
}

# Run app ----
shinyApp(ui = ui, server = server)
