#' Generates QC reports from DIA data from DIA-NN
#'
#' This function takes dia-nn output and generate QC reports that were missing from
#' fragpipe using dia-nn v.2.2.0. All reports are based on each raw file run.
#'
#' @param inputPath Input path where "metadata.csv", "report.tsv", "report.stats.tsv" and the fasta file are.
#' @param jobname Specifies the folder name to save .
#' @param outputpath Output path to save all the outputs.
#' @return Generates the QC reports for each run and saves them in the outputPath.
#' @export

diannQCreport <- function(inputPath, jobname, outputPath){

  metadata <- read_excel(here(inputPath, "sample_metadata.xlsx"))
  
  # General peptide and protein identifications
  DIA_stats <- fread(here(inputPath, "report.stats.tsv")) %>%
    mutate(sample_id = basename(File.Name)) %>%
    mutate(sample_id = str_split_fixed(sample_id, "\\.", n=Inf)[,1]) %>% # to take out .d from the end
    relocate(sample_id, .before = File.Name) %>%
    merge(., metadata, by="sample_id")
  
  # Full report for QC plots
  df1 <- fread(here(inputPath, "report.tsv")) %>%
    merge(., by.x = "Run",
          metadata, by.y = "sample_id")

  # Creating a PDF with all QC plots
  cairo_pdf(here(outputPath, "plots", paste0(jobname, "_QCstats.pdf")),
      width = 8.5,
      height = 11)
  
  p1 <- DIA_stats %>%
    ggplot(aes(x=shortname, y=Precursors.Identified)) +
    geom_bar(stat="identity", fill="black") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) +
    xlab("") +
    ggtitle("Identified Precursors")
  
  p2 <- DIA_stats %>%
    ggplot(aes(x=shortname, y=Proteins.Identified)) +
    geom_bar(stat="identity", fill="black") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) +
    xlab("") +
    ggtitle("Identified Proteins")
  
  gridExtra::grid.arrange(p1, p2, ncol = 2, widths = c(4, 4), heights = c(4, 4))
  
  dev.off()
  
  for(i in unique(df1$shortname)){
    
    cairo_pdf(here(outputPath, "plots", paste0(jobname, "_", i, "_QC.pdf")),
        width = 8.5,
        height = 11)
    
    p3 <- df1 %>%
      filter(shortname == i) %>%
      ggplot(aes(x=RT, y=Ms1.Area)) +
      geom_line() +
      theme_minimal() +
      xlab("Retention Time (RT)") +
      ylab("TIC MS1") +
      ggtitle(paste0(i, " Chromatogram"))
    
    print(paste0("Plotting ", i, " TIC..."))
    
    p4 <- df1 %>%
      filter(shortname == i) %>%
      ggplot(aes(x=RT)) +
      geom_histogram(fill = "black",
                     binwidth = 1) +
      theme_minimal() +
      ggtitle(paste0(i, " RT Distribution"))
    
    print(paste0("Plotting ", i, " RT Distribution..."))
    
    p5 <- df1 %>%
      filter(shortname == i) %>%
      ggplot(aes(x=Precursor.Mz)) +
      geom_histogram(fill = "black") +
      theme_minimal() +
      ggtitle(paste0(i, " m/z Distribution"))
    
    print(paste0("Plotting ", i, " m/z Distribution..."))
    
    p6 <- df1 %>%
      filter(shortname == i) %>%
      ggplot(aes(x=Precursor.Charge)) +
      geom_histogram(fill = "black",
                     binwidth = 1) +
      theme_minimal() +
      ggtitle(paste0(i, " Precursor Charge"))
    
    print(paste0("Plotting ", i, " Precursor Charge..."))
    
    p7 <- df1 %>%
      filter(shortname == i) %>%
      group_by(Protein.Group) %>%
      summarize(precursors = n()) %>%
      ggplot(aes(x=precursors)) +
      geom_histogram(fill = "black") +
      theme_minimal() +
      ggtitle(paste0(i, " Precursor per protein group"))
    
    print(paste0("Plotting ", i, " Precursor per protein group..."))
    
    p8 <- df1 %>%
      filter(shortname == i) %>%
      ggplot(aes(x=Predicted.RT, y=RT)) +
      ggrastr::geom_point_rast() +
      theme_minimal() +
      ggtitle(paste0(i, " RT vs Predicted RT"))
    
    print(paste0("Plotting ", i, " RT vs Predicted RT..."))
    
    p9 <- df1 %>%
      filter(shortname == i) %>%
      ggplot(aes(x=Predicted.IM, y=IM)) +
      ggrastr::geom_point_rast() +
      theme_minimal() +
      ggtitle(paste0(i, " IM vs Predicted IM"))
    
    print(paste0("Plotting ", i, " IM vs Predicted IM..."))
    
    p10 <- df1 %>%
      filter(shortname == i) %>%
      ggplot(aes(x=RT, y=Precursor.Mz)) +
      ggrastr::geom_point_rast(size = 0.1) +
      theme_minimal() +
      ggtitle(paste0(i, " m/z vs RT"))
    
    print(paste0("Plotting ", i, " m/z vs RT..."))
    
    gridExtra::grid.arrange(p3, p4, p5, p6, p7, p8, p9, p10, ncol = 3)
    
    dev.off()
    
    rm(p3, p4, p5, p6, p7, p8, p9, p10)
  }
  
}
  