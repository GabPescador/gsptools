#' Generates post-processing R script
#'
#' This function generates an R script to run post-processing on msdap output.
#' It is a helper function for the processDiannMSdap.R function.
#'
#' @param outputPath Path where results will be saved.
#' @return Creates an R script to run post-processing analysis.
#' @export

generatePostprocessingScript <- function(outputPath) {
  
  msdap_path <- file.path(outputPath, "msdap_results")
  plots_path <- file.path(outputPath, "plots")
  timestamp  <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  glue::glue('
library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)

# Load MSdap output — adjust filename pattern to match your output
df <- data.table::fread(list.files("{msdap_path}", 
                                    pattern = "de_proteins.*\\.tsv", 
                                    full.names = TRUE, 
                                    recursive = TRUE)[1])

p <- df %>%
  ggplot(aes(x=foldchange.log2, y=-log10(pvalue))) +
  geom_point()

ggsave(filename = paste0({plots_path}, "test_plot.pdf"), plot = p)
  ')
}
