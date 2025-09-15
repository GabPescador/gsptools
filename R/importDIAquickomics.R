#' Import DIA tables from ProteoVista quickomics input files
#'
#' This function takes a path to the ProteoVista quickomics input files and cleans the tables so that
#' gene names are nicely after the protein ID, while also saving the output into the outputPath.
#'
#' @param inputPath Path where quickomics files are and should include msqrob_dea_results.csv, quickomics_expression_data.csv, quickomics_gene_protein_table.csv and quickomics_sample_metadata.csv
#' @param outputPath Path to where you want to save the output tables
#' @return Assigns the new dataframes as "dea" and "abundance" in the global environment, while also saving the tables as .csv in the outputPath.
#' @export

importDIAquickomics <- function(inputPath, outputPath){
  # inputPath is the folder with all quickomics files
  # outputPath is the folder where you want to save your resulting tables
  
  dea <- fread(paste0(inputPath, "msqrob_dea_results.csv"))
  abundance <- fread(paste0(inputPath, "quickomics_expression_data.csv"))
  geneNames <- fread(paste0(inputPath, "quickomics_gene_protein_table.csv"))
  metadata <- fread(paste0(inputPath, "quickomics_sample_metadata.csv")) %>%
    select(sampleid, group) %>%
    group_by(group) %>% # this adds back _1, _2, _3, etc for each replicate
    mutate(uniqueSample = paste0(group, "_", row_number())) %>%
    ungroup()
  
  # Making tables nicer
  # DEA
  dea <- merge(dea, geneNames[,c("UniqueID", "gene_symbols", "Gene.Name")],
               by = "UniqueID") %>%
    relocate(c(gene_symbols, Gene.Name), .after = UniqueID)
  assign("dea", dea, envir = .GlobalEnv) # assign the df as a new object
  write_csv(dea, paste0(outputPath, "msqrob_dea_with_names.csv"))
  
  # Abundance
  abundance <- merge(abundance, geneNames[,c("UniqueID", "gene_symbols", "Gene.Name")],
                     by = "UniqueID") %>%
    relocate(c(gene_symbols, Gene.Name), .after = UniqueID)
  names(abundance) <- metadata$uniqueSample[match(names(abundance), metadata$sampleid)] # this renames columns to simpler uniqueIDs from the metadata file
  assign("abundance", abundance, envir = .GlobalEnv)
  write_csv(abundance, paste0(outputPath, "abundance_with_names.csv"))
}