#' Import DIA tables from ProteoVista quickomics input files
#'
#' This function takes a path to the ProteoVista quickomics input files and cleans the tables so that
#' gene names are nicely after the protein ID, while also saving the output into the outputPath.
#'
#' @param inputPath Path where quickomics files are and should include msqrob_dea_results.csv, quickomics_expression_data.csv, quickomics_gene_protein_table.csv and quickomics_sample_metadata.csv
#' @param outputPath Path to where you want to save the output tables
#' @param organism Organism that the data is from to fix capitalization of protien names. Zebrafish = zf, and protiens should start with upper letter but everything else is lower case.
#' @return Assigns the new dataframes as "dea" and "abundance" in the global environment, while also saving the tables as .csv in the outputPath.
#' @export

importDIAquickomics <- function(inputPath, outputPath, organism = "zf"){
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
  if (organism == "zf"){
    
  dea <- merge(dea, geneNames[,c("UniqueID", "gene_symbols", "Gene.Name")],
               by = "UniqueID") %>%
    relocate(c(gene_symbols, Gene.Name), .after = UniqueID) %>%
    rename("adj.P.Val" = "Adj.P.Value",
           "ProteinID" = "UniqueID",
           "ProteinName" = "Gene.Name",
           "contrast" = "test") %>%
    mutate(ProteinName = str_to_sentence(str_to_lower(ProteinName)))
  assign("dea", dea, envir = .GlobalEnv) # assign the df as a new object
  write_csv(dea, paste0(outputPath, "msqrob_dea_with_names.csv"))
  
  # Abundance
  abundance <- merge(abundance, geneNames[,c("UniqueID", "gene_symbols", "Gene.Name")],
                     by = "UniqueID") %>%
    relocate(c(gene_symbols, Gene.Name), .after = UniqueID) %>%
    rename("ProteinID" = "UniqueID",
           "ProteinName" = "Gene.Name") %>%
    mutate(ProteinName = str_to_sentence(str_to_lower(ProteinName)))
  names(abundance) <- metadata$uniqueSample[match(names(abundance), metadata$sampleid)] # this renames columns to simpler uniqueIDs from the metadata file
  assign("abundance", abundance, envir = .GlobalEnv)
  write_csv(abundance, paste0(outputPath, "abundance_with_names.csv"))
  } else {
    
    dea <- merge(dea, geneNames[,c("UniqueID", "gene_symbols", "Gene.Name")],
                 by = "UniqueID") %>%
      relocate(c(gene_symbols, Gene.Name), .after = UniqueID) %>%
      rename("adj.P.Val" = "Adj.P.Value",
             "ProteinID" = "UniqueID",
             "ProteinName" = "Gene.Name",
             "contrast" = "test")
    assign("dea", dea, envir = .GlobalEnv) # assign the df as a new object
    write_csv(dea, paste0(outputPath, "msqrob_dea_with_names.csv"))
    
    # Abundance
    abundance <- merge(abundance, geneNames[,c("UniqueID", "gene_symbols", "Gene.Name")],
                       by = "UniqueID") %>%
      relocate(c(gene_symbols, Gene.Name), .after = UniqueID) %>%
      rename("ProteinID" = "UniqueID",
             "ProteinName" = "Gene.Name")
    names(abundance) <- metadata$uniqueSample[match(names(abundance), metadata$sampleid)] # this renames columns to simpler uniqueIDs from the metadata file
    assign("abundance", abundance, envir = .GlobalEnv)
    write_csv(abundance, paste0(outputPath, "abundance_with_names.csv"))
  }
}