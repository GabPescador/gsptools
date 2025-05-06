#' Transforms a dataframe with protein sequences into FASTA format
#'
#' This function takes a dataframe downloaded from UNIPROT with Entry (Protein ID)
#' and Sequence columns, then outputs a FASTA file.
#'
#' @param df Dataframe input with aminoacid sequence in a columm named 'Sequence' and protein ids in a column names 'Entry'
#' @param organism Specifies the organism the proteins are from. Default is 'DanRer'.
#' @param outputpath Output path to save the final FASTA file.
#' @return Generates a FASTA file based on a dataframe with aminoacid sequences.
#' @export
fastaTransform <- function(df, organism = "DanRer", outputpath){
  # Takes a dataframe downloaded from UNIPROT that has at least:
  # Entry - with uniprot ID
  # Sequence - with aminoacid sequence from UNIPROT
  # THen transforms all sequences into FASTA format
  
  df2 <- data.frame(id = c(df$Entry),
                    Organism = organism,
                    seq = c(df$Sequence))
  df2$fasta_name <- paste0(df2$id, "_", df2$Organism)
  
  zf_fasta <- unlist(
    mapply(function(fasta_name, seq) {
      c(paste0(">", fasta_name), seq)
    }, df2$fasta_name, df2$seq, SIMPLIFY = FALSE)
  )
  writeLines(zf_fasta, paste0(outputpath, organism, "_output.fasta"))
  
}