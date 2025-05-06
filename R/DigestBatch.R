#' Digest multiple protein sequences
#'
#' This function takes a dataframe downloaded from UNIPROT with Entry (Protein ID), Gene Names
#' and Sequence columns, then uses Digest from OrgMassSpecR and loops for multiple sequences
#' in the dataframe. Output is a dataframe with all digested peptides.
#'
#' @param df Dataframe input with aminoacid sequence in a columm named 'Sequence', protein ids in a column names 'Entry' and gene names in a column named 'Gene Names'
#' @param enz Specifies enzyme to be used. Options are "trypsin" (default), "trypsin.strict", "pepsin"
#' @param miss Missed cleavages to account for. Default is 1.
#' @return Generates a dataframe with all predicted digested peptides form all protein sequences.
#' @export
DigestBatch <- function(df, enz = "trypsin", miss = 1){
  # Batch a dataframe of protein sequences and rbind the digestion results
  # df is a dataframe with aminoacid sequence in a column named Sequence
  # df also needs first column to have the protein ID
  # enz specifies enzyme from function Digest
  # miss specifies how many missing cleavages from function Digest
  output <- data.frame()
  for(i in 1:nrow(df)){
    temp <- OrgMassSpecR::Digest(df$Sequence[i],
                   enzyme = enz, 
                   missed = miss)
    temp$Protein <- as.character(df[i,'Entry'])
    temp$Name <- as.character(df[i,'Gene Names'])
    output <- rbind(output, temp)
  }
  return(output)
}