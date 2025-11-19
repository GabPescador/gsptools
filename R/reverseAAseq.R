#' Reverses an aminoacid sequence in FASTA format to use as rev_ decoy.
#'
#' This function takes a single character vector of an aminoacid sequence and reverses
#' the amino acid order to use a rev_ decoy for proteomics searches. It will also take out the
#' M from the end and put it as the first character.
#'
#' @param sequence character vector of your aminoacid sequence.
#' @return Generates a character vector with your reversed sequence.
#' @export

reverseAAseq <- function(sequence){

  # Reverse the sequence
  chars <- strsplit(sequence, "")[[1]]
  reversed_chars <- rev(chars)

  # Move last character to front
  n <- length(reversed_chars)
  final_chars <- c(reversed_chars[n], reversed_chars[1:(n-1)])

  # Collapse back to string
  final_string <- paste(final_chars, collapse = "")
  print(final_string)
}
