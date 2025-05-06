#' Gets names of genes from an upsetR object
#'
#' This function takes an upsetR object and gives you the names for each intersection group
#' This was created by: https://github.com/docmanny
#'
#' @param input upsetR object
#' @return Generates a dataframe with your gene names per group
#' @export
fromList <- function (input) {
  # Make rownames of upsetR to be your list of arguments
  # created by: https://github.com/docmanny
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}