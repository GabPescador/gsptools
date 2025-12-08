#' Create site-ratios based on specific contrast tables for Ube3-APA package
#'
#' This function takes a dataframe of protein abundances and generates site-ratios
#' based on a contrast table for Ube3-APA package.
#'
#' @param df Dataframe with abundances.
#' @param contrastTable Contrast table dataframe to be used.
#' @param diGly Defaults to TRUE. Selects "position" instead of just "id" from main dataframe.
#' @return Generates a dataframe with site-ratios based on the contrast table.
#' @export

subtractContrastsWide <- function(df,
                                  contrastTable,
                                  diGly = TRUE) {
  if (diGly == TRUE){
  result <- df %>%
    select(id, position)  # Start with ID columns
  } else {
    result <- df %>%
      select(id)
  }
  for (i in 1:nrow(contrastTable)) {
    group1 <- contrastTable$group1[i]
    group2 <- contrastTable$group2[i]
    contrast_name <- paste0(group1, "-", group2)
    
    if (group1 %in% names(df) && group2 %in% names(df)) {
      result[[contrast_name]] <- df[[group1]] - df[[group2]]
    }
  }
  
  return(result)
}