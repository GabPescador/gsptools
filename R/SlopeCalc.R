#' Calculates linear regression for proteasome activity assay in batch
#'
#' This function takes a dataframe with different conditions and replicates from an enzymatic assay
#' like proteasome activity assay where you have signal that increases over time.
#' 
#' Note: You should plot the curves before using this to define the cutoff time and use only the
#' exponential time of the curve.
#'
#' @param df Melted dataframe input with "conditions", "type", "replicates", "time" and values.
#' @param cutoff Specifies the time cutoff so that calculation is only done in the exponential curve.
#' @param conditions Character vector with all conditions in the dataframe. Can be used as df$conditions.
#' @param types Character vector with all types in the dataframe. Used to differentiate between normal conditions or conditions where samples were treated with MG132. Default is "normal".
#' @param replicates Character vector with all replicates in the dataframe. Can be used as df$replicates
#' @return Generates a dataframe with all slope values from lm().
#' @export
SlopeCalc <- function(df, cutoff, conditions, types = "normal", replicates){
  
  p <- data.frame()
  for(j in unique(replicates)){
    for(i in unique(conditions)){
      
      df2 <- df %>%
        filter(!.data$type %in% c("mean", "sd")) %>% # in case your table has mean and sd columns that you calculated by hand
        filter(.data$time <= cutoff & .data$type == types) # filters only the time and type of values you want to use
      
      temp <- df2 %>%
        filter(.data$condition == i & .data$replicate == j)
      
      temp2 <- lm(data = temp, formula = value ~ time)
      
      temp3 <- data.frame()
      temp3[1,1] <- i
      temp3[1,2] <- temp2$coefficients[2]
      temp3[1,3] <- j
      temp3[1,4] <- types
      
      p <- rbind(p, temp3)
      
    }
  }
  
  colnames(p) <- c("Condition", "Slope", "Replicate", "Type")
  return(p)
}