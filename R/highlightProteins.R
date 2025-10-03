#' Highlight proteins based on a list of ProteinIDs
#'
#' This function takes a dataframe and a list of ProteinIDs that you want to highlight based on
#' ribosomeProteins or translationProteins and subsets your dataframe for easier plotting. In the
#' case of dea tables, you can specify the contrast you are specifically looking for, e.g.: mat_lys-zyg_lys.
#' This function also assigns a new column "type" where things are classified in "Sig" if adj.p.val <= 0.05
#' and abs(logFC) >= 1, and "Non-Sig" for non-significant proteins.
#'
#' @param df Dataframe you want to subset from
#' @param targets List of ProteinIDs you want to subset, must come from similar format as ribosomeProteins object
#' @param contrasts Name of the contrast you want to specify from the dea dataframe
#' @return A subset dataframe from the input df for plotting
#' @export
highlightProteins <- function(df,
                              targets = ribosomeProteins,
                              contrasts){
  # subsets df with proteins that we want to highlight in the plot and group them
  # based on significance (adj.p.val <= 0.05)
  # targets = df with the target proteins to highlight
  # contrasts = contrast to be used to filter in the case of DEA df

  # First part is for DEA plots
  if ("contrast" %in% names(df)) {
    hl <- filter(df, protein_id %in% targets$protein_id &
                   contrast == contrasts)
    hl$type <- ifelse(hl$adj.P.Val <= 0.05,
                      "Sig",
                      "Non-Sig")

    return(hl)

  } else { # this part is for abundance plots
    hl <- filter(df, protein_id %in% targets$protein_id)

    return(hl)
  }
}
