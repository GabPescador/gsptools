#' Generates match type plots for SILAC output from Fragpipe
#'
#' This function takes a dataframe from the Fragpipe output combine_modified_peptide_label_quant.tsv,
#' calculates the percentage of pairs that were matched either by MS/MS or requantify and generates
#' plots to visualize the percentage of matches.
#'
#' @param df_digly Dataframe corresponding to combine_modified_peptide_label_quant.tsv for the digly data.
#' @param df_input Dataframe corresponding to combine_modified_peptide_label_quant.tsv for the input data. Default is NULL.
#' @return List of dataframes with calculated match type percentages and plots for visualization.
#' @export

correlationMatchType <- function(df_digly, df_input = NULL){
  
  # Identify heavy columns on the dataframe
  heavy_cols <- grep("heavy_match_type", names(df_digly), value = TRUE)
  
  # For each heavy col, find its light pair and compare
  pair_checks <- lapply(heavy_cols, function(h_col) {
    l_col <- gsub("_heavy_", "_light_", h_col)
    if (!l_col %in% names(df_digly)) {
      message("No light pair found for: ", h_col)
      return(NULL)
    }
    
    mismatches <- df_digly[get(h_col) != get(l_col),
                           .(modified_peptide, gene,
                           heavy = get(h_col),
                           light = get(l_col))]
    setnames(mismatches, c("heavy", "light"), c(h_col, l_col))
    list(pair = sub(".*digly_(.*)_heavy_match_type.*", "\\1", h_col),
         n_mismatches = nrow(mismatches),
         mismatches = mismatches)
  })
  
  # Print summary
  summary_list <- setNames(
    lapply(pair_checks, function(check) {
      list(
        pair         = check$pair,
        n_mismatches = check$n_mismatches,
        mismatches   = check$mismatches
      )
    }),
    sapply(pair_checks, `[[`, "pair")  # use the pair name as the list name
  )
  
  # rename columns to be only heavy or light instead of full pair name
  for (i in seq_along(summary_list)) {
    names(summary_list[[i]]$mismatches)[c(3, 4)] <- c("heavy", "light")
  }
  
  mismatches_digly_df <- rbindlist(
    lapply(summary_list, function(x) {
      if (x$n_mismatches > 0) {
        x$mismatches[, pair := x$pair]  # add pair identifier column
      }
    }),
    fill = TRUE  # handles any column differences between pairs
  ) %>%
    mutate(grouping = case_when( # defines any MS/MS and requantify pair as matched
      heavy=="MS/MS" & light=="requantify" ~ "Matched",
      heavy=="requantify" & light=="MS/MS" ~ "Matched",
      heavy=="MS/MS" & light=="MS/MS" ~ "Matched",
      heavy=="requantify" & light=="requantify" ~ "Matched",
      TRUE ~ "Unmatched"
    ))
  
  perc_digly_df <- mismatches_digly_df %>%
    group_by(pair) %>%
    dplyr::count(grouping) %>%
    pivot_wider(id_cols = pair,
                values_from = n,
                names_from = grouping) %>%
    mutate(total = sum(Matched, Unmatched)) %>%
    mutate(pct_matched = Matched/total*100,
           pct_unmatched = Unmatched/total*100)
  
  p1 <- perc_digly_df %>%
    reshape2::melt(id.vars = c('pair', 'Matched', 'Unmatched', 'total')) %>%
    ggplot(aes(x=pair, y=value, fill=variable)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = c("#AA4499", "#BBBBBB")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("digly Match Types") +
    xlab("")
  
  if (!is.null(df_input)){ # same calculations but for input dataframe
    
    # Identify heavy columns
    heavy_cols <- grep("heavy_match_type", names(df_input), value = TRUE)
    
    # For each heavy col, find its light pair and compare
    pair_checks <- lapply(heavy_cols, function(h_col) {
      l_col <- gsub("_heavy_", "_light_", h_col)
      
      if (!l_col %in% names(df_input)) {
        message("No light pair found for: ", h_col)
        return(NULL)
      }
      
      mismatches <- df_input[get(h_col) != get(l_col), .(modified_peptide, gene,
                                                         heavy = get(h_col),
                                                         light = get(l_col))]
      setnames(mismatches, c("heavy", "light"), c(h_col, l_col))
      
      list(pair = sub(".*input_(.*)_heavy_match_type.*", "\\1", h_col),
           n_mismatches = nrow(mismatches),
           mismatches = mismatches)
    })
    
    # Print summary
    summary_list <- setNames(
      lapply(pair_checks, function(check) {
        list(
          pair         = check$pair,
          n_mismatches = check$n_mismatches,
          mismatches   = check$mismatches
        )
      }),
      sapply(pair_checks, `[[`, "pair")  # use the pair name as the list name
    )
    
    # rename columns to be only heavy or light
    for (i in seq_along(summary_list)) {
      names(summary_list[[i]]$mismatches)[c(3, 4)] <- c("heavy", "light")
    }
    
    mismatches_df <- rbindlist(
      lapply(summary_list, function(x) {
        if (x$n_mismatches > 0) {
          x$mismatches[, pair := x$pair]  # add pair identifier column
        }
      }),
      fill = TRUE  # handles any column differences between pairs
    ) %>%
      mutate(grouping = case_when(
        heavy=="MS/MS" & light=="requantify" ~ "Matched",
        heavy=="requantify" & light=="MS/MS" ~ "Matched",
        heavy=="MS/MS" & light=="MS/MS" ~ "Matched",
        heavy=="requantify" & light=="requantify" ~ "Matched",
        TRUE ~ "Unmatched"
      ))
    
    perc_df <- mismatches_df %>%
      group_by(pair) %>%
      dplyr::count(grouping) %>%
      pivot_wider(id_cols = pair,
                  values_from = n,
                  names_from = grouping) %>%
      mutate(total = sum(Matched, Unmatched)) %>%
      mutate(pct_matched = Matched/total*100,
             pct_unmatched = Unmatched/total*100)
    
    p2 <- perc_df %>%
      reshape2::melt(id.vars = c('pair', 'Matched', 'Unmatched', 'total')) %>%
      ggplot(aes(x=pair, y=value, fill=variable)) +
      geom_bar(stat="identity") +
      scale_fill_manual(values = c("#AA4499", "#BBBBBB")) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle("Input Match Types") +
      xlab("")
  }
  
  if (!is.null(df_input)){
    return(list("diGly Mismatches" = mismatches_digly_df,
                "diGly Mismatches Percentage" = perc_digly_df,
                "diGly Plot" = p1,
                "Input Mismatches" = mismatches_df,
                "Input Mismatches Percentage" = perc_df,
                "Input Plot" = p2))
  } else {
    return(list("diGly Mismatches" = mismatches_digly_df,
                "diGly Mismatches Percentage" = perc_digly_df,
                "diGly Plot" = p1))
  }
  
}