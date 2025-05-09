#' Plot stackplots from zebrafish staging
#'
#' This function plots stackplots with zebrafish staging data and performs Fisher's exact test.
#' Sample tables used as input can be found inside "/inst/extdata/Staging_examples/".
#'
#' @param inputDir Input directory where "classification_cluster.csv" and "sample_info.csv" are, with these names in the file name.
#' @param outputDir Output directory where plots and intermediate tables will be saved.
#' @param p Probability cutoff for considering the prediction correct. Default is 0.6
#' @param ctr Character vector with the names of your control conditions. It should match the column "condition" from "sample_info.csv".
#' @param timepoints Numeric vector with which time points you want to plot, and should match the values in column "time" from "sample_info.csv". Defaults to "c(2, 4, 6)"
#' @return Plots and tables saved in the output folder path.
#' @export
#' 
plotStages <- function(inputDir,
                       outputDir,
                       p = 0.6,
                       ctr,
                       timepoints = c(2, 4, 6)){
  
  # Checking that output directory exists, and if not create it
  # Check if the directory exists
  if (!dir.exists(outputDir)) {
    # If it doesn't exist, create the directory
    dir.create(outputDir, recursive = TRUE)  # recursive = TRUE ensures that any missing parent directories are also created
    print(paste("Directory created:", outputDir))
  } else {
    print(paste("Directory already exists:", outputDir))
  }
  
  ################### Importing module ###################
  
  # First imports results files
  raw_results_file <- list.files(path = inputDir,
                                 pattern = "classification_cluster", full.names = TRUE)
  df1 <- read_csv(raw_results_file)
  
  # Read sample info file and fixes possible errors with lower and upper cases
  sample_info_file <- list.files(path = inputDir,
                                 pattern = "sample_info", full.names = TRUE)
  sample_info <- read_csv(sample_info_file)
  sample_info$file <- tolower(sample_info$file)
  sample_info$condition <- tolower(sample_info$condition)
  sample_info$rep <- tolower(sample_info$rep)
  colnames(sample_info) <- c("file", "condition", "time", "rep")
  
  # Finalize by merging all information into one table and adding the stages names
  df2 <- merge(df1, sample_info, by="file")
  df2 <- merge(df2, Stages, by.x="prediction", by.y="class")
  
  ################### Transformation module ###################
  
  # First filter for only predictions above the threshold set
  df2 <- filter(df2, pmax > p)
  
  # Count totals
  df3 <- df2 %>%
    group_by(.data$condition, .data$stage, .data$time, .data$rep) %>%
    count()
  
  # Make sure conditions and stages are ordered correctly
  df3$stage <- factor(df3$stage, c(rev(Stages$stage)[2:5], "Unknown_Dead"))
  df3$condition <- factor(df3$condition, c(ctr, unique(filter(df3, !.data$condition %in% c(ctr))$condition)))
  
  # Calculate totals per replicate
  totals <- df3 %>%
    group_by(.data$condition, .data$time, .data$rep) %>%
    summarise(total = sum(n))
  
  # and totals in general
  totals_rep <- df3 %>%
    group_by(.data$condition, .data$time) %>%
    summarise(total = sum(n))
  
  df4 <- merge(df3, totals, by=c("condition", "time", "rep"))
  
  # Calculate proportion
  df4$proportion <- df4$n/df4$total
  
  # Get how many replicates you have per condition and time point
  df4 <- df4 %>%
    group_by(.data$condition, .data$time) %>%
    mutate(repn = length(unique(rep)))
  
  # Now calculate mean values and std.error per condition and time
  df4_average <- df4 %>%
    group_by(.data$condition, .data$stage, .data$time) %>%
    summarise(mean = (sum(.data$proportion)/unique(.data$repn)), sderr = std.error(.data$proportion), .groups = "keep")
  
  # Now this is a trick to make sure proportions are kept in good order
  # Took me several days and I am still not sure how doing rev() twice works, but it does so we take the small victories
  df4_average <- df4_average %>%
    mutate(stage = fct_relevel(stage, 
                               "Unknown_Dead", "128-256_cells", "High_Sphere", "30_50_Epiboly", "GermRing_Shield")) %>%
    group_by(.data$condition, .data$time) %>%
    mutate(sdpos = rev(cumsum(rev(mean))))
  
  # In the cases where there is only 1 replicates this sets the NA to 0 to prevent errors when plotting
  df4_average$sderr[is.na(df4_average$sderr)] <- 0
  
  ################### Fisher's test module ###################
  
  # Make wide format for the test
  l1 <- list()
  for(i in timepoints){
  df_fisher <- df4 %>%
    filter(.data$time == i) %>%
    select(.data$condition, .data$stage, .data$n) %>%
    group_by(.data$condition, .data$stage) %>%
    summarize(n = sum(n)) %>%
    reshape2::dcast(condition ~ stage)
  
  # Any missing value should be 0
  df_fisher[is.na(df_fisher)] <- 0
  
  # Preparing the table for using test
  df_fisher <- as.data.frame(df_fisher)
  rownames(df_fisher) <- df_fisher$condition
  df_fisher <- df_fisher[,-1]
  
  l1[[paste0("s",i)]] <- df_fisher
  }
  
  results_df_final <- data.frame()
  for(i in names(l1)){
    results <- list()
    conditions <- rownames(l1[[i]])[!rownames(l1[[i]]) %in% ctr[1]]
    df <- as.data.frame(l1[[i]])
    for(j in conditions){
  test_result <- fisher.test(df[c(ctr[1], j),],
              simulate.p.value = TRUE)
  

  results[[paste(ctr[1], "vs", j)]] <- list(
    condition = j,
    time = parse_number(i),
    p.value = test_result$p.value,
    alternative = test_result$alternative)
    
    # Convert results to data.frame
    results_df <- do.call(rbind, lapply(names(results), function(x) {
      data.frame(Comparison = x,
                 Condition = results[[x]]$condition,
                 Time = results[[x]]$time,
                 PValue = results[[x]]$p.value,
                 Alternative = results[[x]]$alternative)
    }))
    
    results_df_final <- rbind(results_df_final, results_df)
    }
  }
  
  results_df_final <- unique(results_df_final)
  
  write_csv(results_df_final, paste0(normalizePath(outputDir), "/Fishers_test.csv"))
  
  ################### Plotting module ###################
  
  # This sets the total number of conditions to annotate inside the plot
  n <- length(unique(df3$condition))
  
  for(i in timepoints){
    df4_average %>%
      filter(.data$time == i) %>%
      mutate(stage = fct_relevel(.data$stage, 
                                 "GermRing_Shield", "30_50_Epiboly", "High_Sphere", "128-256_cells", "Unknown_Dead")) %>%
      ggplot(aes(fill=.data$stage, x=.data$condition, y=.data$mean)) +
      geom_bar(position="fill", stat="identity", colour="grey1", width=0.65) +
      geom_errorbar(aes(ymin = .data$sdpos-.data$sderr, ymax = .data$sdpos), position = "identity", width = 0.3, color="#232023") +
      scale_fill_manual(values = c("white", "grey66", "grey33", "black")) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1),
            axis.text.y = element_text(size=12),
            axis.title = element_text(size=12),
            plot.title = element_text(hjust = 0.5, size=14),
            legend.text = element_text(size=12),
            legend.title = element_blank()) +
      annotate("text",
               y=1.05,
               x=c(1:n),
               label = c(paste0("n = ", filter(totals_rep, .data$time == i)$total[1:n]))
      ) +
      ggtitle(paste0(i, "hpf")) +
      ylab("Decimal Percentage") +
      xlab("")
    
    ggsave(paste0(outputDir, i, "hpf_stage.pdf"), width = 6.5, height = 3)
    
  }
  
}
  
  
  