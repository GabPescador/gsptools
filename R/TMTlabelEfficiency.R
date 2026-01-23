#' Calculates labeling efficiency from TMT files
#'
#' This function takes an inputPath from FragPipe results where several subfolders
#' contain the psm.tsv files generated and calculates labeling efficiency by counting
#' how many PSMs contain the TMT modification divided by the total number of PSMs.
#'
#' @param inputPath Input path where "psm.tsv" files are located.
#' @return Returns a dataframe with calculations per subfolder found.
#' @export

TMTlabelEfficiency <- function(inputPath) {
  
  # List all psm.tsv files in the directory
  files <- list.files(path = inputPath,
                      include.dirs = TRUE,
                      full.names = TRUE,
                      recursive = TRUE,
                      pattern = "psm")
  
  # Read files and name by their parent folder
  data_list <- lapply(files, data.table::fread)
  names(data_list) <- basename(dirname(files))
  
  # Calculate LE for each psm file
  psm_df <- purrr::imap_dfr(data_list, ~ .x %>%
                   dplyr::summarize(label = sum(str_detect(`Assigned Modifications`, "304.2071")),
                                 total = nrow(.),
                                 n_term = sum(str_detect(`Assigned Modifications`, "N-term\\(304.2071")),
                                 lysine = sum(str_detect(`Assigned Modifications`, "K\\(304.2071"))) %>%
                   dplyr::mutate(LE = label/total,
                              LE_nterm = n_term/total,
                              LE_lysine = lysine/total,
                              sample = .y)) %>%
                   dplyr::relocate(sample, .before = label)
  
  return(psm_df)
}