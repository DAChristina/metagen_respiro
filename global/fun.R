combine_assay_sheets <- function(file, sheet) {
  raw <- readxl::read_excel(file,
                            sheet = sheet,
                            col_names = FALSE)
  
  row1 <- as.character(raw[1, ])
  row2 <- as.character(raw[2, ])
  data <- raw[-c(1, 2), ]
  
  row1_filled <- zoo::na.locf(row1,
                              na.rm = FALSE)
  
  combined_names <- dplyr::case_when(
    is.na(row2) | row2 == "" ~ row1_filled,
    row1_filled == row2 ~ row1_filled,
    TRUE ~ paste0(sheet, row1_filled, "_", row2)
  )
  
  names(data) <- combined_names
  data <- data %>%
    janitor::clean_names() %>% 
    dplyr::rename_with(
      ~ gsub("allplex_tm_", "", .x),
      .cols = everything()
    )
  
  return(data)
}

