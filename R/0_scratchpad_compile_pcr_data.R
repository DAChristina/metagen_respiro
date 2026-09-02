library(tidyverse)

source_all <- function(path = "global"){
  files <- list.files(path,
                      pattern = "\\.R$",
                      full.names = TRUE
  )
  invisible(lapply(files, source))
}
source_all()

# detect sheet first
sheet_all <- readxl::excel_sheets(
  "raw_data/Training REspiro_080526.xlsx") %>% 
  glimpse()

test1 <- readxl::read_xlsx("raw_data/Training REspiro_080526.xlsx",
                          sheet = sheet_all[1]) %>% 
  glimpse()

test2 <- readxl::read_xlsx("raw_data/Training REspiro_080526.xlsx",
                           sheet = sheet_all[2],
                           skip = 1) %>% 
  glimpse()


################################################################################
# generate long data format
# annoying R dimension requirement
result <- NULL
for (sheet in sheet_all[-1]) {
  df <- combine_assay_sheets(file = "raw_data/Training REspiro_080526.xlsx",
                             sheet = sheet)
  
  if (is.null(result)) {
    result <- df
  } else {
    result <- dplyr::full_join(result, df,
                               by = "name")
  }
}

final_result <- result %>%
  dplyr::rename(
    # pneumoBacter
    sample_no = sample_no.x,
    patient_id = patient_id.x,
    pneumo_bacter_assay_8_well = well.x,
    type = type.x, # type is consistent across 4 assays
    pneumo_bacter_assay_8_auto_interpretation = auto_interpretation.x,
    pneumo_bacter_assay_8_comment = comment.x,
    
    # Respiro 1A
    respiratory_panel_1a_well = well.y,
    respiratory_panel_1a_auto_interpretation = auto_interpretation.y,
    respiratory_panel_1a_comment = comment.y,
    
    # Respiro 2
    respiratory_panel_2_well = well.x.x,
    respiratory_panel_2_auto_interpretation = auto_interpretation.x.x,
    respiratory_panel_2_comment = comment.x.x,
    
    # Respiro 3
    respiratory_panel_3_well = well.y.y,
    respiratory_panel_3_auto_interpretation = auto_interpretation.y.y,
    respiratory_panel_3_comment = comment.y.y,
    
  ) %>% 
  dplyr::select(-contains(c(".x", ".y"))) %>% 
  glimpse()


# TOBECONTINUED
# add data info from final_results to well_info
well_info <- readxl::read_excel("raw_data/Training REspiro_080526.xlsx",
                                sheet = sheet_all[1]) %>%
  janitor::clean_names() %>% 
  # define auto-interpretation from well only (?)
  # dplyr::mutate(
  #   assay_source = case_when(
  #     well %in% c() ~ "pneumo_bacter",
  #     
  #   )
  # ) %>% 
  glimpse()







