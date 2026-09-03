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

search_taxon <- function(df, pattern) {
  df %>% 
    dplyr::filter(stringr::str_detect(k2_name,
                                      regex(pattern, ignore_case = TRUE)))
}

prepare_batch <- function(df1, df2, batch, batch_address){
  workLab_logbook <- readxl::read_excel(df1,
                                        skip = 1)
  
  # generate targeted barcode list based on labWork data
  # use ls <folder>/* > raw_data/demux_batch<batch>.csv
  demux_list <- read.csv(df2,
                         header = FALSE)
  
  out <- workLab_logbook %>% 
    janitor::clean_names() %>% 
    dplyr::rename_all(~ paste0("workLab_", .)) %>% 
    # adjust with merged columns
    tidyr::fill(everything(),
                .direction = "down") %>%
    # strict to batch & native barcode
    dplyr::mutate(
      workLab_batch = batch,
      workLab_native_barcode = paste0("barcode",
                                      sprintf("%02d", workLab_native_barcode))
    ) %>% 
    dplyr::left_join(
      demux_list %>% 
        dplyr::rename(workDmx_native_barcode = V1) %>% 
        dplyr::mutate(
          workLab_native_barcode = stringr::str_extract(workDmx_native_barcode,
                                                        "barcode\\d{2}"),
          workDmx_address = paste0(
            batch_address,
            workDmx_native_barcode
          ),
          workDmx_filename = gsub(".fastq", "", workDmx_native_barcode),
        )
      ,
      by = "workLab_native_barcode"
    ) %>% 
    dplyr::filter(!grepl("fastq\\.gz", workDmx_address))
  
  # generate wf_metagenomics samplesheet
  wf_metagenomics_batch <- out %>% 
    dplyr::transmute(
      sample_id = paste0(workDmx_filename, ".fastq.gz"),
      long_reads = paste0(batch_address,
                          sample_id),
    ) %>% 
    glimpse()
  
  write.csv(wf_metagenomics_batch, paste0("inputs/wf_metagenomics_samplesheet_batch", batch, ".csv"),
            na = "", quote = FALSE,
            row.names = FALSE)
  
  # generate foldered file list
  write.table(wf_metagenomics_batch %>% 
                dplyr::transmute(
                  long_reads = gsub(".gz", "", long_reads)
                )
              ,
              paste0("inputs/wf_metagenomics_foldered_batch", batch, ".tsv"),
              na = "", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  
  
  return(out)
}

compile_qc_batch <- function(df1, df2, batch){
  
  wf_metagenomics_out <- read.csv(df2)
  
  out <- df1 %>% 
    dplyr::left_join(
      # wf_metagenomics_out_post
      
      # TEMPORARY read.csv data on my local PC
      wf_metagenomics_out %>%
        dplyr::mutate(
          sample_id_post = gsub(".fastq.gz", "", sample_id_post),
        ) %>% 
        dplyr::rename_all(~ paste0("workSeq_", .))
      ,
      by = c("workDmx_filename" = "workSeq_sample_id_post")
    ) %>% 
    # workLab species name correction
    dplyr::mutate(
      workLab_isolate = case_when(
        str_detect(workLab_isolate, "E\\. coli|E\\.coli") ~ "Escherichia coli",
        TRUE ~ workLab_isolate
      )
    )
  
  write.csv(out,
            paste0("inputs/worklab_workSeq_compiled_batch", batch, ".csv"),
            na = "", quote = FALSE,
            row.names = FALSE)
  
  return(out)
}

compile_csvfiles <- function(path = "inputs"){
  files <- list.files(path,
                      pattern = "worklab_workSeq_compiled_batch.*\\.csv$",
                      full.names = TRUE
  )
  
  return(files %>% 
           lapply(function(f){
             read.csv(f) %>% 
               dplyr::mutate(dplyr::across(everything(),
                                           as.character))
           }) %>% 
           dplyr::bind_rows()
  )
}

compile_brackenfiles <- function(path, name_file){
  files <- list.files(path,
                      pattern = name_file,
                      full.names = TRUE,
                      recursive = TRUE
  )
  
  files <- files[!grepl("compiled", basename(files))]
  
  return(files %>% 
           lapply(function(f){
             read.delim(f,
                        sep = "\t",
                        header = TRUE) %>% 
               dplyr::mutate(dplyr::across(everything(),
                                           as.character),
                             sample_id = basename(f),
                             sample_id = gsub(".kraken2_bracken.report", "",
                                              sample_id),
                             )
           }) %>% 
           dplyr::bind_rows()
  )
}

compile_krakenfiles <- function(path, name_file){
  files <- list.files(path,
                      pattern = name_file,
                      full.names = TRUE,
                      recursive = TRUE
  )
  
  files <- files[!grepl("compiled", basename(files))]
  
  # annoying adjustment
  return(files %>% 
           lapply(function(f){
             read.delim(f,
                        sep = "\t",
                        header = FALSE,
                        col.names = c("pct", "reads_clade", "reads_taxon",
                                      "rank", "taxid", "name")
             ) %>% 
               dplyr::mutate(dplyr::across(everything(),
                                           as.character),
                             sample_id = basename(f),
                             sample_id = gsub(".kraken2.report.txt", "",
                                              sample_id),
               )
           }) %>% 
           dplyr::bind_rows()
  )
}


