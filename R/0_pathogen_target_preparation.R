library(tidyverse)
source_all <- function(path = "global"){
  files <- list.files(path,
                      pattern = "\\.R$",
                      full.names = TRUE
  )
  invisible(lapply(files, source))
}
source_all()

# load db
viral_lib_report <- readr::read_tsv(
  "inputs/k2_viral_20260626/library_report.tsv", # instead of inspect.txt
  comment = "#",
  quote = "",
  col_names = c("k2_library", "k2_name", "k2_url"),
  col_types = "ccc"
) %>%
  dplyr::mutate(k2_name = stringr::str_trim(k2_name)) %>%
  glimpse()

viral_inspect <- readr::read_tsv(
  "inputs/k2_viral_20260626/inspect.txt",
  comment = "#",
  col_names = c("k2_pct", "k2_reads_clade", "k2_reads_taxon",
                "k2_rank", "k2_taxid", "k2_name"),
  col_types = "dddccc"
) %>%
  dplyr::mutate(k2_name = stringr::str_trim(k2_name)) %>%
  glimpse()

standard8_lib_report <- readr::read_tsv(
  "inputs/library_report_standard8.tsv",
  comment = "#",
  quote = "",
  col_names = c("k2_library", "k2_name", "k2_url"),
  col_types = "ccc"
) %>%
  dplyr::mutate(k2_name = stringr::str_trim(k2_name)) %>%
  glimpse()

standard8_inspect <- readr::read_tsv(
  "inputs/inspect_standard8.txt",
  comment = "#",
  col_names = c("k2_pct", "k2_reads_clade", "k2_reads_taxon",
                "k2_rank", "k2_taxid", "k2_name"),
  col_types = "dddccc"
) %>%
  dplyr::mutate(k2_name = stringr::str_trim(k2_name)) %>%
  glimpse()

pathogen_list_lib_report <- 
  # panel 1A source: https://www.seegene.com/products/assays/detail/24
  search_taxon(df = viral_lib_report,
               pattern = "Influenza A") %>% 
  dplyr::mutate(panel = "1A") %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_lib_report,
                 pattern = "Influenza B") %>% 
      dplyr::mutate(panel = "1A")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_lib_report,
                 pattern = "Respiratory syncytial virus A") %>% 
      dplyr::mutate(panel = "1A")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_lib_report,
                 pattern = "Respiratory syncytial virus B") %>% 
      dplyr::mutate(panel = "1A")
  ) %>% 
  
  # panel 2 source: https://www.seegene.com/products/assays/detail/30
  dplyr::bind_rows(
    search_taxon(df = viral_lib_report,
                 pattern = "adenovirus") %>% 
      dplyr::mutate(panel = "2")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_lib_report,
                 pattern = "enterovirus") %>% 
      dplyr::mutate(panel = "2")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_lib_report,
                 pattern = "metapneumovirus") %>% 
      dplyr::mutate(panel = "2")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_lib_report,
                 pattern = "parainfluenza") %>% 
      dplyr::mutate(panel = "2")
    # only Human parainfluenza virus 4a available (not 1,2,3)
  ) %>% 
  
  # panel 3 source: https://www.seegene.com/products/assays/detail/31
  dplyr::bind_rows(
    search_taxon(df = viral_lib_report,
                 pattern = "bocavirus") %>% 
      dplyr::mutate(panel = "3")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_lib_report,
                 pattern = "coronavirus 229E") %>% 
      dplyr::mutate(panel = "3")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_lib_report,
                 pattern = "coronavirus NL63") %>% 
      dplyr::mutate(panel = "3")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_lib_report,
                 pattern = "Human coronavirus OC43") %>% 
      dplyr::mutate(panel = "3")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_lib_report,
                 pattern = "rhinovirus") %>% 
      dplyr::mutate(panel = "3")
  ) %>% 
  
  # panel PneumoBacter source https://www.seegene.com/products/assays/detail/28
  dplyr::bind_rows(
    search_taxon(df = standard8_lib_report,
                 pattern = "Bordetella parapertussis") %>% 
      dplyr::mutate(panel = "PneumoBacter")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = standard8_lib_report,
                 pattern = "Bordetella pertussis") %>% 
      dplyr::mutate(panel = "PneumoBacter")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = standard8_lib_report,
                 pattern = "Mycoplasma pneumoniae") %>% 
      dplyr::mutate(panel = "PneumoBacter")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = standard8_lib_report,
                 pattern = "Chlamydophila pneumoniae") %>% 
      dplyr::mutate(panel = "PneumoBacter")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = standard8_lib_report,
                 pattern = "Haemophilus influenzae") %>% 
      dplyr::mutate(panel = "PneumoBacter")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = standard8_lib_report,
                 pattern = "Legionella pneumophila") %>% 
      dplyr::mutate(panel = "PneumoBacter")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = standard8_lib_report,
                 pattern = "Streptococcus pneumoniae") %>% 
      dplyr::mutate(panel = "PneumoBacter")
  ) %>% 
  
  # adjust NCBI name for easier download (~2,200 files!?)
  dplyr::mutate(
    accession = stringr::str_extract(k2_name, "^\\S+"),
    # rest = stringr::str_remove(k2_name, "^\\S+\\s+"),
    description = stringr::str_extract(stringr::str_remove(k2_name, "^\\S+\\s+"),
                                       "^.*(?=,\\s*[^,]+$)"),
    suffix = stringr::str_extract(stringr::str_remove(k2_name, "^\\S+\\s+"),
                                  "(?<=,\\s)[^,]+$")
  ) %>% 
  glimpse()

pathogen_list_inspect <- 
  # panel 1A source: https://www.seegene.com/products/assays/detail/24
  search_taxon(df = viral_inspect,
               pattern = "Influenza A") %>% 
  dplyr::mutate(panel = "1A") %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_inspect,
                 pattern = "Influenza B") %>% 
      dplyr::mutate(panel = "1A")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_inspect,
                 pattern = "Respiratory syncytial virus A") %>% 
      dplyr::mutate(panel = "1A")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_inspect,
                 pattern = "Respiratory syncytial virus B") %>% 
      dplyr::mutate(panel = "1A")
  ) %>% 
  
  # panel 2 source: https://www.seegene.com/products/assays/detail/30
  dplyr::bind_rows(
    search_taxon(df = viral_inspect,
                 pattern = "adenovirus") %>% 
      dplyr::mutate(panel = "2")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_inspect,
                 pattern = "enterovirus") %>% 
      dplyr::mutate(panel = "2")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_inspect,
                 pattern = "metapneumovirus") %>% 
      dplyr::mutate(panel = "2")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_inspect,
                 pattern = "parainfluenza") %>% 
      dplyr::mutate(panel = "2")
    # only Human parainfluenza virus 4a available (not 1,2,3)
  ) %>% 
  
  # panel 3 source: https://www.seegene.com/products/assays/detail/31
  dplyr::bind_rows(
    search_taxon(df = viral_inspect,
                 pattern = "bocavirus") %>% 
      dplyr::mutate(panel = "3")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_inspect,
                 pattern = "coronavirus 229E") %>% 
      dplyr::mutate(panel = "3")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_inspect,
                 pattern = "coronavirus NL63") %>% 
      dplyr::mutate(panel = "3")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_inspect,
                 pattern = "Human coronavirus OC43") %>% 
      dplyr::mutate(panel = "3")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = viral_inspect,
                 pattern = "rhinovirus") %>% 
      dplyr::mutate(panel = "3")
  ) %>% 
  
  # panel PneumoBacter source https://www.seegene.com/products/assays/detail/28
  dplyr::bind_rows(
    search_taxon(df = standard8_inspect,
                 pattern = "Bordetella parapertussis") %>% 
      dplyr::mutate(panel = "PneumoBacter")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = standard8_inspect,
                 pattern = "Bordetella pertussis") %>% 
      dplyr::mutate(panel = "PneumoBacter")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = standard8_inspect,
                 pattern = "Mycoplasma pneumoniae") %>% 
      dplyr::mutate(panel = "PneumoBacter")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = standard8_inspect,
                 pattern = "Chlamydophila pneumoniae") %>% 
      dplyr::mutate(panel = "PneumoBacter")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = standard8_inspect,
                 pattern = "Haemophilus influenzae") %>% 
      dplyr::mutate(panel = "PneumoBacter")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = standard8_inspect,
                 pattern = "Legionella pneumophila") %>% 
      dplyr::mutate(panel = "PneumoBacter")
  ) %>% 
  dplyr::bind_rows(
    search_taxon(df = standard8_inspect,
                 pattern = "Streptococcus pneumoniae") %>% 
      dplyr::mutate(panel = "PneumoBacter")
  ) %>% 
  # taxid can be used to link appearance in our data
  glimpse()

write.csv(pathogen_list_lib_report, "inputs/pathogen_priority_panels_ncbi.csv",
          row.names = F)

write.csv(pathogen_list_inspect, "inputs/pathogen_priority_panels_taxid.csv",
          row.names = F)
