library(tidyverse)
source_all <- function(path = "global"){
  files <- list.files(path,
                      pattern = "\\.R$",
                      full.names = TRUE
  )
  invisible(lapply(files, source))
}
source_all()

pathogen_list_inspect <- read.csv("inputs/pathogen_priority_panels_taxid.csv") %>% 
  dplyr::mutate(dplyr::across(everything(),
                              as.character)
                ) %>% 
  glimpse()

bracken_folder <- Sys.glob("inputs/output_wfmetagenomics_with_host_k2_db_*/bracken")
kraken_folder <- Sys.glob("inputs/output_wfmetagenomics_with_host_k2_db_*/kraken2")

bracken_compiled <- compile_brackenfiles(
  path = bracken_folder,
  name_file = ".*\\_bracken.report$"
) %>% 
  glimpse()

kraken_compiled <- compile_krakenfiles(
  path = kraken_folder,
  name_file = ".*\\kraken2.report.txt$"
) %>% 
  # filter out unclassified & domain-level summary
  # dplyr::filter(!rank %in% c("U", "D")) %>%
  glimpse()

pathogen_hits <- kraken_compiled %>%
  # dplyr::filter(rank == "S") %>%
  dplyr::semi_join(
    pathogen_list_inspect
    ,
    by = c("taxid" = "k2_taxid")
    ) %>% 
  dplyr::mutate(reads_clade = as.numeric(reads_clade),
                # pct is (reads_clade/total_reads_in_sample)*100
                confidence = case_when(
                  reads_clade < 10  ~ "<10, likely noise/background",
                  # reads_clade < 50 ~ "<50, low confidence (verify?)",
                  TRUE ~ "worth investigating"
                )) %>%
  dplyr::relocate(sample_id) %>%
  glimpse()

write.table(bracken_compiled,
            "inputs/bracken_compiled_long.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE
            )

write.table(kraken_compiled,
            "inputs/kraken_compiled_long.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE
            )

write.table(pathogen_hits,
            "inputs/kraken_compiled_priority_pathogen_hits_long.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE
            )
