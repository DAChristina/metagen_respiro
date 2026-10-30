library(tidyverse)

source_all <- function(path = "global"){
  files <- list.files(path,
                      pattern = "\\.R$",
                      full.names = TRUE
  )
  invisible(lapply(files, source))
}
source_all()

current_data <- "raw_data/RESPIRO Summary Lab Results_20261030.xlsx"

# detect sheet first
sheet_all <- readxl::excel_sheets(
  current_data) %>% 
  glimpse()


test <- readxl::read_xlsx(current_data,
                          sheet = sheet_all[1],
                          skip = 3) %>% 
  # janitor::clean_names() %>% 
  glimpse()


# cleanup colnames & prepare df
colnames <- readxl::read_xlsx(current_data,
                              range = readxl::cell_rows(2:4),
                              col_names = FALSE,
                              col_types = "text",
                              .name_repair = "minimal") %>% 
  glimpse()
h <- as.matrix(colnames)

# cleanup stupid excel merged cells
# only top-left filled in csv file type
for (i in 1:2) {
  for (j in 2:ncol(h)) {
    if (is.na(h[i, j])) h[i, j] <- h[i, j - 1]
  }
}

# cleanup
new_names <- apply(h, 2, function(x) paste(na.omit(x), collapse = "_"))
new_names <- janitor::make_clean_names(new_names)

df <- readxl::read_xlsx(current_data,
                        skip = 4,
                        col_names = new_names) %>%
  janitor::remove_empty(c("rows", "cols")) %>%
  dplyr::rename(sample_date = sample_collection_date) %>% 
  glimpse()

long <- df %>%
  dplyr::select(id, starts_with("kit_")) %>%
  tidyr::pivot_longer(starts_with("kit_"),
                      names_to = c("kit", ".value"),
                      names_pattern = "^kit_([^_]+)_(.*)$") %>%
  dplyr::rename(kit_date = date)

# One row per sample & kit & result (split results and Ct)
res_long <- long %>%
  dplyr::mutate(results = stringr::str_squish(results),
                ct = stringr::str_replace_all(ct, "\\s*;\\s*", ",")
                ) %>%
  tidyr::separate_longer_delim(c(results, ct),
                               delim = ",") %>%
  dplyr::mutate(results = stringr::str_squish(results),
                ct = as.numeric(str_squish(ct)),
                status = case_when(
                  is.na(results) ~ "missing",
                  results == "-" ~ "negative",
                  str_to_lower(results) == "rerun" ~ "rerun",
                  TRUE  ~ "positive")) %>% 
  glimpse()

# "-" (negative) and "rerun" are dropped from counts
pos <- dplyr::filter(res_long, status == "positive")

# Recalculate n_total and type PER KIT, then widen
kit_summary <- res_long %>%
  dplyr::group_by(id, kit) %>%
  dplyr::summarise(n_total = sum(status == "positive"),
                   type = if (n_total == 0) "-" else
                     paste(results[status == "positive"], collapse = ","),
                   rerun = any(status == "rerun"),
                   .groups = "drop")

kit_wide <- kit_summary %>%
  tidyr::pivot_wider(names_from = kit,
                     values_from = c(n_total, type, rerun),
                     names_glue = "{.value}_{kit}")

# Ct values in wide format
ct_wide <- pos %>%
  tidyr::pivot_wider(id_cols = id,
                     names_from = c(kit, results),
                     values_from = ct,
                     names_sep = "_") %>%
  janitor::clean_names() %>% 
  glimpse()

# Overall n_total and type per sample (distinct positives across kits)
overall <- pos %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(n_total_new = n_distinct(results),
                   type_new = paste(unique(results), collapse = ","),
                   .groups = "drop") %>% 
  glimpse()

final <- df %>%
  dplyr::select(id, sample_date,
                n_total_old = n_pathogen_found,
                type_old = pathogen_s_found
  ) %>%
  dplyr::left_join(overall, by = "id") %>%
  dplyr::left_join(kit_wide, by = "id") %>%
  dplyr::left_join(ct_wide, by = "id") %>%
  dplyr::mutate(
    # reduce Flu-A in detected type, recalculate Flu A
    n_total_new = ifelse(stringr::str_detect(type_new, "Flu A,Flu A"),
                         n_total_new - 1,
                         n_total_new
    ),
    n_total_rp1a = coalesce(ifelse(stringr::str_detect(type_new, "Flu A,Flu A"),
                                   n_total_new - 1,
                                   n_total_new),
                            0
                            ),
    
    type_new = case_when(
      type_new == "Flu A" ~ "Flu A (unsubtyped)",
      is.na(type_new) & type_old == "Negative" ~ "Negative",
      TRUE ~ type_new
    ),
    type_new = stringr::str_replace(type_new, "Flu A,Flu A", "Flu A"
    ),
    
    # correct n_total_new as 0 to just only negative samples
    n_total_new = ifelse(type_new == "Negative",0, n_total_new
    ),
    
    any_rerun = if_any(starts_with("rerun_"), identity),
    n_matches_old = as.numeric(n_total_old) == n_total_new,
  ) %>% 
  glimpse()


filter(final, !n_matches_old) # test samples where the old n_total disagrees

write.csv(final, "inputs/lab_qPCR_cleaned.csv", row.names = F)
