library(tidyverse)

# load epiData & logbook
df <- readxl::read_excel("raw_data/LOGBOOK SUBJEK RESPIRO_29 JUN 2026.xlsx") %>% 
  janitor::clean_names() %>% 
  glimpse()

# issues ALWAYS in dates!