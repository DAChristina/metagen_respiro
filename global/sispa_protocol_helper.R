parse_nanostats <- function(file_path) {
  lines <- readr::read_lines(file_path)
  
  # Extract numeric/value part after the colon
  extract_val <- function(pattern) {
    line <- grep(pattern, lines, value = TRUE)
    
    if (length(line) == 0) {
      return(NA)
    }
    
    val <- stringr::str_trim(stringr::str_replace(line[1], "^[^:]*:", ""))
    
    # Take only the first value
    val <- stringr::str_split(val, "\\s+")[[1]][1]
    
    # Remove trailing .0
    val <- stringr::str_remove(val, "\\.0$")
    
    return(val)
  }
  
  # Standard NanoStats metrics
  num_reads <- extract_val("^Number of reads:")
  total_bases <- extract_val("^Total bases:")
  mean_len <- extract_val("^Mean read length:")
  median_len <- extract_val("^Median read length:")
  n50   <- extract_val("^Read length N50:")
  mean_q  <- extract_val("^Mean read quality:")
  median_q <- extract_val("^Median read quality:")
  
  # Quality cutoff percentages
  # Actual NanoStats format:
  # >Q10: 1290519 (92.2%) 744.7Mb
  extract_q_pct <- function(q_level) {
    
    line <- grep(
      paste0("^>Q", q_level, ":"),
      lines,
      value = TRUE
    )
    
    if (length(line) == 0) {
      return(NA)
    }
    
    pct <- stringr::str_match(
      line[1],
      "\\(([0-9.]+)%\\)"
    )[, 2]
    
    return(pct)
  }
  
  q10_pct <- extract_q_pct(10)
  q15_pct <- extract_q_pct(15)
  q20_pct <- extract_q_pct(20)
  
  # Longest read
  # Actual NanoStats format:
  # 1: 115508 (12.1)
  longest_line_idx <- grep(
    "(Top 5 longest reads)",
    lines
  )
  
  longest_read_str <- NA
  
  if (length(longest_line_idx) > 0) {
    
    first_longest_line <- lines[longest_line_idx[1] + 1]
    
    match <- stringr::str_match(
      first_longest_line,
      "^1:\\s+(\\d+)\\s+\\(([0-9.]+)\\)"
    )
    
    if (!is.na(match[1, 1])) {
      
      bp <- as.numeric(match[1, 2])
      q_val <- match[1, 3]
      
      longest_read_str <- paste0(
        scales::comma(bp),
        " (Q",
        q_val,
        ")"
      )
    }
  }
  
  # reformat
  fmt_num <- function(x) {
    if (is.na(x)) {
      return(NA)
    }
    
    num <- as.numeric(gsub("[^0-9.]", "", x))
    
    if (!is.na(num)) {
      return(scales::comma(num))
    }
    
    return(x)
  }
  
  # return
  tibble(
    Metric = c(
      "Number of reads",
      "Total bases",
      "Mean read length (bp)",
      "Median read length (bp)",
      "Read length N50 (bp)",
      "Mean read quality (Q)",
      "Median read quality (Q)",
      "Reads >Q10 (%)",
      "Reads >Q15 (%)",
      "Reads >Q20 (%)",
      "Longest read (bp / Q)"
    ),
    
    Value = c(
      fmt_num(num_reads),
      fmt_num(total_bases),
      mean_len,
      median_len,
      n50,
      mean_q,
      median_q,
      q10_pct,
      q15_pct,
      q20_pct,
      longest_read_str
    )
  )
}

parse_kraken_report <- function(file_path, barcode_name, classification_type) {
  # Kraken2 report structure:
  # V1: % clade, V2: clade reads, V3: taxon reads, V4: rank, V5: taxon ID, V6: name
  report <- readr::read_tsv(file_path,
                            col_names = c("pct", "reads", "direct_reads", "rank", "tax_id", "name"),
                            col_types = "ddicic",
                            show_col_types = FALSE)
  
  # Filter for species level (rank S) and pull top 5 hits
  species_data <- report %>%
    dplyr::filter(rank == "S") %>%
    dplyr::mutate(
      name = stringr::str_trim(name),
      barcode = barcode_name,
      classification = classification_type
    ) %>%
    arrange(desc(reads)) %>%
    head(5) %>%
    dplyr::select(barcode, name, pct, reads, classification)
  
  return(species_data)
}

