library(see)

coverage_eval <- function(directory, pattern){
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(data.table)
  
  files <- list.files(directory, pattern=pattern, full.names=TRUE)
  files <- files[!grepl("18-110", files)]
  # Extract sample names from filenames (adjust pattern to your case)
  sample_names <- gsub(paste0(pattern, "$"), "", basename(files))
  
  # Read and combine
  coverage_list <- map2(files, sample_names, ~ {
    df <- fread(.x)   # adjust if header present
    colnames(df) <- c("chr", "start", "end", "exon_id", "coverage")
    df %>% select(exon_id, coverage) %>% mutate(sample = .y)
  })
  
  coverage_long <- bind_rows(coverage_list)
  
  # Reshape into exon × sample matrix
  coverage_wide <- coverage_long %>%
    pivot_wider(names_from = sample, values_from = coverage)
  
  # Baseline albumin data
  alb <- list.files("ALB/", pattern="*ALB*", full.names=TRUE)
  alb <- files[!grepl("18-110", files)]

  # Read and combine
  alb_coverage_list <- map2(alb, sample_names, ~ {
    df <- fread(.x)   # adjust if header present
    colnames(df) <- c("chr", "start", "end", "exon_id", "coverage")
    df %>% select(exon_id, coverage) %>% mutate(sample = .y)
  })
  
  alb_coverage_long <- bind_rows(alb_coverage_list)
  
  # Reshape into exon × sample matrix
  alb_coverage_wide <- alb_coverage_long %>%
    pivot_wider(names_from = sample, values_from = coverage)
  
  # Now, only patient's data
  
  hhl <- list.files(directory, pattern = pattern, full.names = T)
  hhl <- hhl[grepl("18-110", hhl)]
  hhl_names <- gsub(paste0(pattern, "$"), "", basename(hhl))
  coverage_list_hhl <- map2(hhl, hhl_names, ~ {
    df <- fread(.x)   # adjust if header present
    colnames(df) <- c("chr", "start", "end", "exon_id", "coverage")
    df %>% select(exon_id, coverage) %>% mutate(sample = .y)
  })
  
  coverage_wide_hhl <- bind_rows(coverage_list_hhl) %>% 
    pivot_wider(names_from = sample, values_from = coverage)
  
  # Patient's albumin values
  hhl_alb <- list.files("ALB/", pattern = "*ALB*", full.names = T)
  hhl_alb <- hhl_alb[grepl("18-110", hhl_alb)]
  
  coverage_list_hhl_alb <- map2(hhl_alb, hhl_names, ~ {
    df <- fread(.x)   # adjust if header present
    colnames(df) <- c("chr", "start", "end", "exon_id", "coverage")
    df %>% select(exon_id, coverage) %>% mutate(sample = .y)
  })
  
  coverage_wide_hhl_alb <- bind_rows(coverage_list_hhl_alb) %>% 
    pivot_wider(names_from = sample, values_from = coverage)
  
  # Creating Coverage/area ratios for each control sample
  exon_ids <- coverage_wide$exon_id
  coverage_matrix <- coverage_wide %>% select(-exon_id)
  alb_coverage_matrix <- alb_coverage_wide %>% select(-exon_id)
  # compute column sums (total coverage per sample)
  alb_totals <- colSums(alb_coverage_matrix, na.rm = TRUE)
  # divide each exon coverage by total of its sample
  coverage_ratio <- sweep(coverage_matrix, 2, alb_totals, FUN = "/")
  # reattach exon_id
  coverage_ratio_df <- data.frame(exon_id = exon_ids, coverage_ratio)
  
  # Same thing for HHL samples
  
  coverage_matrix_hhl <- coverage_wide_hhl %>% select(-exon_id)
  coverage_matrix_hhl_alb <- coverage_wide_hhl_alb %>% select(-exon_id)
  totals_hhl_alb <- colSums(coverage_matrix_hhl_alb, na.rm = TRUE)
  coverage_ratio_hhl <- sweep(coverage_matrix_hhl, 2, totals_hhl_alb, FUN = "/")
  coverage_ratio_df_hhl <- data.frame(exon_id = exon_ids, coverage_ratio_hhl)
  
  # Calculating the mean ratios
  
  mean_ratio_df <- coverage_ratio_df %>%
    rowwise() %>%
    mutate(mean_ratio = mean(c_across(-exon_id), na.rm = TRUE)) %>%
    ungroup()
  
  final_ratios_HHL <- data.frame(
    exon_id = exon_ids,
    `18-1102` = coverage_ratio_df_hhl$X18.1102/mean_ratio_df$mean_ratio,
    `18-1103` = coverage_ratio_df_hhl$X18.1103/mean_ratio_df$mean_ratio,
    `18-1104` = coverage_ratio_df_hhl$X18.1104/mean_ratio_df$mean_ratio
  )
  
  return(final_ratios_HHL)
}

plotting_coverage <- function(df, exon_number){
  library(ggplot2)
  library(see)
  
  ratio_long <- df %>%
    pivot_longer(-exon_id, names_to = "sample", values_to = "ratio") %>%
    mutate(exon_num = as.numeric(gsub("exon", "", exon_id))) %>%
    arrange(exon_num)
  
  # Ensure exons are ordered 1–23
  ratio_long$exon_id <- factor(ratio_long$exon_id, 
                               levels = paste0("exon", 1:exon_number))
  
  pl <- ggplot(ratio_long, aes(x = exon_id, y = ratio, fill = sample)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme_minimal(base_size = 14) +
    labs(x = "Exon", y = "Coverage ratio", fill = "Sample") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(pl)
}

canonical <- coverage_eval("./", ".regions.bed.gz")
pl_canon <- plotting_coverage(canonical, 23)
pl_final <- pl_canon + scale_fill_okabeito(labels = c("X18.1102" = "III:3", 
                                          "X18.1103" = "II:3", 
                                          "X18.1104" = "II:4"))

ggsave("../WES_coverage_ATP2B2.png", plot = pl_final, dpi = 600, width = 15, height = 8)
