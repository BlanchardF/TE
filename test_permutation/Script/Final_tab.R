library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

# ==============================
# Fonction de parsing des summary.txt
# ==============================
parse_summary <- function(filepath) {
  # Déduire la population (ex: S100)
  sample <- str_extract(filepath, "S[0-9]+")
  
  # Déduire la fréquence
  freq <- case_when(
    str_detect(filepath, "freq_20_40") ~ "20-40",
    str_detect(filepath, "freq_40_60") ~ "40-60",
    str_detect(filepath, "freq_60_80") ~ "60-80",
    str_detect(filepath, "freq_80") ~ "80",
    str_detect(filepath, "filtered_cleaned_freq") ~ "-20",
    TRUE ~ "all"
  )
  
  # Lire le fichier
  df <- read.table(filepath, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(df) <- c("region", "metric", "value")
  
  # Passer les métriques en colonnes
  df <- df %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    mutate(true_count = as.integer(true_count),
           mean_perm = as.numeric(mean_perm),
           pvalue = as.numeric(pvalue))
  
  # Harmoniser les noms de régions
  df <- df %>%
    mutate(region = case_when(
      str_detect(region, "promoters.2000") ~ "Promoter (2000)",
      str_detect(region, "5UTRs") ~ "5UTR",
      str_detect(region, "3UTRs") ~ "3UTR",
      str_detect(region, "exons") ~ "Exon",
      str_detect(region, "introns") ~ "Intron",
      TRUE ~ "Intergenic"
    ),
    Frequencies = freq,
    Population = sample)
  
  return(df)
}

# ==============================
# Lecture de tous les fichiers summary.txt
# ==============================
files <- list.files("~/TE/TE/test_permutation/Results/All",
                    pattern = "summary.txt",
                    recursive = TRUE,
                    full.names = TRUE)

all_data <- map_dfr(files, parse_summary)

# ==============================
# TABLEAU 1 : format large 
# ==============================
tableau_large <- all_data %>%
  select(region, Frequencies, Population, true_count) %>%
  pivot_wider(names_from = Population, values_from = true_count, values_fill = 0) %>%
  mutate(Frequencies = factor(Frequencies, levels = c("all", "80", "60-80", "40-60", "20-40", "-20")),
         region = factor(region, levels = c("Promoter (2000)", "5UTR", "3UTR", "Intron", "Exon", "Intergenic"))) %>%
  arrange(Frequencies, region) %>%
  select(region, Frequencies, S10, S31, S73, S100)

write.table(tableau_large,
            "summary_final.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ==============================
# TABLEAU 2 : format long (toutes les infos)
# ==============================
tableau_long <- all_data %>%
  mutate(Frequencies = factor(Frequencies, levels = c("all", "80", "60-80", "40-60", "20-40", "-20")),
         region = factor(region, levels = c("Promoter (2000)", "5UTR", "3UTR", "Intron", "Exon", "Intergenic"))) %>%
  arrange(Frequencies, region, Population)

write.table(tableau_long,
            "summary_complete.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

