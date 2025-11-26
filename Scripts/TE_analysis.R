#!/beegfs/data/soft/R-4.3.1/bin/Rscript

# ----------------------------
# Chargement des librairies
# ----------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(ggrepel)
  library(ggsci)
  library(viridis)
  library(stringr)
  library(tools)
})

# ----------------------------
# 0. Récupérer et valider les arguments
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("❌ Utilisation : Rscript TE_analysis.R <INS_files...> <TE_count_files...> <out_folder>")
}

# Le dernier argument est le dossier de sortie
out_folder <- tail(args, 1)

# Le reste correspond aux fichiers INS et TE_counts
n_total <- length(args) - 1
if (n_total %% 2 != 0) {
  stop("❌ Le nombre total d'arguments (hors dossier de sortie) doit être pair — même nombre de fichiers INS et TE_counts.")
}

n <- n_total / 2
ins_files <- args[1:n]
count_files <- args[(n + 1):(2 * n)]

# Vérification de l'existence des fichiers
missing_ins <- ins_files[!file.exists(ins_files)]
missing_counts <- count_files[!file.exists(count_files)]

if (length(missing_ins) > 0) {
  warning("⚠️ Fichiers INS manquants :\n", paste(missing_ins, collapse = "\n"))
}
if (length(missing_counts) > 0) {
  warning("⚠️ Fichiers TE_counts manquants :\n", paste(missing_counts, collapse = "\n"))
}

cat("Fichiers INS chargés (", length(ins_files), ") :\n", paste(ins_files, collapse = "\n"), "\n\n")
cat("Fichiers TE_counts chargés (", length(count_files), ") :\n", paste(count_files, collapse = "\n"), "\n\n")

# ----------------------------
# 1. Charger tous les fichiers INS
# ----------------------------
read_INS_bed <- function(file) {
  df <- read.delim(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # Garde le nom complet du genome
  strain <- basename(file) %>% file_path_sans_ext()
  strain <- sub("_INS_100$", "", strain)
  df$Strain <- strain
  
  colnames(df) <- make.names(colnames(df))
  if(!"FREQ" %in% colnames(df)) df$FREQ <- NA
  if(!"SIZE_TE" %in% colnames(df)) df$SIZE_TE <- NA
  if(!"FREQ_WITH_CLIPPED" %in% colnames(df)) df$FREQ_WITH_CLIPPED <- NA
  
  df$FREQ <- as.numeric(df$FREQ)
  df$SIZE_TE <- as.integer(df$SIZE_TE)
  df$FREQ_WITH_CLIPPED <- as.numeric(df$FREQ_WITH_CLIPPED)
  
  return(df)
}

all_INS <- bind_rows(lapply(ins_files, read_INS_bed))
cat("✔️ Fichiers INS chargés :", length(ins_files), "\n")

# ----------------------------
# 2. Charger tous les fichiers TE_counts
# ----------------------------
read_TE_count <- function(file) {
  df <- read.delim(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  strain <- basename(file) %>% file_path_sans_ext()
  strain <- sub("_TE_counts$", "", strain)
  df$Strain <- strain
  return(df)
}

all_TE_counts <- bind_rows(lapply(count_files, read_TE_count))
cat("✔️ Fichiers TE_counts chargés :", length(count_files), "\n")

# ----------------------------
# 3. Définir l'ordre des souches selon l’ordre d’apparition dans les fichiers
# ----------------------------
strain_order <- sapply(ins_files, function(f) basename(f) %>% file_path_sans_ext() %>% sub("_INS_100$", "", .))
all_INS$Strain <- factor(all_INS$Strain, levels = strain_order)
all_TE_counts$Strain <- factor(all_TE_counts$Strain, levels = strain_order)

# ----------------------------
# 4. Préparer les données TE_counts
# ----------------------------
all_TE_counts$Count <- as.numeric(all_TE_counts$Count)

TE_count_sum_family <- all_TE_counts %>%
  group_by(Strain, Family) %>%
  summarize(Count = sum(Count), Class = unique(Class), .groups="drop") %>%
  arrange(Strain, Count)

TE_count_sum_class <- all_TE_counts %>%
  group_by(Strain, Class) %>%
  summarize(Count = sum(Count), .groups="drop") %>%
  arrange(Strain, Count)

TE_count_sum_family$Family <- factor(TE_count_sum_family$Family)
TE_count_sum_family$Class <- factor(TE_count_sum_family$Class)
TE_count_sum_class$Class <- factor(TE_count_sum_class$Class)

# ----------------------------
# 5. Plots TE_family et TE_class
# ----------------------------
A <- ggplot(TE_count_sum_family, aes(x = fct_reorder(Family, Count), y = Count, fill = Family)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), size = 3) +
  coord_flip() +
  labs(x = "TE superfamily", y = "Count") +
  theme_classic() +
  scale_fill_d3(palette = "category20c") +
  theme(legend.position = "none") +
  facet_wrap(~ Strain)
ggsave(file.path(out_folder, "TE_count_sum_family.png"), A , width = 12, height = 8)

B <- ggplot(TE_count_sum_class, aes(x = Strain, y = Count, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 4) +
  labs(x = "Population", y = "Insertion number") +
  theme_classic() +
  scale_fill_d3(palette = "category20c") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12))
ggsave(file.path(out_folder, "TE_count_sum_class.png"), B , width = 12, height = 8)

# ----------------------------
# 6. Préparer les INS pour analyse par chromosome
# ----------------------------
all_INS <- all_INS %>%
  mutate(Family = str_extract(V4, "(?<=/)[^|]+"),
         Class = str_extract(V4, "(?<=#)[^/]+"),
         TE = str_extract(V4, "^[^_]+_[^#]+"))

all_INS$V1 <- gsub("_RagTag", "", all_INS$V1)

All_INS_agg <- all_INS %>%
  group_by(V1, Class, Strain) %>%
  summarise(Count = n(), .groups="drop") %>%
  filter(V1 %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chrX", "chrY", "chr4"))

All_INS_agg$Strain <- factor(All_INS_agg$Strain, levels = strain_order)

C <- ggplot(All_INS_agg, aes(x = Strain, y = Count, fill = Class)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  labs(x = "Chromosome", y = "TE Count") +
  scale_fill_d3(palette = "category20c") +
  theme_classic() +
  facet_grid(~ V1) +
  ylim(0, 140) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.2),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14))
ggsave(file.path(out_folder, "All_INS_agg.png"), C , width = 12, height = 8)

# ----------------------------
# 7. Plots par classe (LTR, LINE, DNA)
# ----------------------------
all_strains <- unique(all_TE_counts$Strain)
all_tes <- unique(all_TE_counts$TE)

TE_count_completed <- expand.grid(TE = all_tes, Strain = all_strains)
merged_data <- TE_count_completed %>%
  left_join(all_TE_counts, by = c("TE", "Strain")) %>%
  mutate(Count = ifelse(is.na(Count), 0, Count)) %>%
  group_by(TE) %>%
  fill(Class, Family, .direction = "downup") %>%
  ungroup()

merged_data$Strain <- factor(merged_data$Strain, levels = strain_order)
my_colors <- c("#a3a500", "#00bf7d", "#00b0f6", "#e76bf3")

plot_by_class <- function(class_name) {
  ggplot(subset(merged_data, Class == class_name), 
         aes(x = Strain, y = Count, fill = Strain)) +
    geom_bar(stat = "identity", color = "black", width = 0.5, linewidth = 0.2) +
    coord_flip() +
    labs(title = class_name, x = "Population", y = "Count") +
    theme_classic() +
    scale_fill_manual(values = my_colors) +
    scale_y_continuous(breaks = function(x) floor(pretty(x))) +
    theme(legend.position = "none",
          axis.text.x = element_text(hjust = 1, size = 7),
          axis.text.y = element_text(angle = 45, size = 6),
          axis.title = element_text(size = 10),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(size = 6, face = "italic"),
          aspect.ratio = 1,
          axis.line = element_line(colour = 'black', size = 0.5)) +
    facet_wrap(~ TE, scales = "free")
}

ggsave(file.path(out_folder, "TE_count_LTR.png"), plot_by_class("LTR"), width = 12, height = 8)
ggsave(file.path(out_folder, "TE_count_LINE.png"), plot_by_class("LINE"), width = 12, height = 8)
ggsave(file.path(out_folder, "TE_count_DNA.png"), plot_by_class("DNA"), width = 12, height = 8)

cat("✅ Analyse terminée. Résultats enregistrés dans :", out_folder, "\n")
