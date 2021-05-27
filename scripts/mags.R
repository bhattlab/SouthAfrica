# Analyze South Africa metagenome-assembled genomes (MAGs)

library(dplyr)
library(ggplot2)
library(cowplot)
library(scales)
library(gtools)

# final list of HIV- samples to include
# samples <- as.character(read.table("/Users/Fiona/Desktop/za_tmp/samples_final.txt")$V1)
setwd("/Users/tamburif/za-microbiome/")

# samples sequenced with nanopore
nmag_samples <- c("C27", "C29", "C33")

# read binning output
# mags
mags <- read.table("input_final/bins_final.txt", sep = '\t', header = T)

# nmags
nmags <- read.table("input_final/binfinal_C27C29C33.tsv", sep = '\t', header = T)
# nmag_comparison <- read.table("/Users/tamburif/scg/za/nanopore_hq.txt", sep = '\t', header = T)

# assign quality
hq_filter <- quo(Completeness > 90 & Contamination < 5 & rna.5S >= 1 & rna.16S >= 1 & rna.23S >= 1 & tRNA >= 18)
mq_filter <- quo(Completeness >= 50 & Contamination < 10)

mags <- mutate(mags, Quality = ifelse(!!hq_filter, "High-quality", ifelse(!!mq_filter, "Medium-quality", "Low-quality")))
mags$Quality <- factor(mags$Quality, levels = c("Low-quality", "Medium-quality", "High-quality"))

nmags <- mutate(nmags, Quality = ifelse(!!hq_filter, "High-quality", ifelse(!!mq_filter, "Medium-quality", "Low-quality")))
nmags$Quality <- factor(nmags$Quality, levels = c("Low-quality", "Medium-quality", "High-quality"))

# medium quality or better bins
m <- mags_filt %>%
  filter(Quality != "Low-quality")
write.table(m, "input_final/MAGs_C27C29C33.txt", sep = "\t", row.names = F, quote = F)

# mags figure
mag_pal <- c("#FFDD8D", "#FF8C35", "#C84019")

mags_count <- mags %>%
  dplyr::count(Quality)

a <- ggplot(mags_count, aes(Quality, n, fill = Quality)) +
  geom_bar(stat = "identity") +
  theme_cowplot(12) +
  scale_fill_manual(values = mag_pal) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(
    x = "",
    y = "Number of genomes"
  )

b <- ggplot(mags, aes(x = Quality, y = N50/1000, fill = Quality)) +
  # geom_violin(trim = F) +
  geom_boxplot() +
  theme_cowplot(12) +
  scale_y_log10() +
  scale_fill_manual(values = mag_pal) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(
    x = "",
    y = "Contig N50 (kilobase)"
  )

plot_grid(a + theme(legend.position = "none"), b, rel_widths = c(0.42, 0.58), labels = c("A", "B"))
ggsave("final_plots/supplementary/mags.png", width = 6, height = 4)

# compare nmags to mags
mags_filt <- filter(mags, Sample %in% nmag_samples)

## summarize mags
mags_filt <- mags_filt %>%
  select(Sample, Bin, Final.Class, Completeness, Contamination, Strain.heterogeneity, Size.Mb, N50, Coverage, Quality) %>%
  arrange(Sample, desc(Quality)) %>%
  filter(Bin != "bin.unbinned")

mags_filt %>% dplyr::count(Quality)

mags_summary <- ggplot(mags_filt, aes(x = Quality, y = N50/1000, fill = Quality)) +
  geom_violin(trim = F) +
  theme_cowplot(12) +
  scale_y_log10() +
  labs(
    y = "Contig N50 (kilobase)"
  )

# ggsave("/Users/Fiona/Desktop/mags_summary.png", mags_summary, device = "png", width = 6, height = 4)


## summarize nmags
nmags <- nmags %>%
  select(Sample, Bin, Final.Class, Completeness, Contamination, Strain.heterogeneity, Size.Mb, N50, Coverage, Quality) %>%
  arrange(Sample, desc(Quality)) %>%
  filter(Bin != "bin.unbinned")

nmags %>% dplyr::count(Quality)

nmags_summary <- ggplot(nmags, aes(x = Quality, y = N50/1000, fill = Quality)) +
  geom_violin(trim = F) +
  theme_cowplot(12) +
  scale_y_log10(label = comma) +
  labs(
    y = "Contig N50 (kilobase)"
  )

# ggsave("/Users/Fiona/Desktop/nmags_summary.png", nmags_summary, device = "png", width = 6, height = 4)

## compare mags and nmags
mags_comp <- mags_filt %>% filter(Sample %in% nmag_samples)
mags_comp$type <- "mag"
nmags$type <- "nmag"

mag_nmag <- rbind(mags_comp, nmags)

labs <- c("Short read", "Long read")
names(labs) <- c("mag", "nmag")

# N50
mag_nmag_N50 <- ggplot(mag_nmag, aes(x = Quality, y = N50/1000, fill = type)) +
  geom_violin(trim = F) +
  theme_cowplot(12) +
  scale_y_log10(label = comma) +
  scale_fill_discrete(labels = c("Short read", "Long read")) +
  labs(
    y = "Contig N50 (kilobase)",
    fill = ""
  )
# ggsave("/Users/Fiona/Desktop/mag_nmag_N50.png", mag_nmag_N50, device = "png", width = 6, height = 4)

# completeness
mag_nmag_completeness <- ggplot(mag_nmag, aes(x = Quality, y = Completeness, fill = type)) +
  geom_violin(trim = F) +
  theme_cowplot(12) +
  # scale_y_log10(label = comma) +
  scale_fill_discrete(labels = c("Short read", "Long read")) +
  labs(
    y = "Completeness (%)",
    fill = ""
  )
# ggsave("/Users/Fiona/Desktop/mag_nmag_completeness.png", mag_nmag_completeness, device = "png", width = 6, height = 4)

# contamination
mag_nmag_contamination <- ggplot(mag_nmag, aes(x = Quality, y = Contamination, fill = type)) +
  geom_violin(trim = T) +
  theme_cowplot(12) +
  # scale_y_log10(label = comma) +
  scale_fill_discrete(labels = c("Short read", "Long read")) +
  scale_y_continuous(trans = "pseudo_log") +
  labs(
    y = "Contamination (%)",
    fill = ""
  )
# ggsave("/Users/Fiona/Desktop/mag_nmag_contamination.png", mag_nmag_contamination, device = "png", width = 6, height = 4)


## taxonomic classifications of dereplicated mags
drep_mags <- read.table("input_final/drep_bins.txt", sep = "\t")
names(drep_mags) <- c("sample", "bin")

drep_mags_info <- mags %>%
  mutate(mag = paste(Sample, Bin)) %>%
  filter(mag %in% paste(drep_mags$sample, drep_mags$bin)) %>%
  select(-mag)

labels <- read.table("input_final/za_labels.tsv", sep = "\t", header = T)

drep_mags_info <- merge(drep_mags_info, labels, by.x = "Sample", by.y = "sample_code")

drep_mags_info <- drep_mags_info %>%
  select(id, Bin, Final.Class, Completeness, Contamination, Strain.heterogeneity, Size.Mb, N50, Coverage, Quality)

drep_mags_info <- drep_mags_info[mixedorder(as.character(drep_mags_info$Bin)), ]
drep_mags_info <- drep_mags_info[mixedorder(as.character(drep_mags_info$id)), ]

names(drep_mags_info) <- c("Sample", "Bin", "Final Class", "Completeness", "Contamination",
                           "Strain heterogeneity", "Size (Mb)", "N50", "Coverage", "Quality")
write.table(drep_mags_info, "final_tables/drep_mags.txt", sep = "\t", row.names = F, quote = F)

