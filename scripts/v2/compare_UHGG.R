library(cowplot)
library(ggplot2)
library(here)
library(tidyverse)

## params ----
global_pal <- c("#E3211C", "#F89897", "#6A3D9A", "#CAB2D6", "#1F78B4", "#A5CEE3")
za_pal <- global_pal[3:4]

## read input data ----
## sample metadata
za_meta <- readRDS(here("rds/za_meta.rds"))

## genomesearch output
gsearch <- read.table(here("input_final/uhgg_ANI_compare/perident_table_uhgg_formatted.tsv"),
                      sep = "\t", header = T,
                      col.names = c("Sample", "Bin", "gsearch.ANI",
                                    "gsearch.Phylum", "gsearch.Species"))

## fastani output
fastani <- read.table(here("input_final/uhgg_ANI_compare/fastani_top.txt"),
                      sep = "\t", header = F,
                      col.names = c("Sample", "Bin","UHGG_genome", "fastani.ANI",
                                    "fragment_mappings", "query_fragments", "AF"))

# take top match per mag
# some mags missing due to low identity
fastani_top <- fastani %>%
  arrange(Sample, Bin, -fastani.ANI) %>%
  group_by(Sample, Bin) %>%
  slice(1)

## gtdbtk classifications
# some mags are missing -- may be archaeal?
gtdbtk <- read.table(here("input_final/mags/gtdbtk.bac120.summary.tsv"),
                     sep = "\t", header = T)
gtdbtk <- gtdbtk %>%
  separate(user_genome, into = c("Sample", "Bin"), sep = "_", extra = "merge")

## dereplication lists
# za 95% ANI
drep_za <- read.table(here("input_final/mags/za_drep95_reps.txt"), sep = "\t",
                      header = F, col.names = c("Sample", "Bin"))

# za plus uhgg reps, 95% ANI
drep_za_uhgg <- read.table(here("input_final/mags/za_uhggreps_drep95_reps.txt"),
                           sep = "\t", header = F, col.names = c("Sample", "Bin"))

## binning dastool pipeline output
mags <- read.table(here("input_final/mags/binning_table_all_simple.tsv"),
                   sep = "\t", header = T, comment.char = "")

## format mag info ----
# add genomesearch and fastani comparisons
# only consider MQ and higher
mags <- mags %>%
  filter(Bin != "unbinned") %>%
  mutate(Quality = ifelse(Completeness > 90 & Contamination < 5 & tRNA >= 18 & rna.16S > 0 & rna.23S > 0 & rna.5S > 0,
                          "High-quality",
                          ifelse(Completeness >= 50 & Contamination <10, "Medium-quality", "Low-quality")),
         Quality = factor(Quality, levels = c("Low-quality", "Medium-quality", "High-quality"))) %>%
  left_join(za_meta, by = c("Sample" = "sample")) %>%
  left_join(gsearch, by = c("Sample", "Bin")) %>%
  left_join(fastani_top, by = c("Sample", "Bin")) %>%
  left_join(gtdbtk[, c("Sample", "Bin", "classification")], by = c("Sample", "Bin"))

mags_mqplus <- mags %>%
  filter(Quality != "Low-quality")

## ANI plots ----
## gsearch
ggplot(mags_mqplus, aes(site, gsearch.ANI, fill = site)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.7), alpha = 0.75, color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  # ylim(c(0, 100)) +
  ylab("Average Amino Acid Identity to UHGG") + 
  scale_fill_manual(values = za_pal) +
  theme_cowplot() +
  theme(legend.position = "none") +
  background_grid(major = "y")

# ggsave("/Users/tamburif/Desktop/gsearch_MAGs.png", width = 5, height = 5)

## fastANI
ggplot(mags_mqplus, aes(site, fastani.ANI, fill = site)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.7), alpha = 0.75, color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  # ylim(c(0, 100)) +
  ylab("Average Amino Acid Identity to UHGG") + 
  scale_fill_manual(values = za_pal) +
  theme_cowplot() +
  theme(legend.position = "none") +
  background_grid(major = "y")

ggsave("/Users/tamburif/Desktop/fastANI_MAGs.png", width = 5, height = 5)

## fastANI/genomesearch correlation
ggplot(mags_mqplus, aes(gsearch.ANI, fastani.ANI, color = site, size = Quality)) +
  geom_point(alpha = 0.75) +
  scale_size_manual(values = c(2, 3)) +
  scale_color_manual(values = za_pal) +
  theme_cowplot() +
  background_grid()

# top 95%
m <- mags_mqplus %>%
  filter(gsearch.ANI > 95,
         fastani.ANI > 95)

ggplot(m, aes(gsearch.ANI, fastani.ANI, color = site, size = Quality)) +
  geom_point(alpha = 0.75) +
  scale_size_manual(values = c(2, 3)) +
  scale_color_manual(values = za_pal) +
  theme_cowplot() +
  background_grid() +
  scale_x_log10() +
  scale_y_log10() +
  ylim(c(95, 100)) +
  xlim(c(95, 100)) +
  geom_abline(slope = 1, intercept = 0) +
  labs(color = "Site") +
  coord_fixed()

ggsave("/Users/tamburif/Desktop/compare_genomesearch_fastANI.png", width = 7, height = 5)

ggplot(mags, aes(site, fastani.ANI, fill = site)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.7), alpha = 0.75, color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  # ylim(c(0, 100)) +
  ylab("Average Amino Acid Identity to UHGG") + 
  scale_fill_manual(values = za_pal) +
  theme_cowplot() +
  theme(legend.position = "none") +
  background_grid(major = "y")

## novel mags ----
## what are the mags with low identity to uhgg?
novel_mags <- mags_mqplus %>%
  filter(is.na(fastani.ANI) | fastani.ANI < 95) %>%
  separate(classification, into = c("D", "P", "C", "O", "F", "G", "S"),
           sep = ";", remove = F) %>%
  mutate(za_rep = paste(Sample, Bin) %in% paste(drep_za$Sample, drep_za$Bin),
         uhgg_rep = paste(Sample, Bin) %in% paste(drep_za_uhgg$Sample, drep_za_uhgg$Bin))

# by phylum
novel_mags %>%
  group_by(P) %>%
  tally() %>%
  arrange(-n)

# by genus
novel_mags %>%
  group_by(C) %>%
  tally() %>%
  arrange(-n)

# za reps
novel_mags %>%
  filter(za_rep)

# za plus uhgg reps
novel_mags %>%
  filter(uhgg_rep)