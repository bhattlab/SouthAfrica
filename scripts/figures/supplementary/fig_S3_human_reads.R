library(cowplot)
library(ggpubr)
library(gtools)
library(here)
library(RColorBrewer)
library(reshape2)
library(tidyverse)

## load data ----
source(here("scripts/load_data.R"))

labs <- read.table(here("input_final/za_labels.tsv"), sep = "\t", header = T)

## readcounts ----

## read PAIRS
za_readcounts <- read.table(here("input_final", "readcounts",
                                 "za_readcounts_preprocessing.txt"),
                            sep = '\t', header = T)

# read PAIRS -- multiply by 2
# filter out extra samples
za_counts <- za_readcounts %>%
  filter(Sample %in% za_meta$sample) %>%
  mutate(raw_reads = raw_reads * 2,
         trimmed_reads = trimmed_reads * 2,
         dedup_reads = dedup_reads * 2,
         host_removed_reads = host_removed_reads * 2,
         orphan_reads = orphan_reads * 2
         ) %>%
  dplyr::select(-starts_with("orphan"), -ends_with("frac"))

# za_counts <- za_readcounts[, c(1:3, 5, 7)]

## stats for intro
za_counts %>%
  mutate(
    Gb_raw = raw_reads * 150 / 1e9,
    Gb_preproc = host_removed_reads * 150 / 1e9
  ) %>%
  summarise(
    median_reads_raw = median(raw_reads) / 1e6,
    median_reads_preproc = median(host_removed_reads) / 1e6,
    range_reads_raw = paste(range(raw_reads) / 1e6, collapse = " "),
    range_reads_preproc = paste(range(host_removed_reads) / 1e6, collapse = " "),
    median_gb_raw = median(Gb_raw),
    median_gb_raw = median(Gb_raw),
    mean_gb_preproc = mean(Gb_preproc)
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

colnames(za_counts) <- c("Sample", "Raw reads", "Deduplicated reads",
                         "Trimmed reads", "Non-human reads")

counts_long <- melt(za_counts, id.vars = "Sample",
                    variable.name = "step", value.name = "reads")
counts_long$reads_m <- (counts_long$reads / 1e6)

# plot readcounts
rc <- ggplot(counts_long, aes(x=reads_m, fill=step)) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  facet_grid(step~., scales = "free_y") +
  labs(
    title = "Preprocessing Readcounts\n",
    x = "\nReads (M)",
    y = "Count\n",
    fill = ""
  ) +
  theme_cowplot(14) +
  background_grid()

# ggsave(here("final_plots/misc/readcounts_preprocessing.png"), rc,
#        height = 12, width = 12)

# save readcounts to supp file
za_counts <- za_counts %>%
  mutate(
    label = labs[match(za_counts$Sample, labs$sample), "label_abbrev"],
    Sample = label
  ) %>%
  dplyr::select(Sample, ends_with("reads"))

za_counts <- za_counts[mixedorder(za_counts$Sample), ]

## do BBR/SWT differ in % host reads?
z <- za_counts %>%
  mutate(
    human_removed = `Deduplicated reads` - `Non-human reads`,
    human_rm_frac = human_removed / `Deduplicated reads`,
    site = gsub("\\d+$", "", Sample)
  )

wilcox.test(human_rm_frac ~ site, z)

z %>%
  group_by(site) %>%
  summarise(mean_human_frac = mean(human_rm_frac))

a <- ggplot(z, aes(human_rm_frac * 100, fill = site)) +
  geom_histogram(color = "white") +
  scale_fill_manual(values = za_pal) +
  theme_cowplot() +
  theme(
    legend.position = "top",
    legend.justification = "center"
  ) +
  labs(
    x = "Human reads removed after de-duplication (%)",
    fill = "Site",
    y = "Sample count"
  ) +
  background_grid()

b <- ggplot(z, aes(site, human_rm_frac * 100, fill = site)) +
  geom_jitter(color = "darkgray", alpha = 0.75, width = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  scale_fill_manual(values = za_pal) +
  theme_cowplot() +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5) +
  labs(
    x = "",
    y = "Human reads removed\nafter de-duplication (%)"
  ) +
  background_grid(major = "y")

plot_grid(a, b, labels = c("A", "B"), ncol = 1)

ggsave(here("final_plots/supplementary/figure_S3_human_reads.png"),
       width = 5, height = 9)

za_counts <- za_counts %>%
  mutate("Human reads removed after de-duplication" =
           `Deduplicated reads` - `Non-human reads`)

write.table(za_counts, here("final_tables/table_s3_readcounts.txt"),
            sep = "\t", row.names = F, quote = F)
