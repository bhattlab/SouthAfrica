library(cowplot)
library(DESeq2)
library(genefilter)
library(ggrepel)
library(ggpubr)
library(here)
library(reshape2)
library(scales)
library(tidyverse)
library(vegan)

## load data ----
source(here("scripts/load_data.R"))

## find differential genera across geography with DESeq2
# count_data <- global_G[genefilter(global_G, pOverA(0.1, 1000)), ]
count_data <- global_G

col_data <- pheno_global %>%
  filter(sample %in% names(count_data)) %>%
  column_to_rownames(var = "sample")

count_data <- count_data[, rownames(col_data)]

stopifnot(all(rownames(col_data) == names(count_data)))

col_data$geography <- ifelse(col_data$site %in% c("Tanzania", "Madagascar"), "nonwestern",
                             ifelse(col_data$site == "South Africa", "za", "western"))
col_data$geography <- factor(col_data$geography, levels = c("za", "nonwestern", "western"))

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~ geography)

# estimate size factors
dds <- estimateSizeFactors(dds, type = "poscounts")

# filter normalized counts
idx <- genefilter(counts(dds, normalized = TRUE), pOverA(0.1, 500))
dds <- dds[idx, ]

# deseq2
dds <- DESeq(dds, parallel = T)

dds <- DESeq(dds, parallel = T)
resultsNames(dds) # lists the coefficients

alpha <- 0.1

res_w <- results(dds, name="geography_western_vs_za", alpha = alpha)
res_filt_w <- data.frame(res_w) %>%
  rownames_to_column("feature") %>%
  filter(padj < alpha) %>%
  arrange(-abs(log2FoldChange)) %>%
  mutate(
    comparison = "western",
    pos = log2FoldChange > 0
  )

res_nw <- results(dds, name="geography_nonwestern_vs_za", alpha = alpha)
res_filt_nw <- data.frame(res_nw) %>%
  rownames_to_column("feature") %>%
  filter(padj < alpha) %>%
  arrange(-abs(log2FoldChange)) %>%
  mutate(
    comparison = "nonwestern",
    pos = log2FoldChange > 0
  )

res_all <- rbind(res_filt_w, res_filt_nw)
res_all <- res_all %>%
  filter(abs(log2FoldChange) > 2) %>%
  select(feature, comparison, pos)
res_all <- dcast(res_all, feature ~ comparison)

plot_features <- res_all %>%
  filter(nonwestern == western) %>%
  pull(feature)

## plot
deseq2_counts <- counts(dds, normalized = TRUE)

# add pseudocount
deseq2_counts[deseq2_counts == 0] <- 0.05

counts_long <- melt(data.matrix(deseq2_counts))
names(counts_long) <- c("genus", "sample", "rel_abundance")

counts_long_filt <- counts_long %>%
  filter(genus %in% plot_features)

counts_meta <- merge(counts_long_filt, col_data %>%
                       rownames_to_column("sample"), by = "sample")

counts_meta$geography <- factor(counts_meta$geography,
                                levels = c("nonwestern", "za", "western"))

ggplot(counts_meta, aes(geography, rel_abundance, fill = geography)) +
  geom_jitter(color = "darkgray", width = 0.25, alpha = 0.75, size = 0.75) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  scale_y_log10(label = comma) +
  theme_cowplot() +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  facet_wrap(~ genus, scales = "free_y", ncol = 4, labeller = label_wrap_gen()) +
  scale_fill_discrete(labels = c("Non-western", "South Africa", "Western")) +
  labs(
    x = "",
    y = "DESeq2 normalized count",
    fill = ""
  ) +
  background_grid(major = "y")

ggsave(here("final_plots/supplementary/figure_S11_za_vs_western_nonwestern.png"),
       width = 11, height = 11)
