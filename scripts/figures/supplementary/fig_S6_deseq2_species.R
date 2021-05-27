library(cowplot)
library(DESeq2)
library(genefilter)
library(ggpubr)
library(here)
library(reshape2)
library(tidyverse)
library(vegan)

## load data ----
source(here("scripts/load_data.R"))

## differentially abundant species with deseq2 ----
count_data <- za_S
count_data <- count_data[, za_meta$sample]

stopifnot(all(za_meta$sample == names(count_data)))

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = za_meta,
                              design = ~ site)

# estimate size factors
dds <- estimateSizeFactors(dds, type = "poscounts")

# filter normalized counts
idx <- genefilter(counts(dds, normalized = TRUE), pOverA(0.2, 500))
dds <- dds[idx, ]

# deseq2
dds <- DESeq(dds, parallel = T)

res <- results(dds, name = "site_Soweto_vs_Bushbuckridge", alpha = 0.05)

resOrdered <- res %>%
  as_tibble(rownames = "species") %>%
  arrange(pvalue) %>%
  mutate(site = ifelse(log2FoldChange < 0, "Bushbuckridge", "Soweto"))

resOrdered %>%
  filter(abs(log2FoldChange) > 1, padj < 0.05) %>%
  arrange(log2FoldChange)

res_df_filt <- resOrdered %>%
  filter(padj < 0.05, abs(log2FoldChange) > 2) %>%
  arrange(log2FoldChange)

res_df_filt$species <- factor(res_df_filt$species,
                            levels = as.character(res_df_filt$species))

ggplot(res_df_filt, aes(log2FoldChange, species, fill = site)) +
  # geom_point(size = 2, color = "black", pch = 21) +
  geom_bar(stat = "identity") +
  theme_cowplot() +
  theme(axis.text.y = element_text(face = "italic")) +
  scale_fill_manual(values = za_pal) +
  scale_x_continuous(breaks = seq(-7.5, 2.5, 2.5)) +
  labs(
    x = "Log2 Fold Change (Soweto/Bushbuckridge)",
    y = "species",
    fill = ""
  ) +
  background_grid() +
  theme(legend.position = "top",
        legend.justification = "center")

ggsave(here("final_plots/supplementary/figure_S6_deseq_species.png"),
       width = 8.5, height = 11)

res_table_S <- resOrdered %>%
  arrange(site, log2FoldChange) %>%
  mutate(rank = "species")

write.table(res_table_S, here("final_tables/deseq_species.txt"), sep = "\t",
            row.names = F, quote = F)
