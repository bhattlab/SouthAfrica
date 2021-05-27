library(cowplot)
library(DESeq2)
library(genefilter)
library(ggplot2)
library(here)
library(RColorBrewer)
library(tidyverse)
library(vegan)

# read metadata ---
za_meta <- readRDS(here("rds/za_meta.rds"))
za_pheno <- readRDS(here("rds/za_pheno.rds"))
za_pheno$site_abbrev <- ifelse(za_pheno$site == "Bushbuckridge", "BBR", "SWT")

# read gtdb data ---
for (lvl in c("S", "G", "P")){
  bracken <- read.table(here(paste0("input_final/gtdb_r95/bracken_", lvl, ".txt")),
                        sep = "\t", header = T)
  rownames(bracken) <- gsub("s__", "", rownames(bracken))
  
  # keep relevant samples
  bracken <- bracken[, za_pheno$sample]
  
  # rm zero features
  bracken <- bracken[rowSums(bracken) > 0, ]
  
  bracken_rel <- sweep(bracken, 2, colSums(bracken), FUN = "/")
  
  assign(paste0("bracken_", lvl), bracken)
  assign(paste0("bracken_", lvl, "_rel"), bracken_rel)
}

## Prevotella plot ---- 
prev <- bracken_S_rel[grep("Prevotella ", rownames(bracken_S_rel)), ]
prev <- prev[genefilter(prev, pOverA(p = 0.3, A = 0)), ]

prev_long <- prev %>%
  head(10) %>%
  rownames_to_column("feature") %>%
  pivot_longer(cols = -feature, names_to = "sample", values_to = "relab") %>%
  left_join(za_meta, by = "sample")

ggplot(prev_long, aes(sample, relab, fill = feature)) +
  geom_bar(stat = "identity") +
  theme_cowplot() +
  facet_grid(~site, space = "free", scales = "free_x")

wc <- prev %>%
  rownames_to_column("feature") %>%
  pivot_longer(cols = -feature, names_to = "sample", values_to = "relab") %>%
  left_join(za_meta, by = "sample") %>%
  group_by(feature) %>%
  summarize(pvalue = wilcox.test(relab ~ site)$p.value,
            mean_swt = mean(relab[site == "Soweto"]),
            mean_bbr = mean(relab[site == "Bushbuckridge"]),
            incBBR = mean_bbr > mean_swt) %>%
  ungroup() %>%
  mutate(padj = p.adjust(pvalue, method = "fdr"))

wc_signif <- wc %>%
  filter(padj < 0.05)

prev_plot <- prev %>%
  rownames_to_column("feature") %>%
  pivot_longer(cols = -feature, names_to = "sample", values_to = "relab") %>%
  left_join(za_meta, by = "sample") %>%
  filter(feature %in% wc_signif$feature) %>%
  mutate(relab = ifelse(relab == 0, 1e-5, relab))

ggplot(prev_plot, aes(site, relab)) +
  geom_boxplot() +
  facet_wrap(~feature) +
  theme_cowplot() +
  background_grid(major = "y") +
  scale_y_log10() +
  ggpubr::stat_compare_means(method = "wilcox")

