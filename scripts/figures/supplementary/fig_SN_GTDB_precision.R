library(cowplot)
library(genefilter)
library(ggpubr)
library(here)
library(MASS)
library(RColorBrewer)
library(reshape2)
library(scales)
library(tidyverse)
library(vegan)

## load data ----
source(here("scripts/load_data.R"))

## read gtdb data ---
for (lvl in c("S", "G", "P")){
  bracken <- read.table(here(paste0("input_final/kraken/gtdb_r95/bracken_", lvl, ".txt")),
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


# cts <- bracken_S_rel[grepl("Prevotella|Faecalibacterium", rownames(bracken_S_rel)), ]
# cts <- cts[genefilter(cts, pOverA(0.1, 0.001)), ]
# 
# cts_long <- cts %>%
#   rownames_to_column("species") %>%
#   pivot_longer(cols = -species, names_to = "sample")

## Prevotella plot ---- 
prev <- bracken_S_rel[grepl("Prevotella |Faecalibacterium", rownames(bracken_S_rel)), ]
prev <- prev[genefilter(prev, pOverA(p = 0.2, A = 0.001)), ]

prev_long <- prev %>%
  # head(10) %>%
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
  filter(padj < 0.01)

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



# bray_dist <- vegdist(t(cts), method = "bray")
# 
# mds <- cmdscale(bray_dist, eig = TRUE, x.ret = TRUE)
# 
# mds_values <- mds$points
# 
# # calculate percentage of variation that each MDS axis accounts for
# mds_var_per <- round(mds$eig/sum(mds$eig) * 100, 1)
# 
# # plot
# mds_data <- data.frame(sample = rownames(mds_values),
#                        x = mds_values[,1],
#                        y = mds_values[,2])
# 
# # merge pheno data
# mds_meta <- merge(mds_data, za_meta, by = "sample")
# 
# mds_plot <- ggplot(mds_meta, aes(x, y, color = site)) +
#   geom_point(size = 2, alpha = 0.85) +
#   scale_color_manual(values = za_pal) +
#   labs(x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
#        y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
#        color = "") +
#   theme_cowplot() +
#   background_grid() +
#   coord_fixed() +
#   theme(legend.position = "top",
#         legend.justification = "center")
