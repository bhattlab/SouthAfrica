library(cowplot)
library(DESeq2)
library(genefilter)
library(ggpubr)
library(here)
library(reshape2)
library(tidyverse)

## load data ----
source(here("scripts/load_data.R"))

## plot top taxa by mean relative abundance ----
top_plot <- function(counts, n = 10){
  counts_long <- reshape2::melt(counts[1:n, ] %>% rownames_to_column(var = "taxon"), variable.name = "sample", value.name = "relab")
  
  counts_long$taxon <- factor(counts_long$taxon, levels = rev(unique(counts_long$taxon)))
  counts_long <- merge(counts_long, za_meta, by = "sample")
  
  counts_long$site <- factor(counts_long$site, levels = c("Soweto", "Bushbuckridge"))
  
  ggplot(counts_long, aes(taxon, relab * 100, fill = site)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.75, color = "darkgray", size = 1.25) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    theme_cowplot(12) +
    theme(legend.position = "top",
      legend.justification = "center") +
    scale_fill_manual(values = za_pal, breaks = c("Bushbuckridge", "Soweto")) +
    labs(x = "",
         y = "Relative abundance (%)",
         color = "",
         fill = "Site") +
    coord_flip() +
    background_grid(major = "x", minor = "x")
}

## top species plot
top_S_bracken <- top_plot(za_S_rel)

## top genus plot
g <- za_G_rel
# rownames(g) <- gsub("miscellaneous", "misc", rownames(g))
top_G_bracken <- top_plot(g)

plot_grid(
  top_S_bracken,
  top_G_bracken + theme(legend.position = "none"),
  ncol = 1,
  labels = c("A", "B"),
  align = "hv",
  axis = "l"
)

ggsave(here("final_plots/supplementary/figure_S1_top_taxa.png"),
       width = 8, height = 10)
