library(cowplot)
library(DESeq2)
library(here)
library(reshape2)
library(RColorBrewer)
library(tidyverse)
library(vegan)

## load data ----
source(here("scripts/load_data.R"))

## mds ----
bray_dist <- vegdist(t(za_S_css), method = "bray")
mds <- cmdscale(bray_dist, eig = T)
mds_scores <- as.data.frame(scores(mds))
names(mds_scores) <- c("x", "y")
mds_var_perc <- round(mds$eig/sum(mds$eig) * 100, 1)

# merge taxonomic data
relab_P <- data.frame(t(za_P_rel)) * 100
relab_G <- data.frame(t(za_G_rel[1:500, ])) * 100
relab_S <- data.frame(t(za_S_rel[1:500, ])) * 100
mds_taxa <- merge(mds_scores, relab_P, by = "row.names")
# mds_taxa <- merge(mds_scores, relab_G, by = "row.names")
mds_taxa <- merge(mds_taxa, relab_G, by.x = "Row.names", by.y = "row.names")
mds_taxa <- merge(mds_taxa, relab_S, by.x = "Row.names", by.y = "row.names")
mds_taxa <- merge(mds_taxa, za_meta, by.x = "Row.names", by.y = "sample")

# Bacteroides/Prevotella ratio
mds_taxa <- mds_taxa %>%
  mutate(
    Bact_Prev_ratio = log2(Bacteroides / Prevotella),
    Bacteroidetes_Firm_ratio = log2(Bacteroidetes / Firmicutes)
  )

gradient_pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

plot_by <- function(color_by){
  
  ggplot(mds_taxa, aes(x, y, color = !!sym(color_by), shape = site)) +
    geom_point(size = 3, alpha = 0.85) +
    # scale_color_gradientn(colours = gradient_pal(100),
    #                       limits = c(min(mds_taxa[, color_by]),
    #                                  max(mds_taxa[, color_by]))) +
    viridis::scale_color_viridis() +
    labs(
      x = paste("MDS1 (", mds_var_perc[1], "%)",sep=""),
      y = paste("MDS2 (", mds_var_perc[2], "%)",sep=""),
      shape = "Site"
    ) +
    theme_cowplot()
}

mds_Bact_Prev_ratio <- plot_by("Bact_Prev_ratio")

mds_Bact_Prev_ratio +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(margin = margin(r = 5, unit = "pt"))
  ) +
  labs(
    color = "log2(Bacteroides/\nPrevotella)",
    shape = "      Site"
  ) +
  guides(color = guide_colorbar(title.position = "left", title.vjust = 0.8),
         shape = guide_legend(title.position = "left",
                              title.vjust = 0.8, ncol = 1)) +
  coord_fixed() +
  background_grid()

ggsave(here("final_plots/supplementary/figure_S4_bacteroides_prevotella.png"),
       width = 6, height = 6)
