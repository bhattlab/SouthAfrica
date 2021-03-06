library(cowplot)
library(genefilter)
library(ggpubr)
library(harrietr)
library(here)
library(reshape2)
library(scales)
library(tidyverse)

# load data ----
load(here("RData/metadata.RData"))
load(here("RData/za_data.RData"))
load(here("RData/global_data.RData"))

# plot distributions
vanish_tax <- data.frame(
  "F" = gsub(".+f__|\\|g__.+", "", vanish_G),
  "G" = gsub(".+g__", "", vanish_G)
)

za_G_VANISH <- za_G_rel[which(rownames(za_G_rel) %in% as.character(vanish_tax$G)), ]
za_F_VANISH <- za_F_rel[which(rownames(za_F_rel) %in% as.character(vanish_tax$F)), ]

za_VANISH_G_long <- reshape2::melt(za_G_VANISH %>% rownames_to_column(var = "G"),
                                   id.vars = "G",
                                   variable.name = "sample",
                                   value.name = "relab")

za_VANISH_G_long <- merge(za_VANISH_G_long, za_meta, by = "sample")

succin_long <- za_VANISH_G_long %>%
  filter(G %in% c("Succinatimonas", "Succinivibrio", "Treponema") &
           site == "Bushbuckridge")

density_plot <- ggplot(succin_long, aes(x = relab * 100)) +
  geom_density(fill = za_pal[1]) +
  scale_x_log10(label = comma) +
  facet_wrap(~ G, scales = "free_y", ncol = 1) +
  theme_cowplot() +
  labs(
    x = "Relative abundance (%)",
    y = "Density"
  ) +
  theme(strip.text = element_text(face = "italic")) +
  background_grid()

# correlate VANISH taxa ----
vanish_cor <- function(counts){
  # correlation
  vanish_cor_G <- cor(t(counts), method = "spearman")
  
  # cluster to order rows
  hc <- hclust(dist(vanish_cor_G))
  row_order <- hc$order
  vanish_cor_G <- vanish_cor_G[row_order, row_order]
  
  # plot
  vanish_cor_G_long <- melt(vanish_cor_G, value.name = "cor")
  
  ggplot(vanish_cor_G_long, aes(Var2, Var1, fill = cor)) +
    geom_tile() +
    geom_text(aes(label = round(cor, 2))) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1)) +
    theme_cowplot() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title = element_blank()
    ) +
    labs(
      fill = "Spearman's\nrho"
    )
}

plot_data <- za_G_VANISH
rownames(plot_data) <- paste0(rownames(plot_data), "\n(",
                              vanish_tax[match(rownames(plot_data),
                                               vanish_tax$G), "F"], ")")

plot_G <- vanish_cor(plot_data)
plot_F <- vanish_cor(za_F_VANISH)

plot_grid(
  plot_grid(density_plot, plot_F, labels = c("a", "b")),
  plot_G,
  ncol = 1,
  labels = c("", "c"),
  rel_heights = c(0.27, 0.63)
)

ggsave(here("final_plots/supplementary/figure_S2_VANISH_cor.png"),
       width = 9, height = 11)
ggsave(here("final_plots/pdf/supp/figure_S2_VANISH_cor.pdf"),
       width = 9, height = 11)
