library(cowplot)
library(DESeq2)
library(here)
library(reshape2)
library(scales)
library(tidyverse)
library(vegan)

## load data ----
source(here("scripts/load_data.R"))

## plot global abundance of vanish taxa ----

# pull vanish taxa relative abundance
vanish_tax <- data.frame(
  "F" = gsub(".+f__|\\|g__.+", "", vanish_G),
  "G" = gsub(".+g__", "", vanish_G)
)

global_G_VANISH <- global_G_rel[which(rownames(global_G_rel) %in% 
                                               vanish_tax$G), ]

counts_long <- reshape2::melt(global_G_VANISH %>% rownames_to_column(var = "G"),
                              id.vars = "G",
                              variable.name = "sample",
                              value.name = "relab")

counts_long$G <- factor(counts_long$G, levels = rev(unique(counts_long$G)))
counts_long <- merge(counts_long, pheno_global, by = "sample")

# facet by family
counts_long$F <- vanish_tax[match(counts_long$G, vanish_tax$G), "F"]

# plot in decreasing order
levels <- counts_long %>%
  group_by(G) %>%
  summarise(mean_relab = mean(relab)) %>%
  arrange(-mean_relab) %>%
  pull(G)

counts_long$G <- factor(counts_long$G, levels = levels)

counts_long$relab[counts_long$relab == 0] <- min(counts_long$relab[counts_long$relab > 0]) / 2

counts_long$site2 <- factor(counts_long$site2, levels = sites)

ggplot(counts_long, aes(G, relab * 100, fill = site2)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25,
                                             dodge.width = 0.85),
             alpha = 0.75, size = 0.75, color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.75,
               position = position_dodge(width = 0.85)) +
  facet_wrap(~ `F`, scales = "free", nrow = 3) +
  scale_y_log10(labels = label_comma(accuracy = 0.00001)) +
  theme_cowplot(14) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    axis.text.x = element_text(face = "italic"),
    strip.text = element_text(face = "italic"),
    axis.title.y = element_text(size = 12),
    legend.spacing.x = unit(0.1, "cm"),
    plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"),
  ) +
  scale_fill_manual(values = global_pal) +
  labs(
    x = "",
    y = "Relative Abundance (%)",
    fill = ""
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  background_grid(major = "y",
                  minor = "y")

ggsave(here("final_plots/supplementary/figure_S11_VANISH_global.png"),
       width = 9, height = 11, bg = "white")
