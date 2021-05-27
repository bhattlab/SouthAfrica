library(cowplot)
library(harrietr)
library(here)
library(reshape2)
library(tidyverse)

## load data ----
source(here("scripts_v2/load_data.R"))

## estimate crAssphage prevalence in population ----
# bracken reports read PAIRS -- divide threshold by 2
crass <- "crAss-like viruses"

n <- 650 # ~1X coverage of a 97kb genome by 150 bp reads
length(which(za_G[crass, ] >= n / 2))

za_meta %>%
  filter(sample %in% names(za_G[, which(za_G[crass, ] >= n / 2)])) %>%
  group_by(site) %>%
  tally()

tbl <- za_meta %>%
  mutate(crass_TRUE = sample %in% names(za_G[, which(za_G[crass, ] >= n / 2)])) %>%
  group_by(site) %>%
  summarize(crass_n = sum(crass_TRUE),
            crass_perc = sum(crass_TRUE)/n())

t <- table(tbl$site, tbl$crass_TRUE)
fisher.test(t)

# cor with other taxa
za_S_rel_filt <- za_S_rel[which(rowSums(za_S_rel) > 0), ]
za_S_rel_filt <- za_S_rel_filt[genefilter(za_S_rel_filt, pOverA(0.1, 0.0001)), ]

cor_long <- melt_dist(cor(t(za_S_rel_filt), method = "spearman"))

cor_long <- cor_long %>%
  filter(grepl("Guerin_crAss", iso1) | grepl("Guerin_crAss", iso1)) %>%
  arrange(-abs(dist)) %>%
  filter(abs(dist) > 0.4)

# plot heatmap
keep <- c(cor_long$iso1, cor_long$iso2)

cor_long_heatmap <- melt(cor(t(za_S_rel_filt), method = "spearman"), value.name = "cor")
cor_long_heatmap <- cor_long_heatmap %>%
  filter(Var1 %in% keep & Var2 %in% keep)

ggplot(cor_long_heatmap, aes(Var2, Var1, fill = cor)) +
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