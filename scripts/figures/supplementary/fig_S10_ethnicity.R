library(cowplot)
library(ggpubr)
library(here)
library(MASS)
library(RColorBrewer)
library(reshape2)
library(tidyverse)
library(vegan)

# load data ----
load(here("RData/metadata.RData"))
load(here("RData/za_data.RData"))

# mds plot ----
# ethnicity info
ethnicity <- read.table(here("input_final/pheno/ethnicity.txt"), sep = "\t",
                        header = T)

# nmds clustering by ethincity
pheno <- merge(za_pheno, awi_eth, by.x = "study_id", by.y = "Microbiome_link")

# mds
vare_dis <- vegdist(t(za_S_rare_css), method = "bray")
mds <- cmdscale(vare_dis, eig = TRUE, x.ret = TRUE)
mds_values <- mds$points

# variation per axis
mds_var_per <- round(mds$eig/sum(mds$eig) * 100, 1)

# Plot
mds_data <- data.frame(sample = rownames(mds_values),
                       x = mds_values[,1],
                       y = mds_values[,2])

# merge pheno data
mds_meta <- merge(mds_data, pheno, by = "sample", all.x = T)
mds_meta$ethnicity[is.na(mds_meta$ethnicity)] <- "Unknown"

mds_meta$ethnicity <- factor(
  mds_meta$ethnicity,
  levels = c(sort(unique(mds_meta$ethnicity[-which(mds_meta$ethnicity
                                                        %in% c("Other", "Unknown"))])),
             "Other", "Unknown"))

ggplot(mds_meta, aes(x, y, color = ethnicity)) +
  geom_point(size = 2, alpha = 0.75) +
  scale_color_manual(values = c(brewer.pal(9, "Set1"), "#4A4A4A")) +
  labs(x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
       y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
       color="Ethnicity") +
  theme_cowplot() +
  coord_fixed() +
  background_grid()

ggsave(here("final_plots/supplementary/figure_S10_za_ethnicity_mds.png"),
       height = 4, width = 6, bg = "white")
ggsave(here("final_plots/pdf/supp/figure_S10_za_ethnicity_mds.pdf"),
       height = 4, width = 6, bg = "white")
