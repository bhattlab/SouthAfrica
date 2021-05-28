library(cowplot)
library(ggpubr)
library(here)
library(MASS)
library(RColorBrewer)
library(reshape2)
library(tidyverse)
library(vegan)

## load data ----
source(here("scripts/load_data.R"))

## mds plot ----
## ethnicity info
awi_data_f <- here("input_final/pheno/awi-phase1.csv")
awi_data <- read.csv(awi_data_f)

awi_ethnicity <- data.frame(
  "ethnicity_name" = c("Zulu", "Xhosa", "Ndebele", "Sotho", "Venda", "Tsonga",
                       "Tswana", "BaPedi", "Zimbabwean", "Other", "Unknown",
                       "Swati", "Other"),
  "ethnicity" = c(1:12, 40))

awi_eth <- awi_data[, c("Microbiome_link", "study_id", "site", "ethnicity")]
awi_eth$ethnicity_name <- awi_ethnicity[match(awi_eth$ethnicity,
                                              awi_ethnicity$ethnicity), "ethnicity_name"]

awi_eth$site_name <- ifelse(awi_eth$site == "6", "Soweto", "Bushbuckridge")

eth_counts <- plyr::count(awi_eth, c("site_name", "ethnicity_name"))

## nmds clustering by ethincity
pheno <- merge(za_pheno, awi_eth, by.x = "study_id", by.y = "Microbiome_link")

# mds
vare_dis <- vegdist(t(za_S_css), method = "bray")
mds <- cmdscale(vare_dis, eig = TRUE, x.ret = TRUE)
mds_values <- mds$points

# variation per axis
mds_var_per <- round(mds$eig/sum(mds$eig) * 100, 1)

## Plot
mds_data <- data.frame(sample = rownames(mds_values),
                       x = mds_values[,1],
                       y = mds_values[,2])

# merge pheno data
mds_meta <- merge(mds_data, pheno, by = "sample", all.x = T)
mds_meta$ethnicity_name[is.na(mds_meta$ethnicity_name)] <- "Unknown"

mds_meta$ethnicity_name <- factor(
  mds_meta$ethnicity_name,
  levels = c(sort(unique(mds_meta$ethnicity_name[-which(mds_meta$ethnicity_name
                                                        %in% c("Other", "Unknown"))])),
             "Other", "Unknown"))

ggplot(mds_meta, aes(x, y, color = ethnicity_name)) +
  geom_point(size = 2, alpha = 0.75) +
  scale_color_manual(values = c(brewer.pal(9, "Set1"), "#4A4A4A")) +
  labs(x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
       y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
       color="Ethnicity") +
  theme_cowplot() +
  coord_fixed() +
  background_grid()

ggsave(here("final_plots/supplementary/figure_S9_za_ethnicity_mds.png"),
       height = 4, width = 5)
