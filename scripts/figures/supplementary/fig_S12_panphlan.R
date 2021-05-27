library(cowplot)
library(ggpubr)
library(here)
library(MASS)
library(reshape2)
library(tidyverse)
library(vegan)

## load data ----
source(here("scripts/load_data.R"))

### panphlan mds function ----

set.seed(1)

# function to ordinate and plot panphlan results table
plot_panphlan <- function(panphlan_f, plot_title = "", legend = F){
  
  # read panphlan pangenome matrix
  panphlan <- read.table(panphlan_f, sep = "\t", header = T, row.names = 1, quote = "")
  colnames(panphlan) <- gsub('.tsv', '', colnames(panphlan))
  
  # filter for relevant samples
  panphlan <- panphlan[, which(names(panphlan) %in% za_pheno$sample |
                                 grepl("REF", names(panphlan)))]
  
  # calculate jaccard distance
  panphlan_t <- t(panphlan)
  dist_matrix <- vegdist(panphlan_t, method = "jaccard", binary = T)
  # colnames(dist_matrix) <- colnames(panphlan)
  
  # calculate pcoa/mds
  mds <- cmdscale(dist_matrix, eig = TRUE, x.ret = TRUE)
  scores <- scores(mds)
  mds_df <- data.frame(sample = rownames(scores),
                       X1 = scores[,1],
                       X2 = scores[,2])
  
  # axis percentage of variation
  mds_var_per <- round(mds$eig/sum(mds$eig) * 100, 1)
  
  # merge mds and pheno
  mds_pheno <- merge(mds_df, za_pheno, by = "sample", all.x = T)
  mds_pheno$site[is.na(mds_pheno$site)] <- "Reference"
  
  mds_pheno$site <- factor(mds_pheno$site,
                           levels = c("Bushbuckridge", "Soweto", "Reference"))
  
  # no labels
  g <- ggplot(mds_pheno, aes(X1, X2, color = site)) +
    geom_point(size = 1.5, alpha = 0.85) +
    labs(x = paste0("MDS 1 (", mds_var_per[1], "%)"),
         y = paste0("MDS 2 (", mds_var_per[2], "%)"),
         color = "",
         title = plot_title) +
    theme_cowplot() +
    background_grid() +
    scale_color_manual(values = c(za_pal, "black")) +
    theme(legend.position = "top",
          legend.justification = "center",
          plot.title = element_text(face = "italic"))
    # coord_fixed()
  
  if (!legend) g <- g + theme(legend.position = "none")
  
  return(g)
}

## panphlan plots ----
sp <- c("Prevotella_copri",
        "Faecalibacterium_prausnitzii",
        "Bacteroides_vulgatus",
        "Eubacterium_siraeum",
        "Alistipes_putredinis",
        "Butyrivibrio_crossotus")

plots <- list()

for (org in sp){
  infile <- here(paste0("input_final/panphlan_v2/", org,
                        "_profile_sensitive.txt"))
  
  p <- plot_panphlan(infile, plot_title = gsub("_", " ", org), legend = F)
  
  plots[[org]] <- p
}

legend <- get_legend(p + theme(legend.position = "top"))

a1 <- plot_grid(plotlist = plots, ncol = 3, align = "hv", axis = "lrbt")

a <- plot_grid(a1, legend, ncol = 1, rel_heights = c(0.95, 0.05))

## adonis
adonis_res <- data.frame()

# sp <- c("Prevotella_copri",
#         "Faecalibacterium_prausnitzii",
#         "Bacteroides_vulgatus",
#         "Eubacterium_siraeum",
#         "Alistipes_putredinis",
#         "Butyrivibrio_crossotus")

for (org in sp){
  panphlan_f <- here(paste0("input_final/panphlan_v2/", org,
                            "_profile_sensitive.txt"))
  
  # read panphlan pangenome matrix
  panphlan <- read.table(panphlan_f, sep = "\t", header = T, row.names = 1, quote = "")
  colnames(panphlan) <- gsub('.tsv', '', colnames(panphlan))
  
  # metadata
  meta <- za_pheno %>%
    filter(sample %in% names(panphlan))
  
  # filter for relevant samples
  panphlan <- panphlan[, meta$sample]
  
  # calculate jaccard distance
  panphlan_t <- t(panphlan)
  dist_matrix <- vegdist(panphlan_t, method = "jaccard", binary = T)
  
  ad <- adonis(dist_matrix~site, meta)
  
  adonis_res <- rbind(adonis_res,
                      data.frame("Species" = org,
                                 "R2" = ad$aov.tab$R2[1],
                                 "pvalue" = ad$aov.tab$`Pr(>F)`[1]))
}

b <- adonis_res %>%
  mutate(R2 = round(R2, 3),
         FDR = p.adjust(pvalue, method = "fdr"),
         Species = gsub("_", " ", Species)) %>%
  ggtexttable(rows = NULL, cols = c("Species", "R2", "Pr(>F)", "FDR"), theme =
                ttheme("classic", base_size = 12, padding = unit(c(4, 4), "mm")))

plot_grid(a, b, ncol = 1, labels = c("A", "B"), rel_heights = c(0.7, 0.3))

ggsave(here("final_plots/supplementary/figure_S12_panphlan.png"),
       width = 10, height = 8)
