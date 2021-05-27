library(ggplot2)
library(vegan)
library(dplyr)
library(ggrepel)
library(reshape2)
library(grid)
library(MASS)
library(gridExtra)
library(cowplot)

# read metadata ---
za_meta <- readRDS("rds/za_meta.rds")
za_pheno <- readRDS("rds/za_pheno.rds")
za_pheno$site_abbrev <- ifelse(za_pheno$site == "Bushbuckridge", "BBR", "SWT")

### panphlan nmds ----

panphlan_f <- "input_final/panphlan_v2/Bacteroides_vulgatus_profile_sensitive.txt"

# function to ordinate and plot panphlan results table
plot_panphlan <- function(panphlan_f){
  
  # read panphlan pangenome matrix
  panphlan <- read.table(panphlan_f, sep = "\t", header = T, row.names = 1, quote = "")
  colnames(panphlan) <- gsub('.tsv', '', colnames(panphlan))
  
  # filter for relevant samples
  panphlan <- panphlan[, which(names(panphlan) %in% za_pheno$sample | grepl("REF", names(panphlan)))]
  
  # calculate jaccard distance
  panphlan_t <- t(panphlan)
  dist_matrix <- vegdist(panphlan_t, method = "jaccard", binary = T)
  # colnames(dist_matrix) <- colnames(panphlan)
  
  # nmds ordinate
  vare_mds0 <- isoMDS(dist_matrix)
  mds_df <- data.frame(vare_mds0$points)
  mds_df$sample <- row.names(mds_df)
  
  ## calculate pcoa/mds
  # mds <- cmdscale(dist_matrix)
  # 
  # mds_df <- data.frame(sample = rownames(mds),
  #                      X1 = mds[,1],
  #                      X2 = mds[,2])
  
  # merge mds and pheno
  mds_pheno <- merge(mds_df, za_pheno, by = "sample", all.x = T)
  mds_pheno$site[is.na(mds_pheno$site)] <- "Reference"
  
  # no labels
  ggplot(mds_pheno, aes(X1, X2, color = site)) +
    geom_point(size = 2, alpha = 0.75) +
    labs(x = "NMDS 1",
         y = "NMDS 2",
         color = "") +
    theme_cowplot() +
    background_grid() +
    theme(
      plot.title = element_text(face = "italic", size = 14),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      legend.position = "bottom"
    )
}

# npcopri + stat_ellipse(type='t', size = 1, show.legend = F)
# plot_panphlan("/Users/Fiona/scg4_fiona/za/06_panphlan/panphlan_npcopri_pangenome_060419.tsv")

# b uniformis
# b_uniformis <- plot_panphlan("input_final/panphlan/panphlan_buniformis16_pangenome.tsv") + ggtitle("Bacteroides uniformis")
b_uniformis <- plot_panphlan("input_final/panphlan/buniformis16_pangenome.tsv") + ggtitle("Bacteroides uniformis")

# b vulgatus
# b_vulgatus <- plot_panphlan("input_final/panphlan/panphlan_bvulgatus16_pangenome.tsv") + ggtitle("Bacteroides vulgatus")
b_vulgatus <- plot_panphlan("input_final/panphlan/bvulgatus16_pangenome.tsv") + ggtitle("Bacteroides vulgatus")

# npcopri
# npcopri <- plot_panphlan("input_final/panphlan/global_panphlan_npcopri_pangenome.tsv") + ggtitle("Prevotella copri")
npcopri <- plot_panphlan("input_final/panphlan/npcopri_pangenome.tsv") + ggtitle("Prevotella copri")

# p copri
p_copri <- plot_panphlan("input_final/panphlan/pcopri16_pangenome.tsv") + ggtitle("Prevotella copri")

# fprausnitzii
# f_prausnitzii <- plot_panphlan("input_final/panphlan/panphlan_fprausnitzii16_pangenome.tsv") + ggtitle("Faecalibacterium prausnitzii")
f_prausnitzii <- plot_panphlan("input_final/panphlan/fprausnitzii16_pangenome.tsv") + ggtitle("Faecalibacterium prausnitzii")

# b ovatus
b_ovatus <- plot_panphlan("input_final/panphlan/bovatus16_pangenome.tsv") + ggtitle("Bacteroides ovatus")

# b stercoris
b_stercoris <- plot_panphlan("input_final/panphlan/bstercoris16_pangenome.tsv") + ggtitle("Bacteroides stercoris")

# b dorei
b_dorei <- plot_panphlan("input_final/panphlan/bdorei16_pangenome.tsv") + ggtitle("Bacteroides dorei")

# e coli
e_coli <- plot_panphlan("input_final/panphlan/ecoli16_pangenome.tsv") + ggtitle("Escherichia coli")

# row1 <- plot_grid(b_uniformis, b_vulgatus, labels = c("A", "B"), align = "h", label_size = 14)
# row2 <- plot_grid(npcopri, NULL, labels = c("C", ""), align = "h", label_size = 14)

# plot_grid(
#   b_uniformis + theme(legend.position = "none"),
#   b_vulgatus + theme(legend.position = "none"),
#   npcopri + theme(legend.position = "none"),
#   get_legend(npcopri),
#   nrow = 2
#   )

# [1] "Bacteroides vulgatus"         "Bacteroides ovatus"           "Prevotella copri"             "Bacteroides uniformis"        "Bacteroides dorei"            "Faecalibacterium prausnitzii"
# [7] "Escherichia coli"             "Bacteroides stercoris" 

main <- plot_grid(
  b_vulgatus + theme(legend.position = "none"),
  b_ovatus + theme(legend.position = "none"),
  p_copri + theme(legend.position = "none"),
  # npcopri + theme(legend.position = "none"),
  b_uniformis + theme(legend.position = "none"),
  b_dorei + theme(legend.position = "none"),
  f_prausnitzii + theme(legend.position = "none"),
  e_coli + theme(legend.position = "none"),
  b_stercoris + theme(legend.position = "none"),
  # get_legend(p_copri),
  ncol = 3,
  align = "hv",
  axis = "lrtb",
  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "")
)

plot_grid(main, get_legend(p_copri), rel_heights = c(0.95, 0.05), ncol = 1)
ggsave("final_plots/supplementary/panphlan_v2.png", width = 10, height = 10)
