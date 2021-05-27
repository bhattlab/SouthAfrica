library(ggplot2)
library(vegan)
library(dplyr)
library(ggrepel)
library(reshape2)
# library(ggdendro)
library(grid)
library(MASS)
library(gridExtra)
library(cowplot)

################################################################################
### Panphlan nmds ##############################################################
################################################################################

setwd("/Users/tamburif/za-microbiome")

# color palettes
global_pal <- c("#E3211C", "#F89897", "#6A3D9A", "#CAB2D6", "#1F78B4", "#A5CEE3")
za_pal <- global_pal[3:4]

# sites to include
sites <- c("Tanzania", "Madagascar", "Bushbuckridge", "Soweto", "Sweden", "United States")

site_pal <- data.frame(color = global_pal, site2 = sites)

# function to ordinate and plot panphlan results table
plot_panphlan <- function(panphlan_f){
  
  # read panphlan pangenome matrix
  panphlan <- read.table(panphlan_f, sep = "\t", header = T, row.names = 1, quote = "")
  colnames(panphlan) <- gsub('^X|_panphlan', '', colnames(panphlan))
  
  # tmp remove ecoli sample A82_02_1FE
  if (grepl("ecoli", panphlan_f)){
    panphlan <- panphlan[, -which(names(panphlan) == "A82_02_1FE")]
  }
  
  # read pheno table
  pheno_global <- readRDS("rds/pheno_global.rds")
  
  # filter for relevant samples
  # filt <- filter(pheno_global, site2 %in% sites)$sample
  # panphlan <- panphlan[, colnames(panphlan) %in% samples$V1]
  panphlan <- panphlan[, which(names(panphlan) %in% pheno_global$sample)]
  
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
  mds_pheno <- merge(mds_df, pheno_global, by = "sample", all.x = T)
  mds_pheno$site2 <- factor(mds_pheno$site2, levels = sites)
  
  # get color palette
  my_pal <- site_pal %>%
    filter(site2 %in% mds_pheno$site2) %>%
    pull(color) %>%
    as.character()
  
  # no labels
  my_plot <- ggplot(mds_pheno, aes(X1, X2, color = site2)) +
    # geom_point(size = 2, color = "black", pch = 21) +
    geom_point(size = 2, alpha = 0.75) +
    scale_color_manual(values = my_pal) +
    labs(x = "NMDS 1",
         y = "NMDS 2",
         color = "") +
    theme_cowplot(14) +
    theme(
      plot.title = element_text(face = "italic", size = 14),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      legend.position = "bottom"
    ) +
    guides(color = guide_legend(nrow = 1))
  
  my_plot
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

# ggsave("rplots/supplementary_figures/figure_3.pdf", fig_3, device = "pdf", width = 11, height = 8)
# ggsave("final_plots/supplementary/panphlan.png", width = 8, height = 8) 

# ### heatmap
# # read panphlan pangenome matrix
# panphlan <- read.table("/Users/Fiona/scg4_fiona/za/06_panphlan/panphlan_buniformis16_pangenome.tsv", sep = "\t", header = T, row.names = 1, quote = "")
# colnames(panphlan) <- gsub('^X|_panphlan', '', colnames(panphlan))
# 
# # read pheno table
# pheno_global <- read.table("/Users/Fiona/Desktop/pheno_global.tsv", sep = '\t', header = T)
# 
# # filter for relevant samples
# filt <- filter(pheno_global, site2 %in% sites)$sample
# panphlan <- panphlan[, colnames(panphlan) %in% as.character(filt)]
# 
# # filter very low and very high prevalence genes
# n_samples <- ncol(panphlan)
# panphlan[(rowSums(panphlan) < n_samples * 0.95) & (rowSums(panphlan) > n_samples * 0.05),]
# 
# panphlan$gene <- row.names(panphlan)
# pan_long <- melt(panphlan)
# 
# ggplot(pan_long, aes(variable, gene)) +
#   geom_tile(colour = "white") +
#   scale_fill_gradient(low = "white", high = "steelblue")
# 
# 
# # contingency table
# panphlan$gene <- row.names(panphlan)
# pan_long <- melt(panphlan)
# pan_long <- filter(pan_long, value == 1)
# pan_long <- merge(pan_long, pheno_global, by.x = "variable", by.y = "sample", all.x = T, all.y = F)
# pan_long <- pan_long[, c(2, 5)]
# 
# ctable <- table(pan_long)
# 
# ################################################################################
# ### P. copri clade abundance ###################################################
# ################################################################################
# 
# # p copri clade-specific markers (manually formatted from file from ade tett)
# clade_markers_long <- read.table("/Users/Fiona/scg4_fiona/za/za-microbiome/copri_markers_long.txt", sep = "\t", header = F, quote = "", na.strings = "")
# colnames(clade_markers_long) <- c("clade", "marker")
# clade_markers_long <- filter(clade_markers_long, !is.na(marker))
# 
# # panphlan marker coverage output file
# marker_coverage <- read.table("/Users/Fiona/scg4_fiona/za/06_panphlan/panphlan_npcopri_coverage_060419.tsv", sep = "\t", header = T, quote = "")
# colnames(marker_coverage)[1] <- "marker"
# 
# # melt coverage
# marker_long <- melt(marker_coverage, id.vars = "marker", variable.name = "sample", value.name = "coverage")
# 
# # remove markers less than desired coverage
# cvg <- 0.5
# marker_filt <- filter(marker_long, coverage > cvg)
# 
# # merge with clade-specific markers
# clades_merge <- merge(marker_filt, clade_markers_long, by = "marker", all.x = T)
# clades_merge_filt <- filter(clades_merge, !is.na(clade))
# 
# # count n markers per clade, find minimum number of markers to consider clade present
# clade_marker_counts <- plyr::count(clade_markers_long, "clade")
# clade_marker_counts$min <- round(clade_marker_counts$freq / 2)
# 
# # count n clade specific markers per sample
# clade_freq <- plyr::count(clades_merge_filt, c("sample", "clade"))
# 
# # filter samples with min number of clade markers present
# clades_present <- merge(clade_freq, clade_marker_counts, by = "clade")
# clades_present_filt <- filter(clades_present, freq.x >= min)
# clades_present_filt$sample <- gsub("^X|_panphlan", "", clades_present_filt$sample)
# 
# # get size of metagenome from readcount * 150 bp reads
# readcounts <- read.table("/Users/Fiona/scg4_fiona/za/za-microbiome/global_readcounts.tsv", sep = "\t", header = T, quote = "", na.strings = "") # most of the data
# colnames(readcounts)[1] <- "sample"
# 
# readcounts_tz <- read.table("/Users/Fiona/scg4_fiona/za/za-microbiome/smits_readcounts.tsv", sep = "\t", header = T, quote = "", na.strings = "") # smits 2017 data
# readcounts_tz$host_removed_reads <- (readcounts_tz$host_removed_reads_linecount / 4) # divide by 4 since this was a wc -l count
# 
# readcounts_hmp <- read.table("/Users/Fiona/scg4_fiona/za/za-microbiome/hmp_readcounts.tsv", sep = "\t", header = T, quote = "", na.strings = "") # hmp data
# readcounts_hmp$host_removed_reads <- (readcounts_hmp$host_removed_reads_linecount_1 / 2) # divide by 2 since this was a wc -l count on forward reads only
# 
# rc_all <- rbind(readcounts[, c("sample", "host_removed_reads")], readcounts_tz[, c(1,3)], readcounts_hmp[, c(1,3)])
# 
# meta_sizes <- data.frame("sample" = rc_all$sample, "bp" = 150 * rc_all$host_removed_reads)
# 
# # get mean coverage of markers per sample
# clades_merge_unfilt <- merge(marker_long, clade_markers_long, by = "marker", all.x = T)
# clades_merge_unfilt <- filter(clades_merge_unfilt, !is.na(clade))
# 
# avg_clade_cvg <- aggregate(coverage ~ sample + clade, clades_merge_unfilt, "mean")
# avg_clade_cvg$sample <- gsub("^X|_panphlan", "", avg_clade_cvg$sample)
# 
# cvg_merge <- merge(avg_clade_cvg, meta_sizes, by = "sample", all.x = T, all.y = F)
# approx_genome_size <- 3.5e6
# cvg_merge$abundance <- cvg_merge$coverage * approx_genome_size / cvg_merge$bp
# 
# # keep only samples with the min number of clade markers present for at least one clade
# cvg_merge_filt <- filter(cvg_merge, sample %in% unique(clades_present_filt$sample))
# 
# # add in global info
# pheno_global <- read.table("/Users/Fiona/Desktop/pheno_global.tsv", sep = '\t', header = T)
# cvg_merge_global <- merge(cvg_merge_filt, pheno_global, by = "sample", all.x = T, all.y = F)
# cvg_merge_global <- filter(cvg_merge_global, !is.na(site2))
# 
# cvg_merge_global$clade <- gsub("_", " ", cvg_merge_global$clade)
# cvg_merge_global$clade <- factor(cvg_merge_global$clade, levels = rev(c("Clade A", "Clade B", "Clade C", "Clade D")))
# 
# # function to draw heatmap by site
# draw_heatmap <- function(df, this_site){
#   
#   plot_df <- filter(df, site2 == this_site)
#   
#   ggplot(plot_df, aes(x = sample, y = clade)) +
#     geom_tile(aes(fill = abundance)) +
#     # scale_fill_gradient2() +
#     # scale_fill_gradient(low = "white", high = "red3", limits = c(0, 1), breaks = c(0, 1)) +
#     scale_fill_gradient(low = "white", high = "red3") +
#     # facet_wrap(~ site2, scales = "free") +
#     labs(
#       title = this_site,
#       x = "Sample",
#       fill = "Relative abundance"
#     ) +
#     theme_classic() +
#     theme(
#       axis.text.y = element_text(size = 12),
#       # legend.position = "none",
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       axis.title.y = element_blank(),
#       title = element_text(size = 14),
#       legend.title = element_text(size = 12)
#     )
# }
# 
# # heatmap for each site
# 
# plot_list <- list()
# 
# for (s in 1:length(sites)) {
#   plot_list[[s]] <- draw_heatmap(cvg_merge_global, sites[s])
# }
# 
# # dimensions for plot
# counts <- plyr::count(cvg_merge_global, "site2")
# counts$prop <- counts$freq / max(counts$freq)
# 
# lay <- rbind(c(rep(1, 10)),
#              c(rep(2, 4), rep(NA, 6)),
#              c(rep(3, 10)),
#              c(rep(4, 5), rep(NA, 5)),
#              c(rep(5, 4), rep(NA, 6)),
#              c(rep(6, 4), rep(NA, 6)),
#              c(rep(7, 5), rep(NA, 5))
# )
# 
# # arrangle all plots
# grid.arrange(grobs = plot_list, layout_matrix = lay)
# # ggsave("/Users/Fiona/Desktop/npcopri_heatmap.png", g, device = png, height = 10, width = 8)
# 
# 
# ################################################################################
# ### Heatmap/hierarchical clustering ############################################
# ################################################################################
# 
# pan_dendro <- as.dendrogram(hclust(dist_matrix))
# 
# # Create dendro
# dendro_plot <- ggdendrogram(data = pan_dendro, rotate = TRUE)
# 
# pan_plot <- panphlan
# pan_plot$gene <- rownames(pan_plot)
# pan_long <- melt(pan_plot, variable.name = "sample")
# 
# order <- order.dendrogram(pan_dendro)
# pan_long$sample <- factor(pan_long$sample, levels = colnames(panphlan)[order])
# 
# # heatmap
# heatmap <- ggplot(pan_long, aes(x = gene, y = sample)) +
#   geom_tile(aes(fill = value)) +
#   # scale_fill_gradient2() +
#   theme(
#     axis.text.y = element_text(size = 6),
#     axis.text.x = element_blank(),
#     legend.position = "none"
#     )
# 
# grid.newpage()
# print(heatmap, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 0.91))
# print(dendro_plot, vp = viewport(x = 0.90, y = 0.5, width = 0.2, height = 1))