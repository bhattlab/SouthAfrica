# Analyses and plots for South African microbiome pilot data
# Fiona Tamburini

library(ggplot2)
library(vegan)
library(reshape2)
library(RColorBrewer)
library(MASS)
library(gtools)
library(genefilter)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(ggsignif)
library(harrietr)
library(rstatix)
library(genefilter)
library(Hmisc)
library(tibble)
library(DESeq2)
library(dplyr)
library(formattable)
library(scales)
library(metagenomeSeq)

setwd("/Users/tamburif/za-microbiome")

################################################################################
########## Setup/input files ###################################################
################################################################################
# # color palettes
# global_pal <- c("#E3211C", "#F89897", "#6A3D9A", "#CAB2D6", "#1F78B4", "#A5CEE3")
# za_pal <- global_pal[3:4]
# 
# # study sites
# sites <- c("Tanzania", "Madagascar", "Bushbuckridge", "Soweto", "Sweden", "United States")
# 
# # vanish taxa
# vanish_F <- c("Prevotellaceae", "Succinivibrionaceae", "Paraprevotellaceae", "Spirochaetaceae")
# 
# taxonomy <- read.table("input_final/kraken_feb2019_inspect_mpa.out", sep = "\t", quote = "", comment.char = "")
# vanish_G <- taxonomy %>%
#   filter(grepl("g__", V1) & grepl(paste(vanish_F, collapse = "|"), V1)) %>%
#   mutate(feature = gsub("\\|s__.+", "", V1)) %>%
#   pull(feature) %>%
#   unique()
# 
# ## za
# za_meta <- readRDS("rds/za_meta.rds")
# za_pheno <- readRDS("rds/za_pheno.rds")
# 
# ## global
# # exclude rampelli participants under the age of 18
# exclude_rampelli <- c("SRR1929485",
#                       "SRR1929574",
#                       "SRR1930128",
#                       "SRR1930132",
#                       "SRR1930134")
# pheno_global <- readRDS("rds/pheno_global.rds")
# pheno_global <- pheno_global %>%
#   filter(!(sample %in% exclude_rampelli))
# 
# 
# ### genbank 2020 files
# for (rank in c("S", "G", "F", "O", "C", "P")){
#   global_reads <- read.table(paste0("input_final/genbank_jan2020/global/bracken_", rank, "_reads.txt"), sep = "\t", check.names = F, quote = "", comment.char = "")
#   rampelli_reads <- read.table(paste0("input_final/genbank_jan2020/rampelli/bracken_", rank, "_reads.txt"), sep = "\t", check.names = F, quote = "", comment.char = "")
#   # exclude non-adult samples
#   rampelli_reads <- rampelli_reads[, which(!(names(rampelli_reads) %in% exclude_rampelli))]
#   
#   global_reads <- merge(global_reads, rampelli_reads, by = "row.names", all = T) %>%
#     column_to_rownames(var = "Row.names")
#   global_reads[is.na(global_reads)] <- 0
#   
#   global_reads <- global_reads[, pheno_global$sample]
#   global_reads <- global_reads[!(rownames(global_reads) %in% c("Homo", "Homo sapiens")), ]
#   global_reads <- global_reads[rev(order(rowMeans(global_reads))), ]
#   
#   za_reads <- global_reads[, za_meta$sample]
#   za_reads <- za_reads[rev(order(rowMeans(za_reads))), ]
#   
#   global_rel <- sweep(global_reads, 2, colSums(global_reads), FUN = "/")
#   za_rel <- sweep(za_reads, 2, colSums(za_reads), FUN = "/")
#   # rel <- rel[rev(order(rowMeans(rel))), ]
#   
#   global_pseudo <- global_reads + 1
#   global_pseudo_rel <- sweep(global_pseudo, 2, colSums(global_pseudo), FUN = "/")
#   
#   za_pseudo <- za_reads + 1
#   za_pseudo_rel <- sweep(za_pseudo, 2, colSums(za_pseudo), FUN = "/")
#   
#   assign(paste0("global_", rank), global_reads)
#   assign(paste0("global_", rank, "_rel"), global_rel)
#   assign(paste0("global_", rank, "_pseudo_rel"), global_pseudo_rel)
#   
#   assign(paste0("za_", rank), za_reads)
#   assign(paste0("za_", rank, "_rel"), za_rel)
#   assign(paste0("za_", rank, "_pseudo_rel"), za_pseudo_rel)
# }
# 
# 
# # vangay data
# vangay_meta <- read.table("input_final/vangay/vangay_metadata.txt", sep = "\t", header = F)
# bioproj <- read.table("input_final/vangay/PRJEB28687_wgs.txt", sep = "\t", header = F)
# 
# vangay_S <- read.table("input_final/vangay/bracken_species_reads.txt", sep = "\t", quote = "", comment.char = "")
# names(vangay_S) <- bioproj[match(names(vangay_S), bioproj$V4), "V9"]
# names(vangay_S) <- gsub("wgs\\.", "", names(vangay_S))
# names(vangay_S) <- gsub("_", "\\.", names(vangay_S))
# 
# vangay_S <- vangay_S[rev(order(rowMeans(vangay_S))), ]
# vangay_S_rel <- sweep(vangay_S, 2, colSums(vangay_S), "/")
# 
# 
# ## load sourmash data
# for (k in c("21", "31", "51")){
#   for (suffix in c("", "_track_abund")){
#     # sourmash <- read.csv(paste0("input_final/sourmash/compare_k", k, ".csv"), header = T)
#     sourmash <- read.csv(paste0("input_final/sourmash/compare_k", k, suffix, ".csv"), header = T)
#     names(sourmash) <- gsub("X.oak.stanford.scg.lab_asbhatt.fiona.za.05_global_sourmash.sourmash.02_trim_kmers.|_concat.fq.abundtrim", "", names(sourmash))
#     rownames(sourmash) <- names(sourmash)
#     
#     # exclude rampelli samples and other extraneous samples
#     keep <- which(names(sourmash) %in% pheno_global$sample)
#     sourmash <- sourmash[keep, keep]
#     assign(paste0("sourmash_k", k, suffix), sourmash)
#   }
# }

# ## label samples something helpful for publishing
# labels <- za_meta %>% arrange(site, sample) %>% group_by(site) %>% mutate(id = paste(tolower(site), row_number(), sep = "_"))
# labels <- as.data.frame(labels)
# write.table(labels, "input_final/za_labels.tsv", sep = "\t", row.names = F, quote = F)

# # supp tables
# # only include features present at at least 0.01% relab in at least one sample
# genus_supp <- za_G_rel * 100
# names(genus_supp) <- labels[match(names(genus_supp), labels$sample), "id"]
# genus_supp <- genus_supp[, mixedorder(names(genus_supp))]
# genus_supp <- genus_supp[rowSums(genus_supp >= 0.01) >= 1, ]
# genus_supp <- genus_supp %>%
#   rownames_to_column(var = "Genus")
# # write.table(genus_supp, "final_tables/table_s3_genus_classification_rel.txt", sep = "\t", row.names = F, quote = F)
# 
# species_supp <- za_S_rel * 100
# names(species_supp) <- labels[match(names(species_supp), labels$sample), "id"]
# species_supp <- species_supp[, mixedorder(names(species_supp))]
# species_supp <- species_supp[rowSums(species_supp >= 0.01) >= 1, ]
# species_supp <- species_supp %>%
#   rownames_to_column(var = "Species")
# # write.table(species_supp, "final_tables/table_s3_species_classification_rel.txt", sep = "\t", row.names = F, quote = F)


## css normalize data
# za
# for (rank in c("S", "G")){
#   mr <- newMRexperiment(get(paste0("za_", rank)))
#   p <- cumNormStatFast(mr)
#   mr_css <- cumNorm(mr, p = p)
#   counts_css <- MRcounts(mr_css, norm = T, log = T)
#   saveRDS(counts_css, paste0("input_final/za_", rank, "_css.rds"))
# }
# za_S_css <- readRDS("input_final/za_S_css.rds")
# za_G_css <- readRDS("input_final/za_G_css.rds")
# 
# # global
# # for (rank in c("S", "G")){
# #   mr <- newMRexperiment(get(paste0("global_", rank)))
# #   p <- cumNormStatFast(mr)
# #   mr_css <- cumNorm(mr, p = p)
# #   counts_css <- MRcounts(mr_css, norm = T, log = T)
# #   saveRDS(counts_css, paste0("input_final/global_", rank, "_css.rds"))
# # }
# global_S_css <- readRDS("input_final/global_S_css.rds")
# global_G_css <- readRDS("input_final/global_G_css.rds")

# ## clr norm
# x <- aldex.clr(global_S, mc.samples = 128, denom = "all", verbose = F, useMC = T)
# global_S_clr <- sapply(getMonteCarloInstances(x), function(x) rowMeans(x))

################################################################################
########## Summary info ########################################################
################################################################################

# number of global metagenomes
pheno_global %>% tally()

# number of global metagenomes by site
pheno_global %>% group_by(site2) %>% tally()

# za samples only, by site
za_meta %>% group_by(site) %>% tally()

## For Scott
# agt - top 5 most variable taxa, at least 1% relab in 10% of population
# agt <- za_species_rel[, filter(za_meta, site == "Bushbuckridge")$sample]
# names(agt) <- za_meta[match(names(agt), za_meta$sample), "study_id"]
# 
# keep <- data.frame(genefilter(agt, pOverA(p=0.1, A=1)))
# colnames(keep) <- "taxon"
# keep$tax <- row.names(keep)
# keep <- filter(keep, taxon == T)$tax
# agt_filt <- agt[keep, ]
# agt_sd <- data.frame("sd" = apply(agt_filt, 1, sd))
# agt_sd$species <- row.names(agt_sd)
# agt_sd <- agt_sd[rev(order(agt_sd$sd)), ]
# 
# snp_test <- agt_filt[agt_sd[1:10,]$species, ]
# write.table(snp_test, "/Users/Fiona/Desktop/snp_assoc_species.tsv", sep = "\t", quote = F)

################################################################################
########## Crassphage ##########################################################
################################################################################
# 
# # estimate crAssphage prevalence in population
# # bracken reports read PAIRS -- divide threshold by 2
# crass <- "crAss-like viruses"
# 
# n <- 650 # ~1X coverage of a 97kb genome by 150 bp reads
# length(which(za_G[crass, ] >= n / 2))
# 
# za_meta %>%
#   filter(sample %in% names(za_G[, which(za_G[crass, ] >= n / 2)])) %>%
#   group_by(site) %>%
#   tally()
# 
# tbl <- za_meta %>%
#   mutate(crass_TRUE = sample %in% names(za_G[, which(za_G[crass, ] >= n / 2)]))
# 
# t <- table(tbl$site, tbl$crass_TRUE)
# fisher.test(t)
# 
# # cor with other taxa
# za_S_rel_filt <- za_S_rel[which(rowSums(za_S_rel) > 0), ]
# za_S_rel_filt <- za_S_rel_filt[genefilter(za_S_rel_filt, pOverA(0.1, 0.0001)), ]
# 
# cor_long <- melt_dist(cor(t(za_S_rel_filt), method = "spearman"))
# 
# cor_long <- cor_long %>%
#   filter(grepl("Guerin_crAss", iso1) | grepl("Guerin_crAss", iso1)) %>%
#   arrange(-abs(dist)) %>%
#   filter(abs(dist) > 0.4)
# 
# # plot heatmap
# keep <- c(cor_long$iso1, cor_long$iso2)
# 
# cor_long_heatmap <- melt(cor(t(za_S_rel_filt), method = "spearman"), value.name = "cor")
# cor_long_heatmap <- cor_long_heatmap %>%
#   filter(Var1 %in% keep & Var2 %in% keep)
# 
# ggplot(cor_long_heatmap, aes(Var2, Var1, fill = cor)) +
#   geom_tile() +
#   geom_text(aes(label = round(cor, 2))) +
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white",
#                        midpoint = 0, limit = c(-1,1)) +
#   theme_cowplot() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#     axis.title = element_blank()
#   ) +
#   labs(
#     fill = "Spearman's\nrho"
#   )


################################################################################
########## Figure 1: taxonomy barplots and VANISH boxplots #####################
################################################################################
# 
# ### Figure 1A: taxonomy barplots
# sort_by <- "Prevotella"
# n_taxa <- 20
# 
# # color palette for n taxa
# # myCols <- colorRampPalette(brewer.pal(9, "Set1"))
# # global_pal <- myCols(n_taxa)
# # global_pal <- sample(global_pal)
# # global_pal[n_taxa + 1] <- "gray"
# # saveRDS(global_pal, "rds/global_pal.rds")
# barplot_pal <- readRDS("rds/global_pal.rds")
# 
# # find top n taxa
# abundance_threshold <- sort(rowSums(za_G_rel), decreasing = T)[n_taxa]
# bracken_plot <- za_G_rel[rowSums(za_G_rel) >= abundance_threshold,]
# 
# # add "other" column
# bracken_plot <- rbind(bracken_plot, t(data.frame("Other" =  1 - colSums(bracken_plot))))
# 
# # melt data frame
# bracken_plot$Genus <- row.names(bracken_plot)
# 
# # remove annoying "miscellaneous" labels that take up space
# bracken_plot$Genus <- gsub("\\(miscellaneous\\)", "", bracken_plot$Genus)
# bracken_long <- melt(bracken_plot, id.vars = "Genus", variable.name = "sample", value.name = "rel_abundance")
# 
# # merge in pheno date
# bracken_pheno <- merge(bracken_long, za_meta, by = "sample")
# bracken_pheno$label <- labels[match(bracken_pheno$sample, labels$sample), "id"]
# 
# # set factor level for correct plotting order
# bracken_pheno$Genus <- factor(bracken_pheno$Genus, levels = bracken_plot$Genus)
# 
# # shorten labels
# bracken_pheno$label <- gsub("bushbuckridge_", "B", bracken_pheno$label)
# bracken_pheno$label <- gsub("soweto_", "S", bracken_pheno$label)
# 
# # plot in order of decreasing relative abundance of desired taxon
# sorted <- bracken_pheno %>% filter(Genus == sort_by) %>% arrange(desc(rel_abundance)) %>% pull(label)
# bracken_pheno$label <- factor(bracken_pheno$label, levels = sorted)
# 
# plot_bracken <- function(counts, title){
#   g <- ggplot(counts, aes(x=label, y=rel_abundance * 100, fill=Genus)) +
#     geom_bar(stat="identity") +
#     labs(
#       title = title,
#       x = "Sample",
#       y = "Relative Abundance (%)"
#     ) +
#     scale_fill_manual(values = barplot_pal) +
#     guides(fill = guide_legend(ncol=1, keywidth = 0.125, keyheight = 0.1, default.unit = "inch")) +
#     theme_cowplot(12) +
#     theme(
#       plot.title = element_text(face = "plain", size = 14),
#       axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
#       legend.text = element_text(size = 10)
#     ) +
#     scale_y_continuous(limits = c(0, 100.1), expand = c(0, 0))
#   
#   return(g)
# }
# 
# agin <- plot_bracken(filter(bracken_pheno, site == "Bushbuckridge"), "Bushbuckridge") + theme(legend.position = "none")
# soweto <- plot_bracken(filter(bracken_pheno, site == "Soweto"), "Soweto")
# 
# a <- plot_grid(agin, soweto, align = "v", ncol = 1, axis = "l")
# 
# 
# ### Figure 1B: VANISH taxa
# vanish_tax <- data.frame(
#   "F" = gsub(".+f__|\\|g__.+", "", vanish_G),
#   "G" = gsub(".+g__", "", vanish_G)
# )
# 
# za_G_VANISH <- za_G_pseudo_rel[which(rownames(za_G_pseudo_rel) %in% as.character(vanish_tax$G)), ]
# za_F_VANISH <- za_F_pseudo_rel[which(rownames(za_F_pseudo_rel) %in% as.character(vanish_tax$F)), ]
# 
# counts_long <- reshape2::melt(za_G_VANISH %>% rownames_to_column(var = "G"),
#                               id.vars = "G",
#                               variable.name = "sample",
#                               value.name = "relab")
# 
# counts_long$G <- factor(counts_long$G, levels = rev(unique(counts_long$G)))
# counts_long <- merge(counts_long, za_meta, by.x = "sample", by.y = "sample")
# 
# # facet by family
# counts_long$F <- vanish_tax[match(counts_long$G, vanish_tax$G), "F"]
# 
# # compare BBR vs SWT means
# counts_long %>%
#   group_by(G, site) %>%
#   summarise(mean = mean(relab)) %>%
#   group_by(G) %>%
#   summarise("inc_BBR" <- mean[site == "Bushbuckridge"] > mean[site == "Soweto"])
# 
# # plot in decreasing order
# levels <- counts_long %>%
#   group_by(G) %>%
#   summarise(mean_relab = mean(relab)) %>%
#   arrange(-mean_relab) %>%
#   pull(G)
# counts_long$G <- factor(counts_long$G, levels = levels)
# 
# # counts_long$relab[counts_long$relab == 0] <- 1e-8
# 
# b <- ggplot(counts_long, aes(G, relab * 100, fill = site)) +
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.75, aes(group = site), size = 1, color = "darkgray") +
#   geom_boxplot(outlier.shape = NA, alpha = 0.75) +
#   facet_grid(~ F, scales = "free", space = "free") +
#   scale_y_log10(labels = c("0.0001", "0.01", "1", "100"), breaks = c(0.0001, 0.01, 1, 100)) +
#   theme_cowplot() +
#   theme(
#     legend.position = "bottom",
#     legend.justification = "center",
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
#     axis.title.y = element_text(size = 12),
#     strip.text = element_text(face = "italic"),
#   ) +
#   scale_fill_manual(values = za_pal) +
#   scale_color_manual(values = za_pal) +
#   labs(
#     x = "",
#     y = "Relative Abundance (%)",
#     fill = ""
#   ) +
#   stat_compare_means(label = "p.signif", label.y = 2.5)
# 
# 
# ### Figure 1
# plot_grid(a, b, labels = c("A", "B"), ncol = 1, rel_heights = c(0.55, 0.45))
# ggsave("final_plots/figure_1.png", width = 10, height = 12)


################################################################################
########## Figure S2: Bimodal distribution of VANISH taxa ######################
################################################################################

# ## Succinatimonas correlation with other genera
# za_G_cor <- cor(t(za_G_rel[genefilter(za_G_rel, pOverA(0.1, 0.0001)), ]), method = "spearman")
# 
# za_G_cor_long <- melt_dist(za_G_cor)
# za_G_cor_long %>%
#   filter(iso1 == "Succinatimonas" | iso2 == "Succinatimonas") %>%
#   arrange(-abs(dist)) %>%
#   head(20)
# 
# ## plot distributions
# za_VANISH_G_long <- reshape2::melt(za_G_VANISH %>% rownames_to_column(var = "G"),
#                id.vars = "G",
#                variable.name = "sample",
#                value.name = "relab")
# za_VANISH_G_long <- merge(za_VANISH_G_long, za_meta, by = "sample")
# 
# succin_long <- za_VANISH_G_long %>%
#   filter(G %in% c("Succinatimonas", "Succinivibrio", "Treponema") & site == "Bushbuckridge")
# 
# density_plot <- ggplot(succin_long, aes(x = relab * 100)) +
#   geom_density(fill = za_pal[1]) +
#   scale_x_log10(label = comma) +
#   facet_wrap(~ G, scales = "free_y", ncol = 1) +
#   theme_cowplot() +
#   labs(
#     x = "Relative abundance (%)",
#     y = "Density"
#   )
# 
# ## correlate VANISH taxa
# vanish_cor <- function(counts){
#   # correlation
#   vanish_cor_G <- cor(t(counts), method = "spearman")
#   
#   # cluster to order rows
#   hc <- hclust(dist(vanish_cor_G))
#   row_order <- hc$order
#   vanish_cor_G <- vanish_cor_G[row_order, row_order]
#   
#   # plot
#   vanish_cor_G_long <- melt(vanish_cor_G, value.name = "cor")
#   
#   ggplot(vanish_cor_G_long, aes(Var2, Var1, fill = cor)) +
#     geom_tile() +
#     geom_text(aes(label = round(cor, 2))) +
#     scale_fill_gradient2(low = "blue", high = "red", mid = "white",
#                          midpoint = 0, limit = c(-1,1)) +
#     theme_cowplot() +
#     theme(
#       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#       axis.title = element_blank()
#     ) +
#     labs(
#       fill = "Spearman's\nrho"
#     )  
# }
# 
# plot_data <- za_G_VANISH
# rownames(plot_data) <- paste0(rownames(plot_data), "\n(", vanish_tax[match(rownames(plot_data), vanish_tax$G), "F"], ")")
# 
# plot_G <- vanish_cor(plot_data)
# plot_F <- vanish_cor(za_F_VANISH)
# 
# plot_grid(
#   plot_grid(density_plot, plot_F, labels = c("A", "B")),
#   plot_G,
#   ncol = 1,
#   labels = c("", "C"),
#   rel_heights = c(0.27, 0.63)
#   )
# 
# ggsave("final_plots/supplementary/figure_S2_VANISH_cor.png", width = 9, height = 11)

################################################################################
########## Figure 2: Bushbuckridge, Soweto comparison ##########################
################################################################################
# 
# ### A) MDS ordination
# # find bray-curtis distance with vegdist
# # vare_dis <- vegdist(t(za_S_rel), method = "bray")
# vare_dis <- vegdist(t(za_S_css), method = "bray")
# 
# # adonis
# meta <- za_meta %>% filter(sample %in% names(za_S_rel))
# meta <- meta[match(names(za_S_rel), meta$sample), ]
# adonis(vare_dis ~ site, data = meta)
# dispersion <- betadisper(vare_dis, group = meta$site)
# permutest(dispersion)
# 
# # calculate mds
# mds <- cmdscale(vare_dis, eig = TRUE, x.ret = TRUE)
# 
# mds_values <- mds$points
# wa_scores <- wascores(mds_values, t(za_S_rel))
# wa_scores <- data.frame(sample = rownames(wa_scores),
#                         x = wa_scores[,1],
#                         y = wa_scores[,2])
# 
# # isolate taxa with strongest contribution to principal coordinate axes
# n_taxa <- 10
# wa_scores_1<- head(arrange(wa_scores, desc(abs(wa_scores$x))), n = n_taxa)
# wa_scores_2<- head(arrange(wa_scores, desc(abs(wa_scores$y))), n = n_taxa)
# wa_scores_final <- rbind(wa_scores_1, wa_scores_2)
# 
# # calculate percentage of variation that each MDS axis accounts for
# mds_var_per <- round(mds$eig/sum(mds$eig) * 100, 1)
# 
# # plot
# mds_data <- data.frame(sample = rownames(mds_values),
#                        x = mds_values[,1],
#                        y = mds_values[,2])
# 
# # merge pheno data
# mds_meta <- merge(mds_data, za_meta, by = "sample")
# 
# mds_plot <- ggplot(mds_meta, aes(x, y, color = site)) +
#   geom_point(size = 2, alpha = 0.75) +
#   scale_color_manual(values = za_pal) +
#   labs(x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
#        y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
#        color = "") +
#   theme_cowplot()
# 
# mds_plot_ci <- mds_plot + stat_ellipse(aes(color = site), type = 't', size = 1, show.legend = F)
# 
# 
# ## plot mds plot highlighting nanopore samples
# nmag_samples <- c("C27", "C29", "C33")
# 
# mds_meta_nmag <- mds_meta %>%
#   mutate(highlight = ifelse(sample %in% nmag_samples, T, F))
# 
# ggplot(mds_meta_nmag, aes(x, y, color = highlight, size = highlight, shape = site)) +
#   geom_point(alpha = 0.75) +
#   scale_color_manual(values = c("darkgrey", "red")) +
#   scale_size_manual(values = c(2, 3)) +
#   labs(
#     x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
#     y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
#     shape = "Site"
#     ) +
#   theme_cowplot() +
#   guides(size = F, color = F)
# ggsave("final_plots/misc/nanopore_mds.png", width = 7, height = 5)
# 
# 
# ### A) Shannon diversity
# # find shannon diversity with vegdist
# shannon_div <- diversity(t(za_S_rel), index = "shannon")
# div <- data.frame("shannon_div" = shannon_div, "sample" = names(shannon_div))
# div_meta <- merge(div, za_meta, by = "sample")
# 
# # pval
# # pval <- compare_means(shannon_div ~ site2, data = div_meta, method = "wilcox.test", p.adjust.method = "fdr")
# 
# div_plot <- ggplot(div_meta, aes(site, shannon_div)) + 
#   # geom_boxplot(outlier.shape = NA, aes(fill = site)) +
#   # geom_jitter(position = position_jitterdodge(jitter.width = 0.75), alpha = 0.5, aes(fill = site), pch = 21) +
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.6), alpha = 0.75, aes(fill = site), color = "darkgray") +
#   geom_boxplot(outlier.shape = NA, aes(fill = site), alpha = 0.5) +
#   labs(x = "",
#        y = "Shannon Diversity",
#        fill="") +
#   scale_color_manual(values = za_pal) +
#   scale_fill_manual(values = za_pal) +
#   scale_y_continuous(limits = c(2, 6.6)) +
#   geom_signif(comparisons = list(c("Soweto", "Bushbuckridge")), map_signif_level = TRUE, y_position = 6.5, test = "wilcox.test") +
#   theme_cowplot(14) +
#   theme(
#     axis.text.x = element_text(size = 12),
#     legend.position = "none"
#   )
# 
# bbr <- div_meta %>% filter(site == "Bushbuckridge") %>% pull(shannon_div)
# swt <- div_meta %>% filter(site == "Soweto") %>% pull(shannon_div)
# wilcox.test(bbr, swt, exact = T)
# 
# 
# ## C) Differential features - DESeq2
# count_data <- za_G
# count_data <- count_data[, order(names(count_data))]
# 
# # col_data <- za_pheno %>% filter(sample %in% names(count_data)) %>% column_to_rownames(var = "sample") 
# col_data <- za_pheno %>% filter(sample %in% names(count_data))
# col_data <- col_data[order(row.names(col_data)), ]
# col_data$site <- factor(col_data$site, levels = c("Bushbuckridge", "Soweto"))
# 
# count_data_filt <- count_data[genefilter(count_data, pOverA(p = 0.20, A = 1000)), ]
# 
# dds <- DESeqDataSetFromMatrix(countData = count_data_filt,
#                               colData = col_data,
#                               design= ~ site)
# # dds <- estimateSizeFactors(dds, type = "poscounts")
# dds <- DESeq(dds)
# resultsNames(dds) # lists the coefficients
# res <- results(dds, name="site_Soweto_vs_Bushbuckridge", alpha = 0.05)
# resOrdered <- res[order(res$pvalue),]
# resOrdered$site <- ifelse(resOrdered$log2FoldChange < 0, "Bushbuckridge", "Soweto")
# 
# res_df <- data.frame(resOrdered, row.names = row.names(resOrdered))
# res_df_filt <- res_df %>% rownames_to_column(var = "genus") %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% arrange(log2FoldChange)
# res_df_filt$genus <- factor(res_df_filt$genus, levels = as.character(res_df_filt$genus))
# # res_df_filt$site <- ifelse(res_df_filt$log2FoldChange < 0, "Bushbuckridge", "Soweto")
# 
# deseq_genera <- ggplot(res_df_filt, aes(log2FoldChange, genus, fill = site)) +
#   # geom_point(size = 2, color = "black", pch = 21) +
#   geom_bar(stat = "identity") +
#   theme_cowplot() +
#   theme(axis.text.y = element_text(face = "italic")) +
#   scale_fill_manual(values = za_pal) +
#   scale_x_continuous(breaks = seq(-2, 2, 1)) +
#   labs(
#     x = "Log2 Fold Change (Soweto/Bushbuckridge)",
#     y = "Genus",
#     fill = ""
#   ) +
#   background_grid()
# 
# res_table_G <- data.frame(resOrdered) %>%
#   rownames_to_column(var = "Feature") %>%
#   arrange(site, log2FoldChange) %>%
#   mutate(rank = "genus")
# # write.table(res_table_G, "final_tables/deseq_genera.txt", sep = "\t", row.names = F, quote = F)
# 
# ## Figure 1
# row1 <- plot_grid(mds_plot_ci, div_plot, labels = c('A', 'B'), align = "h", label_size = 14, rel_widths = c(0.65, 0.35))
# plot_grid(row1, deseq_genera, labels = c("", "C"), ncol = 1, label_size = 14, rel_heights = c(0.4, 0.6))
# 
# ggsave("final_plots/figure_2.png", width = 9, height = 10)


################################################################################
########## Figure S4: Bushbuckridge/Soweto DESeq2 species ######################
################################################################################
# 
# count_data <- za_S
# count_data <- count_data[, order(names(count_data))]
# 
# col_data <- za_pheno %>% filter(sample %in% names(count_data)) %>% column_to_rownames(var = "sample") 
# col_data <- col_data[order(row.names(col_data)), ]
# col_data$site <- factor(col_data$site, levels = c("Bushbuckridge", "Soweto"))
# 
# count_data_filt <- count_data[genefilter(count_data, pOverA(p = 0.20, A = 1000)), ]
# 
# dds <- DESeqDataSetFromMatrix(countData = count_data_filt,
#                               colData = col_data,
#                               design= ~ site)
# dds <- DESeq(dds)
# resultsNames(dds) # lists the coefficients
# res <- results(dds, name="site_Soweto_vs_Bushbuckridge", alpha = 0.05)
# resOrdered <- res[order(res$pvalue),]
# resOrdered$site <- ifelse(resOrdered$log2FoldChange < 0, "Bushbuckridge", "Soweto")
# 
# res_df <- data.frame(resOrdered, row.names = row.names(resOrdered))
# res_df_filt <- res_df %>% rownames_to_column(var = "species") %>% filter(padj < 0.05 & abs(log2FoldChange) > 2) %>% arrange(log2FoldChange)
# res_df_filt$species <- factor(res_df_filt$species, levels = as.character(res_df_filt$species))
# # res_df_filt$site <- ifelse(res_df_filt$log2FoldChange < 0, "Bushbuckridge", "Soweto")
# 
# deseq_species <- ggplot(res_df_filt, aes(log2FoldChange, species, fill = site)) +
#   geom_bar(stat = "identity") +
#   theme_cowplot() +
#   theme(
#     axis.text.y = element_text(face = "italic"),
#     legend.position = "bottom",
#     legend.justification = "center"
#     ) +
#   scale_fill_manual(values = za_pal) +
#   scale_x_continuous(breaks = seq(-5, 5, 2.5)) +
#   labs(
#     x = "Log2 Fold Change (Soweto/Bushbuckridge)",
#     y = "Species",
#     fill = ""
#     ) +
#   background_grid()
# ggsave("final_plots/supplementary/figure_S4_deseq_site_species.png", deseq_species, width = 8.5, height = 10)
# 
# res_table_S <- data.frame(resOrdered) %>%
#   rownames_to_column(var = "Feature") %>%
#   arrange(site, log2FoldChange) %>%
#   mutate(rank = "species")
# 
# res_table <- rbind(res_table_G, res_table_S) %>%
#   select(Feature, rank, log2FoldChange, pvalue, padj)
# 
# names(res_table) <- c("Feature", "Rank", "Log2 fold change (Soweto/Bushbuckridge)", "P-value", "Adjusted p-value")
# 
# write.table(res_table, "final_tables/deseq_supp.txt", sep = "\t", row.names = F, quote = F)


################################################################################
########## Figure S1: Most abundant species and genera in ZA data ##############
################################################################################
# 
# # plot top taxa by mean relative abundance
# top_plot <- function(counts, n = 10){
#   counts_long <- reshape2::melt(counts[1:n, ] %>% rownames_to_column(var = "taxon"), variable.name = "sample", value.name = "relab")
#   
#   counts_long$taxon <- factor(counts_long$taxon, levels = rev(unique(counts_long$taxon)))
#   counts_long <- merge(counts_long, za_meta, by = "sample")
#   
#   counts_long$site <- factor(counts_long$site, levels = c("Soweto", "Bushbuckridge"))
#   
#   ggplot(counts_long, aes(taxon, relab * 100, fill = site)) +
#     geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.75, color = "darkgray", size = 1.25) +
#     geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#     theme_cowplot(12) +
#     theme(
#       legend.position = "top"
#     ) +
#     scale_fill_manual(values = za_pal, breaks = c("Bushbuckridge", "Soweto")) +
#     labs(
#       x = "",
#       y = "Relative abundance (%)",
#       color = "",
#       fill = "Site"
#     ) +
#     coord_flip() +
#     background_grid()
# }
# 
# ## top species plot
# top_S_bracken <- top_plot(za_S_rel)
# 
# ## top genus plot
# g <- za_G_rel
# # rownames(g) <- gsub("miscellaneous", "misc", rownames(g))
# top_G_bracken <- top_plot(g)
# 
# plot_grid(
#   top_S_bracken,
#   top_G_bracken + theme(legend.position = "none"),
#   ncol = 1,
#   labels = c("A", "B"),
#   align = "hv",
#   axis = "l"
#   )
# ggsave("final_plots/supplementary/figure_S1_top_taxa.png", width = 8, height = 10)
# 
# ## top phyla plot
# # top_P_bracken <- top_plot(bracken_P_rel, n = 8)

################################################################################
########## Figure S3: Bacteroides/Prevotella MDS ###############################
################################################################################
# bray_dist <- vegdist(t(za_S_css), method = "bray")
# mds <- cmdscale(bray_dist, eig = T)
# mds_scores <- as.data.frame(scores(mds))
# names(mds_scores) <- c("x", "y")
# mds_var_perc <- round(mds$eig/sum(mds$eig) * 100, 1)
# 
# # merge taxonomic data
# relab_P <- data.frame(t(za_P_rel)) * 100
# relab_G <- data.frame(t(za_G_rel[1:500, ])) * 100
# relab_S <- data.frame(t(za_S_rel[1:500, ])) * 100
# mds_taxa <- merge(mds_scores, relab_P, by = "row.names")
# # mds_taxa <- merge(mds_scores, relab_G, by = "row.names")
# mds_taxa <- merge(mds_taxa, relab_G, by.x = "Row.names", by.y = "row.names")
# mds_taxa <- merge(mds_taxa, relab_S, by.x = "Row.names", by.y = "row.names")
# mds_taxa <- merge(mds_taxa, za_meta, by.x = "Row.names", by.y = "sample")
# 
# # Bacteroides/Prevotella ratio
# mds_taxa <- mds_taxa %>%
#   mutate(
#     Bact_Prev_ratio = log2(Bacteroides / Prevotella),
#     Bacteroidetes_Firm_ratio = log2(Bacteroidetes / Firmicutes)
#     )
# 
# gradient_pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
# 
# plot_by <- function(color_by){
#   
#   ggplot(mds_taxa, aes(x, y, color = !!sym(color_by), shape = site)) +
#     geom_point(size = 3, alpha = 0.75) +
#     scale_color_gradientn(colours = gradient_pal(100), limits = c(min(mds_taxa[, color_by]), max(mds_taxa[, color_by]))) +
#     labs(
#       x = paste("MDS1 (", mds_var_perc[1], "%)",sep=""),
#       y = paste("MDS2 (", mds_var_perc[2], "%)",sep=""),
#       shape = "Site"
#     ) +
#     theme_cowplot()
# }
# 
# # mds_Bacteroides <- plot_by("Bacteroides")
# # mds_Prevotella <- plot_by("Prevotella")
# # mds_Faecalibacterium <- plot_by("Faecalibacterium")
# 
# mds_Bact_Prev_ratio <- plot_by("Bact_Prev_ratio")
# # mds_Bact_Firm_ratio <- plot_by("Bacteroidetes_Firm_ratio")
# 
# ## plot features contributing to separation
# wa_scores <- wascores(mds_scores, t(za_G_css[genefilter(za_G_css, pOverA(p = 0.10, A = 0.01)), ]), expand = T)
# wa_scores <- data.frame(taxonomy = rownames(wa_scores),
#                         x = wa_scores[,1],
#                         y = wa_scores[,2])
# # wa_scores$taxonomy <- gsub(".+__", "", wa_scores$taxonomy)
# # taxa with strongest contribution to axes
# n_taxa <- 7
# wa_scores_1 <- head(arrange(wa_scores, desc(abs(wa_scores$x))), n = n_taxa)
# wa_scores_2 <- head(arrange(wa_scores, desc(abs(wa_scores$y))), n = n_taxa)
# wa_scores_final <- unique(rbind(wa_scores_1, wa_scores_2))
# wa_scores_final <- wa_scores_final %>%
#   filter(taxonomy != "Bacteria: environmental samples")
#   
# 
# # a <- mds_Bacteroides +
# #   geom_point(data = wa_scores_final, aes(x, y), size = 1, shape = 4, inherit.aes = F) +
# #   geom_text_repel(data = wa_scores_final, color = "black", size = 3, aes(x, y, label = taxonomy), inherit.aes = F)
# # b <- mds_Prevotella +
# #   geom_point(data = wa_scores_final, aes(x, y), size = 1, shape = 4, inherit.aes = F) +
# #   geom_text_repel(data = wa_scores_final, color = "black", size = 3, aes(x, y, label = taxonomy), inherit.aes = F)
# # 
# # plot_grid(a, b, ncol = 1, labels = c("A", "B"))
# 
# mds_Bact_Prev_ratio +
#   # geom_point(data = wa_scores_final, aes(x, y), size = 1, shape = 4, inherit.aes = F) +
#   # geom_text_repel(data = wa_scores_final, color = "black", size = 4, aes(x, y, label = taxonomy), inherit.aes = F) +
#   theme(
#     legend.position = "bottom",
#     legend.justification = "center",
#     legend.text = element_text(margin = margin(r = 5, unit = "pt"))
#     ) +
#   labs(
#     color = "log2(Bacteroides/\nPrevotella)",
#     shape = "      Site"
#   ) +
#   guides(color = guide_colorbar(title.position = "left", title.vjust = 0.8),
#          shape = guide_legend(title.position = "left", title.vjust = 0.8, ncol = 1))
#   
# ggsave("final_plots/supplementary/figure_S3_bacteroides_prevotella.png", width = 6, height = 6)
# 
# ggplot(mds_taxa %>% filter(Row.names != "A18"), aes(log2(Bacteroides), log2(Prevotella))) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   theme_cowplot()
# 
# # pull factors most strongly correlated with MDS axes
# # cor <- cor(mds_taxa[, -1], mds_taxa[, -1])
# # cor_long <- melt(cor)
# # cor_long %>%
# #   filter(Var1 == "x" & Var2 != "x") %>%
# #   arrange(-value) %>%
# #   head()


################################################################################
########## Figure S6: VANISH taxa in global cohort #############################
################################################################################
# 
# vanish_tax <- data.frame(
#   "F" = gsub(".+f__|\\|g__.+", "", vanish_G),
#   "G" = gsub(".+g__", "", vanish_G)
# )
# 
# global_G_VANISH <- global_G_pseudo_rel[which(rownames(global_G_pseudo_rel) %in% vanish_tax$G), ]
# counts_long <- reshape2::melt(global_G_VANISH %>% rownames_to_column(var = "G"),
#                               id.vars = "G",
#                               variable.name = "sample",
#                               value.name = "relab")
# 
# counts_long$G <- factor(counts_long$G, levels = rev(unique(counts_long$G)))
# counts_long <- merge(counts_long, pheno_global, by = "sample")
# 
# # facet by family
# counts_long$F <- vanish_tax[match(counts_long$G, vanish_tax$G), "F"]
# 
# # plot in decreasing order
# levels <- counts_long %>%
#   group_by(G) %>%
#   summarise(mean_relab = mean(relab)) %>%
#   arrange(-mean_relab) %>%
#   pull(G)
# counts_long$G <- factor(counts_long$G, levels = levels)
# 
# counts_long$relab[counts_long$relab == 0] <- 1e-3
# 
# ggplot(counts_long, aes(G, relab * 100, fill = site2)) +
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9), alpha = 0.5, aes(group = site2), size = 1, color = "darkgray") +
#   geom_boxplot(outlier.shape = NA, width = 0.9, alpha = 0.75) +
#   facet_wrap(~ F, scales = "free", nrow = 3) +
#   scale_y_log10(labels = comma) +
#   theme_cowplot(14) +
#   theme(
#     legend.position = "bottom",
#     legend.justification = "center",
#     axis.text.x = element_text(face = "italic"),
#     strip.text = element_text(face = "italic"),
#     axis.title.y = element_text(size = 12)
#   ) +
#   scale_fill_manual(values = global_pal) +
#   labs(
#     x = "",
#     y = "Relative Abundance (%)",
#     fill = ""
#   ) +
#   guides(fill = guide_legend(nrow = 1))
# 
# ggsave("final_plots/supplementary/figure_S6_VANISH_global.png", width = 8.5, height = 11)


################################################################################
########## Figure 3: Global comparison #########################################
################################################################################
# 
# # find bray-curtis distance with vegdist
# # vare_dis <- vegdist(t(global_S_rel), method = "bray")
# vare_dis <- vegdist(t(global_S_css), method = "bray")
# 
# # adonis
# meta <- pheno_global %>% filter(sample %in% names(global_S_rel))
# meta <- meta[match(names(global_S_rel), meta$sample), ]
# adonis(vare_dis ~ site2, data = meta)
# dispersion <- betadisper(vare_dis, group = meta$site2)
# permutest(dispersion)
# 
# ## MDS Calculations (Eigen value decomposition )
# mds <- cmdscale(vare_dis, eig = TRUE, x.ret = TRUE)
# 
# ## Calculate weighted species score
# mds_values <- mds$points
# 
# ## Calculate percentage of variation that each MDS axis accounts for
# mds_var_per <- round(mds$eig/sum(mds$eig) * 100, 1)
# 
# ## Plot
# mds_data <- data.frame(sample = rownames(mds_values),
#                        x = mds_values[,1],
#                        y = mds_values[,2])
# 
# # merge pheno data
# mds_meta <- merge(mds_data, pheno_global, by = "sample")
# 
# # flip axes for consistency
# mds_meta$x <- mds_meta$x * -1
# mds_meta$y <- mds_meta$y * -1
# 
# global_mds <- ggplot(mds_meta, aes(x, y, color = site2)) +
#   geom_point(size = 2, alpha = 0.75) +
#   scale_color_manual(values = global_pal) +
#   labs(x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
#        y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
#        color= " ") +
#   theme_cowplot() +
#   theme(
#     legend.position = "top"
#   )
# 
# global_mds_ci <- global_mds + stat_ellipse(aes(color = site2), type = 't', size = 0.5, alpha = 0.75, show.legend = F)
# 
# 
# ## plot features contributing to separation
# wa_scores <- wascores(-1 * mds_values, t(global_G_css[genefilter(global_G_css, pOverA(p = 0.1, A = 5)), ]), expand = T)
# # flip axes here too
# wa_scores <- data.frame(taxonomy = rownames(wa_scores),
#                         x = wa_scores[,1],
#                         y = wa_scores[,2])
# # wa_scores$taxonomy <- gsub(".+__", "", wa_scores$taxonomy)
# # taxa with strongest contribution to axes
# n_taxa <- 7
# wa_scores_1 <- head(arrange(wa_scores, desc(abs(wa_scores$x))), n = n_taxa)
# wa_scores_2 <- head(arrange(wa_scores, desc(abs(wa_scores$y))), n = n_taxa)
# wa_scores_final <- unique(rbind(wa_scores_1, wa_scores_2))
# wa_scores_final <- wa_scores_final %>%
#   filter(taxonomy != "environmental samples")
# 
# global_mds_ci +
#   geom_point(data = wa_scores_final, aes(x, y), size = 1, shape = 4, inherit.aes = F) +
#   geom_text_repel(data = wa_scores_final, color = "black", size = 4, aes(x, y, label = taxonomy), inherit.aes = F)
# 
# 
# ## taxonomy plots a la Smits 2017
# # which taxa correlate most with MDS 1 and 2?
# global_F_rel_filt <- global_F_rel[genefilter(global_F_rel, pOverA(0.1, 0.001)), ]
# 
# global_F_t <- t(global_F_rel_filt)
# mds1_values <- merge(mds_meta[, c("sample", "x", "y")], global_F_t, by.x = "sample", by.y = "row.names", all = T) %>%
#   column_to_rownames(var = "sample")
# 
# mds1_cor <- cor(mds1_values, method = "spearman")
# mds1_cor_long <- melt(mds1_cor, value.name = "cor")
# 
# n_features <- 4
# 
# features_x <- mds1_cor_long %>%
#   filter(Var1 == "x" & Var2 != "x") %>%
#   arrange(-abs(cor)) %>%
#   head(n_features) %>%
#   pull(Var2) %>%
#   as.character()
# 
# features_y <- mds1_cor_long %>%
#   filter(Var1 == "y" & Var2 != "y") %>%
#   arrange(-abs(cor)) %>%
#   head(n_features) %>%
#   pull(Var2) %>%
#   as.character()
# 
# # plot taxa vs mds1
# global_F_rel_filt <- global_F_rel[genefilter(global_F_rel, pOverA(0.1, 0.001)), ]
# global_F_long <- melt(data.matrix(global_F_rel_filt))
# names(global_F_long) <- c("feature", "sample", "rel_abundance")
# mds_taxa <- merge(mds_meta, global_F_long, by = "sample")
# 
# # features <- c("Bacteroides", "Prevotella", "Succinatimonas", "Treponema")
# mds_taxa_plot <- mds_taxa %>%
#   filter(feature %in% features_x)
# mds_taxa_plot$feature <- factor(mds_taxa_plot$feature, levels = features_x)
# 
# scatter_F <- ggplot(mds_taxa_plot, aes(x = x, y = rel_abundance, color = site2)) +
#   geom_point(size = 1) +
#   # geom_smooth(method = "loess", color = "black") +
#   # scale_y_log10() +
#   facet_wrap(feature ~ ., scales = "free", ncol = 1, strip.position = "left") +
#   scale_color_manual(values = global_pal) +
#   theme_cowplot(11) +
#   theme(
#     axis.title = element_blank(),
#     strip.background = element_blank(),
#     strip.placement = "outside",
#     strip.text.y.left = element_text(angle = 90, face = "italic"),
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.line.y = element_blank()
#   )
# 
# # plot mds axes
# mds1 <- ggplot(mds_meta, aes(site2, x, fill = site2)) +
#   geom_jitter(alpha = 0.75, color = "darkgray", width = 0.25) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.75) +
#   scale_fill_manual(values = global_pal) +
#   theme_cowplot(12) +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     legend.position = "none"
#   ) +
#   labs(
#     x = "",
#     y = "MDS 1"
#   )
# 
# mds2 <- ggplot(mds_meta, aes(site2, y, fill = site2)) +
#   geom_jitter(alpha = 0.75, color = "darkgray", width = 0.25) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.75) +
#   scale_fill_manual(values = global_pal) +
#   theme_cowplot(12) +
#   theme(
#     # axis.text.x = element_text(angle = 90, hjust = 1),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     legend.position = "none"
#   ) +
#   labs(
#     x = "",
#     y = "MDS 2"
#   )
# 
# # cohort summary
# cohort_size <- pheno_global %>%
#   group_by(site2) %>%
#   tally()
# 
# cohort_plot <- ggplot(cohort_size, aes(site2, n, fill = site2)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = global_pal) +
#   scale_y_continuous(breaks = seq(0, 150, 25)) +
#   # scale_fill_brewer(palette = "Paired") +
#   theme_cowplot(12) +
#   theme(
#     legend.position = "none",
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#   ) + 
#   labs(
#     x = "",
#     y = "Participants"
#   )
#   # background_grid()
# 
# 
# ### shannon diversity global
# # find shannon diversity with vegdist
# global_S_rare <- rrarefy(t(global_S), min(colSums(global_S)))
# 
# shannon_div <- diversity(global_S_rare, index = "shannon")
# div <- data.frame("shannon_div" = shannon_div, "sample" = names(shannon_div))
# div_meta <- merge(div, pheno_global, by = "sample")
# 
# div_meta <- div_meta[complete.cases(div_meta), ]
# 
# # pval
# pvals <- compare_means(shannon_div ~ site2, data = div_meta, method = "wilcox.test", p.adjust.method = "fdr")
# 
# div_meta <- div_meta[complete.cases(div_meta), ]
# 
# shannon_global <- ggplot(div_meta, aes(site2, shannon_div, fill = site2)) + 
#   geom_jitter(alpha = 0.75, color = "darkgray", width = 0.25) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.75) +
#   # geom_boxplot(outlier.shape = NA) +
#   # geom_jitter(width = 0.25, pch = 21, alpha = 0.5) +
#   labs(x = "",
#        y = "Shannon Diversity",
#        fill="") +
#   scale_fill_manual(values = global_pal) +
#   scale_y_continuous(breaks = seq(1, 6, 1)) +
#   theme_cowplot(12) +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#     legend.position = "none"
#   )
# 
# ### plot multipanel
# col1 <- plot_grid(
#   cohort_plot + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")),
#   mds1 + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")),
#   mds2 + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")),
#   shannon_global + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")),
#   ncol = 1,
#   align = "v",
#   axis = "lbr",
#   labels = c("A", "C", "", "D"),
#   rel_heights = c(0.3, 0.2, 0.2, 0.3)
#   )
# 
# mds_scatter <- plot_grid(
#   global_mds + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")),
#   scatter_F + theme(legend.position = "none", plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")),
#   align = "v",
#   axis = "lbr",
#   ncol = 1,
#   labels = c("B", ""),
#   rel_heights = c(0.45, 0.55)
# )
# 
# plot_grid(
#   col1,
#   mds_scatter,
#   nrow = 1,
#   rel_widths = c(0.4, 0.6)
# )
# # ggsave("final_plots/supplementary/global_mds_shannon_v2.png", width = 9, height = 9)
# ggsave("final_plots/figure_3.png", width = 9, height = 10)
# 
# 
# ## find differential genera across geography with DESeq2
# count_data <- global_G[genefilter(global_G, pOverA(0.1, 1000)), ]
# 
# col_data <- pheno_global %>%
#   filter(sample %in% names(count_data)) %>%
#   column_to_rownames(var = "sample")
# 
# count_data <- count_data[, rownames(col_data)]
# 
# # all(rownames(col_data) == names(count_data))
# 
# col_data$geography <- ifelse(col_data$site %in% c("Tanzania", "Madagascar"), "nonwestern",
#                              ifelse(col_data$site == "South Africa", "za", "western"))
# col_data$geography <- factor(col_data$geography, levels = c("za", "nonwestern", "western"))
# 
# dds <- DESeqDataSetFromMatrix(countData = count_data,
#                               colData = col_data,
#                               design = ~ geography)
# 
# dds <- DESeq(dds, parallel = T)
# resultsNames(dds) # lists the coefficients
# 
# alpha <- 0.1
# 
# res_w <- results(dds, name="geography_western_vs_za", alpha = alpha)
# res_filt_w <- data.frame(res_w) %>%
#   rownames_to_column("feature") %>%
#   filter(padj < alpha) %>%
#   arrange(-abs(log2FoldChange)) %>%
#   mutate(
#     comparison = "western",
#     pos = log2FoldChange > 0
#     )
# 
# res_nw <- results(dds, name="geography_nonwestern_vs_za", alpha = alpha)
# res_filt_nw <- data.frame(res_nw) %>%
#   rownames_to_column("feature") %>%
#   filter(padj < alpha) %>%
#   arrange(-abs(log2FoldChange)) %>%
#   mutate(
#     comparison = "nonwestern",
#     pos = log2FoldChange > 0
#     )
# 
# res_all <- rbind(res_filt_w, res_filt_nw)
# res_all <- res_all %>%
#   filter(abs(log2FoldChange) > 2) %>%
#   select(feature, comparison, pos)
# res_all <- dcast(res_all, feature ~ comparison)
# 
# plot_features <- res_all %>%
#   filter(nonwestern == western) %>%
#   pull(feature)
# 
# ## plot
# deseq2_counts <- counts(dds, normalized=TRUE)
# 
# # add pseudocount
# deseq2_counts[deseq2_counts == 0] <- 0.01
# 
# counts_long <- melt(data.matrix(deseq2_counts))
# names(counts_long) <- c("genus", "sample", "rel_abundance")
# 
# counts_long_filt <- counts_long %>%
#   filter(genus %in% plot_features)
# 
# counts_meta <- merge(counts_long_filt, col_data %>% rownames_to_column("sample"), by = "sample")
# 
# counts_meta$geography <- factor(counts_meta$geography, levels = c("nonwestern", "za", "western"))
# 
# ggplot(counts_meta, aes(geography, rel_abundance, fill = geography)) +
#   geom_jitter(color = "darkgray", width = 0.25, alpha = 0.5) +
#   geom_boxplot(alpha = 0.5) +
#   scale_y_log10(label = comma) +
#   theme_cowplot() +
#   theme(
#     legend.position = "bottom",
#     legend.justification = "center",
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank()
#     ) +
#   facet_wrap(~ genus, scales = "free_y", ncol = 3) +
#   scale_fill_discrete(labels = c("Non-western", "South Africa", "Western")) +
#   labs(
#     x = "",
#     y = "DESeq2 normalized counts",
#     fill = ""
#   )
# ggsave("final_plots/supplementary/za_vs_western_nonwestern.png", width = 10, height = 12)

################################################################################
########## Figure 4: Reference databases are incomplete  #######################
################################################################################
# 
# ### A) classification plot
# # read data with unclassified reads
# unclass <- read.table("input_final/global_unclassified.txt", sep = "\t")
# names(unclass) <- c("sample", "unclassified_reads")
# 
# # rampelli
# unclass_rampelli <- read.table("input_final/global_unclassified_rampelli.txt", sep = "\t")
# names(unclass_rampelli) <- c("sample", "unclassified_reads")
# 
# unclass <- rbind(unclass, unclass_rampelli)
# 
# readcounts <- read.table("input_final/readcounts_global.txt", sep = "\t")
# names(readcounts) <- c("sample", "total_reads")
# 
# readcounts_rampelli <- read.table("input_final/readcounts_global_rampelli.txt", sep = "\t")
# names(readcounts_rampelli) <- c("sample", "total_reads")
# 
# readcounts <- rbind(readcounts, readcounts_rampelli)
# 
# # kracken/bracken results are in read pairs
# unclass$unclassified_reads <- unclass$unclassified_reads * 2
# 
# unclass$total_reads <- readcounts[match(unclass$sample, readcounts$sample), "total_reads"]
# 
# unclass <- unclass %>%
#   mutate(
#     classified_reads = total_reads - unclassified_reads,
#     classified_perc = classified_reads / total_reads
#     )
# 
# # merge pheno data
# unclass_pheno <- merge(pheno_global, unclass, by = "sample")
# 
# ## Plot by site with statistical tests
# # https://github.com/kassambara/ggpubr/wiki/Adding-Adjusted-P-values-to-a-GGPlot
# 
# my_comparisons <- list(
#   c("Bushbuckridge", "Soweto"),
#   # c("Bushbuckridge", "United States"),
#   c("Soweto", "United States")
# )
# 
# pvals <- unclass_pheno %>% pairwise_wilcox_test(classified_perc ~ site2, comparisons = my_comparisons)
# pvals <- pvals %>% add_y_position()
# pvals$y.position <- pvals$y.position + 2
# 
# read_class_plot <- ggplot(unclass_pheno, aes(site2, classified_perc * 100)) + 
#   # geom_boxplot(outlier.shape = NA, aes(fill = site2)) +
#   # geom_jitter(position = position_jitterdodge(jitter.width = 2), alpha = 0.5, pch = 21, aes(fill = site2)) +
#   geom_jitter(alpha = 0.75, color = "darkgray", width = 0.3) +
#   geom_boxplot(outlier.shape = NA, aes(fill = site2), alpha = 0.75) +
#   labs(x = "",
#        y = "Classified reads (%)",
#        fill="Site") +
#   scale_fill_manual(values = global_pal) +
#   scale_y_continuous(breaks = seq(0, 100, 10)) +
#   theme_cowplot() +
#   theme(
#     legend.position = "none",
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#     plot.margin = unit(c(0, 0.5, 0, 0.5), unit = "cm")
#   ) +
#   stat_pvalue_manual(pvals, label = "p.adj.signif", tip.length = 0.01, vjust = 0.5, hide.ns = T, y.position = c(100, 104))
# 
# 
# ##### B) k-mer comparisons
# # perform kmer nmds/mds
# kmer_mds <- function(k){
#   k_matrix_filt <- get(paste0("sourmash_k", k, "_track_abund"))
#   
#   k_dist <- as.dist(1 - k_matrix_filt)
#   
#   cmdscale(k_dist, eig = TRUE, x.ret = TRUE)
#   
#   isoMDS(k_dist)
# }
# 
# # plot kmer mds
# plot_kmer_mds <- function(mds, k){
#   
#   mds <- data.frame(mds$points)
#   mds$sample <- row.names(mds)
#   mds_pheno <- merge(mds, pheno_global, by = "sample", all.x = T, all.y = F)
#   
#   # randomize sample order for more even plotting
#   set.seed(42)
#   rows <- sample(nrow(mds_pheno))
#   mds_pheno <- mds_pheno[rows, ]
#   
#   kmer_mds <- ggplot(mds_pheno, aes(X1, X2, color = site2)) +
#     geom_point(size = 1.5, alpha = 0.75) +
#     scale_color_manual(values = global_pal, guide = guide_legend(nrow = 1)) +
#     labs(
#       # title = paste("k=", k, sep = ""),
#       x = "NMDS 1",
#       y = "NMDS 2",
#       color="",
#       fill="Site"
#       ) +
#     theme_cowplot() +
#     theme(
#       legend.position = "top",
#       plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), unit = "cm")
#       ) +
#     annotate("text", label = paste0("italic(k) ==", k), x = min(mds_pheno$X1) + 0.1, y = max(mds_pheno$X2), parse = T)
#   
#   kmer_mds + stat_ellipse(aes(color = site2), type = 't', size = 0.5, alpha = 0.75, show.legend = F)
# }
# 
# k21_mds <- kmer_mds("21")
# k31_mds <- kmer_mds("31")
# k51_mds <- kmer_mds("51")
# 
# k21 <- plot_kmer_mds(k21_mds, 21)
# k31 <- plot_kmer_mds(k31_mds, 31)
# k51 <- plot_kmer_mds(k51_mds, 51)
# 
# 
# ##### C) compare k-mer and bray pairwise distances
# # species bray-curtis vs kmer angular similarity
# for (method in c("bray", "jaccard")){
#   
#   suffix <- ifelse(method == "bray", "_track_abund", "")
#   
#   veg_dist <- vegdist(t(global_S_css), method = method)
#   dist_long <- melt_dist(as.matrix(veg_dist))
#   dist_long <- merge(dist_long, pheno_global[, -2], by.x = "iso1", by.y = "sample", all.x = T, all.y = F)
#   dist_long <- merge(dist_long, pheno_global[, -2], by.x = "iso2", by.y = "sample", all.x = T, all.y = F)
#   pair_dis_sp <- filter(dist_long, site2.x == site2.y)
#   pair_dis_sp$method <- "Species"
#   
#   k_long <- melt_dist(1 - get(paste0("sourmash_k31", suffix)))
#   k_long <- merge(k_long, pheno_global[, -2], by.x = "iso1", by.y = "sample", all.x = T, all.y = F)
#   k_long <- merge(k_long, pheno_global[, -2], by.x = "iso2", by.y = "sample", all.x = T, all.y = F)
#   pair_dis_k <- filter(k_long, site2.x == site2.y)
#   pair_dis_k$method <- "K-mer"
#   
#   pair_dis <- rbind(pair_dis_sp, pair_dis_k)
#   pair_dis$method <- factor(pair_dis$method, levels = c("Species", "K-mer"))
#   
#   assign(paste0("pair_dis", suffix), pair_dis)
# }
# 
# 
# ## plot all species / k-mer comparisons
# my_comparisons <- list(
#   c("Soweto", "Tanzania"),
#   c("Soweto", "Madagascar"),
#   c("Soweto", "Bushbuckridge"),
#   c("Soweto", "Sweden"),
#   c("Soweto", "United States")
# )
# dist_plot_all <- ggplot(pair_dis_track_abund, aes(x = site2.y, y = dist, fill = site2.y)) +
#   geom_boxplot(outlier.size = 0.75) +
#   facet_wrap(~method, scales = "free", labeller = labeller(groupwrap = label_wrap_gen(10))) +
#   scale_fill_manual(values = global_pal) +
#   labs(
#     x = "",
#     y = "Pairwise Distance"
#   ) +
#   theme_cowplot(12) +
#   theme(
#     legend.position = "none",
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#   ) +
#   stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", vjust = 0.75) +
#   background_grid()
# 
# 
# ## only za/Sweden
# my_comparisons <- list(
#   c("Bushbuckridge", "Sweden"),
#   c("Soweto", "Sweden")
# )
# 
# dat_sp <- pair_dis_track_abund %>%
#   filter(site2.y %in% c("Bushbuckridge", "Sweden"), method == "Species")
# wilcox.test(dist ~ site2.y, dat_sp, alternative = "l")
# # wilcox.test(dat[which(dat$site2.y == "Bushbuckridge"), "dist"], dat[which(dat$site2.y == "Sweden"), "dist"], alternative = "t")
# 
# dat_k <- pair_dis_track_abund %>%
#   filter(site2.y %in% c("Bushbuckridge", "Sweden"), method == "K-mer")
# wilcox.test(dist ~ site2.y, dat_k, alternative = "g")
# 
# ## are dist distributions normal?
# d <- dat %>% filter(site2.y == "Bushbuckridge") %>% pull(dist)
# ggdensity(d)
# ggqqplot(d)
# shapiro.test(sample(d, 5000))
# 
# dist_plot_za_swed <- ggplot(filter(pair_dis_track_abund, site2.y %in% c("Bushbuckridge", "Sweden")), aes(x = site2.y, y = dist, fill = site2.y)) +
#   geom_boxplot(outlier.size = 0.75) +
#   facet_wrap(~method, scales = "free", labeller = labeller(groupwrap = label_wrap_gen(10))) +
#   scale_fill_manual(values = global_pal[c(3, 5)]) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.075))) +
#   # expand_limits(y = 1.1) +
#   labs(
#     x = "",
#     y = "Pairwise Distance"
#   ) +
#   theme_cowplot(12) +
#   theme(
#     legend.position = "none",
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#   ) +
#   # stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test") +
#   stat_compare_means(comparisons = list(c("Bushbuckridge", "Sweden")), label = "p.signif", method = "wilcox.test", vjust = 0.5) +
#   background_grid()
# 
# 
# ## Madagascar/USA
# dist_plot_mad_usa <- ggplot(filter(pair_dis_track_abund, site2.y %in% c("Madagascar", "United States")), aes(x = site2.y, y = dist, fill = site2.y)) +
#   geom_boxplot(outlier.size = 0.75) +
#   facet_wrap(~ method, scales = "free", labeller = labeller(method = label_wrap_gen(10))) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.075))) +
#   scale_fill_manual(values = global_pal[c(2, 6)]) +
#   labs(
#     x = "",
#     y = "Pairwise Distance"
#   ) +
#   theme_cowplot(12) +
#   theme(
#     legend.position = "none",
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#   ) +
#   stat_compare_means(comparisons = list(c("Madagascar", "United States")), label = "p.signif", method = "wilcox.test", vjust = 0.5) +
#   # geom_signif(comparisons = list(c("Madagascar", "United States")), map_signif_level = TRUE, test = "wilcox.test") +
#   background_grid()
# 
# row1 <- plot_grid(get_legend(k21))
# 
# row2 <- plot_grid(
#   read_class_plot,
#   # k21 + theme(legend.position = "top", legend.text = element_text(size = 10), legend.spacing.x = unit(0.1, "cm")),
#   # k21 + theme(legend.position = "none"),
#   k31 + theme(legend.position = "none"),
#   # k51 + theme(legend.position = "none"),
#   # labels = c("A", "B", "C", "D"),
#   labels = c("A", "B"),
#   # nrow = 2
#   nrow = 1
#   )
# 
# row4 <- plot_grid(
#   dist_plot_za_swed,
#   dist_plot_mad_usa,
#   labels = c("D", "E"),
#   # rel_widths = c(3, 2),
#   align = "hv",
#   axis = "bt"
#   )
# 
# plot_grid(
#   row1,
#   row2,
#   dist_plot_all,
#   row4,
#   labels = c("", "", "C", ""),
#   ncol = 1,
#   rel_heights = c(0.05, 0.3, 0.35, 0.3)
#   )
# 
# ggsave("final_plots/figure_4.png", width = 8, height = 12.5)
# 
# 
# ### supplementary: Jaccard
# ## plot all species / k-mer comparisons
# ggplot(pair_dis, aes(x = site2.y, y = dist, fill = site2.y)) +
#   geom_boxplot(outlier.size = 0.75) +
#   facet_wrap(~method, scales = "free", labeller = labeller(groupwrap = label_wrap_gen(10))) +
#   scale_fill_manual(values = global_pal) +
#   labs(
#     x = "",
#     y = "Pairwise Jaccard Distance"
#   ) +
#   theme_cowplot(12) +
#   theme(
#     legend.position = "none",
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#   ) +
#   background_grid()
# ggsave("final_plots/misc/jaccard_dist_supp.png", width = 8, height = 5)

# ################################################################################
# ########## Figure S11: Taxonomy barplots of Nanopore samples ###################
# ################################################################################
# 
# nanopore_bins <- read.table("input_final/bins_by_sample.tsv", sep = "\t", header = T)
# illumina_bins <- read.table("input_final/MAGs_C27C29C33.txt", sep = "\t", header = T)
# 
# plot_nanopore <- function(count_data, sample, n_taxa, rank_label){
#   
#   # color palette for n taxa
#   myCols <- colorRampPalette(brewer.pal(12, "Paired"))
#   global_pal <- myCols(n_taxa)
#   global_pal <- sample(global_pal)
#   global_pal[n_taxa + 1] <- "gray"
#   
#   # find top n taxa
#   bracken_plot <- count_data %>%
#     select(sample) %>%
#     mutate(taxon = rownames(count_data)) %>%
#     melt(variable.name = "sample", value.name = "rel_abundance") %>%
#     arrange(desc(rel_abundance)) %>%
#     head(n_taxa)
#   
#   # add "other" column
#   bracken_plot <- rbind(bracken_plot, data.frame(taxon = "Other", sample = sample, rel_abundance =  1 - sum(bracken_plot$rel_abundance)))
#   
#   # remove annoying "miscellaneous" labels that take up space
#   bracken_plot$taxon <- gsub(" \\(miscellaneous\\)", "", bracken_plot$taxon)
#   # bracken_long <- melt(bracken_plot, id.vars = "taxon", variable.name = "sample", value.name = "rel_abundance")
#   
#   # merge in pheno date
#   # bracken_pheno <- merge(bracken_long, pheno, by = "sample")
#   # bracken_pheno <- merge(bracken_long, za_meta, by = "sample")
#   bracken_plot$label <- labels[match(bracken_plot$sample, labels$sample), "id"]
#   
#   # designate whether taxon was assembled into a MAG or nMAG
#   nmags <- nanopore_bins %>%
#     filter(Sample == sample) %>%
#     mutate(Final.Class = gsub("_", " ", Final.Class))
#   
#   mags <- illumina_bins %>%
#     filter(Sample == sample) %>%
#     mutate(Final.Class = gsub("_", " ", Final.Class))
#   
#   if (rank_label == "Genus"){
#     nmags <- nmags %>%
#       mutate(Final.Class = gsub("\\sbacterium|\\ssp.|\\s.+", "", Final.Class))
#       
#     mags <- mags %>%
#       mutate(Final.Class = gsub("\\sbacterium|\\ssp.", "", Final.Class))
#   }
#   
#   bracken_plot$Final.Class <- gsub("unclassified ", "", bracken_plot$taxon)
#   # bracken_plot$Final.Class <- gsub("_", "\\s", bracken_plot$Final.Class)
#   
#   bracken_plot$taxon <- ifelse(bracken_plot$Final.Class %in% mags$Final.Class, paste0(bracken_plot$taxon, "*"), as.character(bracken_plot$taxon))
#   bracken_plot$taxon <- ifelse(bracken_plot$Final.Class %in% nmags$Final.Class, paste0(bracken_plot$taxon, "†"), as.character(bracken_plot$taxon))
#   
#   # set factor level for correct plotting order
#   bracken_plot$taxon <- factor(bracken_plot$taxon, levels = unique(bracken_plot$taxon))
#   
#   # plot in order of decreasing relative abundance of desired taxon
#   # sorted <- bracken_pheno %>% filter(taxon == sort_by) %>% arrange(desc(rel_abundance)) %>% pull(label)
#   # bracken_pheno$label <- factor(bracken_pheno$label, levels = mixedsort(unique(bracken_pheno$label)))
#   # bracken_pheno$label <- factor(bracken_pheno$label, levels = sorted)
#   
#   ggplot(bracken_plot, aes(label, rel_abundance * 100, fill = taxon)) +
#     geom_bar(stat="identity") +
#     labs(
#       x = "",
#       y = "Relative Abundance (%)",
#       fill = rank_label
#     ) +
#     scale_fill_manual(values = global_pal) +
#     guides(fill = guide_legend(ncol=1, keywidth = 0.125, keyheight = 0.1, default.unit = "inch")) +
#     theme_cowplot() +
#     theme(
#       legend.text = element_text(size = 10),
#       axis.ticks.x = element_blank(),
#       axis.text.x = element_blank()
#     ) +
#     scale_y_continuous(limits = c(0, 100), expand = c(0, 0))
# }
# 
# nmag_samples <- c("C27", "C29", "C33")
# nmag_species <- lapply(nmag_samples, function(x){plot_nanopore(za_S_rel, x, 30, "Species")})
# nmag_genus <- lapply(nmag_samples, function(x){plot_nanopore(za_G_rel, x, 30, "Genus")})
# 
# row1 <- plot_grid(
#   nmag_species[[1]] + ggtitle("Bushbuckridge 105"),
#   nmag_genus[[1]],
#   nrow = 1,
#   align = "hv",
#   axis = "lrtb"
#   )
# 
# row2 <- plot_grid(
#   nmag_species[[2]] + ggtitle("Bushbuckridge 107"),
#   nmag_genus[[2]],
#   nrow = 1,
#   align = "hv",
#   axis = "lrtb"
# )
# 
# row3 <- plot_grid(
#   nmag_species[[3]] + ggtitle("Bushbuckridge 112"),
#   nmag_genus[[3]],
#   nrow = 1,
#   align = "hv",
#   axis = "lrtb"
# )
# 
# plot_grid(row1, row2, row3, labels = c("A", "B", "C"), align = "v", ncol = 1)
# ggsave("final_plots/supplementary/figure_S11_za_nanopore_taxa.png", width = 10, height = 15)
# 
# # check labels
# s <- "C33"
# nanopore_bins %>%
#   filter(Sample == s) %>%
#   arrange(Final.Class) %>%
#   pull(Final.Class) %>%
#   unique()
# 
# illumina_bins %>%
#   filter(Sample == s) %>%
#   arrange(Final.Class)
# 

################################################################################
########## Figure S5: MDS by ethnicity #########################################
################################################################################
# 
# ## ethnicity
# awi_data_f <- "input_final/awigen_data/awi-phase1.csv"
# awi_data <- read.csv(awi_data_f)
# 
# awi_ethnicity <- data.frame("ethnicity_name" = c("Zulu", "Xhosa", "Ndebele", "Sotho", "Venda", "Tsonga", 
#                                                  "Tswana", "BaPedi", "Zimbabwean", "Other", "Unknown", "Swati", "Other"),
#                             "ethnicity" = c(1:12, 40))
# 
# awi_eth <- awi_data[, c("Microbiome_link", "study_id", "site", "ethnicity")]
# awi_eth$ethnicity_name <- awi_ethnicity[match(awi_eth$ethnicity, awi_ethnicity$ethnicity), "ethnicity_name"]
# 
# awi_eth$site_name <- ifelse(awi_eth$site == "6", "Soweto", "Bushbuckridge")
# 
# eth_counts <- plyr::count(awi_eth, c("site_name", "ethnicity_name"))
# 
# ## nmds clustering by ethincity
# pheno <- merge(za_pheno, awi_eth, by.x = "study_id", by.y = "Microbiome_link")
# 
# # mds
# vare_dis <- vegdist(t(za_S_css), method = "bray")
# mds <- cmdscale(vare_dis, eig = TRUE, x.ret = TRUE)
# mds_values <- mds$points
# 
# # variation per axis
# mds_var_per <- round(mds$eig/sum(mds$eig) * 100, 1)
# 
# ## Plot
# mds_data <- data.frame(sample = rownames(mds_values),
#                        x = mds_values[,1],
#                        y = mds_values[,2])
# 
# # merge pheno data
# mds_meta <- merge(mds_data, pheno, by = "sample", all.x = T)
# mds_meta$ethnicity_name[is.na(mds_meta$ethnicity_name)] <- "Unknown"
# 
# mds_meta$ethnicity_name <- factor(
#   mds_meta$ethnicity_name,
#   levels = c(sort(unique(mds_meta$ethnicity_name[-which(mds_meta$ethnicity_name %in% c("Other", "Unknown"))])), "Other", "Unknown"))
# 
# ggplot(mds_meta, aes(x, y, color = ethnicity_name)) +
#   geom_point(size = 2, alpha = 0.75) +
#   scale_color_manual(values = c(brewer.pal(9, "Set1"), "#4A4A4A")) +
#   labs(x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
#        y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
#        color="Ethnicity") +
#   theme_cowplot()
# 
# ggsave("final_plots/supplementary/figure_S5_za_ethnicity_mds.png", height = 4, width = 5)

################################################################################
########## Misc supplementary analysis #########################################
################################################################################

########## Vangay et al data ##########
za_vangay <- merge(za_S, vangay_S, by = "row.names", all = T)
za_vangay[is.na(za_vangay)] <- 0
za_vangay <- za_vangay %>% column_to_rownames("Row.names")
za_vangay <- za_vangay[genefilter(za_vangay, pOverA(0.1, 10)), ]

za_vangay_meta <- rbind(
  data.frame("sample" = za_meta$sample, "site" = za_meta$site, "cohort" = "South Africa"),
  data.frame("sample" = vangay_meta$V1, "site" = vangay_meta$V44, "cohort" = "Vangay") %>% filter(sample %in% names(vangay_S))
)
rownames(za_vangay_meta) <- za_vangay_meta$sample

za_vangay <- za_vangay[, za_vangay_meta$sample]

# css norm
mr <- newMRexperiment(za_vangay, phenoData = AnnotatedDataFrame(za_vangay_meta))
p <- cumNormStatFast(mr)
mr_css <- cumNorm(mr, p = p)
za_vangay_css <- MRcounts(mr_css, norm = T, log = T)

# clr norm
x <- aldex.clr(za_vangay, mc.samples = 128, denom = "all", verbose = F, useMC = T)
za_vangay_clr <- sapply(getMonteCarloInstances(x), function(x) rowMeans(x))

### mds
vare_dis <- vegdist(t(za_vangay_css), method = "bray")

# # adonis
# meta <- za_meta %>% filter(sample %in% names(za_S_rel))
# meta <- meta[match(names(za_S_rel), meta$sample), ]
# adonis(vare_dis ~ site, data = meta)
# dispersion <- betadisper(vare_dis, group = meta$site)
# permutest(dispersion)

# calculate mds
mds <- cmdscale(vare_dis, eig = TRUE, x.ret = TRUE)

mds_values <- mds$points
# wa_scores <- wascores(mds_values, t(za_vangay_css[genefilter(za_vangay_css, pOverA(0.1, 1)), ]))
# wa_scores <- data.frame(sample = rownames(wa_scores),
#                         x = wa_scores[,1],
#                         y = wa_scores[,2])
# 
# # isolate taxa with strongest contribution to principal coordinate axes
# n_taxa <- 10
# wa_scores_1<- head(arrange(wa_scores, desc(abs(wa_scores$x))), n = n_taxa)
# wa_scores_2<- head(arrange(wa_scores, desc(abs(wa_scores$y))), n = n_taxa)
# wa_scores_final <- rbind(wa_scores_1, wa_scores_2)

# calculate percentage of variation that each MDS axis accounts for
mds_var_per <- round(mds$eig/sum(mds$eig) * 100, 1)

# plot
mds_data <- data.frame(sample = rownames(mds_values),
                       x = mds_values[,1],
                       y = mds_values[,2])

# merge pheno data
mds_meta <- merge(mds_data, za_vangay_meta, by = "sample")
mds_meta$site <- factor(mds_meta$site, levels = c("Bushbuckridge", "Soweto", "HmongThai", "Hmong1st", "Karen1st", "Control"))

mds_plot <- ggplot(mds_meta, aes(x, y, color = site, shape = cohort)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c(za_pal, brewer.pal(5, "Dark2")[c(1,2,4,5)])) +
  # scale_color_brewer(palette = "Dark2") +
  labs(x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
       y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
       color = "Population",
       shape = "Study"
       ) +
  theme_cowplot()

mds_plot + stat_ellipse(alpha = 0.75, show.legend = F)
ggsave("final_plots/misc/vangay_supp.png", width = 8, height = 6)


########## Microbes that correlate with age ##########
awi_pheno <- read.csv("input_final/pheno/awi-phase1.csv", header = T)

awi_pheno$sample <- za_meta[match(awi_pheno$Microbiome_link, za_meta$study_id), "sample"]

awi_pheno <- awi_pheno %>%
  filter(!is.na(sample)) %>%
  select(sample, age) %>%
  column_to_rownames(var = "sample")

awi_pheno <- data.frame(t(awi_pheno))

za_S_filt <- za_S_rel[, which(names(za_S_rel) %in% names(awi_pheno))]
za_S_filt <- za_S_filt[genefilter(za_S_filt, pOverA(p = 0.1, A = 0.001)), ]

za_S_age <- rbind(awi_pheno, za_S_filt)

cor_S <- cor(t(za_S_age), method = "spearman")

cor_S_long <- melt(cor_S)

cor_S_long %>%
  filter(Var1 == "age" & Var2 != "age") %>%
  arrange(-value) %>%
  head(20)

########## Correlate pheno data ##########
# convert character columns to factor
vars <- za_pheno %>% select(bp_sys_1, bp_dia_1, glucose_reading, bmi)
vars <- as.data.frame(unclass(vars))
vars <- vars[complete.cases(vars), ]
vars <- data.matrix(vars)

# calculate spearman cor
pheno_cor <- cor(vars, method = "spearman")

# cluster to order rows
hc <- hclust(dist(pheno_cor))
row_order <- hc$order
pheno_cor <- pheno_cor[row_order, row_order]

# plot
pheno_cor_long <- melt(pheno_cor, value.name = "cor")

cor_plot <- ggplot(pheno_cor_long, aes(Var2, Var1, fill = cor)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1)) +
  theme_cowplot(10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title = element_blank()
  )
ggsave("final_plots/misc/pheno_corr.png", width = 9, height = 8.5)

### cor pheno data with species
vars <- za_pheno %>% select(sample, bp_sys_1, bp_dia_1, glucose_reading, bmi)
# vars <- as.data.frame(unclass(vars))
# vars <- vars[complete.cases(vars), ]
# vars <- data.matrix(vars)

za_species_filt <- za_species_rel[genefilter(za_species_rel, pOverA(p = 0.10, A = 0.001)), ]

relab_S <- data.frame(t(za_species_filt)) * 100
cor_data <- merge(vars, relab_S, by.x = "sample", by.y = "row.names")
cor_data$sample <- NULL

cors <- cor(cor_data)
# hc <- hclust(dist(species_cor))
# row_order <- hc$order
# species_cor <- species_cor[row_order, row_order]

cor_data_long <- melt(cors, value.name = "cor")

filter(cor_data_long, Var1 == "glucose_reading" & abs(cor) > 0.3)
filter(cor_data_long, Var1 == "bmi" & abs(cor) > 0.3)
filter(cor_data_long, Var1 == "bp_sys_1" & abs(cor) > 0.3)


########## Global readcounts after preprocessing ##########

# read PAIRS
# this is outdated -- doesn't contain Rampelli data
global_rc <- read.table("input_final/readcounts/global_readcounts.tsv", sep = '\t', header = T)

# filter out extra samples
global_rc <- filter(global_rc, Sample %in% unique(pheno_global$sample))

counts <- readcounts[, c(1:3, 5, 7)]
colnames(counts) <- c("Sample", "Raw reads", "Trimmed reads", "Deduplicated reads", "Non-human reads")

counts_long <- melt(counts, id.vars = "Sample", variable.name = "step", value.name = "reads")
counts_long$reads_m <- (counts_long$reads / 1e6)


# ########## Supp figure: ZA read counts plot ##########
# ## read PAIRS
# za_readcounts <- read.table("input_final/readcounts/readcounts_preproc_za.tsv", sep = '\t', header = T)
# 
# # read PAIRS -- multiply by 2
# # filter out extra samples
# za_readcounts <- filter(za_readcounts, Sample %in% za_meta$sample) %>%
#   mutate(
#     raw_reads = raw_reads * 2,
#     trimmed_reads = trimmed_reads * 2,
#     deduplicated_reads = deduplicated_reads * 2,
#     host_removed_reads = host_removed_reads * 2,
#     orphan_reads = orphan_reads * 2,
#   )
# 
# za_counts <- za_readcounts[, c(1:3, 5, 7)]
# 
# ## stats for intro
# za_counts %>%
#   mutate(
#     Gb_raw = raw_reads * 150 / 1e9,
#     Gb_preproc = host_removed_reads * 150 / 1e9
#     ) %>%
#   summarise(
#     median_gb_raw = median(Gb_raw),
#     median_gb_raw = median(Gb_raw),
#     mean_gb_preproc = mean(Gb_preproc)
#   )
# 
# colnames(za_counts) <- c("Sample", "Raw reads", "Trimmed reads", "Deduplicated reads", "Non-human reads")
# 
# counts_long <- melt(za_counts, id.vars = "Sample", variable.name = "step", value.name = "reads")
# counts_long$reads_m <- (counts_long$reads / 1e6)
# 
# # plot readcounts
# rc <- ggplot(counts_long, aes(x=reads_m, fill=step)) +
#   geom_histogram(binwidth = 1) +
#   scale_x_continuous(breaks = seq(0, 100, 10)) +
#   facet_grid(step~., scales = "free_y") +
#   labs(
#     title = "Preprocessing Readcounts\n",
#     x = "\nReads (M)",
#     y = "Count\n",
#     fill = ""
#   ) +
#   theme_cowplot(14)
# 
# # ggsave("final_plots/misc/readcounts_preprocessing.png", rc, device = "png", height = 12, width = 12)
# 
# # save readcounts to supp file
# za_counts <- za_counts %>%
#   mutate(
#     label = labels[match(za_counts$Sample, labels$sample), "id"],
#     Sample = label
#     ) %>%
#   select(Sample, ends_with("reads"))
# 
# za_counts <- za_counts[mixedorder(za_counts$Sample), ]
# 
# ## do BBR/SWT differ in % host reads?
# z <- za_counts %>%
#   mutate(
#     human_removed = `Deduplicated reads` - `Non-human reads`,
#     human_rm_frac = human_removed / `Deduplicated reads`,
#     site = Hmisc::capitalize(gsub("_.+", "", Sample))
#     )
# 
# wilcox.test(human_rm_frac ~ site, z)
# 
# z %>%
#   group_by(site) %>%
#   summarise(mean_human_frac = mean(human_rm_frac))
# 
# a <- ggplot(z, aes(human_rm_frac * 100, fill = site)) +
#   geom_histogram(color = "white") +
#   scale_fill_manual(values = za_pal) +
#   theme_cowplot() +
#   theme(
#     legend.position = "top",
#     legend.justification = "center"
#   ) +
#   labs(
#     x = "Human reads removed after de-duplication (%)",
#     fill = "Site",
#     y = "Sample count"
#   ) +
#   background_grid()
# 
# b <- ggplot(z, aes(site, human_rm_frac * 100, fill = site)) +
#   geom_jitter(color = "darkgray", alpha = 0.75, width = 0.3) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#   scale_fill_manual(values = za_pal) +
#   theme_cowplot() +
#   theme(legend.position = "none") +
#   stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5) +
#   labs(
#     x = "",
#     y = "Human reads removed\nafter de-duplication (%)"
#   ) +
#   background_grid()
# 
# plot_grid(a, b, labels = c("A", "B"), ncol = 1)
# 
# ggsave("final_plots/supplementary/figure_SX_human_reads.png", width = 5, height = 9)
# 
# za_counts <- za_counts %>%
#   mutate("Human reads removed after de-duplication" = `Deduplicated reads` - `Non-human reads`)
# 
# write.table(za_counts, "final_tables/table_s1_readcounts.txt", sep = "\t", row.names = F, quote = F)

## stats
# raw reads
mean(za_counts$`Raw reads`) / 1e6
range(za_counts$`Raw reads`) / 1e6

# preprocessed
mean(za_counts$`Non-human reads`) / 1e6
range(za_counts$`Non-human reads`) / 1e6


########## Supp figure: correlation of species in ZA data ##########
za_species_filt <- za_species_rel[genefilter(za_species_rel, pOverA(p = 0.10, A = 0.0001)), ]
species_cor <- cor(t(za_species_filt))

hc <- hclust(dist(species_cor))
row_order <- hc$order
species_cor <- species_cor[row_order, row_order]

species_cor_long <- melt(species_cor, value.name = "cor")
ggplot(species_cor_long, aes(Var2, Var1, fill = cor)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue",
                       high = "red",
                       mid = "white",
                       midpoint = 0,
                       limit = c(-1,1)
                       ) +
  theme_cowplot(6) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )


########## Ben SGB mapping ##########

mapping_f <- "/Users/Fiona/Desktop/za_sgb_mapping.txt"
mapping <- read.table(mapping_f, sep = "\t", header = T)
mapping_long <- melt(mapping, measure.vars = c("Oct2018.1", "Oct2018.Segata"), 
                     id.vars = "name", value.name = "percent", variable.name = "database")
mapping_long$db_name <- ifelse(mapping_long$database == "Oct2018.1", "Reference genomes", "Reference genomes + SGBs")

gb <- ggplot(mapping_long, aes(x = db_name, y = percent, fill = db_name)) +
  geom_boxplot() +
  labs(title="",
       x = "Database",
       y = "Mapping Reads (%)\n",
       fill="") +
  scale_fill_manual(values = c("#5775A3", "#6A6A6A")) +
  theme_minimal() +
  my_thm +
  theme(legend.position = "none")


########## Supplementary figures: analysis by BMI ##########
shannon_div <- diversity(t(za_species_rel), index = "shannon")
div <- data.frame("shannon_div" = shannon_div, "sample" = names(shannon_div))
div_meta <- merge(div, za_pheno, by = "sample")
# div_meta <- div_meta[complete.cases(div_meta), ]
plot_data <- div_meta[, c(1:4, 9)]

# remove C9 -- incorrect BMI
plot_data <- filter(plot_data, sample != "C9")
plot_data <- plot_data[complete.cases(plot_data), ]
plot_data$group <- ifelse(plot_data$bmi >= 30, "Obese", "Overweight")
plot_data$group <- ifelse(plot_data$bmi < 25, "Healthy", plot_data$group)
plot_data$group <- factor(plot_data$group, levels = c("Healthy", "Overweight", "Obese"))

# statistical tests
pvals <- compare_means(shannon_div ~ group, data = plot_data, group.by = "site", method = "wilcox.test", p.adjust.method = "BH")
pvals$p.signif <- ifelse(pvals$p.adj < 0.05, "*", "ns")
pvals$p.signif <- ifelse(pvals$p.adj < 0.01 & pvals$p.adj >= 0.001, "**", pvals$p.signif)
pvals$p.signif <- ifelse(pvals$p.adj < 0.001, "***", pvals$p.signif)

# plot shannon diversity by site
bmi_div <- ggplot(plot_data, aes(x=group, y=shannon_div)) + 
  geom_violin(aes(fill = group), trim = F) +
  stat_summary(fun.data=mean_sdl, aes(group=group), position=position_dodge(.9), geom="pointrange", color="black") +
  facet_grid(. ~ site, scales = "free") +
  labs(x = "\nBody Mass Index Group",
       y = "Shannon Diversity\n",
       fill="") +
  theme_bw() +
  scale_fill_brewer(palette = "Blues") +
  my_thm +
  stat_pvalue_manual(data = pvals, label = "p.signif", y.position = c(6.6, 7.1, 7.6), size = 5) +
  stat_compare_means(method = "kruskal.test", label.y = 8, size = 5)

ggsave("rplots/za_bmi_div.png", bmi_div, device = "png", width = 12, height = 8)


########## ZA NMDS by bmi ##########
# remove C9 -- incorrect BMI
plot_data <- filter(mds_pheno, sample != "C9")
plot_data <- plot_data[, c(1:5, 10)]
plot_data <- plot_data[complete.cases(plot_data), ]
plot_data$group <- ifelse(plot_data$bmi >= 30, "Obese", "Overweight")
plot_data$group <- ifelse(plot_data$bmi < 25, "Healthy", plot_data$group)
plot_data$group <- factor(plot_data$group, levels = c("Healthy", "Overweight", "Obese"))

nmds_bmi <- ggplot(plot_data, aes(x=X1, y=X2, fill=group)) +
  geom_point(size=4, color = "black", pch = 21) +
  # scale_shape_manual(values = c(21, 24)) +
  theme_bw() +
  facet_wrap(~ site, nrow = 2) +
  scale_fill_brewer(palette = "Blues") +
  scale_color_brewer(palette = "Blues") +
  labs(x = "NMDS1",
       y = "NMDS2",
       fill="") +
  xlim(-0.6, 0.6) +
  ylim(-0.6, 1) +
  my_thm

# nmds_bmi + stat_ellipse(aes(color = site), type='t', size = 1, show.legend = F)

ggsave("rplots/za_bmi_nmds_.png", nmds_bmi, device = "png", height = 12, width = 11)


########## World map ##########
# library(maps)
# 
# # https://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html
# world <- map_data("world")
# 
# # map_pal <- c("#A50026", "#D73027", "#FDAE61", "#74ADD1", "#4575B4", "#313695")
# # map_sites <- c("Madagascar", "Tanzania", "South Africa", "Finland", "Sweden", "USA")
# 
# map_pal <- c("#E3211C", "#F89897", "#1F78B4", "#A5CEE3")
# map_sites <- c("Tanzania", "Madagascar", "Sweden", "USA")
# 
# world$region_color <- ifelse(world$region %in% map_sites, world$region, NA)
# world$region_color <- factor(world$region_color, levels = map_sites)
# 
# ggplot(data = world) + 
#   geom_polygon(aes(x = long, y = lat, fill = region_color, group = group)) + 
#   coord_fixed(1.3) +
#   scale_fill_manual(values = map_pal, na.value = "#DDDDDD") +
#   guides(fill=FALSE) +  # do this to leave off the color legend
#   theme_void()
