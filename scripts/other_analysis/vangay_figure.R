library(cowplot)
library(here)
library(metagenomeSeq)
library(tidyverse)
library(vegan)

source(here("scripts/load_data.R"))

## load vangay data ----
vangay_meta <- read.table(here("input_final/vangay/vangay_metadata.txt"),
                          sep = "\t", header = F)
bioproj <- read.table(here("input_final/vangay/PRJEB28687_wgs.txt"),
                      sep = "\t", header = F)

vangay_S <- read.table(here("input_final/vangay/bracken_species_reads.txt"),
                       sep = "\t", quote = "", comment.char = "")

names(vangay_S) <- bioproj[match(names(vangay_S), bioproj$V4), "V9"]
names(vangay_S) <- gsub("wgs\\.", "", names(vangay_S))
names(vangay_S) <- gsub("_", "\\.", names(vangay_S))

vangay_S <- vangay_S[rev(order(rowMeans(vangay_S))), ]
vangay_S_rel <- sweep(vangay_S, 2, colSums(vangay_S), "/")

## za + vangay figure ----
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
za_vangay_css <- MRcounts(mr_css, norm = T, log = F)

# clr norm
# x <- aldex.clr(za_vangay, mc.samples = 128, denom = "all", verbose = F, useMC = T)
# za_vangay_clr <- sapply(getMonteCarloInstances(x), function(x) rowMeans(x))

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
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(values = c(za_pal, brewer.pal(5, "Dark2")[c(1,2,4,5)])) +
  # scale_color_brewer(palette = "Dark2") +
  labs(x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
       y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
       color = "Population",
       shape = "Study"
  ) +
  theme_cowplot()

mds_plot +
  stat_ellipse(alpha = 0.75, show.legend = F) +
  background_grid()

ggsave("final_plots/misc/vangay_supp.png", width = 8, height = 6)

## global + za + vangay figure ----
za_vangay <- merge(global_S, vangay_S, by = "row.names", all = T)
za_vangay[is.na(za_vangay)] <- 0
za_vangay <- za_vangay %>% column_to_rownames("Row.names")
za_vangay <- za_vangay[genefilter(za_vangay, pOverA(0.1, 10)), ]

za_vangay_meta <- rbind(
  data.frame("sample" = pheno_global$sample, "site" = pheno_global$site2,
             "cohort" = pheno_global$site),
  data.frame("sample" = vangay_meta$V1, "site" = vangay_meta$V44,
             "cohort" = "Vangay") %>% filter(sample %in% names(vangay_S))
)
rownames(za_vangay_meta) <- za_vangay_meta$sample

za_vangay <- za_vangay[, za_vangay_meta$sample]

# css norm
mr <- newMRexperiment(za_vangay, phenoData = AnnotatedDataFrame(za_vangay_meta))
p <- cumNormStatFast(mr)
mr_css <- cumNorm(mr, p = p)
za_vangay_css <- MRcounts(mr_css, norm = T, log = F)

# clr norm
# x <- aldex.clr(za_vangay, mc.samples = 128, denom = "all", verbose = F, useMC = T)
# za_vangay_clr <- sapply(getMonteCarloInstances(x), function(x) rowMeans(x))

### mds
vare_dis <- vegdist(t(za_vangay_css), method = "bray")

# calculate mds
mds <- cmdscale(vare_dis, eig = TRUE, x.ret = TRUE)

mds_values <- mds$points

# calculate percentage of variation that each MDS axis accounts for
mds_var_per <- round(mds$eig/sum(mds$eig) * 100, 1)

# plot
mds_data <- data.frame(sample = rownames(mds_values),
                       x = mds_values[,1],
                       y = mds_values[,2])

# merge pheno data
mds_meta <- merge(mds_data, za_vangay_meta, by = "sample")
mds_meta$site <- factor(mds_meta$site,
                        levels = c("Tanzania", "Madagascar", "Bushbuckridge",
                                   "Soweto", "Sweden", "United States", "HmongThai",
                                   "Hmong1st", "Karen1st", "Control"))

mds_meta$alpha_lvl <- ifelse(mds_meta$cohort == "Vangay", "Vangay et al.",
                             "Present study")

mds_plot <- ggplot(mds_meta, aes(x, y, color = site, alpha = alpha_lvl,
                                 shape = alpha_lvl)) +
  geom_point(size = 2) +
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_color_manual(values = c(global_pal, brewer.pal(6, "Dark2")[c(1,2,4,6)])) +
  labs(x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
       y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
       color = "Population",
       shape = "Study"
  ) +
  theme_cowplot() +
  guides(alpha = F)

mds_plot +
  # stat_ellipse(alpha = 0.75, show.legend = F) +
  background_grid()

ggsave("final_plots/misc/vangay_supp.png", width = 8, height = 6)
