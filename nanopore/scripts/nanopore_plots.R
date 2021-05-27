library(ggplot2)
library(dplyr)
library(gplots)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(scales)

setwd("/Users/dylanmaghini/scg4/projects/southafrica/nanopore/publicationanalyses/000.cleananalysis")

#### LOAD BINNING DATA ####
# Read in nanopore MAGs
nmags <- read.table("input_final/binfinal_C27C29C33.tsv", sep = '\t', header = T)
nmags$Type = factor("nMAG")
nmags <- nmags %>% select("Sample", "Bin", "Completeness", "Contamination", "tRNA", "rna.16S", "rna.5S", "rna.23S", "N50", "Type", "Total.length", "Coverage")

# Read in Illumina MAGs
mags <- read.csv("input_final/binning_table_all_full.tsv", sep = '\t', header = T)
mags$Type = factor("MAG")
mags <- mags %>% select("Sample", "Bin", "Completeness", "Contamination", "tRNA", "rna.16S", "rna.5S", "rna.23S", "N50", "Type", "Total.length", "Coverage")

# group nanopore and Illumina MAGs into single dataframe
all <- rbind(nmags, mags)

# add a Quality column based on completeness, contamination, number of tRNAs, presence of 16S, 23S, 5S sequences
all <- mutate(all, Quality = ifelse(Completeness > 90 & Contamination < 5 & tRNA >= 18 & rna.16S > 0 & rna.23S > 0 & rna.5S > 0, "High-quality", ifelse(Completeness >= 50 & Contamination <10, "Medium-quality", "Low-quality")))
all$Quality <- factor(all$Quality, levels = c("Low-quality", "Medium-quality", "High-quality"))
all <- mutate(all, ContaminationLog = ifelse(Contamination > 0, Contamination, 0.0001))
all <- mutate(all, N50overSize = N50/Total.length)

# read in table matching sample IDs with study site
sampletosite <- read.table("input_final/za_sample_codes.csv", sep = ",", header = T)
all <- merge(all, sampletosite, by="Sample")

# load color palettes
greens <- c("#99d8c9", "#41ae76", "#006d2c")
oranges <- c("#fee391", "#fe9929", "#cc4c02")

point <- format_format(big.mark = ",", decimal.mark = ".", scientific = FALSE)

#### BINNING STATISTICS ########

# Illumina versus nanopore completeness boxplot
a <- ggplot(all %>% filter(Bin!="bin.unbinned") %>% filter(Bin!="unbinned"), aes(x=reorder(Type), y=Completeness, fill=interaction(Quality, Type))) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.75, size = 1, color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  #geom_boxplot() +
  ylab("Completeness (%)") +
  ylim(0,100) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank()) +
  scale_fill_manual(values=c(greens,oranges)) +
  theme(text = element_text(size = 15))

# Illumina versus nanopore contamination boxplot
b <- ggplot(all %>% filter(Bin!="bin.unbinned") %>% filter(Bin!="unbinned"), aes(x=reorder(Type), y=ContaminationLog, fill=interaction(Quality, Type))) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.75, size = 1, color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  #geom_boxplot() +
  ylab("Contamination (%)") +
  scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100), labels = c("0", "", "0.01", "", "1", "", "100")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank()) +
  scale_fill_manual(values=c(greens,oranges)) +
  annotation_logticks(sides="l") +
  theme(text = element_text(size = 15))

# Illumina versus nanopore N50 boxplot
c <- ggplot(all %>% filter(Bin!="bin.unbinned") %>% filter(Bin!="unbinned"), aes(x=reorder(Type), y=N50, fill=interaction(Quality, Type))) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.75, size = 1, color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  ylab("N50 (bp)") +
  scale_y_log10(labels=point) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank()) +
  scale_fill_manual(values=c(greens,oranges), labels=c("Low Quality nMAG", "Medium Quality nMAG", "High Quality nMAG", "Low Quality MAG", "Medium Quality MAG", "High Quality MAG")) +
  annotation_logticks(sides="l") +
  theme(text = element_text(size = 15)) + 
  theme(legend.title=element_blank())

d <- ggplot(all %>% filter(Bin!="bin.unbinned") %>% filter(Bin!="unbinned"), aes(x=reorder(Type), fill=interaction(Quality, Type))) +
  geom_bar(width=0.6, position = position_dodge(width=0.8), color="black", alpha=0.75) +
  theme_bw() + 
  ylab("Number of Genomes") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank()) +
  scale_fill_manual(values=c(greens,oranges), labels=c("Low Quality nMAG", "Medium Quality nMAG", "High Quality nMAG", "Low Quality MAG", "Medium Quality MAG", "High Quality MAG")) +
  theme(text = element_text(size = 15), legend.title=element_blank())

prow <- plot_grid(d + theme(legend.position = "none"), 
                  a + theme(legend.position = "none"), 
                  b + theme(legend.position = "none"), 
                  c + theme(legend.position = "none"), rel_widths = c(1, 1, 1, 1), nrow = 2, ncol=2, labels = c("A", "B", "C", "D"))
abcdlegend <- get_legend(c + theme(legend.title = element_blank(), legend.box.margin = margin(0,0,0,0)))

plot_grid(prow, abcdlegend, nrow = 1, rel_widths = c(2.8, .8))
ggsave("plots/supp_binning_plots_1.jpg", units="in", width=10, height =7, dpi=500)

# Illumina versus nanopore, coverage versus N50
e <- ggplot(all %>% filter(Sample %in% c("C27", "C29", "C33"), Bin != "bin.unbinned", Bin != "unbinned"), aes(x=Coverage, y=N50, color=interaction(Quality, Type))) +
  geom_point() +
  ylab("N50 (bp)") +
  xlab("Coverage") +
  #scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
  scale_y_log10(labels=point) +
  scale_x_log10(labels=point) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line(colour = "black"), legend.title = element_blank()) +
  scale_color_manual(values=c(greens,oranges), labels=c("Low Quality nMAG", "Medium Quality nMAG", "High Quality nMAG", "Low Quality MAG", "Medium Quality MAG", "High Quality MAG")) +
  annotation_logticks(sides="lb") +
  theme(text = element_text(size = 15))

f <- ggplot(all %>% filter(Sample %in% c("C27", "C29", "C33"), Bin != "bin.unbinned", Bin != "unbinned"), aes(x=Total.length, y=N50, color=interaction(Quality, Type))) +
  geom_abline(intercept=0, slope = 1, color = 'black', alpha = 0.3) +
  geom_point() +
  ylab("N50 (bp)") +
  xlab("Genome Size (bp)") +
  #scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
  scale_y_log10(labels=point) +
  #scale_x_log10(labels=point) +
  scale_x_log10(breaks = c(1000000, 10000000), labels = point) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line(colour = "black"), legend.title = element_blank()) +
  scale_color_manual(values=c(greens,oranges), labels=c("Low Quality nMAG", "Medium Quality nMAG", "High Quality nMAG", "Low Quality MAG", "Medium Quality MAG", "High Quality MAG")) +
  annotation_logticks(sides="lb") +
  theme(text = element_text(size = 15))
f

prow <- plot_grid(e + theme(legend.position = "none"), 
                  f + theme(legend.position = "none"), rel_widths = c(1, 1), nrow = 1, labels = c("A", "B"))
leg <- get_legend(f + theme(legend.title = element_blank(), legend.box.margin = margin(0,0,0,0)))

plot_grid(prow, leg, nrow = 1, rel_widths = c(2, .7))
ggsave("plots/supp_binning_plots_2.jpg", units="in", width=9, height = 4, dpi=500)


#### MAG GC CONTENT COMPARISON ####

# load GC content for nanopore MAGs from each sample
gc33 <- read.table("input_final/C33_all_gc.txt", sep = '\t', header = F)
gc27 <- read.table("input_final/C27_all_gc.txt", sep = '\t', header = F)
gc29 <- read.table("input_final/C29_all_gc.txt", sep = '\t', header = F)
names(gc33) = c("Bin", "GC")
gc33 <- mutate(gc33, Sample = "C33")
names(gc27) = c("Bin", "GC")
gc27<- mutate(gc27, Sample = "C27")
names(gc29) = c("Bin", "GC")
gc29 <- mutate(gc29, Sample = "C29")
gc <- rbind(gc27, gc29, gc33)
gc <- mutate(gc, Type = "nMAG")

# load GC content for in Illumina MAGs for these three samples
igc33 <- read.table("input_final/C33.tsv", sep = '\t', header = F)
igc27 <- read.table("input_final/C27.tsv", sep = '\t', header = F)
igc29 <- read.table("input_final/C29.tsv", sep = '\t', header = F)
names(igc33) = c("Bin", "GC")
igc33 <- mutate(igc33, Sample = "C33")
names(igc27) = c("Bin", "GC")
igc27<- mutate(igc27, Sample = "C27")
names(igc29) = c("Bin", "GC")
igc29 <- mutate(igc29, Sample = "C29")
igc <- rbind(igc27, igc29, igc33)
igc <- mutate(igc, Type = "MAG")

# Combine long and short read data
allgc <- rbind(igc, gc)
all <- merge(all, allgc, by=c("Bin", "Sample", "Type"))

# Plot nanopore MAG GC content (MAGs with high contiguity, facet by quality)
g <- ggplot(all %>% filter(N50 > 1000000, Type=="nMAG"), aes(x=Quality, y=GC, fill = Quality)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.75, size = 1, color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  theme_bw() +
  ylim(c(0.1, 0.70)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank()) +
  scale_fill_manual(values=c(greens)) +
  theme(text = element_text(size = 15), legend.title=element_blank()) +
  stat_compare_means(comparisons = list(c("Low-quality", "High-quality")), method="wilcox.test", label="p.signif")

# Plot Illumina and nanopore bin GC content, faceted by quality
h <- ggplot(all, aes(x=Quality, y=GC, fill = interaction(Quality, reorder(Type)))) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.75, size = 1, color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  theme_bw() +
  ylim(c(0.1, 0.70)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank()) +
  scale_fill_manual(values=c(oranges, greens), labels=c("Low Quality MAG", "Medium Quality MAG", "High Quality MAG", "Low Quality nMAG", "Medium Quality nMAG", "High Quality nMAG")) +
  theme(text = element_text(size = 15),  legend.title=element_blank())

prow <- plot_grid(h + theme(legend.position = "none"), 
                  g + theme(legend.position = "none"), rel_widths = c(1, 1), nrow = 1, ncol=2, labels = c("A", "B"))
newleg <- get_legend(h + theme(legend.title = element_blank(), legend.box.margin = margin(0,0,0,0)))

plot_grid(prow, newleg, nrow = 1, rel_widths = c(2.8, .7))
ggsave("plots/supp_binning_plots_gcall.jpg", units="in", width=11, height =4, dpi=500)

#### READ GC CONTENT COMPARISON ####

# load subsampled Illumina read GC content
reads_ic27 <- read.table('input_final/C27_gc_subsample.tsv', sep = '\t', header = 0, row.names = NULL)
reads_ic29 <- read.table('input_final/C29_gc_subsample.tsv', sep = '\t', header = 0, row.names = NULL)
reads_ic33 <- read.table('input_final/C33_gc_subsample.tsv', sep = '\t', header = 0, row.names = NULL)

# load subsampled nanopore read GC content
reads_nc27<- read.table('input_final/C27_gc_subsample.tsv', sep = '\t', header = 0, row.names = NULL)
reads_nc29<- read.table('input_final/C29_gc_subsample.tsv', sep = '\t', header = 0, row.names = NULL)
reads_nc33<- read.table('input_final/C33_gc_subsample.tsv', sep = '\t', header = 0, row.names = NULL)

# format dataframes
reads_ic27 <- cbind(reads_ic27, "iC27")
reads_ic29 <- cbind(reads_ic29, "iC29")
reads_ic33 <- cbind(reads_ic33, "iC33")
reads_nc27 <- cbind(reads_nc27, "nC27")
reads_nc29 <- cbind(reads_nc29, "nC29")
reads_nc33 <- cbind(reads_nc33, "nC33")
names(reads_ic27)[2] = "GC"
names(reads_ic27)[3] = "Sample"
names(reads_ic29)[2] = "GC"
names(reads_ic29)[3] = "Sample"
names(reads_ic33)[2] = "GC"
names(reads_ic33)[3] = "Sample"
names(reads_nc27)[1] = "GC"
names(reads_nc27)[2] = "Sample"
names(reads_nc29)[1] = "GC"
names(reads_nc29)[2] = "Sample"
names(reads_nc33)[1] = "GC"
names(reads_nc33)[2] = "Sample"
reads_ic27 <- reads_ic27 %>% select(GC, Sample)
reads_ic29 <- reads_ic29 %>% select(GC, Sample)
reads_ic33 <- reads_ic33 %>% select(GC, Sample)

# combine short read and long read dataframes
all_comparos <- rbind(reads_ic27, reads_ic29, reads_ic33, reads_nc27, reads_nc29, reads_nc33)
all_comparos$Sample = factor(all_comparos$Sample, levels = c("iC27", "iC29", "iC33", "nC27", "nC29", "nC33"))


# C27 GC plot
a <- ggplot(all_comparos[all_comparos$Sample == "iC27" | all_comparos$Sample == "nC27",], aes(x=GC, fill = Sample, alpha = 0.3)) + 
  guides(alpha = "none") +
  geom_density() + 
  theme_bw() + 
  ggtitle("Bushbuckridge 105") + 
  xlab("GC Content") + 
  xlim(c(0,1)) +
  scale_fill_manual(values = c("grey", "black"), labels = c("Illumina", "Nanopore")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size=17), text = element_text(size = 15), legend.title = element_blank())
a

# C29 GC plot
b <- ggplot(all_comparos[all_comparos$Sample == "iC29" | all_comparos$Sample == "nC29",], aes(x=GC, fill = Sample, alpha = 0.3)) + 
  guides(alpha = "none") +
  geom_density() + 
  theme_bw() + 
  ggtitle("Bushbuckridge 107") +
  xlab("GC Content") + 
  xlim(c(0,1)) +
  scale_fill_manual(values = c("grey", "black"), labels = c("Illumina", "Nanopore")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size=17), text = element_text(size = 15), legend.title = element_blank())
b

# C33 GC plot
c <- ggplot(all_comparos[all_comparos$Sample == "iC33" | all_comparos$Sample == "nC33",], aes(x=GC, fill = Sample, alpha = 0.3)) + 
  guides(alpha = "none") +
  geom_density() + 
  theme_bw() + 
  ggtitle("Bushbuckridge 112") +
  xlab("GC Content") + 
  xlim(c(0,1)) +
  scale_fill_manual(values = c("grey", "black"), labels = c("Illumina", "Nanopore")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(size=17), text = element_text(size = 15), legend.title = element_blank())
c

#
prow <- plot_grid(a + theme(legend.position = "none"),
                  b + theme(legend.position = "none"),
                  c + theme(legend.position = "none"), rel_widths = c(1, 1, 1), nrow = 1)
abclegend <- get_legend(c + theme(legend.title = element_blank(), legend.box.margin = margin(0,0,0,0)))

plot_grid(prow, abclegend, nrow = 1, rel_widths = c(3, 0.4))
ggsave("plots/supp_gcreads.jpg", units="in", width=11, height = 3, dpi=500)

#### PAIRED nMAGs AND iMAGs ####

# read in and format mobile element data tables
phage <- read.table("input_final/vibrant.tsv", sep = '\t', header = F)
mobilegenes <- read.table("input_final/mobilegenes.tsv", sep = '\t', header = F)
abx <- read.table("input_final/abx.tsv", sep = '\t', header = F)
names(mobilegenes) <- c("LBin", "Sbin", "Transposases", "Recombinases", "Integrases")
names(phage) <- c("LBin", "Sbin", "Phages")
names(abx) <- c("LBin", "Sbin",  "Abx")
total <- merge(mobilegenes, phage, by="LBin")
total <- merge(total, abx, by="LBin")
total <- mutate(total, filler="na") %>% filter(LBin != "nC33_bin.233.fa") # filter out redundant bin

# plot additional transposases in each nMAG
a <- ggplot(total, aes(x=filler, y=Transposases)) +
  geom_jitter(alpha=0.5, height=0) +
  geom_boxplot(outlier.shape = NA, fill = "grey", alpha=0.5) +
  ggtitle("Transposases") +
  ylab("Additional Elements Per\nNanopore MAG") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(), axis.text.x=element_blank()) +
  theme(text = element_text(size = 14), legend.title=element_blank(), plot.title = element_text(hjust = 0.5, size =10))

# plot additional recombinases in each nMAG
b <- ggplot(total, aes(x=filler, y=Recombinases)) +
  geom_jitter(alpha=0.5, height=0) +
  geom_boxplot(outlier.shape = NA, fill = "grey", alpha=0.5) +
  ggtitle("Recombinases") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x=element_blank()) +
  theme(text = element_text(size = 14), legend.title=element_blank(),  plot.title = element_text(hjust = 0.5, size= 10))

# plot additional phages in each nMAG
c <- ggplot(total, aes(x=filler, y=Phages)) +
  geom_jitter(alpha=0.5, height=0) +
  geom_boxplot(outlier.shape = NA, fill="grey", alpha=0.5) +
  ggtitle("Phages") +
  theme_bw() +
  scale_y_continuous(breaks=c(0, 5, 10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x=element_blank()) +
  theme(text = element_text(size = 14), legend.title=element_blank(),  plot.title = element_text(hjust = 0.5, size = 10))

# plot additional Abx resistance genes in each nMAG
d <- ggplot(total, aes(x=filler, y=Abx)) +
  geom_jitter(alpha=0.5, height=0) +
  geom_boxplot(outlier.shape = NA, fill="grey", alpha=0.5) +
  ggtitle("Abx Resistance Genes") +
  scale_y_continuous(breaks = c(0,1,2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x=element_blank()) +
  theme(text = element_text(size = 14), legend.title=element_blank(),  plot.title = element_text(hjust = 0.5, size = 10))

prow <- plot_grid(a + theme(legend.position = "none"),
                  NULL,
                  b + theme(legend.position = "none"), 
                  NULL,
                  c + theme(legend.position = "none"), 
                  NULL,
                  d + theme(legend.position = "none"), rel_widths = c(1.3,0.1,1,0.1,1,0.1, 1), nrow = 1, ncol=7)
ggsave("plots/direct_plots.jpg", units="in", width=8, height = 3, dpi=800)
