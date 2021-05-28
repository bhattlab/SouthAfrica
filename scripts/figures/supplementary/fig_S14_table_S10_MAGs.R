library(cowplot)
library(ggplot2)
library(ggpubr)
library(here)
library(tidyverse)

## params ----
global_pal <- c("#E3211C", "#F89897", "#6A3D9A", "#CAB2D6", "#1F78B4", "#A5CEE3")
za_pal <- global_pal[3:4]

## read input data ----
labels <- read.table(here("input_final/za_labels.tsv"), sep = "\t", header = T)

## sample metadata
za_meta <- readRDS(here("rds/za_meta.rds"))

## genomesearch output
gsearch <- read.table(here("input_final/uhgg_ANI_compare/perident_table_uhgg_formatted.tsv"),
                      sep = "\t", header = T,
                      col.names = c("Sample", "Bin", "gsearch.ANI",
                                    "gsearch.Phylum", "gsearch.Species"))

## fastani output
fastani <- read.table(here("input_final/uhgg_ANI_compare/fastani_top.txt"),
                      sep = "\t", header = F,
                      col.names = c("Sample", "Bin","UHGG_genome", "fastani.ANI",
                                    "fragment_mappings", "query_fragments", "AF"))

# take top match per mag
# some mags missing due to low identity
fastani_top <- fastani %>%
  arrange(Sample, Bin, -fastani.ANI) %>%
  group_by(Sample, Bin) %>%
  slice(1)

## gtdbtk classifications
gtdbtk_bac <- read.table(here("input_final/mags/gtdbtk.bac120.summary.tsv"),
                         sep = "\t", header = T)
gtdbtk_arc <- read.table(here("input_final/mags/gtdbtk.ar122.summary.tsv"),
                         sep = "\t", header = T)
gtdbtk <- rbind(gtdbtk_bac, gtdbtk_arc)
gtdbtk <- gtdbtk %>%
  separate(user_genome, into = c("Sample", "Bin"), sep = "_", extra = "merge")

## dereplication lists
# za 95% ANI
drep_za_95 <- read.table(here("input_final/mags/za_drep95_reps.txt"), sep = "\t",
                         header = F, col.names = c("Sample", "Bin"))

# za 95% ANI
drep_za_99 <- read.table(here("input_final/mags/za_drep99_reps.txt"), sep = "\t",
                         header = F, col.names = c("Sample", "Bin"))

# za plus uhgg reps, 95% ANI
drep_za_uhgg <- read.table(here("input_final/mags/za_uhggreps_drep95_reps.txt"),
                           sep = "\t", header = F, col.names = c("Sample", "Bin"))

## binning dastool pipeline output
mags <- read.table(here("input_final/mags/binning_table_all_simple.tsv"),
                   sep = "\t", header = T, comment.char = "")

## format mag info ----
# add genomesearch and fastani comparisons
# only consider MQ and higher
mags <- mags %>%
  filter(Bin != "unbinned") %>%
  mutate(Quality = ifelse(Completeness > 90 & Contamination < 5 & tRNA >= 18 &
                            rna.16S > 0 & rna.23S > 0 & rna.5S > 0,
                          "High-quality",
                          ifelse(Completeness >= 50 & Contamination <10,
                                 "Medium-quality",
                                 "Low-quality")),
         Quality = factor(Quality, levels = c("Low-quality", "Medium-quality",
                                              "High-quality")),
         near_complete = (Completeness >= 90 & Contamination <= 5 & N50 >= 10e3 &
                          X..contigs.....0.bp. <= 500 & Coverage >= 5 &
                          (Size.Mb*1e3)/X..contigs.....0.bp. >= 5)) %>%
  left_join(za_meta, by = c("Sample" = "sample")) %>%
  left_join(gsearch, by = c("Sample", "Bin")) %>%
  left_join(fastani_top, by = c("Sample", "Bin")) %>%
  left_join(gtdbtk[, c("Sample", "Bin", "classification")], by = c("Sample", "Bin")) %>%
  left_join(labels[, c("sample", "label_abbrev")], by = c("Sample" = "sample")) %>%
  mutate(za_rep95 = paste(Sample, Bin) %in% paste(drep_za_95$Sample, drep_za_95$Bin),
         za_rep99 = paste(Sample, Bin) %in% paste(drep_za_99$Sample, drep_za_99$Bin),
         uhgg_rep = paste(Sample, Bin) %in% paste(drep_za_uhgg$Sample, drep_za_uhgg$Bin),
         novel = (is.na(fastani.ANI) | fastani.ANI < 95) & uhgg_rep)

# "near-complete" mags
mags %>%
  filter(near_complete) %>%
  tally()

mags %>%
  count(Quality)

mags_mqplus <- mags %>%
  filter(Quality != "Low-quality")

mags_mqplus %>%
  filter(za_rep99) %>%
  separate(classification, into = c("D", "P", "C", "O", "F", "G", "S"),
           sep = ";", remove = F) %>%
  count(G, site) %>%
  arrange(-n)

# save to supplementary table
mags_supp <- mags %>%
  left_join(labels[, c("sample", "label")], by = c("Sample" = "sample")) %>%
  filter(za_rep99) %>%
  select(label, Bin, site, best_species, classification, Completeness, Contamination,
         Strain.heterogeneity, Size.Mb, N50, Genes, tRNA, starts_with("rna"),
         za_rep95, uhgg_rep, fastani.ANI, novel) %>%
  mutate(N50 = N50 / 1000)

names(mags_supp) <- c("Sample", "Bin", "Site", "NCBI species", "GTDB taxonomy",
                      "Completeness", "Contamination", "Strain heterogeneity",
                      "Size (Mb)", "N50 (Kb)", "# genes", "# tRNA",
                      "# 16S rRNA", "# 23S rRNA", "# 5S rRNA",
                      "dRep 95% sANI rep (cohort)",
                      "dRep 95% sANI rep (cohort plus UHGG species reps)",
                      "FastANI to closest UHGG genome", "Novel compared to UHGG")

write.table(mags_supp, here("final_tables/table_S10_drep_mags.txt"), sep = "\t",
            quote = F, row.names = F)

## novel mags ----
## what are the mags with low identity to uhgg?
novel_mags <- mags_mqplus %>%
  filter(novel == T) %>%
  separate(classification, into = c("D", "P", "C", "O", "F", "G", "S"),
           sep = ";", remove = F)

write.table(novel_mags, here("final_tables/novel_mags.txt"), sep = "\t",
            quote = F, row.names = F)

## figure ---
data_a <- mags_mqplus %>%
  arrange(Completeness, -Contamination) %>%
  mutate(plot_quality = ifelse(Quality == "Medium-quality" & near_complete,
                               "Near-complete", as.character(Quality)),
         plot_quality = factor(plot_quality,
                               levels = c("High-quality", "Near-complete",
                                          "Medium-quality")))

a <- ggplot(data_a, aes(Completeness, Contamination, color = plot_quality,
                        size = plot_quality)) +
  geom_point(alpha = 0.75) +
  scale_size_manual(values = c(1.5, 0.74, 0.75)) +
  scale_color_brewer(palette = "Set1") +
  theme_cowplot() +
  background_grid() +
  theme(legend.position = "top",
        legend.justification = "center",
        legend.spacing.x = unit(0.1, 'cm')) +
  labs(color = "",
       size = "")

data_b <- mags_mqplus %>%
  mutate(site_abbrev = ifelse(site == "Soweto", "SWT", "BBR"))

b <- ggplot(data_b, aes(site_abbrev, fastani.ANI, fill = site_abbrev)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.7), alpha = 0.75,
              size = 1, color = "darkgray") +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  # ylim(c(0, 100)) +
  labs(x = "", 
       y = "FastANI to UHGG") + 
  scale_fill_manual(values = za_pal) +
  theme_cowplot() +
  theme(legend.position = "none") +
  background_grid(major = "y") +
  geom_hline(yintercept = 95, linetype = "dashed")

c <- novel_mags %>%
  group_by(P) %>%
  mutate(P = gsub("p__", "", P)) %>%
  tally() %>%
  arrange(-n) %>%
  ggtexttable(rows = NULL, cols = c("GTDB phylum", "No. MAGs"),
              theme = ttheme("classic", base_size = 9,
                             padding = unit(c(4, 4), "mm")))

d <- novel_mags %>%
  select(label_abbrev, P, C, O, `F`, G, S) %>%
  arrange(P, C, O, `F`, G, S) %>%
  mutate(across(where(is.character), ~gsub("[a-z]__", "", .))) %>%
  ggtexttable(rows = NULL, cols = c("Participant", "Phylum", "Class", "Order",
                                    "Family", "Genus", "Species"),
              theme = ttheme("classic", base_size = 8,
                             padding = unit(c(4, 4), "mm")))

plot_grid(plot_grid(a, b, labels = c("A", "B"), rel_widths = c(0.6, 0.4)),
          d, ncol = 1, labels = c("", "C"))

# plot_grid(a,
#           plot_grid(b, c, rel_widths = c(0.55, 0.45), labels = c("B", "C")),
#           ncol = 1,
#           rel_heights = c(0.55, 0.45),
#           labels = c("A", ""))

ggsave(here("final_plots/supplementary/figure_S14_MAGs.png"),
       width = 8, height = 9)
