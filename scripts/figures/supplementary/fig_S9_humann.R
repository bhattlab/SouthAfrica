library(cowplot)
library(here)
library(Maaslin2)
library(tidyverse)
library(RColorBrewer)

# metadata ----
load(here("RData/metadata.RData"))
# load(here("RData/za_data.RData"))

rownames(za_meta) <- za_meta$sample

# functions ----

# load humann output
load_humann <- function(fp, remove_unclass = T){
  humann_out <- data.table::fread(fp, sep = "\t", header = T)
  
  names(humann_out) <- gsub("_concat_Abundance-CPM", "", names(humann_out))
  
  if (remove_unclass){
    humann_out <- humann_out %>%
      filter(!grepl("[a-z]__|\\|unclassified", `# Pathway`))
  }
  
  humann_out <- humann_out %>%
    column_to_rownames("# Pathway")
  
  return(humann_out)
}

# maaslin function
run_maaslin_humann <- function(humann_filt, meta_filt, dirname, model,
                               norm = "NONE", trform = "NONE"){
  
  fit_data <- Maaslin2(humann_filt, meta_filt,
                       dirname,
                       transform = trform,
                       analysis_method = model,
                       fixed_effects = c("site"),
                       normalization = norm,
                       standardize = FALSE,
                       plot_heatmap = F,
                       plot_scatter = F)
  
  res <- fit_data$results %>%
    as_tibble() %>%
    filter(metadata == "site") %>%
    arrange(pval) %>%
    mutate(model = model,
           norm = norm)
  
  return(res)
}

# run maaslin on humann data ----

# maaslin parameters
model <- "CPLM"
norm <- "TSS"
transf <- "NONE"

prefix <- paste(model, norm, sep = "_")
outfile <- here(paste0("rds/res_humann_", prefix, ".rds"))

# humann input data
fp <- here("input_final/humann/join_metacyc-pwy-CPM.tsv")

humann_out <- load_humann(fp)

humann_filt <- humann_out[, za_meta$sample]
# humann_filt <- filter_humann(humann_filt)

if (!file.exists(outfile)){
  
  # run maaslin
  dirname <- here("output/maaslin2", prefix)
  
  res <- run_maaslin_humann(humann_filt, za_meta, dirname,
                            model = model, norm = norm, trform = transf)
  
  # save results
  saveRDS(res, outfile)
}

# volcano plot ----
# prefix <- paste(model, norm, sep = "_")
# outfile <- here(paste0("rds/res_humann_", prefix, ".rds"))
# res <- readRDS(outfile)
# 
# labs <- data.frame("label" = rownames(humann_filt),
#                    "feature" = make.names(rownames(humann_filt)))
# 
# res_plot <- res %>%
#   left_join(labs, by = "feature") %>%
#   mutate(label = ifelse(qval < 0.05, label, ""),
#          signif = ifelse(qval < 0.05, T, F),
#          effect = ifelse(pval < 0.05, coef < 0, NA))
# 
# g <- ggplot(res_plot, aes(coef, -log10(pval), label = label, size = signif,
#                           color = effect)) +
#   geom_point() +
#   scale_size_manual(values = c(1.5, 2.5)) +
#   scale_color_manual(values = rev(za_pal), na.value = "#000000") +
#   ggrepel::geom_text_repel(size = 2.5) +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
#   theme_cowplot() +
#   background_grid() +
#   labs(x = "Coefficient") +
#   theme(legend.position = "none")

# coef plot ----
prefix <- paste(model, norm, sep = "_")
outfile <- here(paste0("rds/res_humann_", prefix, ".rds"))
res <- readRDS(outfile)

# feature labels
labs <- data.frame("label" = rownames(humann_filt),
                   "feature" = make.names(rownames(humann_filt)))

# fix formatting
labs$label <- gsub("&beta;", "Beta", labs$label)

res_plot <- res %>%
  left_join(labs, by = "feature") %>%
  filter(qval < 0.05,
         abs(coef) > 0.5) %>%
  mutate(effect = coef < 0,
         label = fct_reorder(label, coef))

a <- ggplot(res_plot, aes(coef, label, fill = effect)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = rev(za_pal), na.value = "#000000",
                    labels = c("Enriched in SWT", "Enriched in BBR")) +
  theme_cowplot() +
  background_grid() +
  labs(x = "Coefficient",
       y = "MetaCyc pathway",
       fill = "") +
  theme(legend.position = "bottom",
        legend.justification = "center",
        plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm"))

# stratified taxon plot of significant pathways ----
humann_strat <- load_humann(fp, remove_unclass = F)
humann_strat <- humann_strat[, za_meta$sample]

display_n <- 11

humann_plot <- humann_strat %>%
  rownames_to_column("feature") %>%
  separate(feature, into = c("pwy", "taxon"), sep = "\\|", fill = "right") %>%
  mutate(pwy = gsub("&beta;", "Beta", pwy),
         taxon = gsub("g__(.+)\\..+", "\\1", taxon)) %>%
  filter(pwy %in% res_plot$label, !is.na(taxon)) %>%
  pivot_longer(cols = -c(pwy, taxon), names_to = "sample",
               values_to = "value") %>%
  left_join(za_meta[, c("sample", "site")], by = "sample") %>%
  group_by(pwy, taxon, site) %>%
  summarise(tot = sum(value)) %>%
  group_by(pwy) %>%
  mutate(tot_rel = tot / sum(tot)) %>%
  group_by(taxon, site) %>%
  mutate(mean_rel = mean(tot_rel)) %>%
  ungroup() %>%
  mutate(taxon = fct_reorder(taxon, mean_rel, .desc = T),
         # taxon = fct_expand(taxon, "Other"),
         group_no = as.integer(factor(taxon))) %>%
  filter(group_no <= display_n | taxon == "unclassified") %>%
  select(-tot, -mean_rel, -group_no)

# add "other" column
other <- humann_plot %>%
  group_by(pwy, site) %>%
  summarise(tot_rel = 1 - sum(tot_rel)) %>%
  mutate(taxon = "Other")

fct_lvls <- c(levels(humann_plot$taxon), "Other")

humann_plot <- bind_rows(humann_plot, other) %>%
  mutate(taxon = factor(taxon, levels = fct_lvls),
         pwy = factor(pwy, levels = levels(res_plot$label)))

# color palette
# myCols <- colorRampPalette(brewer.pal(9, "Set1"))
# pal <- myCols(display_n)
# pal <- sample(pal)
# pal[display_n + 1] <- "lightgray"

pal <- brewer.pal(12, "Set3")
pal <- pal[!(pal == "#D9D9D9")]
pal <- c(pal, "#D9D9D9")

b <- ggplot(humann_plot, aes(pwy, tot_rel * 100, fill = taxon)) +
  geom_bar(stat = "identity", color = "gray60") +
  theme_cowplot() +
  scale_fill_manual(values = pal) +
  scale_y_continuous(expand = c(0, 0)) +
  # scale_fill_brewer(palette = "Set3") +
  facet_wrap(~site, nrow = 1) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        strip.background = element_rect(fill = "white"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0), "cm")) +
  labs(x = "",
       y = "Stratified relative abundance (%)",
       fill = "") +
  coord_flip() +
  guides(fill = guide_legend(ncol = 3))

plot_grid(a, b, ncol = 2, align = "h", axis = "bt", rel_widths = c(0.65, 0.35),
          labels = c("A", "B"))

ggsave(here("final_plots/supplementary/figure_S9_humann_metacyc.png"),
       width = 14, height = 8, bg = "white")

# write output table ----
humann_tbl <- humann_filt

labs <- data.frame(za_pheno)

names(humann_tbl) <- labs[match(names(humann_tbl), labs$sample), "label"]

humann_tbl <- humann_tbl[, gtools::mixedorder(names(humann_tbl))]

humann_tbl <- humann_tbl %>%
  rownames_to_column("MetaCyc Pathway")

write.table(humann_tbl, here("final_tables/table_S8_humann_pwys.txt"),
            sep = "\t", row.names = F, quote = F)
