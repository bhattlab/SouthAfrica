# plot table of conmeds and categories

library(cowplot)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(here)
library(vegan)

# metadata ----
za_meta <- readRDS(here("rds/za_meta.rds"))

# conmeds data table
conmeds <- read.table(here("input_final/conmeds.tsv"), sep = "\t",
                      header = T, quote = "")

conmeds <- conmeds %>%
  filter(!(medicine_CURATED %in% c("", "??", "NONE", "UNKNOWN")))

# set "other" as last level
meds <- unique(conmeds$medicine_CATEGORY)

conmeds$medicine_CATEGORY <- factor(conmeds$medicine_CATEGORY,
                                    levels = c(sort(meds[meds != "OTHER"]), "OTHER"))

conmeds_tbl <- conmeds %>%
  group_by(medicine_CATEGORY, medicine_CURATED) %>%
  tally()

write.table(conmeds_tbl, here("final_tables/table_s2_conmeds.txt"),
            sep = "\t", quote = F, row.names = F)

# stats for text ----
# n participants anti-infective/antibiotic
conmeds %>%
  filter(medicine_CATEGORY == "ANTIBIOTIC") %>%
  dplyr::select(sample, site) %>%
  distinct() %>%
  group_by(site) %>%
  tally()

# n participants anti-hypertensive
conmeds %>%
  filter(medicine_CATEGORY == "ANTI-HYPERTENSIVE") %>%
  dplyr::select(sample, site) %>%
  distinct() %>%
  group_by(site) %>%
  tally()

# conmeds mds plots ----

# mds plot with drug trt colored
global_pal <- c("#E3211C", "#F89897", "#6A3D9A", "#CAB2D6", "#1F78B4", "#A5CEE3")
za_pal <- global_pal[3:4]
za_meta <- readRDS(here("rds/za_meta.rds"))
za_S_css <- readRDS(here("rds/za_S_css.rds"))

vare_dis <- vegdist(t(za_S_css), method = "bray")

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
mds_meta <- merge(mds_data, za_meta, by = "sample")

# plot medicine categories with at least n patients
categories <- conmeds %>%
  dplyr::select(sample, medicine_CATEGORY) %>%
  distinct() %>%
  group_by(medicine_CATEGORY) %>%
  tally() %>%
  filter(n >= 2, !medicine_CATEGORY %in% c("SUPPLEMENT", "OTHER")) %>%
  pull(medicine_CATEGORY) %>%
  as.character()

# categories <- c(categories[categories != "ANTIBIOTIC"], "ANTIBIOTIC")

plot_conmeds <- function(categories, mds_plot){
  
  plot_data <- data.frame()
  
  for (c in categories){
    p <- mds_plot
    
    # which patients
    pts <- conmeds %>%
      filter(medicine_CATEGORY == c) %>%
      pull(sample) %>%
      unique()
    
    # p$is_taking <- ifelse(p$sample %in% pts, c, NA)
    p$taking_med <- ifelse(p$sample %in% pts, T, NA)
    p$category <- c
    p$alpha <- ifelse(p$sample %in% pts, T, F)
    
    plot_data <- rbind(plot_data, p)
  }
  
  # facet labels
  plot_data$category <- Hmisc::capitalize(tolower(as.character(plot_data$category)))
  plot_data$category[plot_data$category == "Nsaid"] <- "NSAID"
  
  # set plotting order so that colored points are plotted on top
  plot_data <- plot_data %>%
    arrange(alpha)
  
  ggplot(plot_data, aes(x, y, shape = site, color = taking_med, alpha = alpha)) +
    geom_point(size = 2) +
    scale_alpha_manual(values = c(0.5, 0.85)) +
    scale_color_manual(values = c("red3"), na.value = "gray") +
    facet_wrap(~ category, ncol = 3) +
    labs(x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
         y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
         color = "") +
    theme_cowplot() +
    labs(
      shape = "Site",
      color = "Medication"
    ) +
    guides(color = F, alpha = F) +
    theme(
      legend.position = "top",
      legend.justification = "center"
    )
}

# plot by category and color by medicine
plot_conmeds_by_drug <- function(category, mds_plot){
  
  # which patients
  pts <- conmeds %>%
    filter(medicine_CATEGORY == category) %>%
    pull(sample) %>%
    unique()
  
  # unique drug combos
  c <- conmeds %>%
    filter(medicine_CATEGORY == category) %>%
    dplyr::select(sample, medicine_CURATED, medicine_CATEGORY) %>%
    distinct() %>%
    mutate(
      value = 1
    )
  
  c_wide <- reshape2::dcast(c, sample ~ medicine_CURATED,
                            value.var = "medicine_CURATED")
  c_wide[is.na(c_wide)] <- ""
  c_wide$label <- apply(c_wide[, 2:ncol(c_wide)], 1, paste, collapse=",")
  c_wide$label <- gsub("^,+|,+$", "", c_wide$label)
  c_wide$label <- gsub(",+", ",", c_wide$label)
  c_wide$label <- Hmisc::capitalize(tolower(c_wide$label))
  
  plot_data <- merge(mds_plot, c_wide, by = "sample", all.x = T)
  
  # n for each med / combo
  n_med <- plot_data %>%
    group_by(label) %>%
    tally() %>%
    mutate(
      label = ifelse(is.na(label), "None", label),
      lab = paste0(label, " (", n, ")")
    ) %>%
    pull(lab)
  
  plot_data$category <- "Antibiotic"
  plot_data$alpha <- ifelse(is.na(plot_data$label), F, T)
  
  # fit legend labels
  plot_data$label <- gsub(",", ", ", plot_data$label)
  
  # plot colored points on top
  plot_data <- plot_data %>%
    arrange(alpha)
  
  ggplot(plot_data, aes(x, y, shape = site, color = label, alpha = alpha)) +
    geom_point(size = 2) +
    scale_alpha_manual(values = c(0.5, 0.85)) +
    scale_color_discrete(na.value = "gray",
                         breaks = sort(unique(plot_data$label))) +
    labs(x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
         y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
         color = "") +
    facet_wrap(~ category) +
    theme_cowplot() +
    theme(legend.position = "bottom") +
    labs(
      shape = "Site",
      color = ""
    ) +
    guides(
      alpha = F,
      shape = F,
      color = guide_legend(ncol = 1)
      )
}

## adonis table
vare_dis <- vegdist(t(za_S_css), method = "bray")

meta <- za_meta %>%
  filter(sample %in% colnames(za_S_css))

meta <- meta[match(colnames(za_S_css), meta$sample), ]

# meds per patient
meds_pp <- conmeds %>%
  dplyr::select(sample, medicine_CATEGORY) %>%
  distinct()

meds_pp_wide <- reshape2::dcast(meds_pp, sample ~ medicine_CATEGORY,
                                value.var = "medicine_CATEGORY")
samples <- meds_pp_wide$sample

meds_pp_wide[!is.na(meds_pp_wide)] <- T
meds_pp_wide[is.na(meds_pp_wide)] <- F

meds_pp_wide$sample <- samples

meta <- merge(meta, meds_pp_wide, by = "sample", all.x = T)
meta[is.na(meta)] <- F

res_tbl <- data.frame()
for (c in categories){
  a_res <- adonis(formula(paste0("vare_dis ~`", c, "`")), data = meta)
  a_tbl <- data.frame(a_res$aov.tab[1, 5:6])
  
  res_tbl <- rbind(res_tbl, a_tbl)
}

# format adonis results table
res_tbl <- res_tbl %>%
  tibble::rownames_to_column("Category") %>%
  mutate(
    R2 = round(R2, 3)
  )
names(res_tbl) <- c("Category", "R2", "Pr(>F)")

res_tbl$Category <- gsub("`", "", res_tbl$Category)
res_tbl$Category <- Hmisc::capitalize(tolower(res_tbl$Category))

res_tbl$Category[res_tbl$Category == "Nsaid"] <- "NSAID"

t <- res_tbl %>%
  mutate(FDR = p.adjust(`Pr(>F)`, method = "fdr"),
         FDR = round(FDR, 3)) %>%
  ggtexttable(rows = NULL, theme =
                ttheme("classic", base_size = 12, padding = unit(c(4, 4), "mm")))


## compile figure
a <- plot_conmeds(categories, mds_meta)
b <- plot_grid(plot_conmeds_by_drug("ANTIBIOTIC", mds_meta), t,
               rel_widths = c(0.395, 0.605), nrow = 1, labels = c("B", "C"))

plot_grid(a, b, ncol = 1, rel_heights = c(0.71, 0.29), labels = c("A", ""))

ggsave(here("final_plots/supplementary/figure_S5_conmeds_mds.png"),
       width = 8, height = 14)
