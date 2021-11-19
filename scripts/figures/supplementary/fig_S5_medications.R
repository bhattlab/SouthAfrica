library(cowplot)
library(tidyverse)
library(ggpubr)
library(here)
library(vegan)

# load data ----
load(here("RData/metadata.RData"))
load(here("RData/za_data.RData"))
load(here("RData/medications.RData"))

# meds data table ----
# set "other" as last level
meds <- meds %>%
  mutate(medicine_category = as.factor(medicine_category),
         medicine_category = fct_relevel(medicine_category, "OTHER", after = Inf))

meds_tbl <- meds %>%
  group_by(medicine_category, medicine_name) %>%
  tally()

write.table(meds_tbl, here("final_tables/table_s2_medications.txt"),
            sep = "\t", quote = F, row.names = F)

# stats for text ----
# n participants anti-infective/antibiotic
meds_meta <- meds %>%
  left_join(za_meta, by = "study_id")

meds_meta %>%
  filter(medicine_category == "ANTIBIOTIC") %>%
  select(sample, site) %>%
  distinct() %>%
  group_by(site) %>%
  tally()

# n participants anti-hypertensive
meds_meta %>%
  filter(medicine_category == "ANTI-HYPERTENSIVE") %>%
  select(sample, site) %>%
  distinct() %>%
  group_by(site) %>%
  tally()

# meds mds plots ----
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
categories <- meds_meta %>%
  select(sample, medicine_category) %>%
  distinct() %>%
  group_by(medicine_category) %>%
  tally() %>%
  filter(n >= 2, !medicine_category %in% c("SUPPLEMENT", "OTHER")) %>%
  pull(medicine_category) %>%
  as.character()

plot_meds <- function(categories, mds_plot){
  
  plot_data <- data.frame()
  
  for (c in categories){
    p <- mds_plot
    
    # which patients
    pts <- meds_meta %>%
      filter(medicine_category == c) %>%
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
    scale_color_manual(values = c("red3"), na.value = "darkgray") +
    facet_wrap(~ category, ncol = 3) +
    coord_equal() +
    labs(x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
         y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
         color = "") +
    theme_cowplot() +
    background_grid() +
    labs(
      shape = "Site",
      color = "Medication"
    ) +
    guides(color = "none", alpha = "none") +
    theme(
      legend.position = "top",
      legend.justification = "center",
      strip.background = element_rect(fill = "gray90"),
      legend.margin = margin(c(-0.2, 0.5, -0.2, 0.5), unit = "cm")
    )
}

# plot by category and color by medicine
plot_meds_by_drug <- function(category, mds_plot){
  
  # which patients
  pts <- meds_meta %>%
    filter(medicine_category == category) %>%
    pull(sample) %>%
    unique()
  
  # unique drug combos
  c <- meds_meta %>%
    filter(medicine_category == category) %>%
    dplyr::select(sample, medicine_name, medicine_category) %>%
    distinct() %>%
    mutate(
      value = 1
    )
  
  c_wide <- reshape2::dcast(c, sample ~ medicine_name,
                            value.var = "medicine_name")
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
    scale_color_discrete(na.value = "darkgray",
                         breaks = sort(unique(plot_data$label))) +
    coord_equal() +
    labs(x = paste("MDS 1 (", mds_var_per[1], "%)",sep=""),
         y = paste("MDS 2 (", mds_var_per[2], "%)",sep=""),
         color = "") +
    facet_wrap(~ category) +
    theme_cowplot() +
    background_grid() +
    theme(legend.position = "bottom",
          legend.justification = "left",
          strip.background = element_rect(fill = "gray90"),
          legend.margin = margin(c(-0.2, 0.5, -0.2, -0.2), unit = "cm"),
          legend.key = element_rect(size = 0.46),
          legend.key.height = unit(0.46, "cm"),
          legend.key.width = unit(0.46, "cm")) +
    labs(
      shape = "Site",
      color = ""
    ) +
    guides(
      alpha = "none",
      shape = "none",
      color = guide_legend(ncol = 1)
      )
}

# adonis table ----
vare_dis <- vegdist(t(za_S_css), method = "bray")

meta <- za_meta %>%
  filter(sample %in% colnames(za_S_css))

meta <- meta[match(colnames(za_S_css), meta$sample), ]

# meds per patient
meds_pp <- meds_meta %>%
  dplyr::select(sample, medicine_category) %>%
  distinct()

meds_pp_wide <- reshape2::dcast(meds_pp, sample ~ medicine_category,
                                value.var = "medicine_category")
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
                ttheme("classic", base_size = 11, padding = unit(c(4, 4), "mm")))


# compile figure ----
a <- plot_meds(categories, mds_meta)
b <- plot_grid(plot_meds_by_drug("ANTIBIOTIC", mds_meta), t,
               rel_widths = c(0.43, 0.57), nrow = 1, labels = c("b", "c"))

plot_grid(a, b, ncol = 1, rel_heights = c(0.71, 0.29), labels = c("a", ""))

ggsave(here("final_plots/supplementary/figure_S5_medications_mds.png"),
       width = 8, height = 12, bg = "white")
ggsave(here("final_plots/pdf/supp/figure_S5_medications_mds.pdf"),
       width = 8, height = 12, bg = "white")
