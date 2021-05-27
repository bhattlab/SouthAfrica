library(here)

## color palettes ----
global_pal <- c("#E3211C", "#F89897", "#6A3D9A", "#CAB2D6", "#1F78B4", "#A5CEE3")
za_pal <- global_pal[3:4]

## study sites ----
sites <- c("Tanzania", "Madagascar", "Bushbuckridge", "Soweto",
           "Sweden", "United States")

## vanish taxa ----
vanish_F <- c("Prevotellaceae", "Succinivibrionaceae", "Paraprevotellaceae",
              "Spirochaetaceae")

taxonomy <- read.table(here("input_final/taxonomy/kraken_feb2019_inspect_mpa.out"),
                       sep = "\t", quote = "", comment.char = "")

vanish_G <- taxonomy %>%
  filter(grepl("g__", V1) & grepl(paste(vanish_F, collapse = "|"), V1)) %>%
  mutate(feature = gsub("\\|s__.+", "", V1)) %>%
  pull(feature) %>%
  unique()

# metadata ----
za_pheno <- readRDS(here("rds/za_pheno.rds"))
za_meta <- readRDS(here("rds/za_meta.rds"))

pheno_global <- readRDS(here("rds/pheno_global.rds"))

# sequence data ----
for (rank in c("S", "G", "F", "O", "C", "P")){
  
  # global
  assign(paste0("global_", rank),
         readRDS(here("rds", paste0("global_", rank, ".rds"))))
  
  assign(paste0("global_", rank, "_rel"),
         readRDS(here("rds", paste0("global_", rank, "_rel.rds"))))
  
  assign(paste0("global_", rank, "_pseudo_rel"),
         readRDS(here("rds", paste0("global_", rank, "_pseudo_rel.rds"))))
  
  # za
  assign(paste0("za_", rank),
         readRDS(here("rds", paste0("za_", rank, ".rds"))))
  
  assign(paste0("za_", rank, "_rel"),
         readRDS(here("rds", paste0("za_", rank, "_rel.rds"))))
  
  assign(paste0("za_", rank, "_pseudo_rel"),
         readRDS(here("rds", paste0("za_", rank, "_pseudo_rel.rds"))))
}

# CSS normalized ----
za_S_css <- readRDS(here("rds/za_S_css.rds"))
za_G_css <- readRDS(here("rds/za_G_css.rds"))
global_S_css <- readRDS(here("rds/global_S_css.rds"))
global_G_css <- readRDS(here("rds/global_G_css.rds"))

# sourmash ----
for (k in c("21", "31", "51")){
  for (suffix in c("", "_track_abund")){
    
    assign(paste0("sourmash_k", k, suffix),
           readRDS(here("rds", paste0("sourmash_k", k, suffix, ".rds"))))
  }
}

