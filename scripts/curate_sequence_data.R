library(genefilter)
library(here)
library(metagenomeSeq)
library(tidyverse)
library(vegan)

set.seed(1)

# metadata ----
load(here("RData/metadata.RData"))

# public global sequence data ----
# exclude rampelli participants under the age of 18
exclude_rampelli <- c("SRR1929485", "SRR1929574", "SRR1930128", "SRR1930132",
                      "SRR1930134")

pheno_global <- read.table(here("input_final/pheno/pheno_global.txt"), sep = "\t",
                           header = T)

pheno_global <- pheno_global %>%
  filter(!(sample %in% exclude_rampelli))

bf_S <- read.table(here("input_final/kraken/genbank_jan2020/burkina_faso/bracken_S_reads.txt"),
                   sep = "\t", header = T, quote = "", comment.char = "")

pheno_bf <- data.frame("sample" = colnames(bf_S), "site" = "Burkina Faso",
                       "site2" = "Burkina Faso")

pheno_global <- rbind(pheno_global, pheno_bf)
pheno_global$site2 <- factor(pheno_global$site2,
                             levels = c("Tanzania", "Madagascar", "Burkina Faso",
                                        "Bushbuckridge", "Soweto", "Sweden",
                                        "United States"))

# kraken/bracken output ----
for (rank in c("S", "G", "F", "O", "C", "P")){
  
  # za reads
  za_path <- here(paste0("input_final/kraken/genbank_jan2020/za/bracken_",
                    rank, "_reads.txt"))
  za_reads <- read.table(za_path, sep = "\t", check.names = F, quote = "",
                         comment.char = "")
  za_reads <- za_reads[, za_meta$sample]
  
  # global reads
  global_path <- here(paste0("input_final/kraken/genbank_jan2020/global/bracken_",
                  rank, "_reads.txt"))
  global_reads <- read.table(global_path, sep = "\t", check.names = F, quote = "",
                             comment.char = "")
  global_reads <- global_reads[, !grepl("[A-C]\\d+$", names(global_reads))]
  
  # rampelli reads
  rampelli_path <- here(paste0("input_final/kraken/genbank_jan2020/rampelli/bracken_",
                    rank, "_reads.txt"))
  rampelli_reads <- read.table(rampelli_path, sep = "\t", check.names = F, quote = "",
                               comment.char = "")
  
  # burkina faso reads
  bf_path <- here(paste0("input_final/kraken/genbank_jan2020/burkina_faso//bracken_",
                         rank, "_reads.txt"))
  bf_reads <- read.table(bf_path, sep = "\t", check.names = F, quote = "",
                         comment.char = "")
  
  # exclude non-adult samples
  rampelli_reads <- rampelli_reads[, which(!(names(rampelli_reads) %in%
                                               exclude_rampelli))]
  
  # merge datasets
  global_reads <- merge(global_reads, za_reads, by = "row.names",
                        all = T) %>%
    column_to_rownames("Row.names")
  
  global_reads <- merge(global_reads, rampelli_reads, by = "row.names",
                        all = T) %>%
    column_to_rownames("Row.names")
  
  global_reads <- merge(global_reads, bf_reads, by = "row.names",
                        all = T) %>%
    column_to_rownames("Row.names")
  
  global_reads[is.na(global_reads)] <- 0
  
  global_reads <- global_reads[, pheno_global$sample]
  global_reads <- global_reads[!(rownames(global_reads) %in% c("Homo", "Homo sapiens")), ]
  global_reads <- global_reads[rev(order(rowMeans(global_reads))), ]
  
  # filter out rare/low-abundance uninformative features
  global_reads <- global_reads[genefilter(global_reads, pOverA(p = 0.01, A = 0)), ]
  
  # za_reads <- global_reads[, za_meta$sample]
  za_reads <- za_reads[genefilter(za_reads, pOverA(p = 0.01, A = 0)), ]
  za_reads <- za_reads[rev(order(rowMeans(za_reads))), ]
  
  global_rel <- sweep(global_reads, 2, colSums(global_reads), FUN = "/")
  za_rel <- sweep(za_reads, 2, colSums(za_reads), FUN = "/")
  
  global_pseudo <- global_reads + 1
  global_pseudo_rel <- sweep(global_pseudo, 2, colSums(global_pseudo), FUN = "/")
  
  za_pseudo <- za_reads + 1
  za_pseudo_rel <- sweep(za_pseudo, 2, colSums(za_pseudo), FUN = "/")
  
  # global
  assign(paste0("global_", rank), global_reads)
  assign(paste0("global_", rank, "_rel"), global_rel)
  assign(paste0("global_", rank, "_pseudo_rel"), global_pseudo_rel)
  
  # za
  assign(paste0("za_", rank), za_reads)
  assign(paste0("za_", rank, "_rel"), za_rel)
  assign(paste0("za_", rank, "_pseudo_rel"), za_pseudo_rel)
}

# CSS normalize ----
for (rank in c("S", "G")){
  cts <- get(paste0("za_", rank))
  mr <- newMRexperiment(cts)
  p <- cumNormStatFast(mr)
  mr_css <- cumNorm(mr, p = p)
  counts_css <- MRcounts(mr_css, norm = T, log = T)
  
  assign(paste0("za_", rank, "_css"), counts_css)
}

# global
for (rank in c("S", "G")){
  cts <- get(paste0("global_", rank))
  mr <- newMRexperiment(cts)
  p <- cumNormStatFast(mr)
  mr_css <- cumNorm(mr, p = p)
  counts_css <- MRcounts(mr_css, norm = T, log = T)
  
  assign(paste0("global_", rank, "_css"), counts_css)
}

# rarefied and normalized ----
za_S_rare <- rrarefy(t(za_S), min(colSums(za_S)))
za_S_rare <- data.frame(t(za_S_rare))

mr <- newMRexperiment(za_S_rare)
p <- cumNormStatFast(mr)
mr_css <- cumNorm(mr, p = p)
za_S_rare_css <- MRcounts(mr_css, norm = T, log = T)
za_S_rare_css <- data.frame(za_S_rare_css)

za_S_rare_rel <- sweep(za_S_rare, 2, colSums(za_S_rare), FUN = "/")

# sourmash data ----
for (k in c("21", "31", "51")){
  for (suffix in c("", "_track_abund")){
    
    sourmash <- read.csv(here(paste0("input_final/sourmash/compare_k", k,
                                     suffix, ".csv")),
                         header = T)
    
    names(sourmash) <- gsub("X.oak.+_trim_kmers.|_concat.fq.abundtrim", "",
                            names(sourmash))
    
    rownames(sourmash) <- names(sourmash)
    
    # exclude rampelli samples and other extraneous samples
    keep <- which(names(sourmash) %in% pheno_global$sample)
    sourmash <- sourmash[keep, keep]
    
    assign(paste0("sourmash_k", k, suffix), sourmash)
  }
}

# study sites ----
sites <- c("Tanzania", "Madagascar", "Burkina Faso", "Bushbuckridge", "Soweto",
           "Sweden", "United States")

# vanish taxa ----
vanish_F <- c("Prevotellaceae", "Succinivibrionaceae", "Paraprevotellaceae",
              "Spirochaetaceae")

taxonomy <- read.table(
  here("input_final/taxonomy/kraken_feb2019_inspect_mpa.out"),
  sep = "\t", quote = "", comment.char = "")

vanish_G <- taxonomy %>%
  filter(grepl("g__", V1) & grepl(paste(vanish_F, collapse = "|"), V1)) %>%
  mutate(feature = gsub("\\|s__.+", "", V1)) %>%
  pull(feature) %>%
  unique()

# color palettes ----
global_pal <- c("#E3211C", "#F89897", "#FF7F00", "#6A3D9A", "#CAB2D6",
                "#1F78B4", "#A5CEE3")
names(global_pal) <- levels(pheno_global$site2)

za_pal <- global_pal[c("Bushbuckridge", "Soweto")]

# save rdata ----

# global sequence data
filelist = ls()[grepl("global_[A-Z].*", ls())]
save(pheno_global, sites, vanish_F, vanish_G, taxonomy, global_pal, za_pal,
     list = filelist, file = here("RData/global_data.RData"))

# za sequence data
filelist = ls()[grepl("za_[A-Z].*", ls())]
save(global_pal, za_pal, list = filelist, file = here("RData/za_data.RData"))

# sourmash data
filelist = ls()[grepl("sourmash_k", ls())]
save(global_pal, za_pal, list = filelist, file = here("RData/sourmash_data.RData"))

# palettes
save(global_pal, za_pal, file = here("RData/palettes.RData"))

# vanish taxa
save(vanish_F, vanish_G, file = here("RData/vanish.RData"))
