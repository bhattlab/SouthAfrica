# Bushbuckridge and Soweto microbiome project

Analysis scripts for Wits-Stanford shotgun microbiome data.

To generate all figures from the manuscript, please follow these steps:

1. Access phenotype data from EGA under dataset ID [EGAD00001006581](https://ega-archive.org/datasets/EGAD00001006581)
2. Edit `scripts/create_metadata.R` to point to the location of the EGA data and run the script to create required RData file
3. install all required packages and run each R script in `scripts/figures` and `scripts/figures/supplementary`

e.g. from the command line:

     find scripts/figures -name "*.R" | xargs -I foo Rscript foo

Required packages:

* DESeq2
* MASS
* Maaslin2
* RColorBrewer
* cowplot
* dplyr
* genefilter
* ggplot2
* ggpubr
* ggrepel
* ggspatial
* gplots
* gtools
* harrietr
* here
* maps
* metagenomeSeq
* reshape2
* rgeos
* rnaturalearth
* rnaturalearthdata
* scales
* sf
* tidyverse
* vegan
