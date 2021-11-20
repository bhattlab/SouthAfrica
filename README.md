# AWI-Gen Microbiome Pilot: Metagenomics

Analysis scripts for Wits-Stanford shotgun microbiome data.

To generate all figures from the manuscript, please follow these steps:

1. Clone this repo and create a new R project (in RStudio, File > New Project > Existing Directory)
2. Access phenotype and medication data from EGA under dataset ID [EGAD00001006581](https://ega-archive.org/datasets/EGAD00001006581). You can submit an access request to the H3Africa Data Access Committee at https://catalog.h3africa.org
3. Edit `scripts/create_metadata.R` to point to the locations of the EGA phenotype and medication data and run this script to create the required RData file
4. Install all required packages and run R scripts in `scripts/figures`, `scripts/figures/supplementary`, and `nanopore/scripts`

     e.g. from the command line:

          find scripts/figures -name "*.R" | xargs -I foo Rscript foo

Required packages for figures:

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
* metagenomeSeq
* reshape2
* scales
* tidyverse
* vegan

Optional packages to reproduce map figure from Supplementary Information (`scripts/figures/supplementary/supp_info_map.R`)
* maps
* rgeos
* rnaturalearth
* rnaturalearthdata
* sf
