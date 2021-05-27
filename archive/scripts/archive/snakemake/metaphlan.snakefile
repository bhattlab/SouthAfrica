import re,os,subprocess
from os.path import join, expanduser, abspath

################################################################################
# specify project directories
DATA_DIR = config["reads_directory"]
READ_SUFFIX = config["read_specification"]
EXTENSION = config["extension"]
BOWTIE2_INDEXES = config["bowtie2"]
PKL = config['pkl']
MARKERS = config['markers']
SPECIES = config['species']

# if gzipped, set this. otherwise not
gz_ext = '.gz' if EXTENSION.endswith('.gz') else ''

# get file names
FILES = [f for f in os.listdir(DATA_DIR) if f.endswith(EXTENSION)]
SAMPLE_PREFIX = list(set([re.split('|'.join(['_' + a + '\.' for a in READ_SUFFIX] + ['_orphans']), i)[0] for i in FILES]))

################################################################################
rule all:
    input:
        "clades.txt",
        "05_strainphlan/RAxML_bestTree.s__" + SPECIES + ".tree"

################################################################################

rule metaphlan2:
    input:
        fwd = DATA_DIR + "{sample}_" + READ_SUFFIX[0] + EXTENSION,
        rev = DATA_DIR + "{sample}_" + READ_SUFFIX[1] + EXTENSION
    output:
        bowtie2 = "02_metaphlan2/{sample}.bowtie2.bz2",
        sam = "02_metaphlan2/{sample}.sam.bz2",
        profile = "02_metaphlan2/{sample}.txt"
    params:
        index = BOWTIE2_INDEXES
    resources:
        time=lambda wildcards, attempt: attempt * 12,
        mem=lambda wildcards, attempt: attempt * 32,
        threads=6
    shell:
        "metaphlan2.py {input.fwd},{input.rev} --bowtie2db /labs/asbhatt/fiona/tools/metaphlan2/db_v20 --bowtie2out {output.bowtie2} --samout {output.sam} --nproc {resources.threads} --input_type fastq > {output.profile}"

rule sample2markers:
    input: "02_metaphlan2/{sample}.sam.bz2",
    output: "03_sample_markers/{sample}.markers"
    resources:
        mem=64,
        time=lambda wildcards, attempt: attempt * 6
    shell:
        "sample2markers.py --ifn_samples {input} --input_type sam "\
        "--output_dir 03_sample_markers/"

rule strainphlan_clades:
    input: sample_markers=expand("03_sample_markers/{sample}.markers", sample=SAMPLE_PREFIX)
    output: "clades.txt"
    resources:
        mem=64,
        time=12,
        threads=8
    shell:
        "strainphlan.py --ifn_samples {input} --output_dir . --nprocs_main {resources.threads} --print_clades_only > clades.txt"

rule extract_markers:
    input:
        mpa_pkl = PKL,
        markers = MARKERS
    output:
        "04_clade_markers/{species}.markers.fasta"
    resources:
        mem=16,
        time=4
    shell:
        "extract_markers.py --mpa_pkl {input.mpa_pkl} "\
        "--ifn_markers {input.markers} --clade s__{wildcards.species} "\
        "--ofn_markers {output}"

rule strainphlan:
    input:
        sample_markers=expand("03_sample_markers/{sample}.markers", sample=SAMPLE_PREFIX),
        mpa_pkl = PKL,
        markers = MARKERS,
        clade_markers=rules.extract_markers.output
    output:
        "05_strainphlan/RAxML_bestTree.s__{species}.tree"
    resources:
        mem=100,
        time=12
    shell:
        "strainphlan.py --mpa_pkl {input.mpa_pkl} "\
        "--ifn_samples {input.sample_markers} "\
        "--ifn_markers {input.clade_markers} "\
        # "--ifn_ref_genomes 0.reference_genomes/foo/ncbi-genomes/*.fna.gz "\
        "--output_dir 05_strainphlan --clades s__{wildcards.species}"
