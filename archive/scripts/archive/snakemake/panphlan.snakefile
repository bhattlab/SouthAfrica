import re,os,subprocess
from os.path import join, expanduser, abspath

################################################################################
# specify project directories
DATA_DIR = config["reads_directory"]
# PROJECT_DIR = config["output_directory"]
READ_SUFFIX = config["read_specification"]
EXTENSION = config["extension"]
# path to pangenome files
PANGENOME_PATH = config["pangenome_path"]
# pangenome names -- must match files in PANGENOME_PATH
PANGENOME_NAMES = config["pangenome_name"]
# if gzipped, set this. otherwise not

gz_ext = '.gz' if EXTENSION.endswith('.gz') else ''
# print(gz_ext)

# convert PROJECT_DIR to absolute path
# if PROJECT_DIR[0] == '~':
# 	PROJECT_DIR = expanduser(PROJECT_DIR)
# PROJECT_DIR = abspath(PROJECT_DIR)

# get file names
FILES = [f for f in os.listdir(DATA_DIR) if f.endswith(EXTENSION)]
SAMPLE_PREFIX = list(set([re.split('|'.join(['_' + a + '\.' for a in READ_SUFFIX] + ['_orphans']), i)[0] for i in FILES]))

# print(expand("{pangenome_name}_pangenome.tsv", pangenome_name=PANGENOME_NAMES))
# print(expand(os.path.join(DATA_DIR, "{sample}_" + READ_SUFFIX[0] + EXTENSION), sample=SAMPLE_PREFIX))
################################################################################
# localrules: get_db

# print(expand("{pangenome_name}_pangenome.tsv", pangenome_name=PANGENOME_NAMES))

rule all:
	input:
		# PANGENOME_PATH + "/" + PANGENOME_NAME + "_pangenome.tsv",
		expand("{pangenome_name}_pangenome.tsv", pangenome_name=PANGENOME_NAMES),
		expand("{pangenome_name}_coverage.tsv", pangenome_name=PANGENOME_NAMES),
		# PANGENOME_NAME + "_pangenome.tsv",
		# PANGENOME_NAME + "_coverage.tsv",
		expand("02_map_results/{sample}_{pangenome_name}.csv.bz2", pangenome_name=PANGENOME_NAMES, sample=SAMPLE_PREFIX)

################################################################################

# troubleshooting panphlan_pangenome_generation.py error -- change to python3

# annotate all genomes with prokka
# rule prokka:
# 	input: "fna/{genome}.fna"
# 	output: "gff/{genome}/{genome}.gff"
# 	shell:
# 		"prokka {input} --force --outdir gff/{wildcards.genome} --prefix {wildcards.genome} --centre X --compliant"
#
# rule make_db:
# 	input: expand("gff/{genome}/{genome}.gff", genome=[os.path.splitext(os.path.basename(g))[0] for g in os.listdir('/labs/asbhatt/fiona/za/Genome_bins/fna')])
# 	output: "panphlan_database/panphlan_pcopri-mag_pangenome.csv"
# 	params:
# 		name="pcopri-mag",
# 		fna_dir="fna",
# 		gff_dir="gff"
# 	resources:
# 		time=72,
# 		mem=32
# 	shell:
# 		"panphlan_pangenome_generation.py -c {params.name} "\
# 		"--i_fna {params.fna_dir} --i_gff {params.gff_dir} -o panphlan_database --verbose"

# add orphan reads here
rule merge:
	input:
		fwd = os.path.join(DATA_DIR, "{sample}_" + READ_SUFFIX[0] + EXTENSION),
		rev = os.path.join(DATA_DIR, "{sample}_" + READ_SUFFIX[1] + EXTENSION)
	output: "01_merged/{sample}.fq.gz"
	resources:
		mem=8,
		time=6
	shell:
		"seqtk mergepe {input} > {01_merged/{wildcards.sample}.fq; gzip 01_merged/{wildcards.sample}.fq}"

rule panphlan_map:
    input: "01_merged/{sample}.fq.gz"
	output: "02_map_results/{sample}_{pangenome_name}.csv.bz2"
	params:
		# pangenome_name = "{pangenome_name}",
		pangenome_path = os.path.join(PANGENOME_PATH, "panphlan_{pangenome_name}")
		# out_folder = "02_map_results/{sample}_{pangenome_name}.csv"
	resources:
		time=8,
		mem=32
	threads: 8
	shell:
		"panphlan_map.py -c {wildcards.pangenome_name} --i_bowtie2_indexes {params.pangenome_path} -i {input} -o {output} --verbose -p {threads} -m {resources.mem}"

rule panphlan_profile:
	input: expand("02_map_results/{sample}_{{pangenome_name}}.csv.bz2", sample=SAMPLE_PREFIX)
	output:
		dna = "{pangenome_name}_pangenome.tsv",
		cov = "{pangenome_name}_coverage.tsv"
	params:
		# pangenome_name = PANGENOME_NAME,
		pangenome_path = os.path.join(PANGENOME_PATH, "panphlan_{pangenome_name}")
	resources:
		time=lambda wildcards, attempt: 12 * attempt,
		mem=250
	shell:
		"panphlan_profile.py -c {wildcards.pangenome_name} --i_bowtie2_indexes {params.pangenome_path} -i 02_map_results --o_dna {output.dna} --o_cov {output.cov} --verbose"
