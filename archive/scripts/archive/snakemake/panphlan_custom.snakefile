import re,os,subprocess
from os.path import join, expanduser, abspath

################################################################################
# specify project directories
DATA_DIR = config["reads_directory"]
READ_SUFFIX = config["read_specification"]
EXTENSION = config["extension"]
PANGENOME_PATH = config["pangenome_path"]
PANGENOME_NAME = config["pangenome_name"]
# if gzipped, set this. otherwise not

gz_ext = '.gz' if EXTENSION.endswith('.gz') else ''
# print(gz_ext)

# get file names
FILES = [f for f in os.listdir(DATA_DIR) if f.endswith(EXTENSION)]
SAMPLE_PREFIX = list(set([re.split('|'.join(['_' + a + '\.' for a in READ_SUFFIX] + ['_orphans']), i)[0] for i in FILES]))

################################################################################
# localrules: get_db

rule all:
	input:
		# PANGENOME_PATH + "/" + PANGENOME_NAME + "_pangenome.csv",
		PANGENOME_NAME + "_pangenome.csv",
		expand("02_map_results/{sample}_{organism}.csv.bz2", organism=PANGENOME_NAME, sample=SAMPLE_PREFIX)

################################################################################

# map to markers with bowtie2
rule bowtie2:
    input: rules.merge.output
	output: "02_map_results/{sample}_" + PANGENOME_NAME + ".csv.bz2"
	params:
		pangenome_name = PANGENOME_NAME,
		pangenome_path = PANGENOME_PATH,
		out_folder = "02_map_results/{sample}_" + PANGENOME_NAME + ".csv"
	resources:
		time=6,
		mem=8
	shell:
		"bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i>} [-S <sam>]"
		"panphlan_map.py -c {params.pangenome_name} --i_bowtie2_indexes {params.pangenome_path} -i {input} -o {params.out_folder} --verbose"

rule panphlan_profile:
	input: expand("02_map_results/{sample}_{organism}.csv.bz2", organism=PANGENOME_NAME, sample=SAMPLE_PREFIX)
	output: PANGENOME_NAME + "_pangenome.csv"
	params:
		pangenome_name = PANGENOME_NAME,
		pangenome_path = PANGENOME_PATH
	resources:
		time=lambda wildcards, attempt: 12 * attempt,
		mem=250
	shell:
		"panphlan_profile.py -c {params.pangenome_name} --i_bowtie2_indexes {params.pangenome_path} -i 02_map_results --o_dna {output} --verbose"
