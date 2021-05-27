import re, os, subprocess
from os.path import join, expanduser, abspath

READ_SUFFIX = ['1', '2', 'orphans'] #['1', '2'] # or ['R1', 'R2']
EXTENSION = ".fq.gz"
DATA_DIR = config['reads']

# get file names
FILES = [f for f in os.listdir(DATA_DIR) if f.endswith(EXTENSION)]
SAMPLE_PREFIX = list(set([re.split('|'.join(['_' + a + '\.' for a in READ_SUFFIX]), i)[0] for i in FILES]))

rule all:
    input:
        # config['db_name'] + ".dmnd",
        expand("02_diamond_filter/{sample}.tsv", sample=SAMPLE_PREFIX)

# make diamond db
rule make_db:
    input: config['database']
    output: config['db_name'] + ".dmnd"
    params: db_name = config['db_name']
    shell:
        "diamond makedb --in {input} -d {params.db_name}"

# rule merge
rule concat:
    input:
        R1 = join(DATA_DIR, "{sample}_" + READ_SUFFIX[0] + EXTENSION),
        R2 = join(DATA_DIR, "{sample}_" + READ_SUFFIX[1] + EXTENSION),
        orp = join(DATA_DIR, "{sample}_orphans" + EXTENSION)
    output: "00_concat/{sample}.fq.gz"
    shell:
        "cat {input} > {output}"

# align reads to db
rule diamond_align:
    input:
        db = rules.make_db.output,
        fq = rules.concat.output
    output: "01_diamond/{sample}.tsv"
    resources:
        mem=lambda wildcards, attempt: 16 * attempt,
        time=lambda wildcards, attempt: 12 * attempt
    shell:
        "diamond blastx -d {input.db} -q {input.fq} --out {output} --outfmt tab"

rule diamond_filter:
    input: rules.diamond_align.output
    output: "02_diamond_filter/{sample}.tsv"
    resources:
        mem=lambda wildcards, attempt: 8 * attempt,
        time=lambda wildcards, attempt: 1 * attempt
    script:
        "cazy_filter.R"

# rule combine:
#     input: expand("01_diamond/{sample}.tsv", sample=SAMPLE_PREFIX)
#     output: "02_final/cazy_all.tsv"
#     shell:
#         "cp {input} 02_final/; "
# # convert to tabular output
# rule diamond_align:
#     input:
#     output:
#     shell:
#         "diamond view -a 8.diamond/output/foo.daa -o 8.diamond/tabular/foo.m8"
