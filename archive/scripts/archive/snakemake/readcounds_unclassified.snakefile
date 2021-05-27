DATA_DIRS = config["raw_reads_directory"]

rule readcounts_unclassified:
	input:
		readcount_files = expand(join(DATA_DIR, "{sample}_") + READ_SUFFIX[0] + EXTENSION, sample=SAMPLES),
	output:
		"readcounts_all.tsv"
	resources:
		time = 1
	run:
		outfile = str(output)
		if (os.path.exists(outfile)):
			os.remove(outfile)
		with open(outfile, 'w') as outf:
			outf.writelines('Sample\traw_reads\ttrimmed_reads\ttrimmed_frac\tdeduplicated_reads\tdeduplicated_frac\thost_removed_reads\thost_removed_frac\torphan_reads\torphan_frac\n')
			for sample in SAMPLE_PREFIX:
				raw_file = join(DATA_DIR, sample + "_") + READ_SUFFIX[0] + EXTENSION
				trimmed_file = join(PROJECT_DIR, "01_processing/01_trimmed/" + sample + "_") + READ_SUFFIX[0] + "_val_1.fq.gz"
				dedup_file = join(PROJECT_DIR, "01_processing/03_sync/" + sample + "_1.fq")
				rmhost_file = join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_1.fq")
				orphans_file = join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_orphans.fq")

				raw_reads = int(file_len(raw_file) / 4)
				trimmed_reads = int(file_len(trimmed_file) / 4)
				dedup_reads = int(file_len(dedup_file) / 4)
				rmhost_reads = int(file_len(rmhost_file) / 4)
				orphans_reads = int(file_len(orphans_file) / 4)

				trimmed_frac = round(trimmed_reads / float(raw_reads), 3)
				dedup_frac = round(dedup_reads / float(raw_reads), 3)
				rmhost_frac = round(rmhost_reads / float(raw_reads), 3)
				orphans_frac = round(orphans_reads / float(raw_reads), 3)

				line = '\t'.join([sample, str(raw_reads),
					str(trimmed_reads), str(trimmed_frac),
					str(dedup_reads), str(dedup_frac),
					str(rmhost_reads), str(rmhost_frac),
					str(orphans_reads), str(orphans_frac)])
				outf.writelines(line+ '\n')
