#!/usr/bin/env nextflow
/*================================================================

								---- MODULE PIPELINE ---------

/*================================================================
The Aguilar Lab presents...

- A pipeline to classify SNPs in microRNA regions and provide an overview of
diseases associated with microRNAs that present SNPs

==================================================================
Version: 0.1
Project repository:
==================================================================
Authors:

- Bioinformatics Design
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)



- Bioinformatics Development
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)


- Nextflow Port
 Jose Eduardo Garcia-Lopez (jeduardogl655@gmail.com)

///////////////////////////////////////////////////////////////

  Define pipeline Name
  This will be used as a name to include in the results and intermediates directory names
*/

pipeline_name = "nf-compare-miRNome"

/*This directories will be automatically created by the pipeline to store files during the run
*/
results_dir = "${params.output_dir}/${pipeline_name}-results/"
intermediates_dir = "${params.output_dir}/${pipeline_name}-intermediate/"

/*================================================================/*

/* MODULE START */


process CAT_TOOLS_TARGETS {

	publishDir "${intermediates_dir}/cat-targets/",mode:"symlink"

	input:
	file TARGETS

	output:
	file "${params.output_name}"

	shell:
  """
 cat *.tsv \
 | grep -P -v "a_Gene_ID\tmiRNA_ID\tUTR_start\tUTR_end\tSite_type\tchrom\ttarget_ID\tprediction_tool" > all.tmp
 echo "a_Gene_ID\tmiRNA_ID\tUTR_start\tUTR_end\tSite_type\tchrom\ttarget_ID\tprediction_tool" > header.txt
 cat header.txt all.tmp > "${params.output_name}"


	"""
}
