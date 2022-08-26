#!/usr/bin/env nextflow
/*================================================================

								---- MODULE PIPELINE ---------

/*================================================================
The Aguilar Lab presents...

- A pipeline to extract and create miRNA and 3'UTR consensus sequences for analysis
   with targetscan and miRmap.

==================================================================
Version: 0.2
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
nextflow.enable.dsl=2

/*This directories will be automatically created by the pipeline to store files during the run
*/
results_dir = "${params.output_dir}/${pipeline_name}-results/"
intermediates_dir = "${params.output_dir}/${pipeline_name}-intermediate/"

/*================================================================/*

/* MODULE START */

process FILTER_RESULTS {
	tag "$CHR"

	publishDir "${results_dir}/filter-results/${params.type}", mode:"copy"

	input:
	tuple val(CHR), file (TARGETS)
	each EMSEMBL_TARGETS
  each Rscript

	output:
	file "*.filtered*"
	tuple val(CHR), path("*.log"), emit: LOG

	shell:
	"""
	  Rscript --vanilla ${Rscript} ${CHR} ${EMSEMBL_TARGETS} ${TARGETS}

 	"""

}
