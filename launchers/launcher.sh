##!/usr/bin/env bash
cd ../
ref_dir="/bodega/projects/miRNome_project/analisis/100GMX/targets/ref_targets/"
alt_dir="/bodega/projects/miRNome_project/analisis/100GMX/targets/alt_targets/"
output_directory="/bodega/projects/miRNome_project/canonical_transcripts/results"

echo -e "======\n Testing NF execution \n======" \
&& rm -rf $output_directory \
&& nextflow run nf-filter_miRnome.nf \
	--ref_dir $ref_dir \
  --alt_dir $alt_dir \
	--output_dir $output_directory \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html \
	-with-timeline $output_directory/`date +%Y%m%d_%H%M%S`_timeline.html \
  && echo -e "======\n  pipeline execution END \n======"
