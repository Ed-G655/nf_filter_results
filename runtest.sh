ref_dir="test/data/ref_targets/"
alt_dir="test/data/alt_targets/"
output_directory="$(dirname $ref_dir)/results"

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
&& echo -e "======\n Basic pipeline TEST SUCCESSFUL \n======"
	#-stub-run \
