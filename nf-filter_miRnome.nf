#!/usr/bin/env nextflow

/*================================================================
The Aguilar Lab presents...

- A pipeline to compare microRNA targets from targetScan and miRmap
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

Pre-processing:

Core-processing:


Pos-processing
Pos1-COMPARE_TARGETS

Analysis:

================================================================*/

/* Define the help message as a function to call when needed *//////////////////////////////
def helpMessage() {
	log.info"""
  ==========================================
	The ${pipeline_name}
  v${version}
  ==========================================

	Usage:

	nextflow run ${pipeline_name}.nf --ref_dir <path to input 1> --alt_dir <path to input 2> [--output_dir path to results ]

		--ref_dir	<- Directory whith the mirmap files;

	  --alt_dir	<- Directory whith the targetScan files;

	  --output_dir     <- directory where results, intermediate and log files will be stored;
	      default: same dir where --ref_dir resides

	  -resume	   <- Use cached results if the executed project has been run before;
	      default: not activated
	      This native NF option checks if anything has changed from a previous pipeline execution.
	      Then, it resumes the run from the last successful stage.
	      i.e. If for some reason your previous run got interrupted,
	      running the -resume option will take it from the last successful pipeline stage
	      instead of starting over
	      Read more here: https://www.nextflow.io/docs/latest/getstarted.html#getstart-resume
	  --help           <- Shows Pipeline Information
	  --version        <- Show version
	""".stripIndent()
}

/*//////////////////////////////
  Define pipeline version
  If you bump the number, remember to bump it in the header description at the begining of this script too
*/
version = "0.1"

/*//////////////////////////////
  Define pipeline Name
  This will be used as a name to include in the results and intermediates directory names
*/
pipeline_name = "nf-filter_miRnome"

/*
  Initiate default values for parameters
  to avoid "WARN: Access to undefined parameter" messages
*/
params.ref_dir = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.alt_dir = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.help = false //default is false to not trigger help message automatically at every run
params.version = false //default is false to not trigger version message automatically at every run

/*//////////////////////////////
  If the user inputs the --help flag
  print the help message and exit pipeline
*/
if (params.help){
	helpMessage()
	exit 0
}

/*//////////////////////////////
  If the user inputs the --version flag
  print the pipeline version
*/
if (params.version){
	println "${pipeline_name} v${version}"
	exit 0
}

/*//////////////////////////////
  Define the Nextflow version under which this pipeline was developed or successfuly tested
  Updated by iaguilar at MAY 2021
*/
nextflow_required_version = '20.01.0'
/*
  Try Catch to verify compatible Nextflow version
  If user Nextflow version is lower than the required version pipeline will continue
  but a message is printed to tell the user maybe it's a good idea to update her/his Nextflow
*/
try {
	if( ! nextflow.version.matches(">= $nextflow_required_version") ){
		throw GroovyException('Your Nextflow version is older than Pipeline required version')
	}
} catch (all) {
	log.error "-----\n" +
			"  This pipeline requires Nextflow version: $nextflow_required_version \n" +
      "  But you are running version: $workflow.nextflow.version \n" +
			"  The pipeline will continue but some things may not work as intended\n" +
			"  You may want to run `nextflow self-update` to update Nextflow\n" +
			"============================================================"
}

/*//////////////////////////////
  INPUT PARAMETER VALIDATION BLOCK
*/

/* Check if the input directory is provided
    if it was not provided, it keeps the 'false' value assigned in the parameter initiation block above
    and this test fails
*/
if ( !params.ref_dir | !params.alt_dir ) {
  log.error " Please provide the --ref_dir AND --alt_dir \n\n" +
  " For more information, execute: nextflow run ${pipeline_name} --help"
  exit 1
}

/*
Output directory definition
Default value to create directory is the parent dir of --input_dir
*/
params.output_dir = file(params.ref_dir).getParent() //!! maybe creates bug, should check

/*
  Results and Intermediate directory definition
  They are always relative to the base Output Directory
  and they always include the pipeline name in the variable pipeline_name defined by this Script

  This directories will be automatically created by the pipeline to store files during the run
*/
resulalt_dir = "${params.output_dir}/${pipeline_name}-results/"
intermediates_dir = "${params.output_dir}/${pipeline_name}-intermediate/"

/*
Useful functions definition
*/

/*//////////////////////////////
  LOG RUN INFORMATION
*/
log.info"""
==========================================
The ${pipeline_name} pipeline
v${version}
==========================================
"""
log.info "--Nextflow metadata--"
/* define function to store nextflow metadata summary info */
def nfsummary = [:]
/* log parameter values beign used into summary */
/* For the following runtime metadata origins, see https://www.nextflow.io/docs/latest/metadata.html */
nfsummary['Resumed run?'] = workflow.resume
nfsummary['Run Name']			= workflow.runName
nfsummary['Current user']		= workflow.userName
/* string transform the time and date of run start; remove : chars and replace spaces by underscores */
nfsummary['Start time']			= workflow.start.toString().replace(":", "").replace(" ", "_")
nfsummary['Script dir']		 = workflow.projectDir
nfsummary['Working dir']		 = workflow.workDir
nfsummary['Current dir']		= workflow.launchDir
nfsummary['Launch command'] = workflow.commandLine
log.info nfsummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "\n\n--Pipeline Parameters--"
/* define function to store nextflow metadata summary info */
def pipelinesummary = [:]
/* log parameter values beign used into summary */
pipelinesummary['Input REF dir']			= params.ref_dir
pipelinesummary['Input ALT dir']			= params.alt_dir
pipelinesummary['Results Dir']		= resulalt_dir
pipelinesummary['Intermediate Dir']		= intermediates_dir
/* print stored summary info */
log.info pipelinesummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================\nPipeline Start"

/*//////////////////////////////
  PIPELINE START
*/

/* enable DSL2*/
nextflow.enable.dsl=2

/*
	READ GENERAL INPUTS
*/

// Load mirmap REF data
Channel
	.fromPath( "${params.ref_dir}*.mirmapout" )
	.set{ mirmap_ref_input }

// Load mirmap alt data
	Channel
		.fromPath( "${params.alt_dir}*.alt.mirmapout" )
		.set{ mirmap_alt_input }

//load targetScan REF data
Channel
	.fromPath( "${params.ref_dir}*.tsout" )
	.set{ ts_ref_input }


	//load targetScan REF data
	Channel
		.fromPath( "${params.alt_dir}*.alt.tsout" )
		.set{ ts_alt_input }


/*
 Load R fileS
*/

/* R_script_1*/
Channel
	.fromPath( "./modules/filter_targets/results_filter.R" )
	.set{ R_script_1}

	/* R_script_2*/
	Channel
		.fromPath( "./modules/resume_report/resume_transcipts.R" )
		.set{ R_script_2}

/* LOAD canonical targets */
			Channel
				.fromPath( "./modules/filter_targets/Emsembl_Canonical.tsv" )
				.set{ EMSEMBL_TARGETS}


	/*	  Import modules */

								/*POS-processing */
include{FILTER_RESULTS as FILTER_RESULTS_REF_MP} from './modules/filter_targets/main.nf' addParams(type: 'REF')
include{FILTER_RESULTS as FILTER_RESULTS_REF_TS} from './modules/filter_targets/main.nf' addParams(type: 'REF')

include{FILTER_RESULTS as FILTER_RESULTS_ALT_TS} from './modules/filter_targets/main.nf' addParams(type: 'ALT')
include{FILTER_RESULTS as FILTER_RESULTS_ALT_MP} from './modules/filter_targets/main.nf' addParams(type: 'ALT')


include{RESUME_REPORT as RESUME_REPORT_REF_MP} from './modules/resume_report/main.nf' addParams(type: 'REF_MP')
include{RESUME_REPORT as RESUME_REPORT_REF_TS} from './modules/resume_report/main.nf' addParams(type: 'REF_TS')
include{RESUME_REPORT as RESUME_REPORT_ALT_TS} from './modules/resume_report/main.nf' addParams(type: 'ALT_MP')
include{RESUME_REPORT as RESUME_REPORT_ALT_MP} from './modules/resume_report/main.nf' addParams(type: 'ALT_TS')


/*  main pipeline logic */
workflow  {
/* pos-processing */
// Define function to get chrom
def get_chrom = { file -> file.baseName.replaceAll(/.alt/,"")}

					// collect targets outputs
					TARGETSCAN_REF = ts_ref_input.map{file -> tuple(file.baseName, file) }
					TARGETSCAN_ALT = ts_alt_input.map{ file -> tuple(get_chrom(file), file) }

						MIRMAP_REF = mirmap_ref_input.map{file -> tuple(file.baseName, file) }
						MIRMAP_ALT = mirmap_alt_input.map{ file -> tuple(get_chrom(file), file) }

						// // join channels
						// REF_TARGETS_TOOLS = TARGETSCAN_REF.join(MIRMAP_REF)
						// ALT_TARGETS_TOOLS = TARGETSCAN_ALT.join(MIRMAP_ALT)


					// 	FILTER_RESULTS TARGETS
					FILTER_RESULTS_REF_MP(MIRMAP_REF, EMSEMBL_TARGETS, R_script_1)
					FILTER_RESULTS_REF_TS(TARGETSCAN_REF, EMSEMBL_TARGETS, R_script_1)

					FILTER_RESULTS_ALT_TS(MIRMAP_ALT, EMSEMBL_TARGETS, R_script_1)
					FILTER_RESULTS_ALT_MP(TARGETSCAN_ALT, EMSEMBL_TARGETS, R_script_1)


					// 	RESUME TARGETS
					RESUME_REPORT_REF_MP(FILTER_RESULTS_REF_MP.out.LOG, R_script_2)
					RESUME_REPORT_REF_TS(FILTER_RESULTS_REF_TS.out.LOG, R_script_2)
					RESUME_REPORT_ALT_TS(FILTER_RESULTS_ALT_TS.out.LOG, R_script_2)
					RESUME_REPORT_ALT_MP(FILTER_RESULTS_ALT_MP.out.LOG, R_script_2)

}
