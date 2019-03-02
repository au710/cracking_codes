#!/usr/bin/env nextflow

params.mode = "flow"
params.runs = "run_list.csv"
params.code = "strong_lim.py"
params.post = "para_rp.py"
params.container = '~/Containers/fen_crack.sif'
params.ncpus = 1

params.maxForks = 50 //max number running concurrently

Channel  // read the "run" column from file
    .fromPath(params.runs)
    .splitCsv(header:true)
    .map{ row-> tuple(row.run, row.active) }
    .set{ r_list }


    def datestamp() {
      now = new Date()
      stamp = now.format("yyyy-MM-dd_HHmmss", TimeZone.getTimeZone('UTC'))
      return stamp
    }

String results_dir = './RESULTS'

runs_file = file(params.runs)
code_file = file(params.code)
common_dir = file("common") // code here for extracting constants
post_dir = file(params.post) // post-processing codes
bin_dir = file('bin')


process simulation {

	maxForks params.maxForks
	errorStrategy 'ignore' // simulation will fail if error encountered

        // make sure edits don't influence main file
	stageInMode "copy"

	tag "simulation run $run"

	input:
	    file "run_list.csv" from runs_file
	    file "strong_lim.py" from code_file
	    file common from common_dir
            file "para_rp.py" from post_dir
            file bin from bin_dir
	    set val(run), val(active)  from r_list

	output:
            file "*.txt"
            file "*.png"
            file "*/*.{pvd,vtu}"
            file "*/*.txt"
	    file "constants.yml"

        String publish_dir = "${results_dir}/EXTRA2"
	publishDir "${publish_dir}/run_${run}", mode: "link"

        validExitStatus 0,1,6,134,139 //check what these are doing

        // this is the heart of the code - note how small it is
        script:
	"""

        # extract required parameters
	echo $run > b.txt
        singularity exec $params.container bin/extract.py --run=$run run_list.csv > constants.yml

        # run the simulation
        singularity exec $params.container python strong_lim.py > a.txt

        pvpython para_rp.py
        

	"""
}
