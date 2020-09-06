#!/usr/bin/env bash
#$ -l h_data=8G -pe shared 1
#$ -cwd
#$ -M theboocock@ucla.edu
#$ -m a
#$ -l time=24:00:00
#$ -e logs/
#$ -o logs/
SIM_NUMBER=${SGE_TASK_ID}
mkdir -p sims_intro3
OUTPUT_SIM=sims_intro3/$SIM_NUMBER

export PATH=/u/home/s/smilefre/project-kruglyak/anaconda3/bin:$PATH
source activate balancing_selection
mkdir -p sim_output_intro3
/u/home/s/smilefre/project-kruglyak/anaconda3/bin/python simulation_template.py --output-sim $OUTPUT_SIM -f 50000 --simulation-type balancing_selection_introgression -r 
Rscript process_summary_stats.R  --input-sim ${OUTPUT_SIM}.out.sim  --parameter-file ${OUTPUT_SIM}.par -f sim_output_intro3
