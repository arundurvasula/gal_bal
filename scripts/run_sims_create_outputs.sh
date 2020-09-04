#!/usr/bin/env bash
#$ -l h_rt=4:00:00,h_data=8G -pe shared 1
#$ -cwd
#$ -M theboocock@ucla.edu
#$ -m a
#$ -e logs/
#$ -o logs/

SIM_NUMBER=${SGE_TASK_ID}
OUTPUT_SIM=sims/$SIM_NUMBER

export PATH=/u/home/s/smilefre/project-kruglyak/anaconda3/bin:$PATH
source activate balancing_selection
/u/home/s/smilefre/project-kruglyak/anaconda3/bin/python simulation_template.py --output-sim $OUTPUT_SIM -f 50000 --simulation-type balancing_selection
Rscript process_summary_stats.R  --input-sim ${OUTPUT_SIM}.out.sim  --parameter-file ${OUTPUT_SIM}.par
