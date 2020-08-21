#!/usr/bin/env bash
#$ -l h_data=20G -pe shared 1
#$ -cwd
#$ -M theboocock@ucla.edu
#$ -m a
#$ -l time=24:00:00
#$ -e logs/
#$ -o logs/
SIM_NUMBER=${SGE_TASK_ID}
OUTPUT_SIM=sims/$SIM_NUMBER
python simulation_template.py --output-sim $OUTPUT_SIM -f 50000
