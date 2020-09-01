#!/usr/bin/env bash

for factor in {5000,10000,20000,30000,40000,50000,100000}
do
   # python simulation_template.py -p parameter_presets/balancing_selection_v1.txt  -f $factor -o sims/factor_$factor --simulation-type balancing_selection
    Rscript process_summary_stats.R  --input-sim sims/factor_${factor}.out.sim
done 
