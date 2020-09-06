# Copyright (c) 2020 Boocock James <james.boocock@otago.ac.nz>
# Author: Boocock James <james.boocock@otago.ac.nz>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import argparse
import glob
import random
import subprocess
import random 
from string import Template
# greg lange
MUTATION_RATE=3.8e-10
## BYXRM
RECOMBINATION_RATE=3.133483e-06
#Kruglyak 
TOTAL_OUTCROSSING=2e-05
### Popsize 10,000,000
POPULATION_SIZE=10e6
### Total generations ###
TOTAL_GENERATIONS=3.2e9

ALT_DOMINANCE = 0.0001
REF_DOMINANCE = 0.98

PULSE_GENERATION_MIN=5e7

def get_parameters(sim_type, out_sim, factor):

    mapping = {}
    mapping["mutation_rate"]  = MUTATION_RATE
    mapping["recombination_rate"] = RECOMBINATION_RATE 
    mapping["population_size"] = int(POPULATION_SIZE)
    mapping["total_generations"] = int(TOTAL_GENERATIONS)
    mapping["sim_type"] = sim_type
    mapping["burnin"] = 8 * int(POPULATION_SIZE/factor)
    mapping["output_sim"] = out_sim + ".sim"
    mapping["output_par"] = out_sim + ".par"
    pulse =int(TOTAL_GENERATIONS/factor -  random.randint(int(1/factor)+ 1,PULSE_GENERATION_MIN/factor))
   # pulse_percentage = random.randint(1,100e7/factor)
   # TODO: Simulate the migration rate also
    pulse_percentage = random.uniform(.5,1.0)
    mapping["pulse_generation"] = pulse + mapping["burnin"]
    mapping["pulse_percentage"] = pulse_percentage
    mapping["pulse_generation_plus_one"] = mapping["pulse_generation"] + 1 
    mapping["max_generations"] = int(TOTAL_GENERATIONS/factor) + mapping["burnin"]
    cloning = random.uniform(1-1/10,1-1/2000) 
    selfing = 1 - TOTAL_OUTCROSSING /(1-cloning) 
    #### 
    mapping["selfing_rate"]   = selfing
    mapping["cloning_rate"] = cloning 

    ### Mitotic recombination
    mapping["alt_advantage"] = random.uniform(2e-4,3e-1)
    mapping["ref_advantage"] = random.uniform(2e-4,3e-1)
    mapping["mitotic_recomb_rate"] = random.uniform(0,1/10.0)
    mapping["factor"] = factor
    mapping["alt_dominance"] = ALT_DOMINANCE
    mapping["ref_dominance"] = REF_DOMINANCE 
    return(mapping)

def read_parameter_file(parameter_in,factor, out_sim):
    mapping = {}
    with open(parameter_in) as in_f:
        for i, line in enumerate(in_f):
            if i > 0:
                l_s = line.split()
                mapping[l_s[0]] = l_s[1]
    mapping["pulse_generation"] = int(mapping["pulse_generation"]) * int(mapping["factor"])
    mapping["factor"] = factor
    mapping["burnin"] = 8 * int(POPULATION_SIZE/factor)
    mapping["max_generations"] = int(TOTAL_GENERATIONS/factor) + mapping["burnin"]
    mapping["pulse_generation"] = int(mapping["pulse_generation"]/factor) 
    mapping["pulse_generation_plus_one"] = mapping["pulse_generation"] + 1 
    mapping["output_sim"] = out_sim + ".sim"
    mapping["output_par"] = out_sim + ".par"
    print(mapping)
    return(mapping)

def main():
    parser = argparse.ArgumentParser(description="Create simulationn template")
    parser.add_argument("-t","--template",dest="template_in",help="Template input file")
    parser.add_argument("--simulation-type",dest="simulation_type",help="Simulation type (burnin,introgression,balancing_selection,intro_balancing)",default="introgression")
    parser.add_argument("-f","--factor",dest="factor",default=10000)
    parser.add_argument("--tree-sequence",dest="tree_sequence",default=False,action="store_true")
    parser.add_argument("-o","--output-sim",dest="output_sim",default="test.sim",required=True)
    parser.add_argument("-p","--parameter-file",dest="parameter_file",help="Runs with a parameter file")
    parser.add_argument("-n","--num-sims",dest="number_of_sims",help="Number of simulations")
    parser.add_argument("-r","--random",dest="random",help="Random input",action="store_true",default=False)
    #parser.add_argument("--pulse-generation-min",dest="
    args = parser.parse_args()
    simulation_type = args.simulation_type
    template = args.template_in
    factor = int(args.factor)
    output_sim= args.output_sim
    tree_sequence = args.tree_sequence

    random_in = args.random
    import random
    if random_in:
        simulation_type = random.choice(["introgression","balancing_selection","balancing_selection_introgression"])
        print(simulation_type)
    if args.parameter_file is not None:
        mapping = read_parameter_file(args.parameter_file, factor, output_sim + ".out")
    else:
        mapping = get_parameters(simulation_type,output_sim + ".out", factor)
    if simulation_type == "introgression":
        slim_template = Template(open("slim/introgression_v1.slim","r").read())
    elif simulation_type == "balancing_selection":
        slim_template = Template(open("slim/multilocus_bal_v1.slim","r").read())
    elif simulation_type == "balancing_selection_introgression":
        slim_template = Template(open("slim/multilocus_bal_intro_v1.slim","r").read())
    with open(output_sim +".par","w") as  out_f:
        out_f.write("parameter\tvalue\n")
        for key, value in mapping.items():
            out_f.write(key + "\t" + str(value) + "\n")
    slim_output = output_sim  +  ".slim"
    with open(output_sim +".slim","w") as out_f:
        out_f.write(slim_template.substitute(mapping))
    import subprocess
    slim_log = output_sim + ".log"
    subprocess.check_call("/u/home/s/smilefre/bin/slim -s 0 {slim_output} > {slim_log} 2>&1".format(slim_output=slim_output,slim_log=slim_log),shell=True) 

if __name__ == "__main__":
    main()

