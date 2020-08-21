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
RECCOMBINATION_RATE=3.133483e-06
#Kruglyak 
TOTAL_OUTCROSSING=1/200000.0
### Popsize 10,000,000
POPULATION_SIZE=10e6
### Total generations ###
TOTAL_GENERATIONS=3.2e9


def main():
    parser = argparse.ArgumentParser(description="Create simulationn template")
    parser.add_argument("-t","--template",dest="template_in",help="Template input file")
    parser.add_argument("--simulation-type",dest="simulation_type",help="Simulation type (burnin,introgression,balancing_selection,intro_balancing)",default="introgression")
    parser.add_argument("-f","--factor",dest="factor",default=10000)
    parser.add_argument("-o","--output-sim",dest="output_sim",default="test.sim")
    args = parser.parse_args()
    simulation_type = args.simulation_type
    template = args.template_in
    factor = int(args.factor)
    output_sim= args.output_sim
    if simulation_type == "introgression":
        slim_template = Template(open("slim/introgression.slim","r").read())
        
        mapping = {}
        mapping["mutation_rate"] = MUTATION_RATE
        mapping["recombination_rate"] = RECCOMBINATION_RATE
        mapping["population_size"] = int(POPULATION_SIZE/factor)
        ### ### 
        ## cloning rate ##
        mapping["max_generations"] = int(TOTAL_GENERATIONS/factor)
        mapping["output_sim"] = output_sim + ".out" 
        ### Generation ###
        pulse = int(TOTAL_GENERATIONS/factor -  random.randint(1,100e7/factor))
        pulse_percentage = random.uniform(.01,.5)
        mapping["pulse_generation"] = pulse
        mapping["pulse_generation_plus_one"] = pulse + 1
        mapping["pulse_percentage"] = pulse_percentage
        mapping["factor" ] = factor

        ## #### 
        cloning = random.uniform(0,1-TOTAL_OUTCROSSING)
        selfing = 1 -TOTAL_OUTCROSSING - cloning 
        mapping["selfing_rate"]   = selfing
        mapping["cloning_rate"] = cloning 
        mapping["burnin"] = 8 * int(POPULATION_SIZE/factor)

        #### 
        with open(output_sim +".par","w") as  out_f:
            out_f.write("parameter\tvalue\n")
            for key, value in mapping.items():
                out_f.write(key + "\t" + str(value) + "\n")
        slim_output = output_sim  +  ".slim"
        with open(output_sim +".slim","w") as out_f:
            out_f.write(slim_template.substitute(mapping))
        import subprocess
        slim_log = output_sim + ".log"
        subprocess.check_call("slim -s 0 {slim_output} > {slim_log} 2>&1".format(slim_output=slim_output,slim_log=slim_log),shell=True) 

    elif simulation_type == "burnin":
        slim_template = Template(open("slim/burnin.slim","r").read())

if __name__ == "__main__":
    main()

