library(tidyverse)
library(optparse)
parser = OptionParser()
parser = add_option(parser,c("-f","--input-folder"),dest="input_folder",help="Input folder")
parser = add_option(parser,c("-o","--output-rds"),dest="output_rds",help="Output RDS")
args = parse_args(parser)
if(is.null(args$input_folder)){
    stop("No input folder given please specify it using the -f or --input-folder parameters")
}
#list(out_slopes=out_df, pi_between_df=tile_df2, mutation_df=mutation_df,site_frequency_spectrum=mac,parameter_in=parameter_in)
rds = list.files(args$input_folder, pattern="*.rds",full.names=T)
base_names = basename(rds)
output_mutation_list = list()
output_parameters = list()
slopes_list=list()
pi_between_dfs=list()
j = 1
base_names_vector = c()
for(i in 1:length(rds)){
	x = readRDS(rds[i])
    if(is.null(x$parameter_in)){
           print(rds[i])
           print("HERE")
           next
    }
	output_mutation_list[[j]] = x$mutation_df
	output_parameters[[j]] =  x$parameter_in
	slopes_list[[j]] = x$out_slopes
	pi_between_dfs[[j]] = x$pi_between_df
    base_names_vector= c(base_names_vector,  base_names[i])
    j = j + 1
}

names(output_mutation_list) =base_names_vector
print(length(output_parameters))
names(output_parameters) = base_names_vector 
names(slopes_list) = base_names_vector
names(pi_between_dfs) = base_names_vector
out_list = list(output_mutation_list=output_mutation_list, output_parameters=output_parameters, slopes_list=slopes_list, pi_between_dfs=pi_between_dfs)
saveRDS(out_list,file=args$output_rds)



