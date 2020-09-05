source("process_summary_stats_fx.R")
library(optparse)
library(tidyverse)
library(plyranges)
parser = OptionParser()
parser = add_option(parser,c("-i","--input-sim"),dest="input_sim")
parser = add_option(parser,c("-w","--window-size"),default=600,dest="window_size")
parser = add_option(parser,c("-s","--window-step"),default=30,dest="window_step")
parser = add_option(parser,c("-g","--genome-width"),default=30000,dest="genome_width", help="Genome width")
parser = add_option(parser,c("-f","--output-folder"),default="sim_output",dest="output_folder", help="Output folder")
parser = add_option(parser,c("-p","--parameter-file"), dest="parameter_file",help="Parameter file")
args = parse_args(parser)
print(args$input_sim)
print(args)
if(is.null(args$input_sim)){
	stop("No sim input given please specific using -i or --input-sim")
}
if(is.na(args$parameter_file)){
    stop("No parameter file give please specify it using -p or --parameters")
}
parameter_in = read.table(args$parameter_file, header=T, sep="\t")
window_width= as.numeric(args$window_size)
window_step = as.numeric(args$window_step)
genome_width = as.numeric(args$genome_width)
output_folder = args$output_folder
input_sim = args$input_sim
tile_genome = GRanges(seqnames="chr1",IRanges(c(seq(1,5000-window_width+ 1,by=window_step),seq(5001,10000-window_width+ 1,by=window_step),
          seq(10001,15000-window_width+1,by=window_step),seq(15001,20000-window_width,by=window_step),
          seq(20001,25000-window_width+1,by=window_step), seq(25001,30000-window_width+1,by=window_step)),
        end=c(seq(1 -1 + window_width,5000,by=window_step),seq(5001 - 1 + window_width,10000,by=window_step),
              seq(10001 -1 + window_width, 15000,by=window_step),seq(15001-1+window_width,20000,by=window_step),
              seq(20000 -1 + window_width,25000,by=window_step), seq(25001-1 + window_width,30000,by=window_step))))
### CHeck which alleles are present. ### 
# m2, m3,m4,m5,m6,m7.
### Must be missing some alleles
genomes_in2 = read_slim(input_sim,recode_recurrent = F,keep_maf = 1)
genomes_in2$Haplotypes = as.matrix(genomes_in2$Haplotypes)
genotypes = t(sapply(1:(nrow(genomes_in2$Haplotypes)/2),function(x) {genomes_in2$Haplotypes[(2*x-1),] +genomes_in2$Haplotypes[(2*x),]}))
ids  = genomes_in2$Mutations[genomes_in2$Mutations$type != "m1",] %>% arrange(type) %>% select(colID)
mac = apply(genomes_in2$Haplotypes,2,function(x){sum(x)})
#genotypes = genomes_in2$Haplotypes[(nrow(genomes_in2$Haplotypes)/2 +1):nrow(genomes_in2$Haplotypes),]
###
###
eld = calculate_eld_gal(genotypes, genomes_in2,mac)
### GAL statistics #### 
ids  = genomes_in2$Mutations[genomes_in2$Mutations$type != "m1",] %>% arrange(type) %>% select(colID)
mutation_df = genomes_in2$Mutations[ids$colID,]
mutation_df$mac = mac[ids$colID]
mutation_df$eld = eld
if(nrow(ids) != 6){
    print("Atleast one allele fixed")
diversity_all = apply(genomes_in2$Haplotypes[,-ids$colID],2,diversity_pi)
    windowed_ds = GRanges(seqnames="chr1",ranges=IRanges(genomes_in2$Mutations$position[-ids$colID],end=(genomes_in2$Mutations$position[-ids$colID])))
    windowed_ds$div = diversity_all
    pi_between_df = tile_genome %>% group_by_overlaps(windowed_ds) %>% summarise(pi_all = sum(div)/window_width) %>% as.data.frame
    input_base = basename(input_sim)
    mutation_list = list(mutation_df=mutation_df,pi_between_df=pi_between_df)    
    out_slope_summary_stats = paste(output_folder,"/", paste(input_base,".rds",sep=""),sep="")
    saveRDS(mutation_list, file=out_slope_summary_stats)
    quit()
}
input_out_mutation_stats = paste(input_sim,".mutation_stats",sep="")
#write.table(mutation_df, file=input_out_mutation_stats,quote=F,row.names=F,col.names=T,sep="\t")
sim_genome = GRanges(seqnames="chr1",ranges = IRanges(1,end=genome_width))
sim_genome = c(chr1=genome_width)
# 200 amino-acid windows 10 codon stop
# 600 dna windows 30 base-pair step
diversity_all = apply(genomes_in2$Haplotypes[,-ids$colID],2,diversity_pi)
#windowed_ds = GRanges(seqnames="chr1",ranges=IRanges(genomes_in2$Mutations$position[-ids$colID],end=(genomes_in2$Mutations$position[-ids$colID])))
#pi = tile_genome %>% group_by_overlaps(windowed_ds) %>% summarise(pi=sum(div)/window_width)

m1  = (genomes_in2$Haplotypes[,ids$colID[1]] == 1)
diversity1 = apply(genomes_in2$Haplotypes[m1,-ids$colID],2,diversity_pi)
diversity2 = apply(genomes_in2$Haplotypes[!m1,-ids$colID],2,diversity_pi)
windowed_ds = GRanges(seqnames="chr1",ranges=IRanges(genomes_in2$Mutations$position[-ids$colID],end=(genomes_in2$Mutations$position[-ids$colID])))
windowed_ds$div = diversity_all
windowed_ds$div1 = diversity1
windowed_ds$div2 =  diversity2
pi_ff = tile_genome %>% group_by_overlaps(windowed_ds) %>% summarise(pi1=sum(div1)/window_width,pi2=sum(div2)/window_width)
x = as.matrix(genomes_in2$Haplotypes[m1,-genomes_in2$Mutations$colID[genomes_in2$Mutations$type != "m1"]])
y = as.matrix(genomes_in2$Haplotypes[!m1,-genomes_in2$Mutations$colID[genomes_in2$Mutations$type != "m1"]])
rm_idx_v1 = genomes_in2$Mutations$colID[genomes_in2$Mutations$type == "m1"]
diversity_between_groups =sapply(1:ncol(y),function(z){
    one = sum(x[,z] == 1) 
    two = sum(y[,z] == 0) 
    total = (nrow(x) * (nrow(y) ))
    tmp_one = ((one *two) / total)
    one = sum(x[,z] == 0) 
    two = sum(y[,z] == 1) 
    total = (nrow(x) * (nrow(y) ))
    tmp_out = tmp_one + ((one *two) / total)
    return(tmp_out)
})
windowed_ds$div_between = diversity_between_groups
pi_between = tile_genome %>% group_by_overlaps(windowed_ds) %>% summarise(pi_between=sum(div_between)/window_width, pi1=sum(div1)/window_width, pi2=sum(div2)/window_width, pi_all = sum(div)/window_width)
out_sim_input = paste(input_sim,".diversity",sep="")
tile_df = as.data.frame(tile_genome)
tile_df$query = as.numeric(rownames(tile_df))
#tile_df %>% full_join(pi_between_df,by=c("query"="query")) %>% write.table(out_sim_input, quote=F,row.names=F,col.names=T,sep="\t")
pi_between_df = as.data.frame(pi_between)
### Get slopes 
tile_g2 = GRanges(seqnames="chr1",IRanges(c(1,5001,10001,15001,20001,25001),end=c(5000,10000,15000,20000,25000,30000)))
### Calculate eLD ###
df_out = tile_g2 %>% group_by_overlaps(tile_genome) %>% as.data.frame
pi_between_df$query = NULL
df_out = cbind(df_out,as.data.frame(pi_between_df))
out_df = data.frame()
for (i in unique(df_out$query)){
    df2 = df_out[df_out$query == i,]
    df2$row_number = 1:nrow(df2)
    m1 = ((lm(log(pi_between)~ row_number,data=df2)))
    intercept= (summary(m1)$coef[1,1])
    slope = (summary(m1)$coef[2,1])
    max_ds = max(df2$pi_between)
    p = data.frame(row_number=max(df2$row_number))
    predict1 = (predict(m1,p))
    #print(exp(predict1))
    #print(max_ds)
    tmp_row = data.frame(max_ds=max_ds,i=i, slope=slope,predict=exp(predict1),intercept=intercept)
    out_df = rbind(out_df, tmp_row)
}
input_base = basename(input_sim)
#nnprint(head(tile_df))
#print(head(pi_between_df))
pi_between_df = as.data.frame(pi_between)
tile_df2 = tile_df %>% full_join(pi_between_df,by=c("query"="query"))
mutation_list = list(out_slopes=out_df, pi_between_df=tile_df2, mutation_df=mutation_df,site_frequency_spectrum=mac,parameter_in=parameter_in)
out_slope_summary_stats = paste(output_folder,"/", paste(input_base,".rds",sep=""),sep="")
saveRDS(mutation_list, file=out_slope_summary_stats)

#plot((pi$pi),ylim=c(0,max((pi$pi)) + 0.05))
#points((pi_ff$pi1),col="red")
#points((pi_ff$pi2),col="blue")
#abline(v=c(10000/window_width*3,20000/window_width*3,60000/window_width*3))
#abline(v=c(5000/window_width *3 ,15000/window_width * 3,25000/window_width*3),col="green")
#plot(log(pi$pi),ylim=c(-8,max(log(pi$pi)) + 0.05))
#points(log(pi_ff$pi1),col="red")
#points(log(pi_ff$pi2),col="blue")
#abline(v=c(10000/window_width*3,20000/window_width*3,60000/window_width*3))
#abline(v=c(5000/window_width *3 ,15000/window_width * 3,25000/window_width*3),col="green")
#j = 0
#for(i in seq(1,300,by=50)){
  #print(i)
#  query = pi$query[1:50] * 100 - 99
#  aa = lm(log(pi$pi[i:(i+49)])~ query)
#  //  aa = glm(log(pi$pi[i:(i+49)]) ~ exp(1:50))
    
#aa = glm(rev(pi$pi[i:(i+50)]) ~ rev(exp(pi$query[i:(i+50)])))
# // }
#  print(summary(aa))
#  j =  j + 1
#  pred = (predict(aa,data.frame(x=query)))
#  names(pred) = i:(i+49)
#  lines(y=pred,x=names(pred))
  #query = pi_ff$query[1:50] * 100 - 99
#  pi_ff$pi1[pi_ff$pi1 == 0] = 1
#  pi_ff$pi2[pi_ff$pi2 == 0] = 1
  
#  aa = lm(log(pi_ff$pi1[i:(i+49)])~ query)
#  pred = (predict(aa,data.frame(x=query)))
#  names(pred) = i:(i+49)
#  lines(y=pred,x=names(pred),col="red")
  
#  aa = lm(log(pi_ff$pi2[i:(i+49)])~ query)
#  pred = (predict(aa,data.frame(x=query)))
#  names(pred) = i:(i+49)
#  lines(y=pred,x=names(pred),col="green")
  
# // names(pred) = 
#}
#aa = glm(pi$pi[1:50] ~ exp(pi$query[1:50]))
#plot(aa)
