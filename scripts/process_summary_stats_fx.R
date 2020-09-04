### Diversity stats
library(GenomicRanges)

diversity_pi = function(x){
    one = sum( x== 1)
    zero = sum(x==0)
    all = one + zero
    number_comparisons = (all * (all-1)/2)
    if(number_comparisons == 0){return(0)}
    div = (one * zero)/number_comparisons 
	#print(div)
    return(div)
}

library(Matrix)
extract_tempIDs <- function(mutString, rarePos){
  tids <- as.numeric(strsplit(mutString, split = " ", fixed = TRUE)[[1]][-c(1:2)]) + 1
  #Subset colID by tids (tempID) and retain the colIDs that are non-zero,
  #i.e. the rare varaints
  rarePos[tids][rarePos[tids] > 0]
}

read_slim <- function(file_path,
                      keep_maf = 0,
                      recomb_map = NULL,
                      pathway_df = NULL,
                      recode_recurrent = TRUE, populations=2){
  #NOTE: Time to read file ~19 secs
  message("Reading Slim File")
  exDat = readLines(file_path)
  
  #The default output for Slim  is a .txt file with the following headings:
  # - OUT: Contains number of generations, type (i.e. "A" for autosome), and
  #   the file name
  #
  # - Populations:
  #   next line is pop description, "p1", pop size, and
  #   type (i.e. "H" hermaphroditic)
  #
  # - Mutations:
  #   each mutation gets a separate line with (1) temp id, (2) permanent id,
  #   (3) mutation type,  (4) base position, (5) selection coefficient, 0 for
  #   neutral model, (6) dominance coefficient (0.5 for neutral), (7) id of
  #   sub-population in which the mutation arose, (8) the generation in which
  #   the mutation arose, and (9) the prevalence of the mutation.
  #
  # - Individuals:
  #   lists which genomes belong to which individual, redundant since this
  #   information is also listed at the beginning of each genome.
  #
  # - Genomes:
  #   each line lists one of the genomes for an individual, i.e. individual 1's
  #   genomes are contained in lines 1 & 2, individual 2's genomes are contained
  #   in lines 3 & 4, etc. Each line begins with the genome id, followed by the
  #   type (i.e. "A" for autosome), and a list of mutations.  The mutations are
  #   identified by the temporary ID (tempID) specified in the Mutations output.
  
  #find heading locations
  PopHead <- which(exDat == "Populations:")
  MutHead <- which(exDat == "Mutations:")
  IndHead <- which(exDat == "Individuals:")
  GenHead <- which(exDat == "Genomes:")
  
  popCount <- as.numeric(unlist(strsplit(exDat[PopHead + 1], split = " "))[2]) * 2
  
  #-----------#
  # Mutations #
  #-----------#
  #NOTE: time to create Mutations data set ~4 secs
  message("Creating Mutations Dataset")
  
  #extract mutation data from slim's Mutation output
  #only retaining the tempID, position, and prevalence of each mutation
  MutOut <- do.call(rbind, strsplit(exDat[(MutHead + 1):(IndHead - 1)], split = " ", fixed = TRUE))
  MutData <- data.frame(tempID = as.numeric(MutOut[, 1]),
                        type = MutOut[, 3],
                        position = as.numeric(MutOut[, 4]),
                        #selCoef = as.numeric(MutOut[, 5]),
                        #domCoef = as.numeric(MutOut[, 6]),
                        #pop = MutOut[, 7],
                        #genNo = as.numeric(MutOut[, 8]),
                        prevalence = as.numeric(MutOut[, 9]),
                        stringsAsFactors = TRUE)
  
  #add 1 to temp ID so that we can easily associate mutations
  #to columns by default slim's first tempID is 0, not 1.
  MutData$tempID <- MutData$tempID + 1
  #First position in slim is 0, not 1
  MutData$position <- MutData$position + 1
  
  #calculate the derived allele frequency
  #If recode option is TRUE we calculate the derived
  #allele frequecy by first summing the prevelance
  #for the identical mutations and then dividing
  #by the population size
  if (recode_recurrent) {
    calc_afreq <- MutData %>%
      group_by(.data$type, .data$position) %>%
      mutate(afreq = sum(.data$prevalence)/(2*popCount))
    MutData$afreq <- calc_afreq$afreq
  } else {
    MutData$afreq <- MutData$prevalence/(2*popCount)
  }
  
  #order Mutation dataset by tempID, so that (later) we can order
  #the mutations in each haplotypes by increasing genomic position
  MutData <- MutData[order(MutData$tempID), ]
  
  #Identify the future sparseMatrix column ID of variants with <= keep_maf.
  #Variants with large maf are assigned column ID 0, i.e. thrown away.
  #Variants with sufficiently small maf are assigned increasing colIDs.
  #These are like new tempIDs, they do not reflect genomic position.
  MutData$colID <- cumsum(MutData$afreq <= keep_maf)*(MutData$afreq <= keep_maf)
  
  #Using the identified colID, create dataframe of rare mutations only
  RareMutData <- MutData[MutData$colID > 0, ]
  
  #-----------#
  # Genotypes #
  #-----------#
  #NOTE: jpos is the most computationally expensive task (uses strsplit)
  #this chuck takes ~2.2 minutes to run (for genome-wide, exon-only data)
  #This is an improvement from the old time: ~ 6 mins
  message("Creating Sparse Haplotypes Matrix")
  
  #determine future row and column position of each mutation listed in genomes
  #row will correspond to person, column will correspond to the tempID of the
  #mutation
  jpos <- lapply(1:(2*popCount), function(x){
    extract_tempIDs(mutString = exDat[GenHead + x],
                    rarePos = MutData$colID)
  })
  
  ipos <- lapply(1:length(jpos), function(x){
    rep(x, length(jpos[[x]]))
  })
  
  #create sparse matrix containing mutations(columns) for each individual(row)
  GenoData <- sparseMatrix(i = unlist(ipos),
                           j = unlist(jpos),
                           x = rep(1, length(unlist(jpos))),
                           dims = c(2*popCount, nrow(RareMutData)))
  
  #order by genomic postion of rare mutation
  GenoData <- GenoData[, order(RareMutData$position)]
  RareMutData <- RareMutData[order(RareMutData$position),]
  
  #Re-format tempID so that it corresponds to the column
  # (in GenoData) that the mutation is stored in
  RareMutData$colID <- 1:nrow(RareMutData)
  RareMutData <- RareMutData[, -1] #remove the old tempID
  
  #-----------------------------#
  # Re-code Identical Mutations #
  #-----------------------------#
  # Assuming that all mutations are of the same type, different mutations at
  # the same site are actually identical mutations from different lineages.
  # For simplicity, we recode these mutations so that they are only cataloged
  # once.
  if (recode_recurrent) {
    message("Recoding Identical Mutations")
    
    #For each mutation type, check to see if there
    #are any identical mutations to re-code
    if (any(sapply(unique(RareMutData$type),
                   function(x){any(duplicated(RareMutData$position[RareMutData$type == x]))}))) {
      
      #determine which mutation types have identical mutations that need
      #to be recoded
      rc_types <- unique(RareMutData$type)[sapply(unique(RareMutData$type),
                                                  function(x){any(duplicated(RareMutData$position[RareMutData$type == x]))})]
      
      #we have to use a loop to recode identical mutations by type
      #because recoding involves combining and removing columns from the
      #Haplotype matrix and then recoding the colIDs in the Mutations data frame.
      for (i in 1:length(rc_types)) {
        com_id_muts <- combine_identicalmutations(mutmap = RareMutData,
                                                  hapmat = GenoData,
                                                  mut_type = rc_types[i])
        RareMutData <- com_id_muts[[1]]
        GenoData <- com_id_muts[[2]]
      }
    }
  }
  
  
  
  
  #------------------#
  # Re-map Mutations #
  #------------------#
  if (!is.null(recomb_map)) {
    message("Remapping Mutations")
    RareMutData <- reMap_mutations(mutationDF = RareMutData,
                                   recomb_map)
  } else {
    RareMutData$chrom <- 1
  }
  
  row.names(RareMutData) = NULL
  
  #Create unique marker names
  RareMutData$marker <- make.unique(paste0(RareMutData$chrom, sep = "_", RareMutData$position))
  
  #reduce RareMutData, to the columns we actually need
  #really should clean this up soon
  RareMutData <- RareMutData[, c("colID", "chrom", "position",
                                 "afreq", "marker", "type")]
  
  # RareMutData <- RareMutData[, c("colID", "chrom", "position",
  #                                "afreq", "marker", "type",
  #                                "selCoef", "domCoef", "pop")]
  
  #----------------------#
  # Identify Pathway RVs #
  #----------------------#
  #if pathway data has been supplied, identify pathway SNVs
  #if (!is.null(pathway_df)) {
  #  message("Identifying Pathway SNVs")
  #  RareMutData <- identify_pathwaySNVs(markerDF = RareMutData,
  #                                      pathwayDF = pathway_df)
  #}
  
  
  return(SimRVSequences::SNVdata(Haplotypes = GenoData, Mutations = RareMutData))
}
library(arrangements)
calculate_eld_gal = function(genotypes, genomes_in,mac){
  all_p = permutations(n=3,k=3,replace = T) - 1
  y = genomes_in$Mutations[genomes_in$Mutations$type %in% c("m2","m3","m4"),]$colID
  if(any(mac[y] == 0)){
    print("LOST")
  }else{
    prob_all = apply(all_p, 1,function(x){sum(genotypes[,y[1]] == x[1] & genotypes[,y[2]] == x[2] & genotypes[,y[3]] == x[3])/nrow(genotypes)})
    prob_null = apply(all_p, 1,function(x){sum(genotypes[,y[1]] == x[1])/nrow(genotypes) * sum(genotypes[,y[2]] == x[2])/nrow(genotypes) * sum(genotypes[,y[3]] ==x[3])/nrow(genotypes)})
    prob_idx = prob_all != 0
    s_ld = -sum(prob_all[prob_idx] * log(prob_all[prob_idx]))
    prob_idx_null = prob_null != 0
    s_le = -sum(prob_null[prob_idx_null] * log(prob_null[prob_idx_null]))
    if(s_le == 0){
      return(1)
    }else{
      return(1-s_ld/s_le)
    }
  }
}

