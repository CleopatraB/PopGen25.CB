

#Load in your code
#ource("C:/Users/cleob/Gitub/PopGen25/scr/Babor")

#Load in test data

test_data=read.table("C:/Users/cleob/Gitub/PopGen25.CB/Test/test_data.txt")
pops = read.table("C:/Users/cleob/Gitub/PopGen25.CB/Test/pop_file.txt")
#These are tests

freq_GT = function(site,ploidy=2){
  clean_site = na.omit(site)
  p = sum(clean_site)/(ploidy*length(clean_site))
  return(p)
}
#Assuming diploid
allele_freqs = apply(test_data,1,freq_GT)
if(all(allele_freqs == c(0,0.5,0.05,0.25,0.3))){
    print("Diploid freqs ok.")
} else {
    print("Diploid freqs broken!")
}
#Assuming haploid
allele_freqs_hap = apply(test_data,1,function(x) freq_GT(x,1))
if(all(allele_freqs_hap == c(0,1,0.1,0.5,0.6))){
    print("Haploid freqs ok.")
} else {
    print("Haploid freqs broken!")
}

#And copy more or less the above for every function you write, 
#with output corresponding to what you _should_ get at each locus.

# test tree to make sure things are looking good
test_tree = your_function(some_data)
stopifnot(is.rooted(test_tree))

#Function to LD-prune a VCF. Makes calls to subfunctions for each chromosome, and each of those runs sub-functions on each window.
#For a variety of reasons, we drop _all_ alleles with no variation (e.g. allele frequency ==0 or 1).

ld_prune = function(vcf,window_size,cutoff=0.1){
  GT = extract_haps(vcf)
  freqs = apply(GT,1,function(x) freq_GT(x,ploidy=1))
  GT = GT[which(freqs*(1-freqs)>0),]
  chr_id = gsub(":.*:.*:.*","",rownames(GT))
  chrs = sort(unique(chr_id))
  sub_GT = lapply(chrs,function(x) ld_prune_chr(GT[chr_id==x,],window_size=window_size,cutoff=cutoff)) 
  return(sub_GT)
}

#Within each chromosome, generate windows and find only the sites to be dropped.
ld_prune_chr = function(GT,window_size,cutoff=0.1){
  pos = as.numeric(gsub(".*:([0-9]*):.*","\\1",rownames(GT)))
  breaks = unique(c(seq(from=min(pos)-1,to=max(pos),by=window_size),max(pos)))
  bins = cut(pos,breaks,right=TRUE)
  #Drop bins with 1 or fewer entries
  bad_bins = names(which(table(bins)<2))
  bins = factor(bins,exclude=bad_bins)
  drop = unlist(sapply(levels(bins),function(x) ld_prune_window(GT[which(bins==x),],cutoff=cutoff)))
  sub_GT = GT[setdiff(rownames(GT),drop),]
  return(sub_GT)
}

#Within each window, calculate correlations, keep only the first allele with high correlations
ld_prune_window = function(GT,cutoff=0.1){
  cor_mat = cor(t(GT))
  to_check = 1:dim(cor_mat)[1]
  drop = c()
  while(length(to_check)>0){
    high_cor = which(cor_mat[to_check[1],]>cutoff)
    drop = c(drop,setdiff(high_cor,to_check[1]))
    to_check = setdiff(to_check,high_cor)
  }
  if(length(drop)>0){
    return(rownames(GT)[drop])
  }
}

plink_path <- "C:/Users/cleob/AppData/Local/Temp/RtmpKynqf5/downloaded_packages"  # full path to your PLINK executable
vcf_file   <- "C:/Users/cleob/Gitub/PopGen25/labs/data/muc19_subsample.vcf.gz"  # your VCF file
output_prefix <- "sub_sample"     
getwd("C:/Users/cleob/Gitub/PopGen25/labs/data/muc19_subsample.vcf.gz")
ld_command <- paste(
  plink_path,
  "--vcf", vcf_file,
  "--double-id",
  "--allow-extra-chr",
  "--maf 0.01",        # remove rare SNPs
  "--geno 0.1",        # remove SNPs missing >10%
  "--mind 0.5",        # remove samples missing >50%
  "--chr 12",          # restrict to chromosome 12
  "--thin 0.1",        # keep 10% of SNPs at random
  "--r2",              # calculate pairwise LD
  "--ld-window 100",   # consider 100 SNPs per window
  "--ld-window-kb 1000", # max distance 1000kb
  "--ld-window-r2 0",  # report all r² values
  "--out", output_prefix
)
system(ld_command)
list.files(pattern = output_prefix)
log_file <- paste0(output_prefix, ".log")
if(file.exists(log_file)){
  log <- readLines(log_file)
  cat(tail(log, 10), sep = "\n")  # last 10 lines
} else {
  message("Log file not found. PLINK may not have run correctly.")
}
ld_file <- paste0(output_prefix, ".ld.gz")
if(file.exists(ld_file)){
  library(data.table)
  ld_data <- fread(cmd = paste("gzip -dc", ld_file))
  head(ld_data)
} else {
  message("LD file not found. Check PLINK output.")
}


require(ggplot2)
require(tidyverse)
admix=read_table("Gitub/PopGen25/labs/data/admixture_est.admix.16.Q",col_names=FALSE) # makes the admix table 

names(admix) = c("X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X12,X13,X14,XQ5,X16") #Give some putative names for the populations, depening on their number

admix= as.data.frame(admix)
admix$ind =  colnames(vcf@gt[,-1]) #How can you name each individual (hint: the vcf has the names of each of these)
#Pivot to long format- graphs the columns to long format
df_long = pivot_longer(admix,1:16,names_to="Pop",values_to="admix")

ggplot(df_long,aes(x=ind,y=admix,fill=Pop))+geom_col(col=NA)