#Load in your code
source("C:/Users/cleob/Gitub/PopGen25/scr/Babor")

#Load in test data

test_data=read.table("C:/Users/cleob/Gitub/PopGen25.CB/Test/test_data.txt")
pops = read.table("C:/Users/cleob/Gitub/PopGen25.CB/Test/pop_file.txt")
#These are tests

freq_GT = function(site,ploidy=2){
  clean_site = na.omit(site)
  p = sum(clean_site)/(ploidy*length(clean_site)) #How will ploidy go in here?
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

#To furthur accomodate, we need to change the d_ij for haploid data
d_ij = function(site_i,site_j){
    p_i = freq_GT(site_i,1)
    p_j = freq_GT(site_j,1)
    p_ij = freq_ij(site_i,site_j)
    D_ij = p_ij - p_i*p_j
    return(D_ij)
}

#R^2 
cor_mat = cor(t(GT))

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
