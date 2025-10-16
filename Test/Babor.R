#Babor.R
#Set working directory to start with this notebook
#setwd("C:/Users/cleob/Gitub/PopGen25/src")
#source("src/Babor.R")

#Opening a VCF 
setwd("C:/Users/cleob/Gitub/PopGen25/labs/data")
require("vcfR")
require("readr")
samples = read_table("C:/Users/cleob/Gitub/PopGen25/labs/data/muc_samples.txt")
vcf = read.vcfR("C:/Users/cleob/Gitub/PopGen25/labs/data/muc19_subsample.vcf.gz")
GT = extract.gt(vcf,element="GT",as.numeric=TRUE)

#Finding an allele frequency
freq_GT = function(site){
  clean_site = na.omit(row)
  num_derived = sum(site)
  total_sites = 2*length(site)
  return(num_derived/total_sites)
}

#Sample Heterozygosity 
sample_het = function(site){
  n = length(site)*2
  p = freq_GT(site)
  h = n/(n-1)*2*p*(1-p)
  return(h)
}

#Calculating dxy 
dxy_site = function(p1,p2){
  dxy = p1*(1-p2)+p2*(1-p1)
  return(dxy)
}

#FST site
fst_site = function(p1,p2,n1,n2){
    f1 = n_ceu/(n_ceu+n_mxl)
    f2 = 
    num = f1*p1*(1-p1)+f2*p2*(1-p2)
    denom = (f1*p1+f2*p2)*(f1*(1-p1)+f2*(1-p2)) 
    fst = 1-num/denom
    return(fst)
}

#Process GT
process_GT = function(GT,samples){
    ids = gsub(":[A,T,C,G]*:[A,T,C,G]*","",rownames(GT))
    chr = gsub(":.*","",ids)
    pos = as.numeric(gsub(".*:","",ids))
    summary_df = data.frame(chr=chr,pos=pos)
    summary_df$tot_het = apply(GT,1,sample_het)
    pops = unique(samples$population)
    #Get indices for each population
    idx_pops = lapply(pops,function(x) {
        which(colnames(GT) %in% samples$sample_id[which(samples$population == x)])
        })
    #Name them to make retrieval easier.
    names(idx_pops) = pops
    #Check for NAs per site per pop
    obs_alleles = lapply(1:length(pops),function(x){
        apply(GT[,idx_pops[[x]]],1,function(x) length(!is.na(x)))
    })
    names(obs_alleles) = pops
    summary_df$obs = obs_alleles
    # Now, let's get allele frequencies         
    ps = lapply(1:length(pops),function(x) {
        apply(GT[,idx_pops[[x]]],1,freq_GT)
    })
    names(ps) = pops
    summary_df$p = ps
    # Heterozygosities
    hets = lapply(1:length(pops),function(x) {
        apply(GT[,idx_pops[[x]]],1,sample_het)
    })
    names(hets) = pops
    summary_df$het = hets
    #All possible pairs of populations:
    pop_pairs = combn(pops,2)
    #Calculate fst for each
    fst =lapply(1:dim(pop_pairs)[2],function(x) {
            pop1 = pop_pairs[1,x]
            pop2 = pop_pairs[2,x]
            n1 = obs_alleles[[pop1]]
            n2 = obs_alleles[[pop2]]
            p1 = ps[[pop1]]
            p2 = ps[[pop2]]
            return(fst_site(p1,p2,n1,n2))
    })
    #Add names to make comparisons clear
    names(fst) = sapply(1:dim(pop_pairs)[2],function(x) paste(pop_pairs[,x],collapse="-"))
    #And dxy
    dxy =lapply(1:dim(pop_pairs)[2],function(x){
        p1 = ps[[pop_pairs[1,x]]]
        p2 = ps[[pop_pairs[2,x]]]
        return(dxy_site(p1,p2))
    })
    #Adding names to dxy as well
    names(dxy) = sapply(1:dim(pop_pairs)[2],function(x) paste(pop_pairs[,x],collapse="-"))
    #Don't forget to add these to the dataframe
    summary_df$fst = fst
    summary_df$dxy = dxy
    return(summary_df)
}

#Calculating LD 
d_ij = function(site_i,site_j){
    p_i = freq_GT(site_i)
    p_j = freq_GT(site_j)
    p_ij = intersect(which(site_i==1),which(site_j==1))/length(site_i)
    D_ij = p_ij - p_i*p_j
    return(D_ij)
}
#Function to run on each row
allele_2_num = function(site,name){
  div_allele = gsub(".*:.*:([A,T,C,G]+):","",name)
  GT_row = as.numeric(site==div_allele)
  return(GT_row)
}

#Structure for whole vcf
extract_haps = function(vcf){
  haps = extract.haps(vcf)
  names = rownames(haps)
  num_mat = t(sapply(1:dim(haps)[1],function(x) allele_2_num(haps[x,],names[x])))
  rownames(num_mat) = names
  colnames(num_mat) = colnames(haps)
  return(num_mat)
}

freq_ij = function(site_i,site_j){
    #How can we make sure that we find every case where both i and j are equal to 1, _and_ we get the number of sites correct?
    n_ij = length(intersect(which(site_i==1),which(site_j==1)))
    n_sites = length(intersect(which(!is.na(site_i)),which(!is.na(site_j))))
    return(n_ij/n_sites)
}
#This freq GT inlcudes a ploidy arguement to handle diploid data 
freq_GT = function(site,ploidy=2){
  clean_site = na.omit(site)
  p = sum(clean_site)/(ploidy*length(clean_site)) #How will ploidy go in here?
  return(p)
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

#multiple ploidy prodiced multiple vectorssum(ploidy{which(!is.na(site_i))])

#Solution 1 Thin Data
random_rows = sample(1:dim(GT)[1],size=round(0.1*dim(GT)[1]),replace=FALSE)
thin_GT = GT[random_rows,]
cor_mat = cor(t(thin_GT))

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
ld_prune_chr = function(GT,window_size=250,cutoff=0.1){
  pos = as.numeric(gsub(".*:([0-9]*):.*","\\1",rownames(GT)))
  pos = pos[!is.na(pos)]
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