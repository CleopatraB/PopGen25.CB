#Finding an allele frequency
freq_GT = function(site, ploidy=1){
  clean_site = na.omit(site)
  sum(clean_site)/(ploidy*length(clean_site))
}
#Sample Heterozygosity 
sample_het = function(site, ploidy=1){
  n = length(site)*ploidy
  p = freq_GT(site, ploidy)
  h = n/(n-1) * 2 * p * (1-p)
}
#FST site
fst_site = function(p1, p2, n1, n2){
  f1 = n1/(n1+n2)
  f2 = n2/(n1+n2)
  num = f1*p1*(1-p1) + f2*p2*(1-p2)
  denom = (f1*p1 + f2*p2)*(f1*(1-p1) + f2*(1-p2))
  fst = 1 - num/denom
}
#dxy for each site
dxy_site = function(p1, p2){
  p1*(1-p2) + p2*(1-p1)
}

#GT
process_GT = function(GT, samples, ploidy=1){
  GT = as.matrix(GT)
  ids = gsub(":[A,T,C,G]*:[A,T,C,G]*", "", rownames(GT))
  chr = gsub(":.*", "", ids)
  pos = as.numeric(gsub(".*:", "", ids))
  
  pops = sort(unique(samples$population))
  pops_char = as.character(pops)
  
   #Get indices for each population
  idx_pops = lapply(pops, function(x){
    which(colnames(GT) %in% samples$samples[which(samples$population==x)])
  })
    #Name them to make retrieval easier.
  names(idx_pops) = pops_char
  
#Check for NAs per site per pop
  obs_alleles = lapply(idx_pops, function(idx){
    apply(GT[, idx, drop=FALSE], 1, function(x) sum(!is.na(x)))
  })
  
    # Now, let's get allele frequencies
  ps = lapply(idx_pops, function(idx){
    apply(GT[, idx, drop=FALSE], 1, function(x) freq_GT(x, ploidy=ploidy))
  })
  
  # Heterozygosities
  hets = lapply(idx_pops, function(idx){
    apply(GT[, idx, drop=FALSE], 1, function(x) sample_het(x, ploidy=ploidy))
  })
  
  #average r2 for each pop
  avg_r2 = sapply(idx_pops, function(idx){
    mat = GT[, idx, drop=FALSE]
    cor_mat = (cor(t(mat))^2)
    mean(cor_mat[upper.tri(cor_mat)], na.rm=TRUE)
  })
  names(avg_r2) = pops_char

    #All possible pairs of populations, addition if 
  pop_pairs = combn(pops_char, 2)
  if(is.vector(pop_pairs)) pop_pairs = matrix(pop_pairs, nrow=2)
  
    #Calculate fst for each
  fst = lapply(1:ncol(pop_pairs), function(i){
    pop1 = pop_pairs[1,i] 
    pop2 = pop_pairs[2,i]
    p1 = ps[[pop1]] 
    p2 = ps[[pop2]]
    n1 = obs_alleles[[pop1]] 
    n2 = obs_alleles[[pop2]]
    mapply(fst_site, p1, p2, n1, n2)
  })
  names(fst) = apply(pop_pairs, 2, paste, collapse="-")
  
    #And dxy
  dxy = lapply(1:ncol(pop_pairs), function(i){
    pop1 = pop_pairs[1,i]
    pop2 = pop_pairs[2,i]
    p1 = ps[[pop1]]
    p2 = ps[[pop2]]
    mapply(dxy_site, p1, p2)
  })

  
  #Adding names to dxy as well
  names(dxy) = apply(pop_pairs, 2, paste, collapse="-")
  
  summary_list = list(
    chr=chr, pos=pos,
    sample_het=apply(GT,1,sample_het),
    obs=obs_alleles, p=ps, het=hets,
    avg_r2=avg_r2, fst=fst, dxy=dxy
  )
  
  return(summary_list)
}

#-----------------------------
# Run the processing function
#-----------------------------
results <- process_GT(test_data, pop_file, ploidy=1)

#-----------------------------
# Check mean values per population
#-----------------------------
for(pop in sort(unique(pop_file$population))){
  pop_chr <- as.character(pop)
  idx <- which(pop_file$population==pop)
  
  # Mean allele frequency
  expected_mean_af <- mean(sapply(1:nrow(test_data), function(i) freq_GT(test_data[i, idx])), na.rm=TRUE)
  stopifnot(all.equal(expected_mean_af, mean(results$p[[pop_chr]], na.rm=TRUE), tolerance=1e-8))
  
  # Mean heterozygosity
  expected_mean_het <- mean(sapply(1:nrow(test_data), function(i) sample_het(test_data[i, idx])), na.rm=TRUE)
  stopifnot(all.equal(expected_mean_het, mean(results$het[[pop_chr]], na.rm=TRUE), tolerance=1e-8))
  
  # Average r²
  mat <- test_data[, idx, drop=FALSE]
  if(ncol(mat)>=2){
    cor_mat <- suppressWarnings(cor(t(mat), use="pairwise.complete.obs")^2)
    expected_avg_r2 <- mean(cor_mat[upper.tri(cor_mat)], na.rm=TRUE)
  } else expected_avg_r2 <- NA
  stopifnot(all.equal(expected_avg_r2, results$avg_r2[pop_chr], tolerance=1e-8))
}

cat("All tests passed! FST and DXY calculated for all population pairs.\n")

#-----------------------------
# Optional: Print summary
#-----------------------------
cat("\nMean allele frequencies per population:\n")
print(sapply(results$p, mean))

cat("\nMean heterozygosity per population:\n")
print(sapply(results$het, mean))

cat("\nAverage r² per population:\n")
print(results$avg_r2)

cat("\nPairwise FST:\n")
print(sapply(results$fst, mean))

cat("\nPairwise DXY:\n")
print(sapply(results$dxy, mean))
