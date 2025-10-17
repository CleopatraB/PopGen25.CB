#-----------------------------
# Basic functions
#-----------------------------

# Allele frequency
freq_GT <- function(site, ploidy=1){
  clean_site <- na.omit(site)
  p <- sum(clean_site) / (ploidy * length(clean_site))
  return(p)
}

# Sample heterozygosity
sample_het <- function(site, ploidy=2){
  n <- length(site) * ploidy
  p <- freq_GT(site, ploidy)
  h <- n / (n - 1) * 2 * p * (1 - p)
  return(h)
}

# FST per site
fst_site <- function(p1, p2, n1, n2){
  f1 <- n1 / (n1 + n2)
  f2 <- n2 / (n1 + n2)
  num <- f1*p1*(1-p1) + f2*p2*(1-p2)
  denom <- (f1*p1 + f2*p2) * (f1*(1-p1) + f2*(1-p2))
  fst <- 1 - num/denom
  return(fst)
}

# DXY per site
dxy_site <- function(p1, p2){
  return(p1*(1-p2) + p2*(1-p1))
}

#-----------------------------
# Main process function
#-----------------------------

process_GT <- function(GT, samples){
  GT <- as.matrix(GT)
  
  # Extract chromosome and position from rownames
  ids <- gsub(":[A,T,C,G]*:[A,T,C,G]*", "", rownames(GT))
  chr <- gsub(":.*", "", ids)
  pos <- as.numeric(gsub(".*:", "", ids))
  
  summary_list <- list("chr"=chr, "pos"=pos)
  summary_list$sample_het <- apply(GT, 1, sample_het)
  
  pops <- sort(unique(samples$population))
  
  # Indices of samples for each population
  idx_pops <- lapply(pops, function(pop_id){
    which(colnames(GT) %in% samples$samples[which(samples$population == pop_id)])
  })
  names(idx_pops) <- pops
  
  # Observed alleles per site per population
  obs_alleles <- lapply(idx_pops, function(idx){
    apply(GT[, idx, drop = FALSE], 1, function(site) sum(!is.na(site)))
  })
  summary_list$obs <- obs_alleles
  
  # Allele frequencies per population
  ps <- lapply(idx_pops, function(idx){
    apply(GT[, idx, drop = FALSE], 1, freq_GT)
  })
  names(ps) <- pops
  summary_list$p <- ps
  
  # Heterozygosity per population
  hets <- lapply(idx_pops, function(idx){
    apply(GT[, idx, drop=FALSE], 1, sample_het)
  })
  names(hets) <- pops
  summary_list$het <- hets
  
  # Population pairs
  pop_pairs <- combn(pops, 2)
  if(is.vector(pop_pairs)) pop_pairs <- matrix(pop_pairs, nrow=2)  # Ensure 2-row matrix
  
  # FST per pair
  fst <- lapply(1:ncol(pop_pairs), function(pair_idx){
    pop1 <- pop_pairs[1, pair_idx]
    pop2 <- pop_pairs[2, pair_idx]
    n1 <- obs_alleles[[pop1]]
    n2 <- obs_alleles[[pop2]]
    p1 <- ps[[pop1]]
    p2 <- ps[[pop2]]
    mapply(fst_site, p1, p2, n1, n2)
  })
  names(fst) <- apply(pop_pairs, 2, paste, collapse="-")
  summary_list$fst <- fst
  
  # DXY per pair
  dxy <- lapply(1:ncol(pop_pairs), function(pair_idx){
    p1 <- ps[[pop_pairs[1, pair_idx]]]
    p2 <- ps[[pop_pairs[2, pair_idx]]]
    mapply(dxy_site, p1, p2)
  })

  names(dxy) <- apply(pop_pairs, 2, paste, collapse="-")
  summary_list$dxy <- dxy
  
  # Average r² (LD) per population
  avg_r2 <- sapply(idx_pops, function(idx){
    mat <- GT[, idx, drop=FALSE]
    if(ncol(mat) < 2){
    cor_mat <- cor(t(mat))^2
    mean(cor_mat[upper.tri(cor_mat)], na.rm=TRUE)
  })
  names(avg_r2) <- pops
  summary_list$avg_r2 <- avg_r2
  
  return(summary_list)
}
