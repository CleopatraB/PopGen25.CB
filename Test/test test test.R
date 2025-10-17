#-----------------------------
# Allele frequency
#-----------------------------
freq_GT <- function(site, ploidy = 1) {
  clean_site <- na.omit(site)
  if (length(clean_site) == 0) return(NA_real_)
  p <- sum(clean_site) / (ploidy * length(clean_site))
  return(p)
}

#-----------------------------
# Sample heterozygosity
#-----------------------------
sample_het <- function(site, ploidy = 1) {
  clean_site <- na.omit(site)
  n <- length(clean_site) * ploidy
  if (n <= 1) return(NA_real_)
  p <- freq_GT(clean_site, ploidy)
  h <- n / (n - 1) * 2 * p * (1 - p)
  return(h)
}

#-----------------------------
# Fst per site (robust)
#-----------------------------
fst_site <- function(p1, p2, n1, n2) {
  f1 <- n1 / (n1 + n2)
  f2 <- n2 / (n1 + n2)
  num <- f1 * p1 * (1 - p1) + f2 * p2 * (1 - p2)
  denom <- (f1 * p1 + f2 * p2) * (f1 * (1 - p1) + f2 * (1 - p2))
  fst <- ifelse(denom == 0, NA, 1 - num / denom)
  return(fst)
}

#-----------------------------
# Dxy per site
#-----------------------------
dxy_site <- function(p1, p2) {
  return(p1 * (1 - p2) + p2 * (1 - p1))
}

#-----------------------------
# Main processing function
#-----------------------------
process_GT <- function(GT, samples, ploidy = 1) {
  
  # Extract chromosome and position info
  ids <- gsub(":[A,T,C,G]*:[A,T,C,G]*", "", rownames(GT))
  chr <- gsub(":.*", "", ids)
  pos <- suppressWarnings(as.numeric(gsub(".*:", "", ids)))
  
  summary_list <- list(
    chr = chr,
    pos = pos,
    sample_het = apply(GT, 1, function(x) sample_het(x, ploidy = ploidy))
  )
  
  pops <- sort(unique(samples$population))
  
  # Get indices for each population
  idx_pops <- lapply(pops, function(x) {
    which(colnames(GT) %in% samples$sample_id[samples$population == x])
  })
  names(idx_pops) <- pops
  
  # Observed alleles per site per population
  obs_alleles <- lapply(idx_pops, function(idx) {
    apply(GT[, idx, drop = FALSE], 1, function(x) sum(!is.na(x)))
  })
  summary_list$obs <- obs_alleles
  
  # Allele frequencies per population
  ps <- lapply(idx_pops, function(idx) {
    apply(GT[, idx, drop = FALSE], 1, function(x) freq_GT(x, ploidy = ploidy))
  })
  summary_list$p <- ps
  
  # Heterozygosity per population
  hets <- lapply(idx_pops, function(idx) {
    apply(GT[, idx, drop = FALSE], 1, function(x) sample_het(x, ploidy = ploidy))
  })
  summary_list$het <- hets
  
  # Pairwise FST
  pop_pairs <- combn(pops, 2)
  fst <- lapply(1:ncol(pop_pairs), function(x) {
    pop1 <- pop_pairs[1, x]
    pop2 <- pop_pairs[2, x]
    n1 <- obs_alleles[[pop1]]
    n2 <- obs_alleles[[pop2]]
    p1 <- ps[[pop1]]
    p2 <- ps[[pop2]]
    fst_site(p1, p2, n1, n2)
  })
  names(fst) <- apply(pop_pairs, 2, paste, collapse = "-")
  summary_list$fst <- fst
  
  # Dxy
  dxy <- lapply(1:ncol(pop_pairs), function(x) {
    p1 <- ps[[pop_pairs[1, x]]]
    p2 <- ps[[pop_pairs[2, x]]]
    dxy_site(p1, p2)
  })
  names(dxy) <- apply(pop_pairs, 2, paste, collapse = "-")
  summary_list$dxy <- dxy
  
  # Average r² per population
  avg_r2 <- sapply(idx_pops, function(idx) {
    mat <- GT[, idx, drop = FALSE]
    if (ncol(mat) < 2) return(NA_real_)
    cor_mat <- suppressWarnings(cor(t(mat), use = "pairwise.complete.obs")^2)
    mean(cor_mat[upper.tri(cor_mat)], na.rm = TRUE)
  })
  names(avg_r2) <- pops
  summary_list$avg_r2 <- avg_r2
  
  return(summary_list)
}

#-----------------------------
# Run the processing function
#-----------------------------
results <- process_GT(GT = test_data, samples = pop_file, ploidy = 1)

#-----------------------------
# Validation test
#-----------------------------
for (pop in sort(unique(pop_file$population))) {
  idx <- which(pop_file$population == pop)
  
  # Mean allele frequency
  expected_mean_af <- mean(sapply(1:nrow(test_data), function(i) freq_GT(test_data[i, idx])), na.rm = TRUE)
  actual_mean_af <- mean(results$p[[as.character(pop)]], na.rm = TRUE)
  if (!is.na(expected_mean_af) && !is.na(actual_mean_af)) {
    stopifnot(all.equal(expected_mean_af, actual_mean_af, tolerance = 1e-8))
  }
  
  # Mean heterozygosity
  expected_mean_het <- mean(sapply(1:nrow(test_data), function(i) sample_het(test_data[i, idx])), na.rm = TRUE)
  actual_mean_het <- mean(results$het[[as.character(pop)]], na.rm = TRUE)
  if (!is.na(expected_mean_het) && !is.na(actual_mean_het)) {
    stopifnot(all.equal(expected_mean_het, actual_mean_het, tolerance = 1e-8))
  }
  
  # Average r²
  mat <- test_data[, idx, drop = FALSE]
  if (ncol(mat) >= 2) {
    cor_mat <- suppressWarnings(cor(t(mat), use = "pairwise.complete.obs")^2)
    expected_avg_r2 <- mean(cor_mat[upper.tri(cor_mat)], na.rm = TRUE)
  } else {
    expected_avg_r2 <- NA
  }
  actual_avg_r2 <- results$avg_r2[as.character(pop)]
  if (!is.na(expected_avg_r2) && !is.na(actual_avg_r2)) {
    stopifnot(all.equal(expected_avg_r2, actual_avg_r2, tolerance = 1e-8))
  }
}

cat("✅ All mean tests passed!/n/n")
cat("Mean allele frequencies per population:/n")
print(sapply(results$p, function(x) mean(x, na.rm = TRUE)))
cat("/nMean heterozygosity per population:/n")
print(sapply(results$het, function(x) mean(x, na.rm = TRUE)))
cat("/nAverage r² per population:/n")
print(results$avg_r2)