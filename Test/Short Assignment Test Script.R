#Load in your code
Source ("C:/Users/cleob/Gitub/PopGen25.CB/Short")

#Load in test data

test_data=read.table("C:/Users/cleob/Gitub/PopGen25.CB/Test/test_data.txt")
pops = read.table("C:/Users/cleob/Gitub/PopGen25.CB/Test/pop_file.txt")

#mean values


results = process_GT(test_data, pop_file)

#mean values per population
for(pop in sort(unique(pop_file$population))){
  pop_chr = as.character(pop)
  idx = which(pop_file$population==pop)
  
  #mean allele frequency
expected_mean_af = mean(sapply(1:nrow(test_data), function(i) freq_GT(test_data[i, idx])))
  stopifnot(all.equal(expected_mean_af, mean(results$p[[pop_chr]])))
  
  #mean heterozygosity
  expected_mean_het = mean(sapply(1:nrow(test_data), function(i) sample_het(test_data[i, idx])))
  stopifnot(all.equal(expected_mean_het, mean(results$het[[pop_chr]])))

  # avg. r²
  mat = test_data[, idx]
  if(ncol(mat)>=2){
    cor_mat = (cor(t(mat))^2)
    expected_avg_r2 = mean(cor_mat[upper.tri(cor_mat)])
  } else expected_avg_r2 = NA
  stopifnot(all.equal(expected_avg_r2, results$avg_r2[pop_chr]))
}


