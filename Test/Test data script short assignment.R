#makes the rows
row_invariant = matrix(rep(0, 30), nrow = 1)
row_invariant2 = matrix(rep(1, 30), nrow = 1)
one_tenth = matrix(c(1, rep(0,9), 1, rep(0,9), 1, rep(0,9)), nrow = 1)
one_half_fst1 = matrix(c(rep(0,5), rep(1,5), rep(0,5), rep(1,5), rep(0,5), rep(1,5)), nrow = 1)
one_half_fst0 = matrix(c(0,1,1,1,0, 0,1,1,1,0, 0,1,1,1,0, 0,1,1,1,0, 0,1,1,1,0, 0,1,1,1,0)[1:30], nrow = 1)

# makes random rows of data with 30 columns for the 3 populations, 10 each. From internet not class
random_rows = matrix(sample(c(0,1), 35 * 30, replace = TRUE), nrow = 35, ncol = 30)

#puts the rows together
test_data = rbind(row_invariant, row_invariant2, one_tenth, one_half_fst1, one_half_fst0, random_rows)
#makes column names for each pop A, B, C
colnames(test_data) = c(paste0("A",1:10), paste0("B",1:10), paste0("C",1:10))
#makes row names 
rownames(test_data) = paste("chr1:", (1:nrow(test_data))*100, ":A:T", sep="")
#output
write.table(test_data, file = "C:/Users/cleob/Gitub/PopGen25.CB/Test/test_data.txt", quote = FALSE)
pop_file = data.frame( samples = colnames(test_data), population = c(rep(1,10), rep(2,10), rep(3,10)))
write.table(pop_file, file = "C:/Users/cleob/Gitub/PopGen25.CB/Test/pop_file.txt", row.names = FALSE, quote = FALSE)
