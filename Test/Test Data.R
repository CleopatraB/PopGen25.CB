#Generating test data
row_invariant = rep(0,10)
row_invariant2 = rep(1,10)
one_tenth = c(1,rep(0,9))
one_half_fst1 = c(rep(0,5),rep(1,5))
one_half_fst0 = c(0,1,1,1,0,0,1,1,1,0)
fake_data = rbind(row_invariant,row_invariant2,one_tenth,one_half_fst1,one_half_fst0)
colnames(fake_data) = c("A","B","C","D","E","a","b","c","d","e")
#We also want the rows to have names that work with our code:
rownames(fake_data) = paste("chr1:",1:5*100,":A:T",sep="")
write.table(fake_data,file="C:/Users/cleob/Gitub/PopGen25.CB/Test/test_data.txt",quote=FALSE)
pop_file = data.frame(samples=colnames(fake_data),population=c(rep(1,5),rep(2,5)))
write.table(pop_file,file="C:/Users/cleob/Gitub/PopGen25.CB/Test/pop_file.txt",row.names=FALSE,quote=FALSE)
