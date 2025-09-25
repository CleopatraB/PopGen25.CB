#Generating test data
row_invariant = rep(0,10)
row_invariant2 = rep(1,10)
one_tenth = c(1,rep(0,9))
one_half_fst1 = c(rep(0,5),rep(1,5))
one_half_fst0 = c(0,1,1,1,0,0,1,1,1,0)
fake_data = rbind(row_invariant,row_invariant2,one_tenth,one_half_fst1,one_half_fst0)
names(fake_data) = c("A","B","C","D","E","a","b","c","d","e")
#We also want the rows to have names that work with our code:
row.names(fake_data) = paste("chr1:",1:5*100,":A:T",sep="")
write.table(fake_data,file="<your test_dir/test_data.txt>",quote=FALSE)
pop_file = data_frame(samples=names(fake_data),population=c(rep(1,5),rep(2,5)))
write.table(pop_file,row.names=FALSE,quote=FALSE)


