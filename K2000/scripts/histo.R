tab = head(read.table("histo.txt"),200)
plot(tab$V2,tab$V1,type="h",log="y",xlab="kmer coverage",ylab="count",main="")
