tab = head(read.table("histo.txt"),20000)
plot(tab$V2,tab$V1,type="h",log="xy",xlab="unitig coverage",ylab="count",main="freq. unitigs.")
