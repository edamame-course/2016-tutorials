data=read.table("data.txt", header=TRUE, sep="\t",row.names=1)
library(vegan)

#Question2

#no.otus
nrow(data)
#58

#no. individuals in community
sum(rowSums(data))
#838

#no. singletons
data.pa=1*(data>0)
sum(rowSums(data.pa)==1)
#11

#Question3A- sampling curve

#sampling curve
data.sa=specaccum(data, method="collector")
plot(data.sa, xlab="No. samples observed", ylab="Cumulative No. Unique OTUs observed", col="red", lwd=2, main="Sampling curve")

#Question3b - rank abundance
o=sort(rowSums(data)/838, decreasing=TRUE)
plot(o, main="Rank abundance distribution", xlab="OTUs ranked by abundance", ylab="Relative abundance of OTU", col="red")

