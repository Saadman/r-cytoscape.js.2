
sortNetwork<-function(x){
#File<- read.csv(file=x,sep="\t")

sorted<-x[order(x$Pvals),] 

topgenes<-head(sorted,100)

}
