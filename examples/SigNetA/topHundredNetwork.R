topHundredNetwork<-function(File=NULL){
library(igraph)
library(BioNet)
library(DLBCL)
data(interactome)
library(networkD3)
source("geneInfoFromPortals.R")
source("sortNetwork.R")
source("rashidplotmodule.R")
if(!is.null(File))
{
  File<-read.csv(file=File,sep='\t')
}
else{
File<-read.csv(file="files/sig_try3.tsv",sep='\t')
}
sortedFile<-sortNetwork(File)
logic<-sortedFile

geninfo<-geneInfoFromPortals(geneList=as.character(logic$GeneID),symbol=T,names=F)
geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))

#load("/opt/raid10/genomics/mario/ucGithub/netlincs/scratch/data/weightedGraphStringPPI_10.rda")
#ppiGW.copy <- delete.edges(ppiGW, which(E(ppiGW)$weight <=0.7))


subnet <- subNetwork(geneLabels, interactome,neighbors = "none") 
#subnet <- subNetwork(geneLabels, ppiGW,neighbors = "none") 

subnet <- rmSelfLoops(subnet)

logFC<-as.numeric(logic$coefficients)
names(logFC)<-geneLabels
module<-subnet
colorNet<-rashidplotmodule(module, diff.expr = logFC)

dev.off();




library(rcytoscapejs)


id <- nodes(module)
name <- id
nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
print(nodeData)

nodeData$color<- rep("#00FF0F",nrow(nodeData))  #changed color of nodes
nodeData$shape <- "ellipse"  #default shape
nodeData$href <- paste0("http://www.ncbi.nlm.nih.gov/gene/",gsub("[\\(\\)]", "", regmatches(nodeData$name, gregexpr("\\(.*?\\)", nodeData$name))))
nodeData$geneID<-gsub("[\\(\\)]", "", regmatches(nodeData$name, gregexpr("\\(.*?\\)", nodeData$name)))
nodeData$name<-sub(" *\\(.*", "", nodeData$name)
nodeData$Diff_Exp="none"
for(i in 1:length(name)){
  nodeData[i,3]<-colorNet$c[i];
  nodeData[i,7]<-colorNet$d[i]
}

statNet<<-nodeData

ltn<-unlist(lapply(edgeL(module),function(x) length(x[[1]])))


source<-unlist(lapply(1:length(ltn),function(x) rep(id[x],ltn[x])))
target<-unlist(lapply(edgeL(module), function(x) id[unlist(x)]))

vect<-c()
for(i in 1:length(target))  #extracting the value from the key value pair
  vect[i]<-target[[i]]




edgeData <- data.frame(source, target, stringsAsFactors=FALSE)


network <- createCytoscapeJsNetwork(nodeData, edgeData)

for(i in 1:length(target)){
  
  network$edges[[i]]$data$edgeTargetShape="none"  #making undirected graphss
  
}
for(i in 1:length(target)){
  for(j in i:length(target)){
    if(network$edges[[i]]$data$source == network$edges[[j]]$data$target)
      network$edges[[j]]$data$target= "none"
    
  }
  
} #to remove doubly linked edges 
#layout=layout_with_fr(colorNet$n) 
#rcytoscapejs(network$nodes, network$edges,layout="layout_with_fr", showPanzoom=FALSE)
source("rashidcytoscapejs.R")
rashidcytoscapejs(network$nodes, network$edges,showPanzoom=TRUE)
#IMPORTANT INFO:
#esultant  module  visualized  in  Cytoscape.    Signicantly  upregu-
#lated genes are coloured in red, genes that show signicant downregulation are
#coloured in green, for the contrast B vs.  T cells.  The score of the nodes is shown
#by the shape of the nodes, circles indicate a positive score, diamonds a negative
#score

}





