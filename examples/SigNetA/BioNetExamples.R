library(igraph)
library(BioNet)
library(DLBCL)
data(dataLym)
data(interactome)
pvals <- cbind(t = dataLym$t.pval, s = dataLym$s.pval)
rownames(pvals) <- dataLym$label
pval <- aggrPvals(pvals, order = 2, plot = FALSE)

#subnet <- subNetwork(dataLym$label, interactome)
subnet <- rmSelfLoops(subnet)
subnet
fb <- fitBumModel(pval, plot = FALSE)
scores <- scoreNodes(subnet, fb, fdr = 0.001)
score<-scores
head(scores)
length(scores)
pdf("/opt/raid10/genomics/rashid/cytoproject/NetworkAnalysis/scripts/new_his.pdf")
hist(scores)
dev.off()
summary(scores)
pdf("/opt/raid10/genomics/rashid/cytoproject/NetworkAnalysis/scripts/new_plot.pdf")
plot(log10(pval[match(names(scores),names(pval))]),scores)
dev.off()
module <- runFastHeinz(subnet, scores)
typeof(module)
logFC <- dataLym$diff
names(logFC) <- dataLym$label
pdf("/opt/raid10/genomics/rashid/cytoproject/NetworkAnalysis/scripts/new_plot3.pdf")
plotModule(module, scores = scores, diff.expr = logFC)#plots the static gene network
dev.off();

source("/opt/raid10/genomics/rashid/cytoproject/NetworkAnalysis/scripts/rashidplotmodule.R")
pdf("/opt/raid10/genomics/rashid/cytoproject/NetworkAnalysis/scripts/new_plotRashid.pdf")
colorNet<-rashidplotmodule(module, scores = scores, diff.expr = logFC)

dev.off();




library(rcytoscapejs)


id <- nodes(module)
name <- id
nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
nodeData$color<- rep("#00FF0F",nrow(nodeData))  #changed color of nodes
nodeData$shape <- "none"  #default shape

nodeData$href <- paste0("http://www.google.com/search?q=Seinfeld%20", nodeData$name)
nodeData$Diff_Exp="none"
nodeData$score="none"
for(i in 1:length(name)){
  nodeData[i,3]<-colorNet$c[i];
  nodeData[i,6]<-colorNet$d[i]
  nodeData[i,7]<-colorNet$sc[i]
  
}
for(i in 1:length(name)){
  if(colorNet$s[i]=="csquare")
    colorNet$s[i]<-"rectangle"
  else
    colorNet$s[i]<-"ellipse"
    
 nodeData[i,4]<-colorNet$s[i];
  
}
ltn<-unlist(lapply(edgeL(module),function(x) length(x[[1]])))

source<-unlist(lapply(1:37,function(x) rep(id[x],ltn[x])))
target<-unlist(lapply(edgeL(module), function(x) id[unlist(x)]))
vect<-c()
for(i in 1:length(target))  #extracting the value from the key value pair
  vect[i]<-target[[i]]




edgeData <- data.frame(source, target, stringsAsFactors=FALSE)
edgeData$color<- rep("red",nrow(edgeData)) 

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

source("/opt/raid10/genomics/rashid/cytoproject/NetworkAnalysis/scripts/rashidcytoscapejs.R")
rashidcytoscapejs(network$nodes, network$edges,showPanzoom=TRUE)
#IMPORTANT INFO:
#esultant  module  visualized  in  Cytoscape.    Signicantly  upregu-
#lated genes are coloured in red, genes that show signicant downregulation are
#coloured in green, for the contrast B vs.  T cells.  The score of the nodes is shown
#by the shape of the nodes, circles indicate a positive score, diamonds a negative
#score




