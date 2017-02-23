###GENERATE FILE UPLOAD NETWORK HERE#######

#logic<-logic[is.finite(logic$P_Value) & is.finite(logic$Diff_Exp), ]
#getting gene information from server##
source("/opt/raid10/genomics/software/RPrograms/source/geneInfoFromPortals.R")
geninfo<-geneInfoFromPortals(geneList=as.character(logic$geneID),symbol=T,names=F)
geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
#end of gene info  and gene label generation#
# for(i in 1:length(logic[2]))
pval<-as.numeric(logic$P_Value)

names(pval)<-geneLabels
logFC<-as.numeric(logic$Diff_Exp)
names(logFC)<-geneLabels
subnet <- subNetwork(geneLabels, interactome)
subnet <- rmSelfLoops(subnet)
pdf("/opt/raid10/genomics/rashid/cytoproject/NetworkAnalysis/scripts/fitbum2.pdf")
system.time( fb <- fitBumModel(pval, plot = TRUE))
dev.off()
system.time(scores <- scoreNodes(subnet, fb, fdr = 0.001))


system.time(module <- runFastHeinz(subnet, scores))
#a problem:some files give "-Inf" values.#

colorNet<-rashidplotmodule(module, scores = scores, diff.expr = logFC)
library(rcytoscapejs)

# id <- c("Jerry", "Elaine", "Kramer", "George")
id <- nodes(module)
name <- id
nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
nodeData$color<- rep("#00FF0F",nrow(nodeData))  #changed color of nodes
nodeData$shape <- "none"  #default shape
nodeData$href <- paste0("http://www.google.com/search?q=Seinfeld%20", nodeData$name)
for(i in 1:length(name)){
  nodeData[i,3]<-colorNet$c[i];
  
}
for(i in 1:length(name)){
  if(colorNet$s[i]=="csquare")
    colorNet$s[i]<-"rectangle"
  else
    colorNet$s[i]<-"ellipse"
  
  nodeData[i,4]<-colorNet$s[i];
  
}
ltn<-unlist(lapply(edgeL(module),function(x) length(x[[1]])))
source<-unlist(lapply(1:length(ltn),function(x) rep(id[x],ltn[x])))
target<-unlist(lapply(edgeL(module), function(x) id[unlist(x)]))
vect<-c()
for(i in 1:length(target))  #extracting the value from the key value pair
  vect[i]<-target[[i]]



# source <- c("Jerry", "Jerry", "Jerry", "Elaine", "Elaine", "Kramer", "Kramer", "Kramer", "George")
# target <- c("Elaine", "Kramer", "George", "Jerry", "Kramer", "Jerry", "Elaine", "George", "Jerry")
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
  
}
source("/opt/raid10/genomics/rashid/cytoproject/NetworkAnalysis/scripts/rashidcytoscapejs.R")
cyNetwork <- network
rashidcytoscapejs(network$nodes, network$edges,showPanzoom=TRUE)
###END FILE GENERATE UPLOAD NETWORK###