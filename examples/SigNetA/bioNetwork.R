##BIONET ALGORITHM###

bioNetwork<-function(File=NULL){
  library(igraph)
  library(BioNet)
  library(DLBCL)
  data(interactome)
  source("geneInfoFromPortals.R")
  source("sortNetwork.R")
  source("rashidplotmodule.R")
  if(!is.null(File))
  {
    logic<-read.csv(file=File,sep='\t')
  }
  else{
    logic<-read.csv(file="files/sig_try3.tsv",sep='\t')
  }
  
  
  source("geneInfoFromPortals.R")
  
  geninfo<-geneInfoFromPortals(geneList=as.character(logic$GeneID),symbol=T,names=F)
  geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
  pval<-as.numeric(logic$Pvals)
  
  names(pval)<-geneLabels
  
  logFC<-as.numeric(logic$coefficients)
  names(logFC)<-geneLabels
  
  
  subnet <- subNetwork(geneLabels, interactome)
  subnet <- rmSelfLoops(subnet)
  
  
  system.time( fb <- fitBumModel(pval, plot = FALSE))
  #err2<<-try(scoreNodes(subnet, fb, fdr = 0.1),silent=TRUE)
  #if(class(err2)=="try-error"){
  # output$input_error=renderText("No significant subnetwork generated.Please upload another Signature.")
  # }
  #else{
  #output$input_error=renderText("")
  system.time(scores <- scoreNodes(subnet, fb, fdr = 0.1))
  
  #err<<-try(runFastHeinz(subnet, scores),silent=TRUE)
  # if(class(err) == "try-error"){
  #   
  #   
  #   output$input_error=renderText("No significant subnetwork generated.Please upload another Signature.")
  #   stopifnot(class(err) == "try-error")
  #   
  # }
  
  
  # else{
  #output$input_error=renderText("")
  system.time(module <- runFastHeinz(subnet, scores))
  
  source("rashidplotmodule.R")
  pdf("wor.pdf")
  colorNet<-rashidplotmodule(module, scores = scores, diff.expr = logFC)
  dev.off()
  library(rcytoscapejs)
  
  
  id <- nodes(module)
  name <- id
  nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
  
  nodeData$color<- rep("#00FF0F",nrow(nodeData))  #changed color of nodes
  nodeData$shape <- "none"  #default shape
  nodeData$href <- paste0("http://www.ncbi.nlm.nih.gov/gene/",gsub("[\\(\\)]", "", regmatches(nodeData$name, gregexpr("\\(.*?\\)", nodeData$name))))
  nodeData$geneID<-gsub("[\\(\\)]", "", regmatches(nodeData$name, gregexpr("\\(.*?\\)", nodeData$name)))
  nodeNameEntrez<-nodeData$name
  nodeData$name<-sub(" *\\(.*", "", nodeData$name)
  
  nodeData$Diff_Exp="none"
  nodeData$score="none"
  for(i in 1:length(name)){
    nodeData[i,3]<-colorNet$c[i];
    nodeData[i,7]<-colorNet$d[i]
    nodeData[i,8]<-colorNet$sc[i]
    
  }
  for(i in 1:length(name)){
    if(colorNet$s[i]=="csquare")
      #colorNet$s[i]<-"rectangle"
      colorNet$s[i]<-"ellipse"
    
    else
      colorNet$s[i]<-"ellipse"
    
    nodeData[i,4]<-colorNet$s[i];
    
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
    
  }
  source("rashidcytoscapejs.R")
  #plotInput(network,network$nodes,network$edges)
  rashidcytoscapejs(network$nodes, network$edges,showPanzoom=TRUE)
  #}
  
  ###END FILE GENERATE UPLOAD NETWORK###
  #}
  
  
} #end of bionet algorithm