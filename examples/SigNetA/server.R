library(rcytoscapejs)
library(DT)
library(shiny)
library(shinyjs)
library(DLBCL)
library(BioNet)
library(igraph)
library(CLEAN)
library(CLEAN.Hs)
library(networkD3)
data(interactome)
statNet<<-NULL

shinyServer(function(input,output,session){
  
  
 
  observe({
    
    if (is.null(input$file1) || input$file1 == "") {
      shinyjs::disable("downloadNetworkImage")
      shinyjs::hide("saveImage")
      shinyjs::disable("downloadEnrich")
      shinyjs::hide("algorithm")
      shinyjs::hide("layout")
    } else {
      
      shinyjs::enable("downloadNetworkImage")
      shinyjs::show("saveImage")
      shinyjs::enable("downloadEnrich")
      shinyjs::show("algorithm")
      shinyjs::show("layout")
      
      
      Sys.sleep(2)
      
      # Hide the loading message when the rest of the server function has executed
      hide(id = "loading-content", anim = TRUE, animType = "fade")
      
    }
    
    toggle(condition = input$file1 , selector = "#tabs li a[data-value=cytonet]")
    toggle(condition = input$file1, selector = "#tabs li a[data-value=datTable]")
    toggle(condition = input$file1, selector = "#tabs li a[data-value=irich]")
    
    
    
    
  })
  observeEvent(input$saveImage, {
    
    session$sendCustomMessage(type="saveImage", message="NULL")
    
  })
  
  output$contents <- DT::renderDataTable({
   
    if (is.null(input$file1) && input$file1 == "") {
      
      paste("No genes analyzed")
    } else {
      
      
    #  if(class(err) == "try-error"){ 
      #  
      #  output$input_error1=renderText("No significant data generated.Please upload a new signature.")
      #  stopifnot(class(err) == "try-error")
    #  }
     # else{
       # output$input_error1=renderText("")
        inFile <- input$file1
        
        if (is.null(inFile)  ) #not necessary to check for current version.
        {
          return(NULL)
        }
        tablechange()
       
      }
   # }
    
  })
  tablechange<-reactive({
    if(input$algorithm=="1" || input$algorithm=="2"  || input$layout=="1" ||input$layout=="2"){
    statNet<-statNet[-c(1,3,4)]
    #  colnames(statNet)<-c("geneName","NCBI Information","geneID","Diff_Exp","score")
    # statNet<-statNet[c("geneName","geneID","Diff_Exp","score","NCBI Information")]
    colnames(statNet)<-c("geneName","NCBI Information","geneID","Diff_Exp")
    statNet<-statNet[c("geneName","geneID","Diff_Exp","NCBI Information")]
    statNet
    }
  })
  
  output$clickedNode = renderPrint({
    input$clickedNode
  })
  
  output$connectedNodes = renderPrint({
    input$connectedNodes
  })
  
  
  output$enrichmentAnalysis<-DT::renderDataTable({
    if (is.null(input$file1) && input$file1 == "" ) {
      paste("No genes analyzed")
    } else {
      
      
     # if(class(err) == "try-error"){ 
        
       # output$input_error2=renderText("No significant enrichment analysis data generated.Please upload another signature.")
      #  stopifnot(class(err) == "try-error")
      #}
     # else{
       # output$input_error2=renderText("")
        
        enrichmentChange()
        
       #this is where enrichment functionality was
        
        
      }
      
      
  #  }
  }
  
  
  )
  
  enrichmentChange<-reactive({
    if((input$algorithm=="1" && input$layout=="1") || (input$algorithm=="1" && input$layout=="2"))
    {
      logic<-read.csv(input$file1$datapath,sep="\t")
      
      
      logic<-sortNetwork(logic)
      source("geneInfoFromPortals.R")
      
      geninfo<-geneInfoFromPortals(geneList=as.character(logic$GeneID),symbol=T,names=F)
      geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
      
      subnetGenes<-gsub("[\\(\\)]", "", regmatches( statNet$id, gregexpr("\\(.*?\\)",  statNet$id)))
      allGenes<-gsub("[\\(\\)]", "", regmatches(geneLabels, gregexpr("\\(.*?\\)", geneLabels)))
      
      typeof(allGenes)
      result<-geneListEnrichment(subnetGenes, allGenes, functionalCategories = "GO", species = "Hs", minGenesInCategory=10, maxGenesInCategory=1000, inBkg=TRUE, sigFDR = 0.1, verbose=TRUE)
      
      
      
      
      options(scipen = 999)
      ID<-result$categories
      Name<-result$Description
      PValue<-format(result$FisherPValue,digits=2,scientific = TRUE)
      FDR<-format(result$FisherFDR,digits=2,scientific = TRUE)
      NetworkGenes<-result$nGenesInCategory
      AllGenes<-result$nAllGenesInCategory
      logOR<-round(result$logOR,2)
      tableEnrichment<<-data.frame(ID,Name,PValue,FDR,NetworkGenes,AllGenes,logOR)
      
    }
    
    else if((input$algorithm=="2" && input$layout=="2")||(input$algorithm=="2" && input$layout=="1")){
      logic<-read.csv(input$file1$datapath,sep="\t")
      
      
      
      source("geneInfoFromPortals.R")
      
      geninfo<-geneInfoFromPortals(geneList=as.character(logic$GeneID),symbol=T,names=F)
      geneLabels<-apply(geninfo,1,function(x) paste(x[2],"(",as.integer(x[1]),")",sep=""))
      
      subnetGenes<-gsub("[\\(\\)]", "", regmatches( statNet$id, gregexpr("\\(.*?\\)",  statNet$id)))
      allGenes<-gsub("[\\(\\)]", "", regmatches(geneLabels, gregexpr("\\(.*?\\)", geneLabels)))
      
      typeof(allGenes)
      result<-geneListEnrichment(subnetGenes, allGenes, functionalCategories = "GO", species = "Hs", minGenesInCategory=10, maxGenesInCategory=1000, inBkg=TRUE, sigFDR = 0.1, verbose=TRUE)
      
      
      
      
      options(scipen = 999)
      ID<-result$categories
      Name<-result$Description
      PValue<-format(result$FisherPValue,digits=2,scientific = TRUE)
      FDR<-format(result$FisherFDR,digits=2,scientific = TRUE)
      NetworkGenes<-result$nGenesInCategory
      AllGenes<-result$nAllGenesInCategory
      logOR<-round(result$logOR,2)
      tableEnrichment<<-data.frame(ID,Name,PValue,FDR,NetworkGenes,AllGenes,logOR)
      
    }
    
    
  })
  
  output$genenodes <- renderText({ 
    
    if (is.null(input$file1) || input$file1 == "") {
      paste("No genes analyzed")
    } else {
      
      
      newlist<-c()
      
      
      for (i in 1:length(statNet$id)){
        append(statNet$id,"\n",i)
        newlist[i]<-statNet$id[i]
      }
      
      newlist
      
      
      
      
    }
    
    
    
    
  })
  
  
 
   
    
    
# #EXAMPLE NETWORK Generation#
#   observeEvent(input$example,{
#     
#     # toggle(condition = input$file1 , selector = "#tabs li a[data-value=cytonet]")
#     # toggle(condition = input$file1, selector = "#tabs li a[data-value=datTable]")
#     # toggle(condition = input$file1, selector = "#tabs li a[data-value=irich]")
#     
#   
#      output$plot1<- 
#      renderRcytoscapejs({
#      print("Example Network")
#        source("NetworkFunctions.R")
#        topHundredNetwork()
#      })
#    
#    
#     
#  
#       
#     })
  
  output$general_ui <- renderUI({
   if(input$algorithm=="1" && input$layout=="1"){
     rcytoscapejsOutput("plot", height="800px",width="900px")
   }
    else if(input$algorithm=="2" &&  input$layout=="1"){
      rcytoscapejsOutput("plot", height="800px",width="900px")
    }
    
    #LAYOUT
     else if(input$algorithm=="1" && input$layout=="2"){
      simpleNetworkOutput("simple")
     }
    else if(input$algorithm=="2" && input$layout=="2"){
      simpleNetworkOutput("simple")
    }
    
    
    
  })
  datasetInput <- reactive({
    if(input$algorithm=="1"){
      source("topHundredNetwork.R")
      
      ret<-topHundredNetwork(input$file1$datapath)
      
    }
    
    else if(input$algorithm=="2"){
      source("bioNetwork.R")
      
      ret<-bioNetwork(input$file1$datapath)
    }
    
    
    
  })
  
  layoutInput<-reactive({
    
    if(input$layout=="2"){
      if(input$algorithm=="1"){
        source("D3layout.R")
        ret<-D3layout(input$file1$datapath)
      }
      else if(input$algorithm=="2"){
      source("D3bioNetwork.R")
      ret<-D3bioNetwork(input$file1$datapath)
      }
    }
    
    
  }
    
  )
  
  output$plot <- renderRcytoscapejs({
  #  source("NetworkFunctions.R")
    
    
    
    Sys.sleep(2);
    #inFile2<-input$file1;
    #source("NetworkFunctions.R")
    #ret<-topHundredNetwork(inFile2$datapath)
    datasetInput()
    
    
    
    
    
    
    
    
    
    
    
  })
  
  output$simple<-renderSimpleNetwork({
    
    layoutInput()
    
  })

  
  
  plotInput <- function(network,nodes,edges){
    
    cyNetwork <- network
    
    rashidcytoscapejs(nodes,edges,showPanzoom=TRUE)
  }
  
  
  
  output$downloadNetworkImage <- downloadHandler(
    filename = "Shinynet.csv",
    
    content = function(file) {
      statNet<-statNet[-c(1,3,4)]
      colnames(statNet)<-c("geneName","NCBI Information","geneID","Diff_Exp","score")
      statNet<-statNet[c("geneName","geneID","Diff_Exp","score","NCBI Information")]
      write.csv(statNet, file)
      
      
    })  
  output$downloadEnrich <- downloadHandler(
    filename = "Shinynet.csv",
    
    content = function(file) {
      
      write.csv(tableEnrichment, file)
      
      
    })  
  
  
  
  
  
}
)
