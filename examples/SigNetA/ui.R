library(rcytoscapejs)
library(shiny)
library(shinyjs)
library(DT)
shinyUI(fluidPage(
  useShinyjs(),
  tags$style(type="text/css","
             img{margin: -2% 53% 3% ;} #change it to 65 for television
             "),
  title="SigNetA",
  
  titlePanel( img(src='finalResize.jpg', align = "right")),
  
  sidebarLayout(
    
    
    sidebarPanel(
      h4("File Upload"),
      fileInput('file1', 'Choose file to upload',
                accept = c(
                  'text/csv',
                  'text/comma-separated-values',
                  'text/tab-separated-values',
                  'text/plain',
                  '.csv',
                  '.tsv',
                  '.xls'
                )),
      
      h4("Download Network"),
      downloadButton('downloadNetworkImage', 'Download Network Data'),
      actionLink("saveImage", "Download as PNG"),
      h4("Download Enrichment Results"),
      downloadButton('downloadEnrich', 'Download Enrichment Data'),
      #actionButton("example", "Run Example"),
     
      selectInput("algorithm",
                  label="select an algorithm",choices=list("topHundredNetwork"=1,"Bionet"=2),                  
                  selected=1),
      selectInput("layout",
                  label="select a layout",choices=list("Fruchterman Reingold"=1,"D3layout"=2),                  
                  selected = 1),
      tags$hr(),
      width = 3,
      tags$style(type="text/css", "
                 #loadmessage {
                 position: fixed;
                 top: 400px;
                 left: 400px;
                 width: 50%;
                 padding: 8px 0px 8px 0px;
                 text-align: center;
                 font-weight: bold;
                 font-size: 100%;
                 color: #FFFFFF;
                 background-color: #337ab7;
                 border-radius:25px;
                 
                 z-index: 105;
                 }
                 
                 "), #tags style for load bar
      
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"
      )
      ),
    mainPanel(
      
      tabsetPanel(id = 'tabs',
                  tabPanel( "About",
                            includeHTML("index.html")
                            
                  ),
                  tabPanel("Network",
                           value="cytonet",
                           
                           h4("Network"),
                           
                           conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                            tags$div("Loading...",id="loadmessage")
                           ),
                           
                           textOutput(outputId="input_error"),
                           uiOutput("general_ui")
                        # rcytoscapejsOutput("plot", height="800px",width="900px")
                          
                           
                           
                           
                           
                  ),
                  tabPanel("Network Data",
                           value="datTable",
                           
                           conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                            tags$div("Loading...",id="loadmessage")
                           ),
                           textOutput(outputId="input_error1"),
                           
                           DT::dataTableOutput('contents')
                           
                  ),
                  
                  tabPanel( "Gene Ontology(GO) Enrichment Analysis",
                            value="irich",
                            
                            
                            conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                             tags$div("Loading...",id="loadmessage")
                            ),
                            textOutput(outputId="input_error2"),
                            
                            DT::dataTableOutput('enrichmentAnalysis')
                            
                            
                  )
      ), 
      
      tags$head(tags$script(src="cyjs.js"))
      
      
    )
    )
  )) 

