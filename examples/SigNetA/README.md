# SigNetA_StandaloneApp

INSTALLATION:
Download this application from github(e.g. download as ZIP)
Open terminal and go to the directory of the folder that was downloaded and open R console.
If packrat(a dependency manager package) is installed, turn it off by typing packrat::off()

Installing packages(if not already installed):

Install devtools first. 
command:
install.packages(devtools)
    
rcytoscapjs:
Commands:
library("devtools");
devtools::install_github("cytoscape/r-cytoscape.js");

For more info:https://github.com/cytoscape/r-cytoscape.js

shiny:
Command:
install.packages("shiny")
    
shinyjs:
Commands:
install.packages("shinyjs")
    
Bionet:
Command:
source("https://bioconductor.org/biocLite.R");
biocLite("BioNet")
    
For more info:https://www.bioconductor.org/packages/release/bioc/html/BioNet.html
    
CLEAN:
Command:
install_github("uc-bd2k/CLEAN")
    
CLEAN.Hs:
Command:
install_github("uc-bd2k/CLEAN.Hs")
    
RMySQL:
Command:
install.packages("RMySQL")

SigNetAcytoscapejs:
Command:
install_github("Saadman/SigNetAcytoscapejs")


    
While running the application there might be a prompt to install new packages if it doesn't already exist in the users       library.
    
RUNNING THE APPLICATION:

Command:
library(shiny)
runApp()


*Tested on Chrome browser
    
    
    
    
    
    

