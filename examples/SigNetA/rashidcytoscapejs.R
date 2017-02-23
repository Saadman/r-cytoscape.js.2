rashidcytoscapejs<-function (nodeEntries, edgeEntries,layout = "cola", width = 1500, 
          height = 1800, showPanzoom = TRUE, highlightConnectedNodes = TRUE, 
          boxSelectionEnabled = TRUE) 
{
  x = list()
  x$nodeEntries <- nodeEntries
  x$edgeEntries <- edgeEntries
  x$layout <- layout
  x$showPanzoom <- showPanzoom
  x$highlightConnectedNodes <- highlightConnectedNodes
  x$boxSelectionEnabled <- boxSelectionEnabled
  htmlwidgets::createWidget(name = "rcytoscapejs", x, width = width, 
                            height = height, package ="rcytoscapejs")
}
