#original plotmodule function

rashidplotmodule<-function (network, layout = layout.fruchterman.reingold, labels = NULL, 
          diff.expr = NULL, scores = NULL, main = NULL, vertex.size = NULL, 
          ...) 
{
  require(igraph)
  if (is(network, "graphNEL")) {
    network <- igraph.from.graphNEL(network) #difference between graphs?
  }
  if (is.null(V(network)$name)) {
    V(network)$name <- as.character(V(network))
  }
  if (is.null(labels)) {
    if ("geneSymbol" %in% list.vertex.attributes(network)) {
      labels <- V(network)$geneSymbol
    }
    else {
      labels <- V(network)$name
    }
  }
  shapes <- rep("circle", length(V(network)))
  names(shapes) <- V(network)$name
  if (!is.null(scores) && !is.null(names(scores))) {
    shapes[intersect(names(which(scores < 0)), V(network)$name)] <- "circle"
  }
  if (is.null(scores) && "score" %in% list.vertex.attributes(network)) {
    scores <- V(network)$score
    names(scores) <- V(network)$name
    shapes[names(which(scores < 0))] <- "circle"
  }
  if (!is.null(diff.expr) && !is.null(names(diff.expr))) {
    #coloring <- .node.color(network, diff.expr)
    print("First Coloring")
    require(gplots)
    cols = greenred(23)
    bins = c(-10,seq(-.25,.25,0.025),10)
    V(network)$color = cols[cut(diff.expr,breaks=bins)]
    coloring<-V(network)$color
    print(coloring)
  }
  else {
    coloring <- "SkyBlue2"
  }
  if (is.null(diff.expr) && "diff.expr" %in% list.vertex.attributes(network)) {
    diff.exprs = V(network)$diff.expr
    names(diff.exprs) <- V(network)$name
    #coloring <- .node.color(network, diff.exprs)
    print("Second Coloring")
    ###mario edits
    require(gplots)
    cols = greenred(53)
    bins = c(-100,seq(-5,5,0.2),100)
    V(network)$color = cols[cut(V(network)$diff.expr,breaks=bins)]
    coloring<-V(network)$color
    ##############
  }
  max.labels <- max(nchar(labels))
  network.size = length(V(network))
  vertex.size2 <- 8
  cex = 0.6
  if (network.size < 50) {
    if (max.labels > 2) {
      labels.dist <- 0.5
    }
    else {
      vertex.size2 <- 15
      labels.dist <- 0
    }
  }
  if (network.size < 100 && network.size >= 50) {
    if (max.labels > 2) {
      labels.dist <- 0.5
    }
    else {
      labels.dist <- 0
    }
  }
  if (network.size >= 100) {
    if (max.labels > 3) {
      labels.dist <- 0.5
    }
    else {
      labels.dist <- 0
    }
  }
  if (!is.null(vertex.size)) {
    vertex.size2 <- vertex.size
    labels.dist <- vertex.size/15
  }
  plot(network, layout = layout, vertex.size = vertex.size2, 
       vertex.label = labels, vertex.label.cex = cex, vertex.label.dist = labels.dist, 
       vertex.color = coloring, vertex.label.family = "sans", 
       vertex.shape = shapes, main = main, ...)
  return(list(c=coloring,n=network,s=shapes,l=layout,d=diff.expr,sc=scores));
}