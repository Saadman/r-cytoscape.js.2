library(igraph)
load("/opt/raid10/genomics/mario/ucGithub/netlincs/scratch/data/weightedGraphStringPPI_10.rda")
ppiGW.copy <- delete.edges(ppiGW, which(E(ppiGW)$weight <=0.7))