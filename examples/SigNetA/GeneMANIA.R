#   GeneMANIA label propagation algorithm
#   	Solve y = (I - L)f where L = D - W is the graph Laplacian matrix
#
#   Input: 
#       the bias vector y, in the paper, y=1 for seed genes, y=-1 for negative example genes, y=k for all unlabeled genes 
#       where k= (n+ - n-)/n where n+ and n- are the numbers of positive and negative genes. 
#   	
#   Citation:
#   	Sara Mostafavi, Debajyoti Ray, David Warde-Farley, Chris Grouios, Quaid Morris
#   	GeneMANIA: a real-time multiple association network integration algorithm for predicting gene function
#   	Genome Biol. 2008; 9(Suppl 1): S4. 
#
###load network (igraph) from genemania, calculate ginv of (I-L) matrix
load(file="/opt/raid10/genomics/data/networks/functionalNetworks/genemania/genemaniaG.rda")
library(igraph)
library(MASS)
require(org.Hs.eg.db)

L = graph.laplacian(as.undirected(genemaniaG), normalized = F)
I= diag(nrow(L))
#I_minus_L_inv=ginv(I-data.matrix(L))     #It took 7 hours to run on gimm4.
#save(I_minus_L_inv, file="/opt/raid10/genomics/huan/network/GeneMania/genemania_inv_matrix.rda")

###Seed genes (using example from GeneMania website, find geneids)
load(file="/opt/raid10/genomics/huan/network/GeneMania/genemania_inv_matrix.rda")
seeds.symbol=c("MRE11A","RAD51", "MLH1", "MSH2", "DMC1", "RAD51AP1", "RAD50", "MSH6", "XRCC3", "PCNA", "XRCC2")

seeds.gid = sapply(mget(seeds.symbol,org.Hs.egSYMBOL2EG,ifnotfound=NA),function(x)x[1])
seeds.gid=seeds.gid[1:5]

###As an example to run the algorithm, we set y=1 to label seed genes, y=0 for all unlabeled genes in the network. 
# y is the input
y=rep(0, nrow(L))
names(y)=colnames(L)
y[which(names(y) %in% seeds.gid)]=1

# apply genemania algorithm on y
f = I_minus_L_inv %*% y 
names(f)=colnames(L)   

#f[order(f, decreasing=T)][1:20]