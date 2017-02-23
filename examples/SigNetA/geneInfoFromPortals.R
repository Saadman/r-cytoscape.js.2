geneInfoFromPortals<-function(geneList,symbol=T,names=F){
library(RMySQL)
mycon <- dbConnect(MySQL(), user='public', dbname="GeneDB",host="10.165.4.231", port = 4040, password='public')
s<-paste("'",paste(geneList,collapse="','"),"'",sep="")
if(names) Fsql<-paste("select GeneID,Symbol,description from GeneInfo where GeneID in (",s,")",sep="")
else Fsql<-paste("select GeneID,Symbol from GeneInfo where GeneID in (",s,")",sep="")
rsSampleData<-dbSendQuery(mycon, Fsql)
geneData<-fetch(rsSampleData,n=-1)
dbDisconnect(mycon)
geneData
}

geneEntrez2Symbols<-function(geneList){
library(RMySQL)
mycon <- dbConnect(MySQL(), user='public', dbname="GeneDB",host="10.165.4.231", port = 4040, password='public')
s<-paste("'",paste(geneList,collapse="','"),"'",sep="")
Fsql<-paste("select Symbol from GeneInfo where GeneID in (",s,")",sep="")
rsSampleData<-dbSendQuery(mycon, Fsql)
geneData<-fetch(rsSampleData,n=-1)
dbDisconnect(mycon)
symbols<-unlist(geneData)
names(symbols)<-NULL
return(symbols)
}

getGeneLists<-function(geneList){
library(RMySQL)
mycon <- dbConnect(MySQL(), user='public', dbname="GeneDB", host="10.165.4.231", port = 4040, password='public')
s<-paste("'",paste(geneList,collapse="','"),"'",sep="")
Fsql<-paste("select GeneID,GeneList,GeneListTable from GeneToListMap where GeneID in (",s,")",sep="")
rsSampleData<-dbSendQuery(mycon, Fsql)
geneListTable<-fetch(rsSampleData,n=-1)
dbDisconnect(mycon)
geneListTable
}
