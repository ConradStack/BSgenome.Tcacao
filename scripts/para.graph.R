
require(igraph)

setwd("~/sa/datasets/Matina_1.6_v1.1/")
load("paralogs.rda")

na.inds = which(!is.na(para.percent))
from = txnames[na.inds]
to = paras[na.inds]
weight = para.percent[na.inds]

gr = graph.edgelist(cbind(from,to),TRUE)
E(gr)$weight <- weight

#######
# Check
# ID=Thecc1EG033914t5;Parent=Thecc1EG033914; ... paralog=Thecc1EG000562t1,28%;
E(gr)[inc("Thecc1EG033914t5")] # -> Thecc1EG000562t1
E(gr)[inc("Thecc1EG033914t5")]$weight # = 28
#
#######

clusts = clusters(gr)

# largest subgraph (cluster):
cat(sprintf("Found %d clusters\n",clusts$no))
largest = which.max(clusts$csize)
linds = which(clusts$membership==largest)
#large = subgraph(gr,linds)

# High paralogy:
subgr = subgraph.edges(gr,E(gr)[ weight >= 40 ])



######  Ortholog graph:  Matina16 -> Criollo (CIRAD)
require(AnnotationDbi)
txdb=loadDb("cacao11genes_pub3i.good.V4.sqlite")
ks = keys(txdb,"TXNAME")
to.cirad = select(txdb,ks, c("TXCIRAD"), "TXNAME")
to.cirad = to.cirad[which(to.cirad$TXCIRAD!=""),]
addsplits = to.cirad[0,]
tosplit = which(grepl(",",to.cirad$TXCIRAD))
for(modrow in tosplit){
	print(modrow)
	newcirad=strsplit(to.cirad$TXCIRAD[modrow],",")[[1]]
	newname = rep(to.cirad$TXNAME[modrow],length(newcirad))
	addsplits <- rbind(addsplits,cbind(TXNAME=newname,TXCIRAD=newcirad))
}
synteny = rbind(to.cirad[-tosplit,],addsplits)




