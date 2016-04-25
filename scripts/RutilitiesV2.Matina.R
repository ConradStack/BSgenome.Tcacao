
autoinc <- function(x,pat="V",delim="\\.",start.at=2,filext=NULL){
	retval = NULL
	tmp = strsplit(x,delim)[[1]]
	pm = grepl(sprintf("^%s[0-9]{1,}",pat),tmp)
	if(!any(pm)){
		extind = length(tmp)
		if(!is.null(filext))
			extind = which(tmp==filext)
		stopifnot(extind!=1)
		news = c(tmp[1:(extind-1)],sprintf("%s%d",pat,start.at),tmp[extind])
		retval = paste(news,collapse=gsub("\\\\","",delim))
	} else {
		change = tmp[which(pm)]
		stopifnot(length(change) == 1)
		oldint = as.integer(sub( sprintf("^%s([0-9]{1,})",pat), "\\1", change))
		#oldint = as.integer(gsub("[^0-9]","",change))
		change = sub( sprintf("^%s([0-9]{1,})",pat) , sprintf("%s%d",pat,oldint+1) , change)
		tmp[which(pm)] <- change
		retval = paste(tmp,collapse=gsub("\\\\","",delim))
	}
	return(retval)
}

fixname <- function(x) {
	if(length(x) > 1) {
		x = paste(x,collapse=" ")
	}
	gsub("'","-prime",x)
}

fixqual <- function(x, key="Express", delim=":", ...) {
	retval = ""
	tmpind = which(grepl(sprintf("^%s%s",key,delim),x, ...)); #stopifnot(length(tmpind)==1)
	if(length(tmpind) == 1) {
		retval = strsplit(x[tmpind],delim)[[1]][2]
	}
	retval
}

fixequiv <- function(x, na.str="na", agg=NULL) {
	retval = NA
	if( !(length(x)==1 && x == "na") ){
		retval = sub("^([_a-zA-Z0-9]+)\\/.+$","\\1",x)
		if(length(retval) > 1 && !is.null(agg) && is.character(agg))
			retval = paste(retval,collapse=agg)
	}
	retval
}



##########################################
# Todo:
# - change the gff file back from Parent=... to gene=... (mRNA lines only)
#
##########################################


setwd("~/sa/datasets/Matina_1.6_v1.1/")
library(GenomicRanges); 
library(rtracklayer)
library(GenomicFeatures)
library(RSQLite)
require(BSgenome.Tcacao.MarsInc.v11)
genome <- BSgenome.Tcacao.MarsInc.v11
sqlite_file = "cacao11genes_pub3i.good.sqlite"
gff_file = "cacao11genes_pub3i.good.ParentMod.tidy.gff" #  gt gff3validator cacao11genes_pub3i.good.ParentMod.tidy.gff
#gff_file = "cacao11genes_pub3i.good.tidy.gff" #  Parent= back to gene=
while(file.exists(sqlite_file)){
	sqlite_file = autoinc(sqlite_file)
}



# Create initial TxDb object+database:
sinfo = as.data.frame(seqinfo(genome))
junk = data.frame(chrom=rownames(sinfo),length=sinfo$seqlengths,is_circular=FALSE)
txdb = makeTranscriptDbFromGFF(gff_file,"gff3",chrominfo=junk,species="Theobroma cacao")
saveDb(txdb, file=sqlite_file)
isnew = TRUE
drv <- dbDriver("SQLite")
con = dbConnect(drv,sqlite_file)
newfieldnames = sprintf("tx_%s", c("modelname","orthoname","estgroup","qualclass","qualexpress","qualhomo","tair","uniref","cirad","paralog") )
sqlstrs = sprintf("ALTER TABLE transcript ADD %s TEXT",newfieldnames)
for(sql in sqlstrs){
	rs = dbSendQuery(con, sql)
	print(dbGetInfo(rs)$statement)
}



gff=NULL
#if(FALSE) { 
if(file.exists("gff.rda")) { 
	load("gff.rda")
} else {
	print("loading gff from file...")
	gff = import.gff3(gff_file)
	save(gff,file="gff_good_tidy.rda")
}
drv <- dbDriver("SQLite")
con = dbConnect(drv,sqlite_file)
minds = which(gff$type == "mRNA")
trans = dbReadTable(con,"transcript")
stopifnot(length(minds)==nrow(trans))


gall = gff[minds]
metas = mcols(gall)
all(metas$type=="mRNA")
txnames = metas$ID
transmap = sapply(trans$tx_name,function(x) which(txnames == x))
stopifnot(length(unique(transmap)) == length(transmap))



mnames = sapply(metas$Name, fixname )
otmp = sapply(metas$oname, fixname )
onames = character(length(mnames))
for(ii in seq(length(otmp))) {
	tmp = otmp[[ii]]
	if(length(tmp)!=0){
		onames[ii] <- ifelse(tolower(tmp)=="same",mnames[ii],tmp)
	} else {
		onames[ii] <-""
	}
}
ests = metas$estgroup
ests[is.na(ests)] <- ""
qualclass = sapply(metas$quality,fixqual,key="class",ignore.case=TRUE)
qualexpress = sapply(metas$quality,fixqual,key="express",ignore.case=TRUE)
qualhomo = sapply(metas$quality,fixqual,key="homology",ignore.case=TRUE)
if(debug){
	table(qualclass,useNA="always")
	table(qualexpress,useNA="always")
	table(qualhomo,useNA="always")
}
tair = sapply(metas$Dbxref,fixqual,key="tair",ignore.case=TRUE)
tair = sub("^([a-zA-Z0-9]+\\.?[0-9]*)\\/.+$","\\1",tair)
uniref = sapply(metas$Dbxref,fixqual,key="uniref[0-9]+?",delim="_",ignore.case=TRUE)
uniref = sub("^([a-zA-Z0-9]+\\.?[0-9]*)\\/.+$","\\1",uniref)
tctx = sapply( metas$equiv2 ,fixequiv, agg="," )
tctx[is.na(tctx)] <- ""
paras = sapply(metas$paralog, function(xx){
	if(length(xx)<2){
		return("")
	} else {
		xx[1]
	}
})
para.percent = sapply(metas$paralog, function(xx){
	if(length(xx)<2){
		return(NA)
	} else {
		return( as.numeric(gsub("[^0-9]","",xx[2])) )
	}
})



trans$tx_modelname <- mnames[transmap]
trans$tx_orthoname <- onames[transmap]
trans$tx_estgroup <- ests[transmap]
trans$tx_qualclass <- qualclass[transmap]
trans$tx_qualexpress <- qualexpress[transmap]
trans$tx_qualhomo <- qualhomo[transmap]
trans$tx_tair <- tair[transmap]
trans$tx_uniref <- uniref[transmap]
trans$tx_cirad <- tctx[transmap]
trans$tx_paralog <- paras[transmap]


#tmod = cbind(trans,DBID=trans$X_tx_id)


upfields = dbListFields(con,"transcript")
upsub = upfields[7:length(upfields)]
tmod = cbind(trans[,upsub],DBID=trans[,1])
ttest = head(tmod)

upcols = paste(sprintf("%s=?",upsub),collapse=", ")
sql = sprintf("UPDATE transcript SET %s WHERE _tx_id = ?",upcols)

dbest = dbSendPreparedQuery(con, sql, bind.data=tmod)

dbDisconnect(con)

# Check if DB can still be loaded:
rm(txdb)
txdb=loadDb(sqlite_file)



# # (TODO) Need to fix collapsed CDS/Exons (i.e., same gene model and different transcripts "share" CDS/Exons - duplicates are NOT retained by makeTranscriptDbfromGFF)
# # select group_concat(_tx_id,",") as txids, count(*) as nshare from splicing group by _exon_id"
# # all together:
# # select * from (select group_concat(_tx_id,",") as txids, count(*) as nshare from splicing group by _exon_id) where nshare > 1
# # test = cdsBy(txdb)
# # vect = logical(length(test)); count = 1
# # for(ii in seq(length(test)) ){
# # 	xx= test[[ii]]
# # 	if(!is.null(xx$exon_rank) && length(xx$exon_rank) !=0){
# # 		vect[count] = ifelse(all(diff(xx$exon_rank)==1),TRUE,FALSE)
# # 	} else {
# # 		vect[count] = NA
# # 	}
# # 	count = count + 1
# # 	if(count%%1000==0) print(count)
# # }
# # collapsed.exons = !vect
# # save(collapsed.exons,file="collapsed.exons.rda")
# load("collapsed.exons.rda")



# # Paralog graph
# # fixequiv( head(metas$equiv2) )

# # 6K chip data
# library(VariantAnnotation)  #
# sixk = readVcf("6kSNPs.vcf",genome=seqinfo(genome))
# target = rowData(sixk)
# #locs = locateVariants(target, txdb, AllVariants())
# locs <- locateVariants(target, txdb, CodingVariants())
# anno = as.data.frame(locs)

# ## See:
# ## VariantAnnotation:::.predictCodingGRangesList

# #for(ii in seq(nrow(sixk))){
# granno = NULL
# for(ii in 1:5){
# 	grs = tryCatch(predictCoding(sixk[ii],subject=txdb, seqSource=matina16), error= function(e) NA )
# 	if(typeof(grs)=="S4" || !is.na(grs)){
# 		print(ii)
# 		if(is.null(granno)){
# 			granno = grs
# 		} else {
# 			granno = append(granno,grs)
# 		}
# 	}
# }


