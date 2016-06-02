
## TODO:  Need to comment the below code better


require(mmutils)

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
#
##########################################

### (defunct 5-26-2016) Previous QA steps (using genometools from command line):
# gt gff3 -sort yes -retainids yes -checkids yes -addids no -tidy yes -typecheck sofa.obo -xrfcheck GO.xrf_abbr -o cacao11genes_pub3i.mod.cleaned.gff -v yes cacao11genes_pub3i.mod.gff 2> gffcleaning.log
# sed 's/^scaffold_10r/scaffold_10/g' cacao11genes_pub3i.mod.cleaned.gff > cacao11genes_pub3i.mod.cleaned.renamed.gff
# perl -pi -e "s/cacao1mito/cacao1rdna_scaffold_175/g;" cacao11genes_pub3i.mod.cleaned.renamed.gff
# gt gff3validator -typecheck sofa.obo -xrfcheck GO.xrf_abbr cacao11genes_pub3i.mod.cleaned.gff  2> gffvalidation.log # validate gff file, make changes as needed (had to remove some poorly supported genemodels)

### QA / pre-processing steps:
# < import annotations in to geneious, then export them.>
# Execute the following commands:
# 
system("sed '/Polymorphism/d' ~/code/BSgenome.Tcacao/scratch/Theobroma_cacao.main_genome.scaffolds.gff > ~/code/BSgenome.Tcacao/scratch/Theobroma_cacao.main_genome.scaffolds.preproc.gff")

system("gt gff3 -force yes -sort yes -retainids yes -checkids yes -addids no -tidy yes -typecheck ~/sa/Matina_1.6_v1.1/sofa.obo -o ~/code/BSgenome.Tcacao/scratch/Theobroma_cacao.main_genome.scaffolds.preproc2.gff -v yes ~/code/BSgenome.Tcacao/scratch/Theobroma_cacao.main_genome.scaffolds.preproc.gff 2> gffpreproc.log")


setwd("~/sa/Matina_1.6_v1.1/")
library(GenomicRanges); 
library(rtracklayer)
library(GenomicFeatures)
library(RSQLite)
require(BSgenome.Tcacao.MarsInc.v11)
genome <- BSgenome.Tcacao.MarsInc.v11
sqlite_file = "cacao11genes_pub3i.sqlite"
gff_file = "~/code/BSgenome.Tcacao/scratch/Theobroma_cacao.main_genome.scaffolds.preproc2.gff"
gff_fixed = sub("preproc2\\.gff$","ParentsFixed.gff",gff_file)
gff_fixed_tmp = sub("preproc2\\.gff$","ParentsFixed.tmp.gff",gff_fixed)
gff_fixed_tidy = sub("gff$","tidy.gff",gff_fixed)
gff_orig = "cacao11genes_pub3i.gff"  														# "original" pub3i annotation, downloaded from cacaogenomedb.org in 2015 
gff_trans = "~/code/BSgenome.Tcacao/scratch/mazhao_CACAO_Transposable_Element_Annotation_JCS.gff3"							# mazhao's transposable element annotation, from v1.1 reference documents (pre-my-arrival)


# # Get output file:
# while(file.exists(sqlite_file)){
# 	sqlite_file = autoinc(sqlite_file)
# }

########
{
	print(sprintf("Reading from:  %s",basename(gff_file)))
	print(sprintf("Modified gff:  %s",basename(gff_fixed) ))
	print(sprintf("Outputing to:  %s",sqlite_file))
}
########

# Create initial TxDb object+database:
# sinfo = as.data.frame(seqinfo(genome))
# #junk = data.frame(chrom=rownames(sinfo),length=sinfo$seqlengths,is_circular=FALSE)
# txdb = makeTxDbFromGFF(gff_file,"gff3",chrominfo=seqinfo(matina16), organism="Theobroma cacao")
# saveDb(txdb, file=sqlite_file)


# gff=NULL
# #if(FALSE) { 
# if(file.exists("gff.rda")) { 
# 	load("gff.rda")
# } else {

print("loading gff from file...")
gff = import.gff3(gff_file)
save(gff,file="gff.rda")
# gff.orig <- gff
# if(any(gff$source=="Geneious"))
# gff = gff[-which(gff$source == "Geneious"),]
# 

#}

## Fix 'Parent' attribute for mRNA entries, and other small issues
#
minds = which(gff$type == "mRNA")  					# mRNA indices
ginds = which(gff$type == "gene")
gene.ids = sub("(t[0-9]{1,2})$","",gff$ID[minds])	# strip off transcript indication from transcript.id, making it a valid gene.id
gene.ids.len = nchar(gene.ids)
stopifnot( length(table(gene.ids.len))==1 )  		# all gene.ids should have the same length
bad = setdiff(gene.ids, gff$ID[gff$type=="gene"])	# all gene.ids should have a corresponding gene ID in the gff file
stopifnot(length(bad)==0)
stopifnot(gff$gene[minds] == gene.ids)
gff$Parent[minds] <- gene.ids 						# Set the Parent attribute of each mRNA so that they are linked to their corresponding gene.id
# mpars = as.character(gff$Parent[minds])
# mmap = match(gff$ID[ginds],mpars)
# stopifnot(gff$type[minds[mmap]] == "mRNA") 	# sanity check
# gnames = gff$Name[minds[mmap]]
# stopifnot(length(gnames)==length(ginds)) 	# another sanity check
# gff$Name[ginds] <- gnames
# gff$Name[minds] <- NA
mnames.list <- gff$Name
#mnames = sapply(mnames.list,function(x) ifelse(length(x)==1,x,paste(x,collapse=",")))
mnames = sapply(mnames.list, fixname)
mids = gff$ID #[minds]
gff$Name <- NULL
export.gff3(gff,gff_fixed)							# write "fixed" gff file


## Read gff file into data table:
#
gfftab = read.table(gff_fixed, sep="\t", header=FALSE, stringsAsFactors=FALSE)

mrna.inds = which(gfftab[,3] == "mRNA" | gfftab[,3] == "transposable_element")
tokens = tokenize(gfftab[mrna.inds,9],";","=")
mrna.ids = sapply( (gfftab[mrna.inds,9]), function(x) tokenize(grep("^ID",tokenize(x,";")[[1]], value=T),"=")[[1]][2]  )
all(mrna.ids %in% mids) # sanity check
mmap = match(mrna.ids,mids)
stopifnot(!any(is.na(mmap))) # sanity check
new.attr = sprintf("model=%s;", mnames[mmap] )
stopifnot(length(new.attr)==length(mrna.inds))
gfftab[mrna.inds,9] <- sprintf("%s%s",gfftab[mrna.inds,9], new.attr)

# fix strand levels 
gfftab[which(gfftab[,2] == "mazhao"),7] <- '.'
gfftab[which(gfftab[,7] == "*"),7] <- '?'

writeLines(c(
	"##gff-version 3",
	"##source-version RutilitiesV2.Matina.R",
	sprintf("##date %s",date()),
	"#program: overbestgenes, selection of best gene set by evidence scores",
	"#reference: http://arthropods.eugenes.org/EvidentialGene/",
	"#version: 2011.09.29",
	"#author: d. g. gilbert, gilbertd at indiana edu",
	"#scoretype: homolog:9,paralog:1,ovpro:1,ovrna:4,est:1,rseq:1,nintron:50,inqual:2,inerr:20,intr:0,terepeat:-2,UTR:5,CDS:1",
	"#dropscore: *homolog:80,*paralog:149,ovpro:50,*nintron:2,inqual:20,+CDS:201",
	"#sources: n=11: AUGepir1,AUGepir1a,AUGepir3,AUGie3,AUGpier6,AUGpier8,AUGpier8a,AUGpiern7,mar1g.mar11f,mar7g.",
	"#"
), gff_fixed_tmp)
write.table(gfftab, gff_fixed_tmp, append=TRUE, sep="\t",row.names=FALSE, col.names=FALSE, quote=FALSE)
#

## Use genometools to do some additional tidying:
#
system(sprintf("gt gff3 -force yes -sort yes -retainids yes -checkids yes -addintrons yes -addids no -tidy yes -typecheck ~/sa/Matina_1.6_v1.1/sofa.obo -o %s -v yes %s 2> gfftidy.log",gff_fixed_tidy,gff_fixed_tmp)) 
# tail gfftidy.log
# NB -> (5.26.2016) The gene model 'Thecc1EG008240' has exons that overlap (by 1 base), this one is NOT included in the database for this reason.

## Create TxDb object from the 'fixed' gff file, and save it.
#
if(file.exists(gff_fixed_tidy)){
	txdb = makeTxDbFromGFF(gff_fixed_tidy,"gff3",chrominfo=seqinfo(matina16), organism="Theobroma cacao")
	saveDb(txdb, file=sqlite_file)
}

isnew = TRUE
drv <- dbDriver("SQLite")
con = dbConnect(drv,sqlite_file)
newfieldnames = sprintf("tx_%s", c("modelname","orthoname","estgroup","qualclass","qualexpress","qualhomo","tair","uniref","cirad","paralog","paraperc") )
newfieldtypes = c(rep("TEXT",10),"INTEGER")
sqlstrs = sprintf("ALTER TABLE transcript ADD %s %s",newfieldnames,newfieldtypes)
for(sql in sqlstrs){
	rs = dbSendQuery(con, sql)
	print(dbGetInfo(rs)$statement)
}
dbDisconnect(con)


drv <- dbDriver("SQLite")
con = dbConnect(drv,sqlite_file)
trans = dbReadTable(con,"transcript")
{  	# DEBUG
	cat(sprintf("Number of mRNA records in gff file:	%d\n",length(minds)))
	cat(sprintf("Number of mRNA records in TxDb file:	%d\n",length(trans$X_tx_id)) )
	# (5-26-2016: The length of the minds object will / should actually be longer than the trans table because 1+ transcripts were dropped from trans when the gff file was being converted to a txdb object)
}

# The mRNA lines from the original gff file
metas = gff[minds,]
stopifnot(sapply(mnames[minds],length)!=0)
metas$ModelName <- unlist(mnames[minds])
stopifnot(metas$type=="mRNA")

transmap = match(trans$tx_name, metas$ID)
stopifnot(!is.na(transmap))
stopifnot(length(unique(transmap)) == length(transmap)) # want a unique match

otmp = sapply(metas$oname, fixname )
onames = character(length(otmp))
for(ii in seq(length(otmp))) {
	tmp = otmp[[ii]]
	if(length(tmp)!=0){
		onames[ii] <- ifelse(tolower(tmp)=="same", metas$ModelName[ii] ,tmp)
	} else {
		onames[ii] <- ""
	}
}

ests = metas$estgroup
ests[is.na(ests)] <- ""
qualclass = sapply(metas$quality,fixqual,key="class",ignore.case=TRUE)
qualexpress = sapply(metas$quality,fixqual,key="express",ignore.case=TRUE)
qualhomo = sapply(metas$quality,fixqual,key="homology",ignore.case=TRUE)
if(F){
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


# Add columns to 'trans'(cript) table:
trans$tx_modelname <- metas$ModelName[transmap]
trans$tx_orthoname <- onames[transmap]
trans$tx_estgroup <- ests[transmap]
trans$tx_qualclass <- qualclass[transmap]
trans$tx_qualexpress <- qualexpress[transmap]
trans$tx_qualhomo <- qualhomo[transmap]
trans$tx_tair <- tair[transmap]
trans$tx_uniref <- uniref[transmap]
trans$tx_cirad <- tctx[transmap]
trans$tx_paralog <- paras[transmap]
trans$tx_paraperc <- para.percent[transmap]

# Update sqlite database and close connection:
upfields = dbListFields(con,"transcript")
upsub = upfields[7:length(upfields)]
tmod = cbind(trans[,upsub],DBID=trans[,1])
ttest = head(tmod); ttest

upcols = paste(sprintf("%s=?",upsub),collapse=", ")
sql = sprintf("UPDATE transcript SET %s WHERE _tx_id = ?",upcols)

dbest = dbSendPreparedQuery(con, sql, bind.data=tmod)
dbDisconnect(con)



## Verify that data in sqlite database matches what is in the original annotation document
#
# Check if DB can still be loaded:

rm(txdb)
txdb = loadDb(sqlite_file)
orig = import.gff3("~/code/BSgenome.Tcacao/scratch/cacao11genes_pub3i.gff.gz")

# Compare some number of entries at random between original gff and txdb object
n.entries = 10000
had.issues = c("Thecc1EG008240t1")
to.compare = sample(intersect(orig$ID,gff$ID[minds]),n.entries , replace=FALSE) 
to.compare = setdiff(to.compare, had.issues)
orig.inds = match(to.compare,orig$ID)
orsel = orig[orig.inds,]

# columns(txdb) - columns that can be returned
txsel = select(txdb, keys=to.compare, keytype="TXNAME", columns=c("GENEID","TXNAME","TXMODELNAME","TXID","TXSTART","TXEND","TXESTGROUP","TXCIRAD"))

for(ii in 1:length(to.compare)) 
{
	if(ii %% 100 == 0) cat(sprintf("%d out of %d finished\n",ii,n.entries))

	stopifnot( orsel$ID[ii] == txsel$TXNAME[ii] )
	stopifnot( orsel$gene[ii] == txsel$GENEID[ii] )
	stopifnot( start(orsel[ii,]) == txsel$TXSTART[ii] )
	stopifnot( end(orsel[ii,]) == txsel$TXEND[ii] )
	if(!is.na(orsel$estgroup[ii]))
		stopifnot( orsel$estgroup[ii] == txsel$TXESTGROUP[ii] )
	if( !any(is.na(orsel$equiv2[[ii]])) && length(orsel$equiv2[[ii]])==1 && orsel$equiv2[[ii]] != "na" ) 
		stopifnot( grepl(txsel$TXCIRAD[ii],orsel$equiv2[ii]) )
}


## Create package skeleton for TxDb objects (only needs to be done once) 
#
if( FALSE ) {
	
	# For reference while package is being tested.  Should remove these objects from the package later to streamline it, but for now they stay.
	# gff.orig = import.gff3(gff_orig)
	# gff.trans = import.gff3(gff_trans)

	package.skeleton(
			name="TxDb.Tcacao.Geneious.pub3i",
			path="/Users/miamimac2/code/BSgenome.Tcacao",
			list=c(pub3i.gff="orig")
	)
	dir.create("/Users/miamimac2/code/BSgenome.Tcacao/TxDb.Tcacao.Geneious.pub3i/inst")
	dir.create("/Users/miamimac2/code/BSgenome.Tcacao/TxDb.Tcacao.Geneious.pub3i/inst/extdata")
	dir.create("/Users/miamimac2/code/BSgenome.Tcacao/TxDb.Tcacao.Geneious.pub3i/R")

	file.copy(sqlite_file,"/Users/miamimac2/code/BSgenome.Tcacao/TxDb.Tcacao.Geneious.pub3i/inst/extdata/TxDb.Tcacao.Geneious.pub3i.sqlite")
	#saveDb(txdb,"/Users/miamimac2/code/BSgenome.Tcacao/TxDb.Tcacao.Geneious.pub3i/inst/extdata/TxDb.Tcacao.Geneious.pub3i.sqlite")

	
}








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


