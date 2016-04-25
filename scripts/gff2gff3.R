

dat = read.table("mazhao_CACAO_Transposable_Element_Annotation_New.gff",sep="\t",header=T, fill=TRUE)
N = nrow(dat)
gff3 = data.frame(
	seqid = dat$Scaffold,
	source = rep("mazhao",N),
	type = rep("transposable_element",N),
	start = dat$Start_Position,
	end = dat$End_Position,
	score =rep(0.05,N),
	strand =rep("?",N),
	phase =rep(".",N),
	attributes =character(N)
)


ids = with(dat,sprintf("TcTE_%d",Element_ID))
gnames=  with(dat,sprintf("%s_%s_%s",Order,Super_Family,Family))
gff3$attributes <- sprintf("ID=%s;Name=%s",ids,gnames)

write.table(gff3,file="mazhao_CACAO_Transposable_Element_Annotation_JCS.gff3",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
