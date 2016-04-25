require(ape)

reps = read.table("~/sa/datasets/Matina_1.6_v1.1/mazhao_CACAO_Transposable_Element_Annotation_New.gff",header=T,sep="\t",fill=T,stringsAsFactors=F)

i=1 
fastout <- function(x){
	header=sprintf(">%s_%d_%d %s_%s_%s_%s_%s",x["Scaffold"],as.integer(x["Start_Position"]),as.integer(x["End_Position"]),x["Class"],x["Sub_Class"],x["Order"],x["Super_Family"],x["Family"])
	start=as.integer(x["Start_Position"])
	end =as.integer(x["End_Position"])
	scaf = x["Scaffold"]
	c(header,as.character(subseq(scaflist$reference[[scaf]],start,end)))
}

shit = apply(reps,1,fastout)
shi2 = as.vector(shit)
writeLines(shi2,"~/sa/datasets/Matina_1.6_v1.1/mazhao_CACAO_Transposable_Element_Annotation_New.TMP")
tmp = read.dna("~/sa/datasets/Matina_1.6_v1.1/mazhao_CACAO_Transposable_Element_Annotation_New.TMP",format="fasta")
write.dna(tmp,file="~/sa/datasets/Matina_1.6_v1.1/mazhao_CACAO_Transposable_Element_Annotation_New.fasta",format="fasta")
file.

# ff=file("~/sa/datasets/Matina_1.6_v1.1/mazhao_CACAO_Transposable_Element_Annotation_New.fasta","w+")
# for(str in shi2){
# 	if(grepl("^>",str)){
# 		cat(str,file=ff)
# 	} else {
# 		cat
# 	}
# 	cat("\n",file=ff)
# #(as.matrix(shi2),file="~/sa/datasets/Matina_1.6_v1.1/mazhao_CACAO_Transposable_Element_Annotation_New.fasta",
# }
# flush(ff)
# close(ff)

