
setwd("~/sa/datasets/Criollo")
# Criollo
library(GenomicRanges); 
library(rtracklayer)
library(GenomicFeatures)
require(AnnotationDbi)
require(BSgenome.Tcacao.CIRAD.v10)
require(VariantAnnotation)
genome <- BSgenome.Tcacao.CIRAD.v10

# Adding CDS regions, or trying to:
# gt gff3 -force yes -sort yes -retainids yes -addids yes -tidy -o Theobroma_cacao_v1.sorted.gff3 -v Theobroma_cacao_v1.gff3
# gt cds -force yes -matchdesc yes -v yes -startcodon yes -finalstopcodon yes -seqfile Theobroma_cacao_v1.0.pseudomolecule.fna -o Theobroma_cacao_v1.CDS.gff3 Theobroma_cacao_v1.sorted.gff3
# (FAILED: gt keeps renaming stuff in the .cds.gff3 output file )

gff_file = "Theobroma_cacao_v1.sorted.gff3"
sqlite_file = "BSgenome.Tcacao.CIRAD.v10.sqlite"

#


#gff=import.gff3(gff_file, genome=seqinfo(genome))
#gff=import.gff3(gff_file)
minds = which(gff$type == "mRNA")
gall = gff[minds]
metas = mcols(gall)


trans = dbReadTable(con,"transcript")
stopifnot(length(minds)==nrow(trans))

sinfo = as.data.frame(seqinfo(genome))
junk = data.frame(chrom=rownames(sinfo),length=sinfo$seqlengths,is_circular=FALSE)
txdb = makeTranscriptDbFromGFF(gff_file,"gff3",chrominfo=junk,species="Theobroma cacao")
saveDb(txdb, file="BSgenome.Tcacao.CIRAD.v10.V2.sqlite")
#

txdb = loadDb(sqlite_file)
cnames = columns(txdb)
gffin = asGFF(txdb)

transcriptsBy(txdb)
