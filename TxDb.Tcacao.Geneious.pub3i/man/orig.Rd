\name{orig}
\alias{orig}
\docType{data}
\title{
An object containing gene model annotation for Matina v1.1 reference genome
}
\description{
An object (from rtracklayer:::import.gff3) that contains the original annotation 
data for version v1.1 of the Matina 1-6 reference genome.  Downloaded from 
www.cacaogenomedb.com (cacao11genes_pub3i.gff.gz)
}
\usage{data("orig")}
\format{
  The format is:
Formal class 'GRanges' [package "GenomicRanges"] with 6 slots
  ..@ seqnames       :Formal class 'Rle' [package "S4Vectors"] with 4 slots
  .. .. ..@ values         : Factor w/ 171 levels "cacao1chloroplast",..: 3 54 94 124 134 147 153 159 166 8 ...
  .. .. ..@ lengths        : int [1:171] 106159 96198 85415 80872 89502 65703 46021 56208 106736 50772 ...
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots
  .. .. ..@ start          : int [1:792382] 1369 1369 1369 1508 1895 1895 2055 2055 2497 2497 ...
  .. .. ..@ width          : int [1:792382] 1468 1468 385 246 73 73 183 183 340 248 ...
  .. .. ..@ NAMES          : NULL
  .. .. ..@ elementType    : chr "integer"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ strand         :Formal class 'Rle' [package "S4Vectors"] with 4 slots
  .. .. ..@ values         : Factor w/ 3 levels "+","-","*": 1 2 1 2 1 2 1 2 1 3 ...
  .. .. ..@ lengths        : int [1:22923] 294 71 14 17 84 56 8 5 38 6 ...
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ elementMetadata:Formal class 'DataFrame' [package "S4Vectors"] with 6 slots
  .. .. ..@ rownames       : NULL
  .. .. ..@ nrows          : int 792382
  .. .. ..@ listData       :List of 33
  .. .. .. ..$ source   : Factor w/ 3 levels "evg3","evg3alt",..: 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..$ type     : Factor w/ 4 levels "gene","mRNA",..: 1 2 3 4 3 4 3 4 3 4 ...
  .. .. .. ..$ score    : num [1:792382] 2630 2630 0 0 16 0 0 0 0 0 ...
  .. .. .. ..$ phase    : int [1:792382] NA NA NA 0 NA 0 NA 2 NA 2 ...
  .. .. .. ..$ ID       : chr [1:792382] "Thecc1EG000001" "Thecc1EG000001t1" NA NA ...
  .. .. .. ..$ gene     : chr [1:792382] NA "Thecc1EG000001" NA NA ...
  .. .. .. ..$ quality  :Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots
  .. .. .. .. .. ..@ elementType    : chr "character"
  .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots
  .. .. .. .. .. .. .. ..@ end            : int [1:792382] 0 5 5 5 5 5 5 5 5 5 ...
  .. .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ unlistData     : chr [1:308193] "Class:Medium" "Express:None" "Homology:ParalogMedium" "Intron:None" ...
  .. .. .. ..$ aaSize   : chr [1:792382] NA "249" NA NA ...
  .. .. .. ..$ cdsSize  :Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots
  .. .. .. .. .. ..@ elementType    : chr "character"
  .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots
  .. .. .. .. .. .. .. ..@ end            : int [1:792382] 0 2 2 2 2 2 2 2 2 2 ...
  .. .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ unlistData     : chr [1:123290] "76%" "750/981" "62%" "618/977" ...
  .. .. .. ..$ Name     :Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots
  .. .. .. .. .. ..@ elementType    : chr "character"
  .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots
  .. .. .. .. .. .. .. ..@ end            : int [1:792382] 0 1 1 1 1 1 1 1 1 1 ...
  .. .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ unlistData     : chr [1:76690] "Uncharacterized protein (39%P)" "Cystathionine beta-synthase (CBS) family protein (82%T)" "Cystathionine beta-synthase (CBS) family protein (80%T)" "Phospholipid/glycerol acyltransferase family protein (64%T)" ...
  .. .. .. ..$ oname    :Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots
  .. .. .. .. .. ..@ elementType    : chr "character"
  .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots
  .. .. .. .. .. .. .. ..@ end            : int [1:792382] 0 1 1 1 1 1 1 1 1 1 ...
  .. .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ unlistData     : chr [1:63075] "Unknown" "Uncharacterized protein (87%U)" "Uncharacterized protein (85%U)" "Uncharacterized protein (76%U)" ...
  .. .. .. ..$ ortholog :Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots
  .. .. .. .. .. ..@ elementType    : chr "character"
  .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots
  .. .. .. .. .. .. .. ..@ end            : int [1:792382] 0 1 1 1 1 1 1 1 1 1 ...
  .. .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ unlistData     : chr [1:87442] "11%" "frave:gene01181" "87%" "poptr:POPTR_0013s15730.2" ...
  .. .. .. ..$ paralog  :Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots
  .. .. .. .. .. ..@ elementType    : chr "character"
  .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots
  .. .. .. .. .. .. .. ..@ end            : int [1:792382] 0 2 2 2 2 2 2 2 2 2 ...
  .. .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ unlistData     : chr [1:84613] "Thecc1EG030493t1" "38%" "Thecc1EG034062t1" "51%" ...
  .. .. .. ..$ genegroup:Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots
  .. .. .. .. .. ..@ elementType    : chr "character"
  .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots
  .. .. .. .. .. .. .. ..@ end            : int [1:792382] 0 2 2 2 2 2 2 2 2 2 ...
  .. .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ unlistData     : chr [1:79312] "PLA9_G290" "31/31/1" "PLA9_G6641" "1/11/9" ...
  .. .. .. ..$ equiv1   :Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots
  .. .. .. .. .. ..@ elementType    : chr "character"
  .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots
  .. .. .. .. .. .. .. ..@ end            : int [1:792382] 0 1 1 1 1 1 1 1 1 1 ...
  .. .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ unlistData     : chr [1:63475] "CGD0000036/78.83" "CGD0000016/C99.77" "CGD0000016/C96.84" "CGD0000024/I100" ...
  .. .. .. ..$ equiv2   :Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots
  .. .. .. .. .. ..@ elementType    : chr "character"
  .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots
  .. .. .. .. .. .. .. ..@ end            : int [1:792382] 0 1 1 1 1 1 1 1 1 1 ...
  .. .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ unlistData     : chr [1:63747] "Tc01_t000030/10.23" "na" "na" "Tc01_t000050/C99.60" ...
  .. .. .. ..$ intron   :Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots
  .. .. .. .. .. ..@ elementType    : chr "character"
  .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots
  .. .. .. .. .. .. .. ..@ end            : int [1:792382] 0 2 2 2 2 2 2 2 2 2 ...
  .. .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ unlistData     : chr [1:111450] "0%" "0/6" "100%" "10/10" ...
  .. .. .. ..$ express  :Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots
  .. .. .. .. .. ..@ elementType    : chr "character"
  .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots
  .. .. .. .. .. .. .. ..@ end            : int [1:792382] 0 2 2 2 2 2 2 2 2 2 ...
  .. .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ unlistData     : chr [1:123290] "2%" "eq:0" "94%" "eq:75" ...
  .. .. .. ..$ oid      : chr [1:792382] NA "AUGpier8a:AUGpier8ap1s_1g1t1" NA NA ...
  .. .. .. ..$ score    : chr [1:792382] NA "2630" NA NA ...
  .. .. .. ..$ scorevec :Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots
  .. .. .. .. .. ..@ elementType    : chr "character"
  .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots
  .. .. .. .. .. .. .. ..@ end            : int [1:792382] 0 13 13 13 13 13 13 13 13 13 ...
  .. .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ unlistData     : chr [1:784931] "57" "196" "0" "0" ...
  .. .. .. ..$ protein  : chr [1:792382] NA "MSLMGKFSHMPKLQEVRIAFKKIDLTRAYEVRWLDYNHVLIHLSNEQDFNRIWIKQVWFIANQKIRVFKWTLEFESDKESTVVCVEHDYRKPPVDQVWIVVKNKETEIDCIVLGNKARPPGTTKSQP"| __truncated__ NA NA ...
  .. .. .. ..$ Parent   :Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots
  .. .. .. .. .. ..@ elementType    : chr "character"
  .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots
  .. .. .. .. .. .. .. ..@ end            : int [1:792382] 0 0 1 2 3 4 5 6 7 8 ...
  .. .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ unlistData     : chr [1:684089] "Thecc1EG000001t1" "Thecc1EG000001t1" "Thecc1EG000001t1" "Thecc1EG000001t1" ...
  .. .. .. ..$ est      : chr [1:792382] NA NA NA NA ...
  .. .. .. ..$ isoform  : chr [1:792382] NA NA NA NA ...
  .. .. .. ..$ Dbxref   :Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots
  .. .. .. .. .. ..@ elementType    : chr "character"
  .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots
  .. .. .. .. .. .. .. ..@ end            : int [1:792382] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ unlistData     : chr [1:87775] "TAIR:AT5G10860.1/82%" "UniRef50_B9I794/87%" "TAIR:AT5G10860.1/80%" "UniRef50_B9I794/85%" ...
  .. .. .. ..$ estgroup : chr [1:792382] NA NA NA NA ...
  .. .. .. ..$ rseq     :Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots
  .. .. .. .. .. ..@ elementType    : chr "character"
  .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots
  .. .. .. .. .. .. .. ..@ end            : int [1:792382] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ unlistData     : chr [1:595535] "141/141" "B_g11013t00001" "66/263" "cacao3v2sc1Loc1222t24" ...
  .. .. .. ..$ intr     : chr [1:792382] NA NA NA NA ...
  .. .. .. ..$ pro      :Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots
  .. .. .. .. .. ..@ elementType    : chr "character"
  .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots
  .. .. .. .. .. .. .. ..@ end            : int [1:792382] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ unlistData     : chr [1:218461] "258/258" "m7ricco:29907.m000629" "82/82" "m7ricco:29907.m000629" ...
  .. .. .. ..$ rx       : chr [1:792382] NA NA NA NA ...
  .. .. .. ..$ groupname:Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots
  .. .. .. .. .. ..@ elementType    : chr "character"
  .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots
  .. .. .. .. .. .. .. ..@ end            : int [1:792382] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. .. ..@ unlistData     : chr [1:2268] "CheY-like two-component responsive regulator family protein" "CheY-like two-component responsive regulator family protein" "Basic helix-loop-helix DNA-binding superfamily protein" "Serine/threonine-protein kinase DDB_G0267514" ...
  .. .. .. ..$ terepeat : chr [1:792382] NA NA NA NA ...
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ seqinfo        :Formal class 'Seqinfo' [package "GenomeInfoDb"] with 4 slots
  .. .. ..@ seqnames   : chr [1:171] "cacao1chloroplast" "cacao1mito" "scaffold_1" "scaffold_100" ...
  .. .. ..@ seqlengths : int [1:171] NA NA NA NA NA NA NA NA NA NA ...
  .. .. ..@ is_circular: logi [1:171] NA NA NA NA NA NA ...
  .. .. ..@ genome     : chr [1:171] NA NA NA NA ...
  ..@ metadata       : list()
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(orig)
## maybe str(orig) ; plot(orig) ...
}
\keyword{datasets}
