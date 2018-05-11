# BSgenome.Tcacao
R packages containing cacao reference genomes in bioconductor (BSgenome) format.


## Installation steps
1. Install prerequisite R packages (devtools, BSgenome)
2. Install BSgenome.Tcacao packages from github

In R, 

        # install devtools (if needed)
        install.packages("devtools")
        
        # install BSgenome from bioconductor
        source("https://bioconductor.org/biocLite.R")
        biocLite(c("BSgenome","VariantAnnotation","AnnotationDbi","GenomicFeatures"))
        
        # install Matina 1-6 reference genome package and accompanying gene model annotations
        devtools:::install_github("ConradStack/BSgenome.Tcacao/BSgenome.Tcacao.MarsInc.v11")
        devtools:::install_github("ConradStack/BSgenome.Tcacao/TxDb.Tcacao.Geneious.pub3i")

        # install Criollo reference genome
        devtools:::install_github("ConradStack/BSgenome.Tcacao/BSgenome.Tcacao.CIRAD.v10")
