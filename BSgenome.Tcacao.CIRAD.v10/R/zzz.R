###
###

.pkgname <- "BSgenome.Tcacao.CIRAD.v10"

.seqnames <- sprintf("Tc%02d",c(0:10))

.circ_seqs <- NULL

.mseqnames <- NULL

.onLoad <- function(libname, pkgname)
{
    if (pkgname != .pkgname)
        stop("package name (", pkgname, ") is not ",
             "the expected name (", .pkgname, ")")
    extdata_dirpath <- system.file("extdata", package=pkgname,
                                   lib.loc=libname, mustWork=TRUE)

    ## Make and export BSgenome object.
    bsgenome <- BSgenome(
        organism="Theobroma cacao",
        species="cacao",
        provider="CIRAD",
        provider_version="v10",
        release_date="Oct. 2012",
        release_name="CriolloCIRAD",
        source_url="http://cocoagendb.cirad.fr",
        seqnames=.seqnames,
        circ_seqs=.circ_seqs,
        mseqnames=.mseqnames,
        seqs_pkgname=pkgname,
        seqs_dirpath=extdata_dirpath
    )

    ns <- asNamespace(pkgname)

    objname <- pkgname
    assign(objname, bsgenome, envir=ns)
    namespaceExport(ns, objname)

    old_objname <- "criollo"
    assign(old_objname, bsgenome, envir=ns)
    namespaceExport(ns, old_objname)
}

