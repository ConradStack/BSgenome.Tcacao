\name{BSgenome.Tcacao.CIRAD.v10}
\docType{package}

\alias{BSgenome.Tcacao.CIRAD.v10-package}
\alias{BSgenome.Tcacao.CIRAD.v10}
\alias{criollo}

\title{Reference genome sequences for Criollo variety of Theobroma cacao (version 1.0)}

\description{
  Full genome sequences for Theobroma cacao
}

\details{
  
}

\note{
  This BSgenome data package was made from the following source data files:
  \preformatted{
-- information not available --

  }

  See \code{?\link[BSgenome]{BSgenomeForge}} and the BSgenomeForge
  vignette (\code{vignette("BSgenomeForge")}) in the \pkg{BSgenome}
  software package for how to make a BSgenome data package.
}

\author{The Bioconductor Dev Team}

\seealso{
  \itemize{
    \item \link[BSgenome]{BSgenome} objects and the
          \code{\link[BSgenome]{available.genomes}} function
          in the \pkg{BSgenome} software package.
    \item \link[Biostrings]{DNAString} objects in the \pkg{Biostrings}
          package.
    \item The BSgenomeForge vignette (\code{vignette("BSgenomeForge")})
          in the \pkg{BSgenome} software package for how to make a BSgenome
          data package.
  }
}

\examples{
BSgenome.Tcacao.CIRAD.v10
genome <- BSgenome.Tcacao.CIRAD.v10
seqlengths(genome)
genome$Tc01  # same as genome[["Tc01"]]

## See the GenomeSearching vignette in the BSgenome software
## package for some examples of genome-wide motif searching using
## Biostrings and the BSgenome data packages:
if (interactive())
    vignette("GenomeSearching", package="BSgenome")
}

\keyword{package}
\keyword{data}
