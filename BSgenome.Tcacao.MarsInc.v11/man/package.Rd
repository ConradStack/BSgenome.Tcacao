\name{BSgenome.Tcacao.MarsInc.v11}
\docType{package}

\alias{BSgenome.Tcacao.MarsInc.v11-package}
\alias{BSgenome.Tcacao.MarsInc.v11}
\alias{matina16}

\title{Reference genome sequences for Matina 1-6 variety of Theobroma cacao (version 1.1)}

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
BSgenome.Tcacao.MarsInc.v11
genome <- BSgenome.Tcacao.MarsInc.v11
seqlengths(genome)
genome$scaffold_1  # same as genome[["scaffold_1"]]

## See the GenomeSearching vignette in the BSgenome software
## package for some examples of genome-wide motif searching using
## Biostrings and the BSgenome data packages:
if (interactive())
    vignette("GenomeSearching", package="BSgenome")
}

\keyword{package}
\keyword{data}
