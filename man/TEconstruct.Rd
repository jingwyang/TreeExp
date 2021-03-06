% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructor.R
\name{TEconstruct}
\alias{TEconstruct}
\title{Construct a taxaExp object.}
\usage{
TEconstruct(ExpValueFP = NULL, taxa = "all", subtaxa = "all",
  rmOut = FALSE, verbose = FALSE)
}
\arguments{
\item{ExpValueFP}{a text file contains processed expression data (RPKM or TPM for example).
Row names correspond with gene names,
and column names correspond with taxon and subtaxon names.}

\item{taxa}{one single string or a vector of strings specifying main taxa selected for
constructing \code{taxaExp} object.
Taxa names are extracted from row names given in gene length file.
If one single string "all" is given,
all the taxa in the row names will be matched and selected ("all" by default).}

\item{subtaxa}{one single string or a vector of strings sepcifying sub taxa selected for
constructing \code{taxaExp} object.
If one single string "all" is given,
all the sub taxa in the row names will be matched and selected ("all" by default).}

\item{rmOut}{a logical sepcifying whether to remove expression outliers
while constructing \code{taxaExp} objects (TRUE by default).}

\item{verbose}{a logical specifying whether to print more information on the screen
while constructing \code{taxaExp} objects (FALSE by default).}
}
\value{
returns an object of class \code{Taxa} (S3 class, a list of \code{taxonExp} objects).
}
\description{
Constructor function for \code{taxaExp} objects.
This function takes in a expression value file from RNA-seq,
and construct a \code{taxaExp} object from which user can extract information
for display or for further analysis.
}
\examples{

taxa.objects = TEconstruct(ExpValueFP = system.file('extdata/primate_brain_expvalues.txt',
 package = 'TreeExp'), taxa = "all", subtaxa = c("HIP", "CB"))

}
