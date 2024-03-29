% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_GeneSymbol.R
\name{get_GeneEnsembl2Symbol}
\alias{get_GeneEnsembl2Symbol}
\title{ID convention (from ENSEMBL to SYMBOL)}
\usage{
get_GeneEnsembl2Symbol(ensemble.arr, dataset = "mmusculus_gene_ensembl")
}
\arguments{
\item{ensemble.arr, }{a character array. The ENSEMBL ID array.}

\item{dataset, }{must be either "mmusculus_gene_ensembl" or "hsapiens_gene_ensembl".}
}
\value{
A data frame contains 3 columns, "Gene", "Symbol" and "description".
}
\description{
This function takes an array of ENSEMBL ID and convert it into GENE SYMBOL.
}
