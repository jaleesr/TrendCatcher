#' ID convention (from ENSEMBL to SYMBOL)
#'
#' This function takes an array of ENSEMBL ID and convert it into GENE SYMBOL.
#'
#' @param ensemble.arr, a character array. The ENSEMBL ID array.
#' @param dataset, must be either "mmusculus_gene_ensembl" or "hsapiens_gene_ensembl".
#' @return A data frame contains 3 columns, "Gene", "Symbol" and "description".
#'
#' @export
#'
#'
get_GeneEnsembl2Symbol<-function(ensemble.arr, dataset = "mmusculus_gene_ensembl"){
  ensembl = biomaRt::useMart("ensembl", dataset =dataset)
  dat = biomaRt::getBM(
    values = ensemble.arr,
    filters = c("ensembl_gene_id"),
    attributes = c("ensembl_gene_id", "external_gene_name", "description"),
    mart = ensembl
  )
  colnames(dat)<-c("Gene", "Symbol", "description")
  return(dat)
}
#' ID convention (from SYMBOL to ENSEMBL)
#'
#' This function takes an array of SYMBOL and convert it into GENE ENSEMBL.
#'
#' @param symbol.arr, a character array. The SYMBOL array.
#' @param dataset, must be either "mmusculus_gene_ensembl" or "hsapiens_gene_ensembl".
#' @return A data frame contains 3 columns, "Gene", original ID, "Symbol" and "description".
#'
#' @export
#'
#'
#'
get_Symbol2GeneEnsembl<-function(symbol.arr, dataset = "mmusculus_gene_ensembl"){
  ensembl = biomaRt::useMart("ensembl", dataset = dataset)
  dat = biomaRt::getBM(
    values = symbol.arr,
    filters = c("external_gene_name"),
    attributes = c("ensembl_gene_id", "external_gene_name", "description"),
    mart = ensembl
  )
  colnames(dat)<-c("Gene", "Symbol", "description")
  return(dat)
}
