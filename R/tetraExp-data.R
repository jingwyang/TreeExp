#' @title taxaExp object created from expression data
#'
#' @name tetraExp
#'
#' @description large \code{taxaExp} object created from
#' six tissues' expression data of nine tetrapod species
#'
#' @docType data
#'
#' @format a \code{taxaExp} object, contains 53 \code{taxonExp} with RPKM value on 5636 genes.
#'
#' @references
#' Brawand,D. et al. (2011) The evolution of gene expression levels in mammalian organs.
#' Nature, 478, 343-348.
#'
#' @examples
#' data(tetraExp)
#' tetraExp
#' tetraExp[[1]]
#' tetraExp[[1]]$exp_value
NULL
