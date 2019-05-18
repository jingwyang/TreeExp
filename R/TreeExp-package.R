
#' @title the TreeExp package for phylogenetic analysis of expression data
#'
#' @description \code{TreeExp} is an R package developed to provides useful
#' phylogenetic tools applicable to RNA-seq data. The package can be applied to
#' comparative expression evolution analysis, which includes but not liminited to:
#'
#' (i) pairwise expression distance estimation;
#'
#' (ii) relative rate test for transcriptome evolution;
#'
#' (iii) the strength of expression conservation estimation;
#'
#' (iv) ancestral transcriptome inference.
#'
#' Statistical methods implemented in
#' the package was based on Ornstein-Uhlenbeck (OU) model of transcriptome
#' evolution which claims that expression changes are constrained by stabilizing selection.
#'
#' @name TreeExp
#'
#' @docType package
#'
#' @import ape
#' @import phytools
#' @import stats
#' @import utils
#' @import parallel
#' @import reshape2
#' @import methods
#' @author Jingwen Yang (jingw.yang.fudan@gmail.com); Hang Ruan (hang.ruan@hotmail.com).
#'
NULL
