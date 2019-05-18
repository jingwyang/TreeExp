#' @title Relative Rate Test for Gene Expression
#' @description statistically compare the relative expression rate
#' of a given set of genes between two taxa (or species)
#' @name RelaRate.test
#' @rdname RelaRate.test
#' @param expTable an exprssion level table: column corresponds to median
#' expression value of all biological samples within one taxa_subtaxa group;
#' row corresponds to othologous genes
#' @param x taxa name
#' @param y taxa name
#' @param outgroup outgroup name
#' @param Pi value of parameter Pi which measures the variance of optima
#'  among genes under sOU model. Pi = 0 by default.
#' @param alternative character string specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less".
#' @return A list containing the following components:
#' \itemize{
#' \item Z_score the z-statistic
#' \item alternative    records the value of the input argument alternative: "greater", "less" or "two.sided".
#' \item p.value    the p-value for the test
#' }
#' @author Jingwen Yang
#' @examples
#' data(tetraExp)
#' exp_table <-exptabTE(tetraExp, taxa = 'all', subtaxa = 'Brain',rowindex = 1:100)
#' ztest <- RelaRate.test(expTable = exp_table, x = 'human', y = 'chimpanzee',
#' outgroup = 'macaque', alternative = 'greater')
#' ztest
#' @export
#'
#'
RelaRate.test <- function( expTable = NULL, x = NULL, y = NULL, outgroup = NULL, Pi = 0,
                          alternative = c("two.sided", "less", "greater")
                          ){
  taxon_names = tolower(unlist(lapply(colnames(expTable), function(x) unlist(strsplit(x, "_"))[1])))
  x = tolower(x)
  y = tolower(y)
  alt = match.arg(alternative)
  geneN = nrow(expTable)
  outgroup = tolower(outgroup)
  if (! (x %in% taxon_names) ){
    stop(paste0(date(), ": species x is not in the expression table"))
  }
  if (! y %in% taxon_names){
    stop(paste0(date(), ": species y is not in the expression table"))
  }
  if (! outgroup %in% taxon_names){
    stop(paste0(date(), ": species outgroup is not in the expression table"))
  }

  x_colN = which (taxon_names  == x)
  y_colN = which (taxon_names ==  y)
  out_colN = which (taxon_names ==  outgroup)
  sub_expT = expTable[,c(x_colN, y_colN, out_colN)]

  dist.mat = dist.sou(sub_expT)
  Dab = dist.mat[2,1]
  Dac = dist.mat[3,1]
  Dbc = dist.mat[3,2]
  delta_D = Dac - Dbc
  Da = (Dab + Dac - Dbc)/2
  Db = Dab - Da
  Dc = (Dac+Dbc-Dab)/2
  if((Dac-Db)<0){
    Dc = 0
    }
  rho_ac = cor(sub_expT[,1], sub_expT[,2])
  rho_bc = cor(sub_expT[,2], sub_expT[,3])
  rho_c  = Pi + (1-Pi)* exp(-Dc)

  Var_Dac = ((1 - rho_ac^2)^2) / (geneN-3) / ((rho_ac - Pi)^3)
  Var_Dbc = ((1 - rho_bc^2)^2) / (geneN-3) / ((rho_bc - Pi)^3)
  Var_Dc =  ((1 - rho_c^2)^2)/ (geneN-3) / (rho_c - Pi)^2

  Var_deltaD = Var_Dac +Var_Dbc - 2*Var_Dc
  Z_score = delta_D / sqrt(Var_deltaD)

  ################
  if (alternative == "less") {
    pval <- pnorm(Z_score)
  }
  else if (alternative == "greater") {
    pval <- 1 - pnorm(Z_score)
  }
  else {
    pval <- 2 * pnorm(-abs(Z_score))
  }

  ######
  rval <- list(Z_score = Z_score, alternative = alternative, p.value = pval)

  rval
}

