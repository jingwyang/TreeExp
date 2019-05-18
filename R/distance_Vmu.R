# Internal function for estimation of parameter pi which measures the variance of optima among genes under sOU model

.est.pi = function (expT){
  species_num <- ncol(expT)
  species_mean_exp_levels <- apply(expT, 1, mean)
  Var_N_genes <- var(species_mean_exp_levels)

  ############calculate inversed correlation matrix ###############
  dis.mat <- dist.pea(expT)
  corr.mat <- as.matrix(1 - as.dist(dis.mat))
  diag(corr.mat) <- 1
  corrMatInv <- solve(corr.mat)

  #############calculate mean of the parameter Q #####

  Q_gene = estParaQ(exptable = expT, corrmatinv = corrMatInv )
  Q_mean = mean(Q_gene)

  ###########calculate W when alpha is big  #####

  W = species_num / Q_mean

  ########## calculate pi ######

  est_pi <- Var_N_genes * W / (1+ W * Var_N_genes)
  est_pi
}

#' @title Estimate sOU expression distance when expression optima vary among genes
#'
#' @description estimate sOU expression distance when expression optima vary among
#'  genes (sOU_v distance)
#'
#' @name dist.sou_v
#' @rdname dist.sou_v
#' @param exptable an exprssion level table: column corresponds to median
#' expression value of all biological samples within one taxa_subtaxa group;
#' row corresponds to othologous genes
#' ('\code{spe}' for Spearman's correlation coefficient; \code{pea}' for
#' Pearson's correlation correlation coefficient).
#'
#' @return returns a list containing estimated expression distance
#' (sOU_v) matrix and parameter \code{pi} which measures the variance
#'  of optima among genes under sOU model
#' @author Jingwen Yang
#' @examples
#' data('tetraExp')
#' exp_table <- exptabTE(tetraExp, taxa = "all", subtaxa = "heart")
#' dismat <- dist.sou_v(exptable = exp_table)
#' dismat
#' @export



dist.sou_v = function (exptable = NULL){

  Distance <- list()
  species_num <- ncol(exptable)
  pi = .est.pi (exptable)

  dis.mat <- matrix(0, nrow = species_num, ncol = species_num)

  for ( i in 1: (species_num-1)) {
    for (j in (i+1):species_num){
      dis.mat[j,i] <- -log(1- (1-cor(exptable[,i], exptable[,j]))/(1-pi))
    }
  }

  colnames(dis.mat) <- colnames(exptable)
  rownames(dis.mat) <- colnames(dis.mat)

  Distance$pi <- pi
  Distance$distance <- dis.mat +t(dis.mat)
  Distance
}
