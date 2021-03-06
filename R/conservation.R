#' @title Generate an inversed correlation matrix between expression profiles of species
#'
#' @name corrMatInv
#'
#' @description Generate an inversed correlation matrix between expression profiles of species
#'
#' @param objects a vector of objects of class \code{taxonExp} or an object of class \code{taxaExp}
#' @param taxa one single character or a vector of characters specifying main taxa to generate
#' an inversed correlation matrix.
#' If one single character "all" is given,
#' all the taxa included in the \code{taxaExp} will be matched and included ("all" by default).
#' @param subtaxa one single character specifying sub taxa to be included in generating
#' an inversed corrlation matrix.
#' @param method specifying which distance method (Spearman or Pearson) to be used
#' to estimate correlations between expression profiles ("spe" by default)
#'
#' @return returns an  inversed correlation matrix
#' @examples
#' data('tetraExp')
#' species.group <- c("Human", "Chimpanzee", "Bonobo", "Gorilla",
#' "Macaque", "Mouse", "Opossum", "Platypus")
#' inv.corr.mat <- corrMatInv(tetraExp, taxa = species.group, subtaxa = "Brain")
#'
#' @export
corrMatInv = function(objects = NULL, taxa = "all", subtaxa = NULL, method = c("spe", "pea")) {

  method <- match.arg(method)

  dis.mat <- expdist(objects, taxa = taxa,subtaxa = subtaxa, method = method)
  corr.mat <- as.matrix(1 - as.dist(dis.mat))
  diag(corr.mat) <- 1

  solve(corr.mat)
}

#' @title Estimation of parameters of selection strength Gamma Distribution
#'
#' @name estParaGamma
#'
#' @description Estimation of parameters of selection strength Gamma Distribution
#'
#' @param exptable an expression level table: column corresponds to median expression value of all biological samples
#' within one taxa_subtaxa group; row corresponds to othologous genes
#' @param corrmatinv an inversed correlation matrix between specie
#'
#' @return returns a vector of parameters estimated from gamma distribution
#' @examples
#' data('tetraExp')
#' species.group <- c("Human", "Chimpanzee", "Bonobo", "Gorilla",
#' "Macaque", "Mouse", "Opossum", "Platypus")
#' inv.corr.mat <- corrMatInv(tetraExp, taxa = species.group, subtaxa = "Brain")
#' brain.exptable <- exptabTE(tetraExp, taxa = species.group, subtaxa = "Brain" ,logrithm = TRUE)
#' gamma.paras <- estParaGamma(exptable =brain.exptable, corrmatinv =inv.corr.mat)
#' @export
estParaGamma = function(exptable = NULL, corrmatinv = NULL) {

  species_mean_exp_levels <- apply(exptable, 1, mean) # mean expression value of orthologous genes

  species_num <- ncol(exptable)
  gene_num <- nrow(exptable)

  if (ncol(corrmatinv) != species_num) {

    stop(paste0(date(),": species number in expression table and matrix do not match!"))

  }

  Q_gene <- vector(mode = "numeric", length = gene_num)

  for (k in 1:gene_num) { # for each orthologous gene

    for (i in 1:species_num) { # for each species

      for (j in 1:species_num) {

        Q_gene[k] = Q_gene[k] + corrmatinv[i,j] *
          (exptable[k,i] - species_mean_exp_levels[k]) *
          (exptable[k,j] - species_mean_exp_levels[k])

      }
    }
  }

  Q_fstm_est <- mean(Q_gene)
  Q_sndm_est <- mean(Q_gene^2)
  H_est <- Q_sndm_est / Q_fstm_est^2

  alpha_est <- 2 + (1 + 2/species_num) / (H_est - (1 + 2/species_num))
  W_est <- alpha_est / (alpha_est - 1) * species_num / Q_fstm_est

  para_vec <- list("alpha" = alpha_est, "W_average" = W_est, "speNum" = species_num, "geneNum" = gene_num)

  para_vec
}

#' @title Estimation of parameter Q from expression table
#'
#' @name estParaQ
#'
#' @description Estimation of parameter Q from expression table
#'
#' @param exptable an exprssion level table: column corresponds to median expression value of all biological samples
#' within one taxa_subtaxa group; row corresponds to othologous genes
#' @param corrmatinv an inversed correlation matrix between species
#'
#' @return returns a vector of parameter Qs from given expression table
#' @examples
#' data('tetraExp')
#' species.group <- c("Human", "Chimpanzee", "Bonobo", "Gorilla",
#' "Macaque", "Mouse", "Opossum", "Platypus")
#' inv.corr.mat <- corrMatInv(tetraExp, taxa = species.group, subtaxa = "Brain")
#' brain.exptable <- exptabTE(tetraExp, taxa = species.group, subtaxa = "Brain" ,logrithm = TRUE)
#' brain.Q <- estParaQ(brain.exptable, corrmatinv = inv.corr.mat)
#' @export
estParaQ = function(exptable = NULL, corrmatinv = NULL) {

  species_mean_exp_levels <- apply(exptable, 1, mean) # mean expression value of orthologous genes

  species_num <- ncol(exptable)
  gene_num <- nrow(exptable)

  if (ncol(corrmatinv) != species_num) {

    stop(paste0(date(),": species number in expression table and matrix do not match!"))

  }

  Q_gene <- vector(mode = "numeric", length = gene_num)

  for (k in 1:gene_num) { # for each orthologous gene

    for (i in 1:species_num) { # for each species

      for (j in 1:species_num) {

        Q_gene[k] = Q_gene[k] + corrmatinv[i,j] *
          (exptable[k,i] - species_mean_exp_levels[k]) *
          (exptable[k,j] - species_mean_exp_levels[k])

      }
    }
  }

  Q_gene

}

#' @title Bayesian estimation of expression conservation from previous data set
#'
#' @name  estParaWBayesian
#'
#' @description Bayesian estimation of expression conservation from previous data set
#'
#' @param Q_gene a vector specifying parameter Q estimated from \code{estParaQ}
#' @param gammaparas a vector specifying parameters estimated from \code{estParaGamma}
#'
#' @return returns a vector of exptations and confidence intervals of parameter Ws estimated from previous Q estimation
#' @examples
#' data('tetraExp')
#' species.group <- c("Human", "Chimpanzee", "Bonobo", "Gorilla",
#' "Macaque", "Mouse", "Opossum", "Platypus")
#' inv.corr.mat <- corrMatInv(tetraExp, taxa = species.group, subtaxa = "Brain")
#' brain.exptable <- exptabTE(tetraExp, taxa = species.group, subtaxa = "Brain" ,logrithm = TRUE)
#' gamma.paras <- estParaGamma(exptable =brain.exptable, corrmatinv =inv.corr.mat)
#' brain.Q <- estParaQ(brain.exptable, corrmatinv = inv.corr.mat)
#' brain.post<- estParaWBayesian(brain.Q, gamma.paras)
#' @export
estParaWBayesian = function(Q_gene = NULL, gammaparas = NULL) {

  W_gene_exp <- numeric(length = length(Q_gene))

  alpha_est <- gammaparas$alpha
  W_est <- gammaparas$W_average
  species_num <- gammaparas$speNum

  for (i in 1:length(Q_gene)) {

    W_gene_exp[i] <- (alpha_est + species_num / 2) /
      (alpha_est + Q_gene[i] * W_est) * W_est

  }

  W_gene_var <- (W_est ^ 2) / (alpha_est + species_num / 2) * W_gene_exp ^ 2

  tmp <- sqrt(W_gene_var) * qgamma(0.025, alpha_est)
  W_gene_ci95 <- cbind(W_gene_exp - tmp, W_gene_exp + tmp)

  list(w=W_gene_exp, ci95=W_gene_ci95)
}

