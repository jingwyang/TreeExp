
# Internal function for estimating Euclidean distance
.euc.dist = function(x,y) sqrt(sum((x-y)^2))

# Internal function for estimating Cosine similarity
.cosine.sim = function(x,y) sum(x*y) / (sqrt(sum(x^2))*sqrt(sum(y^2)))

# Internal function for estimating Kullback-Leibler Divergence
.kld = function(x,y) sum(x * log2(x/y))

# Internal function for estimating Tanimoto distance
.tani.dist = function(x,y) sum(pmax(x,y) - pmin(x,y)) / sum(pmax(x,y))

# Internal function for estimating Jaccard similarity
.jac.sim = function(x,y) sum(x*y) / (sum(x^2) + sum(y^2) - sum(x*y))

#' @title Internal functions for estimating pair-wise expression distances
#'
#' @name distances
#'
#' @rdname distances
#' @param expMat an exprssion level matrix
#' @description Several published methods to estimate pair-wise expression distances
#'
#' @references
#' Chen H, He X. 2016. The Convergent Cancer Evolution toward a Single Cellular Destination.
#' Mol Biol Evol 33:4-12
#'
#' Gu X, Su Z. 2007. Tissue-driven hypothesis of genomic evolution and sequence-expression correlations.
#' Proc Natl Acad Sci USA 104:2779-2784.
#'
#' Pereira V, Waxman D, Eyre-Walker A. 2009. A problem with the correlation coefficient as a measure of gene expression divergence.
#' Genetics 183:1597-1600.
#'
#' Sudmant PH, Alexis MS, Burge CB. 2015. Meta-analysis of RNA-seq expression data across species, tissues and studies.
#' Genome Biol 16:287.
#'
# Jesen-Shannon divergence
#' @rdname distances
#' @return returns expression distance matrix
#' @examples
#' data('tetraExp')
#' expression_table <- exptabTE(tetraExp, taxa = "all",
#' subtaxa = "Brain")
#' dismat <- dist.pea(expression_table)
#'
#' @export
dist.jsd = function (expMat = NULL) {

  object_n <- ncol(expMat)
  gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nrow = object_n, ncol = object_n)

  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      if (any(expMat[,i] == 0) || any(expMat[,j] == 0)) {

        mu <- (expMat[,i] + expMat[,j] + .02) / 2
        dis.mat[j,i] <- sqrt(.5 * .kld(expMat[,i]+.01, mu) + .5 * .kld(expMat[,j]+.01, mu))

      } else {

        mu <- (expMat[,i] + expMat[,j]) / 2
        dis.mat[j,i] <- sqrt(.5 * .kld(expMat[,i], mu) + .5 * .kld(expMat[,j], mu))

      }

    }

  }
  colnames(dis.mat) <- colnames(expMat)
  rownames(dis.mat) <- colnames(dis.mat)
  dis.mat + t(dis.mat)

}

# Pearson distance
#' @rdname distances
#'
#' @export
dist.pea = function (expMat = NULL) {

  object_n <- ncol(expMat)
  gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nrow = object_n, ncol = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      dis.mat[j,i] <- 1 - cor(expMat[,i],expMat[,j])

    }

  }

  #browser()
  colnames(dis.mat) <- colnames(expMat)
  rownames(dis.mat) <- colnames(dis.mat)
  dis.mat  + t(dis.mat)

}

# Spearman distance
#' @rdname distances
#'
#' @export
dist.spe = function (expMat = NULL) {

  object_n <- ncol(expMat)
  gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nrow = object_n, ncol = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      dis.mat[j,i] <- 1 - cor(expMat[,i],expMat[,j],
                              method = "spearman")

    }

  }

  #browser()
  colnames(dis.mat) <- colnames(expMat)
  rownames(dis.mat) <- colnames(dis.mat)
  dis.mat  + t(dis.mat)

}

# Euclidean distance
#' @rdname distances
#'
#' @export
dist.euc = function (expMat = NULL) {


  object_n <- ncol(expMat)
  gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nrow = object_n, ncol = object_n)

  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      dis.mat[j,i] <- .euc.dist(expMat[,i],expMat[,j])
    }

  }

  colnames(dis.mat) <- colnames(expMat)
  rownames(dis.mat) <- colnames(dis.mat)
  dis.mat  + t(dis.mat)

}

# Cosine distance
#' @rdname distances
#'
#' @export dist.cos
dist.cos = function (expMat = NULL) {

  object_n <- ncol(expMat)
  gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nrow = object_n, ncol = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      dis.mat[j,i] <- 1-.cosine.sim(expMat[,i],expMat[,j])

    }

  }
  colnames(dis.mat) <- colnames(expMat)
  rownames(dis.mat) <- colnames(dis.mat)
  dis.mat  + t(dis.mat)

}

# Tanimoto distance
#' @rdname distances
#'
#' @export
dist.tani = function (expMat = NULL) {

  object_n <- ncol(expMat)
  gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nrow = object_n, ncol = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      dis.mat[j,i] <- .tani.dist(expMat[,i],expMat[,j])

    }

  }
  colnames(dis.mat) <- colnames(expMat)
  rownames(dis.mat) <- colnames(dis.mat)
  dis.mat  + t(dis.mat)

}

# Jaccard distance
#' @rdname distances
#'
#' @export
dist.jac = function (expMat = NULL) {

  object_n <- ncol(expMat)
  gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nrow = object_n, ncol = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      dis.mat[j,i] <- 1 - .jac.sim(expMat[,i],expMat[,j])

    }

  }

  colnames(dis.mat) <- colnames(expMat)
  rownames(dis.mat) <- colnames(dis.mat)
  dis.mat  + t(dis.mat)

}

# Converntional expression distance
#' @rdname distances
#'
#' @export
dist.ced = function (expMat = NULL) {

  object_n <- ncol(expMat)
  gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nrow = object_n, ncol = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      #dis.mat[j,i] <- V11+V22-2*V12
      dis.mat[j,i] <- (.euc.dist(expMat[,i],expMat[,j]))^2 / gene_n

    }

  }
  colnames(dis.mat) <- colnames(expMat)
  rownames(dis.mat) <- colnames(dis.mat)
  dis.mat  + t(dis.mat)

}

# Distance based on stationary Ornstein-Uhlenback model
#
#' @rdname distances
#'
#' @export
dist.sou = function (expMat = NULL) {

  object_n <- ncol(expMat)
  #gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nrow = object_n, ncol = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      V11 <- var(expMat[,i])
      V22 <- var(expMat[,j])
      V12 <- cov(expMat[,i], expMat[,j])

      if (V12 > 0) {
        dis.mat[j,i] <- -log(V12/sqrt(V11*V22))
      } else {
        dis.mat[j,i] <- 0
        warning(paste0(date(),
            sprintf(": correlation between %d and %d may be negative or
                    equals zero, replace with 0", j, i)))
      }

    }

  }
  colnames(dis.mat) <- colnames(expMat)
  rownames(dis.mat) <- colnames(dis.mat)
  dis.mat + t(dis.mat)

}

