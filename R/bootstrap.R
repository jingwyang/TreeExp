# Internal function for estimating Euclidean distance
.FUN = function(x, method, ourgroup) {
    mat <- switch (method,
                       sou = {dist.sou(t(x))},
                       pea = {dist.pea(t(x))},
                       spe = {dist.spe(t(x))},
                       euc = {dist.euc(t(x))},
                       cos = {dist.cos(t(x))},
                       jsd = {dist.jsd(t(x))},
                       tani = {dist.tani(t(x))},
                       jac = {dist.jac(t(x))},
                       ced = {dist.ced(t(x))},
                       sou_v = {(dist.sou_v(t(x)))$distance}
    )
    root(NJ(mat), ourgroup, resolve.root = TRUE)
    }


#' @title Bootstrapping expression phylogeny
#'
#' @description  bootstrap by resampling gene (gene, transcript, exon, etc..)
#'
#' @name boot.exphy
#'
#' @param phy an object of class \code{phylo}.
#' @param x an exprssion level table: column corresponds to median expression value of all biological samples
#' within one taxa_subtaxa group; row corresponds to othologous genes
#' @param outgroup a vector of mode numeric or character specifying the outgroup
#' @param method  specifying which distance method to be used
#' to estimate expression phylogeny in bootstrapping.
#' @param B the number of bootstrap replicates.
#' @param block the number of columns in x that will be resampled together (see details).
#' @param trees a logical specifying whether to return the bootstrapped trees (FALSE by default).
#' @param quiet a logical specifying whether to print more information on the screen
#' while performing bootstrapping(FALSE by default).
#' @param rooted if "phy" is a rooted tree, a character of the root node's label when constructing "phy";
#' if "phy" is unrooted tree, NULL (NULL by default).
#' @param jumble a logical value. By default, the rows of x are randomized to avoid artificially
#' too large bootstrap values associated with very short branches.
#' @param mc.cores the number of cores (CPUs) to be used (passed to parallel).
#' @return similar to \code{boot.phylo} in \code{ape}, \code{boot.exphy} returns a numeric vector
#' which \emph{i}-th element is the number associated to the \emph{i}-th node of \code{phy}.
#' If trees = TRUE, \code{boot.exphy} returns a list whose first element (named "BP") is like before,
#' and the second element("trees") is a list with the bootstrapped trees.
#'
#' @rdname boot.exphy
#'
#' @examples
#' library('ape')
#' data(tetraExp)
#' dismat <- expdist(tetraExp, taxa = "all",
#'                  subtaxa = "Brain",
#'                  method = "pea")
#' tr <- root(NJ(dismat), "Chicken_Brain", resolve.root = TRUE)
#' plot(tr)
#' exp_table <- exptabTE(tetraExp, taxa = "all", subtaxa = "Brain")
#' bs <- boot.exphy(phy = tr, x = exp_table, method = "sou",
#'                  B = 100, outgroup = "Chicken_Brain")
#' nodelabels(bs)
#'
#' @export

boot.exphy = function (phy = NULL, x, outgroup = NULL, method = c("sou", "pea", "spe","euc", "cos", "jsd",
                                                 "tani", "jac","sou_v"), B = 100, block = 1, trees = FALSE,
                       quiet = FALSE, rooted = is.rooted(phy), jumble = TRUE, mc.cores = 1) {
    x<-t(x)
    method<-match.arg(method)
    if (is.null(dim(x)) || length(dim(x)) != 2)
      stop("the data 'x' must have two dimensions (e.g., a matrix or a data frame)")
    if (anyDuplicated(rownames(x)))
      stop("some labels are duplicated in the data: you won't be able to analyse tree bipartitions")
    boot.tree <- vector("list", B)
    y <- nc <- ncol(x)
    nr <- nrow(x)
    if (block > 1) {
      a <- seq(1, nc - 1, block)
      b <- seq(block, nc, block)
      y <- mapply(":", a, b, SIMPLIFY = FALSE)
      getBootstrapIndices <- function() unlist(sample(y, replace = TRUE))
    }
    else getBootstrapIndices <- function() sample.int(y, replace = TRUE)
    if (!quiet) {
      prefix <- "\rRunning bootstraps:      "
      suffix <- paste("/", B)
      updateProgress <- function(i) cat(prefix, i, suffix)
    }
    if (mc.cores == 1) {
      for (i in 1:B) {
        boot.samp <- x[, getBootstrapIndices()]
        if (jumble)
          boot.samp <- boot.samp[sample.int(nr), ]
        boot.tree[[i]] <- .FUN(boot.samp, method, outgroup)
        if (!quiet && !(i%%100))
          updateProgress(i)
      }
    }
    else {
      if (!quiet)
        cat("Running parallel bootstraps...")
      foo <- function(i) {
        boot.samp <- x[, getBootstrapIndices()]
        if (jumble)
          boot.samp <- boot.samp[sample.int(nr), ]
        .FUN(boot.samp, method, outgroup)
      }
      boot.tree <- mclapply(1:B, foo, mc.cores = mc.cores)
      if (!quiet)
        cat(" done.")
    }
    if (!quiet)
      cat("\nCalculating bootstrap values...")
    if (jumble) {
      boot.tree <- .compressTipLabel(boot.tree, ref = phy$tip.label)
      boot.tree <- .uncompressTipLabel(boot.tree)
      boot.tree <- unclass(boot.tree)
    }
    if (rooted) {
      pp <- prop.part(boot.tree)
      ans <- prop.clades(phy, part = pp, rooted = rooted)
    }
    else {
      phy <- reorder(phy, "postorder")
      ints <- phy$edge[, 2] > Ntip(phy)
      ans <- countBipartitions(phy, boot.tree)
      ans <- c(B, ans[order(phy$edge[ints, 2])])
    }
    if (!quiet)
      cat(" done.\n")
    if (trees) {
      class(boot.tree) <- "multiPhylo"
      ans <- list(BP = ans, trees = boot.tree)
    }
    ans
  }
