#' @title Expression distance matrix generated from a \code{taxaExp} object
#'
#' @name expdist
#' @description Generate an expression distance matrix from an object of \code{taxaExp} class
#' using a specified distance method
#'
#' @param objects a vector of objects of class \code{taxonExp} or an object of class \code{taxaExp}
#' @param taxa one single character or a vector of characters specifying main taxa selected for
#' calculating expression distance.
#' If one single character "all" is given,
#' all the taxa included in the \code{taxaExp} will be matched and selected ("all" by default).
#' @param subtaxa one single character or a vector of characters sepcifying sub taxa selected for
#' calculating expression distance.
#' If one singke character "all" is given,
#' all the subtaxa included in the \code{taxaExp} will be matched and selected ("all" by default).
#' @param rowindex a vector of numbers corresponded to indices of selecting rows
#' @param method specifying which distance method to be used
#' to estimate expression phylogeny in bootstrapping.
#' @param logrithm a logical specifying whether to apply expression value log2 tranformation (TRUE by default).
#'
#' @return returns an expression distance matrix
#'
#' @examples
#' data(tetraExp)
#' library('ape')
#' dismat <- expdist(tetraExp, taxa = "all",
#'                  subtaxa = "Brain",
#'                  method = "pea")
#' tr <- root(NJ(dismat), "Chicken_Brain")
#' plot(tr)
#'
#' @export
expdist = function (objects = NULL, taxa = "all", subtaxa = "all", rowindex = NULL,
                    method = c( "sou", "sou_v","pea", "spe","euc", "cos", "jsd",
                                "tani", "jac"), logrithm = TRUE)
{

  if (is.null(objects) || !is(objects) == "taxaExp") {
    stop(paste0(date(), ": no valid taxaExp objects input!"))
  }

  flag1 <- TRUE
  flag2 <- TRUE

  if (any(grepl("all",taxa, ignore.case = TRUE))) {flag1 = FALSE}
  else { taxa <- gsub("\\s+","",taxa)}

  if (any(grepl("all",subtaxa, ignore.case = TRUE))) {flag2 = FALSE}
  else { subtaxa <- gsub("\\s+","",subtaxa)}

  objects_n <- length(objects)
  objects_new_n <- 0

  if ( flag1 || flag2)

  {
    #browser()

    for (i in 1:objects_n)

    {
      if (flag1 && flag2) {
        if (any(grepl(objects[[i]]$taxon_name,taxa, ignore.case=TRUE))
            && any(grepl(objects[[i]]$subTaxon_name, subtaxa, ignore.case=TRUE)))
        {objects_new_n <- objects_new_n + 1}

      } else {
        if (any(grepl(objects[[i]]$taxon_name,taxa,ignore.case=TRUE))
            ||  any(grepl(objects[[i]]$subTaxon_name, subtaxa, ignore.case=TRUE)))
        {objects_new_n <- objects_new_n + 1}
      }

    }

    objects_new <- vector("list",length = objects_new_n)

    counter <- 1

    for (i in 1:objects_n)

    {

      if (flag1 && flag2) {
        if (any(grepl(objects[[i]]$taxon_name,taxa,ignore.case=TRUE))
            &&  any(grepl(objects[[i]]$subTaxon_name, subtaxa, ignore.case=TRUE)))
        {
          objects_new[[counter]] <- objects[[i]]
          counter <- counter + 1
        }

      } else {
        if (any(grepl(objects[[i]]$taxon_name,taxa,ignore.case=TRUE))
            ||  any(grepl(objects[[i]]$subTaxon_name, subtaxa, ignore.case=TRUE)))
        {
          objects_new[[counter]] <- objects[[i]]
          counter <- counter + 1
        }
      }

    }

    class(objects_new) <- "taxaExp"

    objects <- objects_new

  } else {

    objects_new <- vector("list", length = objects_new_n)
    counter <- 1

    for (i in 1:objects_n) {
        objects_new[[counter]] <- objects[[i]]
        counter <- counter + 1
    }

  }

  if (length(objects_new) == 0) {

    stop(paste0(date(),": taxa and subtaxa name not found."))

  }
  #browser()

  method<-match.arg(method)

  message(paste0(date(), ": using ", method, " to calculate pair-wise distance"))

  object_n <- length(objects)

  gene_n <- objects[[1]]$gene_num

  message(paste0(date(),": input ",object_n, " taxa"))
  message(paste0(date(),": total ", gene_n, " genes"))

  #initialization

  expVal <- matrix(0, nrow = gene_n, ncol = object_n)


  taxon_names <- vector("character", length = object_n)

  for (i in 1:object_n) {

    taxon_names[i] = paste0(objects[[i]]$taxon_name, "_", objects[[i]]$subTaxon_name)

    expVal[,i] = apply(objects[[i]]$exp_value,1,median)

  }

  if (!is.null(rowindex)) {

    expVal <- expVal[rowindex,]

  }

  if (logrithm) {

    expVal <- apply(expVal, c(1,2), function (x) log2(x+1))

  }
  #browser()
  colnames(expVal) <- taxon_names

  dis.mat <- switch (method,

    sou = {dist.sou(expVal)},

    pea = {dist.pea(expVal)},

    spe = {dist.spe(expVal)},

    euc = {dist.euc(expVal)},

    cos = {dist.cos(expVal)},

    jsd = {dist.jsd(expVal)},

    tani = {dist.tani(expVal)},

    jac = {dist.jac(expVal)},

    ced = {dist.ced(expVal)},
    sou_v = {dist.sou_v(expVal)}

  )

  dis.mat
  #as.dist(dis.mat)

}
