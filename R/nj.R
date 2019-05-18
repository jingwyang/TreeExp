#' @title Neighbor-joining
#'
#' @name NJ
#'
#' @rdname NJ
#'
#' @description The famous neighbor-joining tree estimation function from Saitou and Nei (1987).
#' @param X a distance matrix
#' @return returns an NJ-tree object of class "phylo"
#' @examples
#' library('ape')
#' data(tetraExp)
#' dismat <- expdist(tetraExp, taxa = "all",
#'                  subtaxa = "Brain",
#'                  method = "pea")
#' tr <- root(NJ(dismat), "Chicken_Brain")
#' plot(tr)
#'
#' @export
NJ = function (X) {

  reorder(nj(X), "postorder")

}
