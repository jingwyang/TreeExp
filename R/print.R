#' @title Consice display of \code{taxaExp} or \code{taxonExp} object
#' @param objects an object of class \code{taxaExp} or class \code{taxonExp}.
#' @param details a logical specifying whether to print taxa and subtaxa names.
#' @param printlen the number of biological replicates title to print (6 by default).
#' @param ... further arguments passed to or from other methods.
#' @return returns a consice display of \code{taxaExp} or \code{taxonExp} object
#' @rdname print
#'
#' @examples
#' data(tetraExp)
#' print(tetraExp, details = TRUE)
#' print(tetraExp[[1]], printlen = 6)
#'
#'
#' @export
#'
print <- function(objects, ...) {

  UseMethod("print")
}

#' @rdname print
#' @export

print.taxaExp <- function(objects, details = FALSE, ...) {

  N <- length(objects)
  cat("\n",N, "taxonExp objects", "\n")

  if (details) {
    cat("\n")
      for (i in 1:N) {
        cat ("object", i, ":", objects[[i]]$taxon_name,
             "\t", objects[[i]]$subTaxon_name, "\n")
      }
  }
}

#' @rdname print
#' @export

print.taxonExp <- function(objects, printlen = 6, ...) {

  cat("\nOne taxonExp object\n")

  cat("Taxon name: ", objects$taxon_name, "\n")

  cat("Subtaxon name: ", objects$subTaxon_name, "\n")

  cat("Total gene number: ", objects$gene_num, "\n")

  cat("Total bio replicates number: ", objects$bioRep_num, "\n")

  cat("Bio replicates titles:\n")

  if (objects$bioRep_num > printlen) {
    cat(paste0("\t", paste(objects$bioRep_id[1:printlen]),
               collapse = ", "), ", ...\n")

  } else {
    print(unlist(objects$bioRep_id))
  }
  cat("\n")
}

