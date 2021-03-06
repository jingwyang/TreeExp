
# Internal function for estimating gamma values

.solveAnEquation <- function( left = NULL ) {

  right <- matrix(c(1,0,1,1,0,0,
                    0,1,1,0,1,0,
                    1,1,0,1,1,0,
                    1,1,1,0,0,1,
                    0,1,0,1,0,1,
                    1,0,0,0,1,1),
              nrow=6,ncol=6,byrow=TRUE)


  ret<-solve(right, left)  # gammaD, gammaE, a, b, c, d

  #names(ret) <- c("gammaD","gammaE","a","b","c","d")

  ret

}


#' @title Estimating gamma from taxaExp objects.
#'
#' @name estgamma
#' @rdname gamma
#'
#' @param objects a vector of objects of class \code{taxonExp} or an object of class \code{taxaExp}
#' @param taxa one single character or a vector of characters specifying main taxa selected for
#' calculating gammaD and gammaE.
#' If one single character "all" is given,
#' all the taxa included in the \code{taxaExp} will be matched and selected ("all" by default).
#' @param subtaxa one single character or a vector of characters sepcifying sub taxa selected for
#' calculating gammaD and gammaE.
#'#' If one single character "all" is given,
#' all the taxa included in the \code{taxaExp} will be matched and selected ("all" by default).
#' @param method specifying which distance method to be used
#' to estimate expression phylogeny in bootstrapping.
#'
#' @return returns a data frame of six columns, refer to details for more information.
#'
#' @examples
#' data(tetraExp)
#' gamma.df <- estgamma(tetraExp, taxa = "all",
#'                      subtaxa = c("Brain", "Cerebellum"),
#'                      method = "pea")
#' gamma.df
#'
#' @references
#' Gu,X. 2016. Understanding tissue expression evolution: from expression phylogeny to phylogenetic network.
#' Brief. Bioinformatics, 17, 249-254.
#' @export
estgamma = function (objects = NULL, taxa = NULL, subtaxa = NULL,
                    method = c("sou", "ced", "pea", "souln", "nbdln", "euc", "cos", "jsd"))
{
  #if(verbose) message(date())

  if (is.null(objects) || !is(objects) == "taxaExp") {
    stop(paste0(date(), ": no valid taxaexp objects input"))
  }

  taxa_n_in_objects <- length(attr(objects,"taxa"))
  subtaxa_n_in_objects <- length(attr(objects,"subtaxa"))

  if (any(taxa_n_in_objects < 2, subtaxa_n_in_objects < 2)) {
    stop(paste0(date(), ": TE objects need contain at least two taxa and two subtaxa to calculate gammas"))
  }


  if (any(grepl("all", taxa, ignore.case = TRUE))) {
    taxa <- unique(unlist(lapply(objects, function(x) x$taxon_name)))
  }

  if (any(grepl("all", subtaxa, ignore.case = TRUE))) {
    subtaxa <- unique(unlist(lapply(objects, function(x) x$subTaxon_name)))
  }

  disMat <- expdist(objects = objects, taxa = taxa, subtaxa = subtaxa, method = method)

  #browser()

  taxon_names_arr <- NULL
  subtaxon_names_arr <- NULL

  gammaD_arr <- NULL
  gammaE_arr <- NULL
  ibA_arr <- NULL # internal branch
  ibB_arr <- NULL
  ibC_arr <- NULL
  ibD_arr <- NULL

  # for all taxa subtaxa combinations

  for (i in 1:(length(taxa)-1)) {

    for (j in (i+1):length(taxa)) {

      for (k in 1:(length(subtaxa)-1)) {

        for (l in (k+1):length(subtaxa)) {

          #browser()

          four_names <- c(paste0(taxa[i], "_", subtaxa[k]), paste0(taxa[i], "_", subtaxa[l]),
                          paste0(taxa[j], "_", subtaxa[k]), paste0(taxa[j], "_", subtaxa[l]))

          left_arr <- c(disMat[four_names[2], four_names[1]], disMat[four_names[3], four_names[1]],
                          disMat[four_names[3], four_names[2]], disMat[four_names[4], four_names[1]],
                          disMat[four_names[4], four_names[2]], disMat[four_names[4], four_names[3]])

          taxon_names <- paste(sort(c(taxa[i],taxa[j])), collapse = "-")
          subtaxon_names <- paste(sort(c(subtaxa[k],subtaxa[l])), collapse = "-")

          taxon_names_arr <- c(taxon_names_arr, taxon_names)
          subtaxon_names_arr <- c(subtaxon_names_arr, subtaxon_names)

          gamma_arr <- .solveAnEquation(left_arr)

          gammaD_arr <- c(gammaD_arr, gamma_arr[1])
          gammaE_arr <- c(gammaE_arr, gamma_arr[2])
          ibA_arr <- c(ibA_arr, gamma_arr[3])
          ibB_arr <- c(ibB_arr, gamma_arr[4])
          ibC_arr <- c(ibC_arr, gamma_arr[5])
          ibD_arr <- c(ibD_arr, gamma_arr[6])


        }
      }
    }
  }

  gamma.df <- data.frame(species = taxon_names_arr, tissues = subtaxon_names_arr,
                        gammaD = gammaD_arr, gammaE = gammaE_arr,
                        inBranchA = ibA_arr, inBranchB = ibB_arr,
                        inBranchC = ibC_arr, inBranchD = ibD_arr)

  gamma.df

}

