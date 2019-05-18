
#' @title Construct a taxaExp object.
#'
#' @description  Constructor function for \code{taxaExp} objects.
#' This function takes in a expression value file from RNA-seq,
#' and construct a \code{taxaExp} object from which user can extract information
#' for display or for further analysis.
#'
#' @name TEconstruct
#'
#' @rdname TEconstruct
#'
#' @param ExpValueFP a text file contains processed expression data (RPKM or TPM for example).
#' Row names correspond with gene names,
#' and column names correspond with taxon and subtaxon names.
#' @param taxa one single string or a vector of strings specifying main taxa selected for
#' constructing \code{taxaExp} object.
#' Taxa names are extracted from row names given in gene length file.
#' If one single string "all" is given,
#' all the taxa in the row names will be matched and selected ("all" by default).
#' @param subtaxa one single string or a vector of strings sepcifying sub taxa selected for
#' constructing \code{taxaExp} object.
#' If one single string "all" is given,
#' all the sub taxa in the row names will be matched and selected ("all" by default).
#' @param rmOut a logical sepcifying whether to remove expression outliers
#' while constructing \code{taxaExp} objects (TRUE by default).
#' @param verbose a logical specifying whether to print more information on the screen
#' while constructing \code{taxaExp} objects (FALSE by default).
#'
#' @return returns an object of class \code{Taxa} (S3 class, a list of \code{taxonExp} objects).
#'
#' @examples
#'
#' taxa.objects = TEconstruct(ExpValueFP = system.file('extdata/primate_brain_expvalues.txt',
#'  package = 'TreeExp'), taxa = "all", subtaxa = c("HIP", "CB"))
#'
#' @export
TEconstruct = function(ExpValueFP=NULL, taxa="all", subtaxa="all",
                       rmOut=FALSE, verbose=FALSE) {

  # check file handle
  if((is.null(ExpValueFP))){
    stop(paste0(date(),": must provide expression values file path and gene length file path"))
  }

  # check file existance
  if(!file.exists(ExpValueFP)){
    stop(paste0(date(),": fail to open file, check your filename or path"))
  }
  #browser()

  # input
  exp_value_df <- read.table(ExpValueFP, header = TRUE)
  row.names(exp_value_df) <- exp_value_df[,1]
  exp_value_df <- exp_value_df[,-1]

# gene.info.df <- read.table(geneInfoFP,header=T)
# remove sample with low read counts

  invalid_arr <- NULL

  for (i in 2:ncol(exp_value_df)) {
    if (mean(exp_value_df[,i]) < 1) {
      invalid_arr <- c(invalid_arr,i)
    }
  }

  message(paste0(date(),": removing ", length(invalid_arr), " sample(s) with ultra-low read counts"))

  if (length(invalid_arr) != 0) {
    invalid_arr = 0 - invalid_arr
    exp_value_df <- exp_value_df[,invalid_arr]
  }

  # gene number and taxon number
  gene_n <- nrow(exp_value_df)

  # get taxon names from read counts file
  taxon_names <- unique(lapply(colnames(exp_value_df), function(x) unlist(strsplit(x, "_"))[1]))
  taxon_n <- length(taxon_names)

  # normalize<-match.arg(normalize)
  # message(paste0(date(), ": using ", normalize, " to normalize raw read counts"))

  # get taxon names
  #browser()
  cat("\n")
  message(paste0(date(),": start constructiong TE objects"))


  if (!any(grepl("all", taxa, ignore.case = TRUE))) {

    taxon_names <- gsub("\\s+", "", taxa)
    taxon_n <- length(taxon_names)

  }

  message(paste0(date(),": total Taxon number ", taxon_n))

  #browser()

  # get subtaxon number
  subtaxon_names <- unique(lapply(colnames(exp_value_df), function(x) unlist(strsplit(x, "_"))[2]))
  subtaxon_n <- length(subtaxon_names)


  if (!any(grepl("all", subtaxa, ignore.case = TRUE))) {

    subtaxon_names <- gsub("\\s+", "", subtaxa)
    subtaxon_n <- length(subtaxon_names)

  }

  message(paste0(date(),": total sub taxon number ", subtaxon_n))

  title <- lapply(colnames(exp_value_df), function(x) unlist(strsplit(x, "_"))[1]) # taxon names
  subtitle <- lapply(colnames(exp_value_df), function(x) unlist(strsplit(x,"_"))[2]) # subtaxon names
  first_two_names <- unique(paste(title,subtitle,sep="_"))

  index <- intersect(unlist(lapply(taxon_names, function(x) grep(x, first_two_names, ignore.case = TRUE))),
          unlist(lapply(subtaxon_names,function(x) grep(x, first_two_names, ignore.case = TRUE))))

  objects_names <- first_two_names[index]

  objects_number <- length(objects_names)

  #browser()

  cat("\n")
  # get gene names
  #gene.names <- read.counts.df[,1]

  message(paste0(date(),": now constructing ",objects_number, " TE objects..."))

  if (!verbose) progbar <- txtProgressBar(style = 3)

  # initialization

  taxonExp.objects <- vector("list",length = objects_number)
  # the number of TE objects constructed is based on seleted numnber

  # for each taxon

  objects_counter <- 0

  for (i in 1:objects_number) {

    #browser()
    if (verbose) message(paste0(date(),": proceeding taxon ", objects_names[i]))

    # get all the sample names matching objects names
    # bundle all the biological replicates into one TE object

    #browser()

    ttl <- unlist(strsplit(objects_names[i], "_"))[1] #taxon title
    subttl <- unlist(strsplit(objects_names[i], "_"))[2] # subtaxon title
    #ttl <- lapply(names, function(x) unlist(strsplit(x, "_"))[1]) # taxon names
    #subttl <- lapply(names, function(x) unlist(strsplit(x,"_"))[2]) # subtaxon names

    idx <- grep(objects_names[i],colnames(exp_value_df), ignore.case = TRUE)
    names <- strsplit(colnames(exp_value_df)[idx],"_")
    repttl <- unlist(lapply(names, function(x) unlist(strsplit(x,"_"))[3])) # biological replicates title names
    gene_names <- rownames(exp_value_df)  # gene names


    # foreach subtaxon
    bio_rep_n <- length(repttl) # biological replicates number
#    omega <- NULL # omega estimated overdispersion parameter

    exp_val <- apply(exp_value_df[idx], c(1,2), as.numeric)

    objects_counter = objects_counter + 1

    if (verbose) message(paste0(date(),": wrapping up into objects"))

    #browser()
    oneObject <- list(exp_value=exp_val,
                      taxon_name = ttl,subTaxon_name = subttl,
                      gene_num = gene_n, gene_name = gene_names,
                      bioRep_num = bio_rep_n, bioRep_id = repttl)

    class(oneObject) <- "taxonExp"

    taxonExp.objects[[objects_counter]] <- oneObject

    #browser()

    if (verbose) message(paste0(date(), ": ", objects_counter, " TE objects constructed"))

    if (verbose) cat("\n")

    if (!verbose) setTxtProgressBar(progbar, objects_counter/objects_number)


  }

  class(taxonExp.objects) <- "taxaExp"

  attr(taxonExp.objects, "taxa") <- unlist(taxon_names)
  attr(taxonExp.objects, "subtaxa") <- unlist(subtaxon_names)

  cat("\n")
  message(date(),": construction complete.")

  taxonExp.objects

}
