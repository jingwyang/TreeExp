#' @title No negative or zero branch length
#'
#' @name no0br
#' @rdname no0br
#'
#' @description This function does a small tweak on the tree to
#' replace negative-length branch with zero-length and/or assign
#' zero-length branch between root and MRCA of the ingroup a small value
#'
#' @param phy an tree of 'phylo' class
#' @param remove.0length a logical specifying whether to replace
#' negtative branch length with 0 (TRUE by default)
#' @param root.length a logical specifying whether to assign a
#' small value to the branch between root and MRCA of
#' the ingroup (FALSE by default). If TRUE, the tree should be rooted first
#' @return returns the tree without zero-length or negative-length branch
#'
#' @export
no0br = function(phy, remove.0length = TRUE, root.length= FALSE) {

    if (!inherits(phy, "phylo"))
        stop(paste0(date(),"tree input is not of class \"phylo\""))

    if (is.null(phy$edge.length))
        stop(paste0(date(),": tree has no branch lengths"))

    tmp <- sort(phy$edge.length)
    tmp1 <- tmp[tmp>0]

    if (remove.0length){
        if (any(tmp < 0)) {
            warning(paste0(date(),": there are negative-length branches in
                       the tree, replacing them with a small length"))
            phy$edge.length[phy$edge.length < 0] <- tmp1[1] / 1e3
          }
          else
            cat(paste0(date(),": the tree has no zero-length branch, nothing
                   to do"))
        }

    if (root.length){
      phy$edge.length[1] <- tmp1[1]
      phy$edge.length[length(phy$edge.length)] <-
      phy$edge.length[length(phy$edge.length)] - phy$edge.length[1]
    }
    phy

}

#' @title Generate an inversed variance matrix from expression profile
#' across species
#'
#' @name varMatInv
#' @rdname varMatInv
#'
#' @description This function generate an inversed variance matrix
#' from expression profiles of one-to-one orthologous genes across species
#'
#' @param objects a vector of objects of class \code{taxonExp} or
#' an object of class \code{taxaExp}
#' @param phy an expression character tree
#' @param taxa one single character or a vector of characters sepcifying
#' taxa to generate an inversed variance matrix. If one single character
#' "all" is given, all the taxa included in the \code{taxaExp}
#' will be matched and included ("all" by default).
#' @param subtaxa one single character specifying sub taxa to be
#' included in generating an inversed variance matrix
#'
#' @return returns an inversed variance matrix
#'
#' @export
varMatInv = function(objects , phy, taxa = "all", subtaxa) {

    if (!inherits(phy, "phylo"))
        stop(paste0(date(),"tree input is not of class \"phylo\""))

    if (is.null(phy$edge.length))
        stop(paste0(date(),": tree has no branch lengths which is a
                    necessity for \"varMatInv\""))

    if (length(subtaxa) > 1 || subtaxa == "all")
        stop(paste0(date(),": only one subtaxon are allowed here"))

    ### using -ln(rho) to estimate pairwise expression distance
    dismat <- expdist(objects, taxa = taxa, subtaxa = subtaxa, method = "sou")

    if (!all(row.names(dismat) %in% phy$tip.label ))
        stop(paste0(date(),": taxa or subtaxa names do not match
                    perfectly with tree tip labels, please check them."))

    n_tip <- Ntip(phy)
    n_node <- Nnode(phy)

    ### extract distances from the tree
    nodes_dist <- dist.nodes(phy)
    corrmat <- apply(nodes_dist, c(1,2), function(x) exp(-x))

    ### stationary variance
    exp_table <- exptabTE(objects, taxa = "all", subtaxa = subtaxa,
                          logrithm = TRUE)
    stat_var <- mean(apply(exp_table, 2, var))
    var_corrmat <- stat_var * corrmat

    solve(var_corrmat)

}

#' @title Ancestral Expression Estimation
#'
#' @name aee
#' @rdname aee
#'
#' @description This function esitmates ancestral expression profile and
#' related statistical uncertainty
#'
#' @param x a vector of known expression profile, preferably log-transformed
#' expression levels (e.g. log RPKM)
#' @param phy an unrooted phylogenetic tree in the form of class "phylo"
#' @param mat a matrix generated from "varMatInv" function
#' @param CI a logical specifying whether to return the 95% confidence
#' intervals of the estimated ancestral expression levels
#'
#' @return returns a list containing estimated ancestral expression profile
#' as well as other requested parameters
#'
#' @examples
#' library('ape')
#' data('tetraExp')
#' dismat <- expdist(tetraExp, taxa = "all", subtaxa = "Brain", method = "sou")
#' exp_tree <- NJ(dismat)
#' exp_tree <- no0br(exp_tree)
#' var_mat <- varMatInv(objects = tetraExp,phy = exp_tree,taxa = "all",
#' subtaxa = "Brain")
#' exp_table <- exptabTE(tetraExp, taxa = "all", subtaxa = "Brain")
#' exp_one <- aee(exp_table[1,], exp_tree, var_mat)
#' exp_tree$node.label <- exp_one$est
#' exp_tree <- root(exp_tree, outgroup = "Chicken_Brain",
#' resolve.root = TRUE)
#' plot(exp_tree,show.node.label = TRUE)
#'
#' @export
aee = function(x, phy, mat, CI = TRUE) {

    ### checking input formats
    if (!inherits(phy,"phylo"))
        stop(paste0(date(),": \"phy\" input is not of class \"phylo\""))

    if (is.null(phy$edge.length))
        stop(paste0(date(),": tree has no branch lengths which is a
                    necessity for \"aee\""))

    if (!is.null(names(x))) {
        if (all(names(x) %in% phy$tip.label))
            x <- x[phy$tip.label]
        else warning(paste0(date(),
            "characters do not match perfectly between expression profile
            vector names and tree tip labels,
            only using tree tip labels in the following analysis"))
    }

    ### checking if the tree is rooted
    n_tip <- Ntip(phy)
    n_node <- Nnode(phy)

#    if (n_tip != n_node + 1)
#        stop(paste0(date(),"tree is not rooted, please make sure tree is properly rooted. "))

    ancestral <- list()

    expr <- numeric(length = n_tip + n_node) ### expression values vector initiate
    ### given expression values of tips
    expr[1:n_tip] <- if (is.null(names(x))) x else as.numeric(x[phy$tip.label])
    ### ancestral nodes expression values to be estimated
    expr[(n_tip+1):(n_tip+n_node)] <- NA

    tr_edges <- phy$edge

    all_tips <- 1:n_tip

    for (i in (n_tip+1):(n_tip+n_node)) {

        mu <- mean(expr[all_tips])

        beta <- unlist(lapply(all_tips, function(x) - mat[x,i] / mat[i,i]))
        beta0 <- mu * (1 - sum(beta))

        expr[i] <- beta0 + sum(beta * expr[all_tips])

    }


    ancestral$est <- expr[(n_tip+1):(n_tip+n_node)]
    ### if calculating confidence interval

    if (CI) {

        ci95 <- matrix(nrow = n_node, ncol = 2)

        for (i in (n_tip+1):(n_tip+n_node)) { ### for every node

        tmp = sqrt(1 / mat[i,i]) * qnorm(0.025)
        ci95[(i-n_tip),] = c(expr[i] + tmp, expr[i] - tmp)

        }

        ancestral$ci95 <- ci95
    }

    ancestral
}
