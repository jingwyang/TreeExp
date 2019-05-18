## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE---------------------------------------------------------
#  install.packages('devtools')
#  devtools::install_github("jingwyang/TreeExp")

## ---- eval=FALSE---------------------------------------------------------
#  library('TreeExp')

## ---- eval=FALSE---------------------------------------------------------
#  ?RelaRate.test()

## ---- eval=FALSE---------------------------------------------------------
#  help(RelaRate.test)

## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github('jingwyang/TreeExp',build_opts = c("--no-resave-data", "--no-manual"))

## ---- eval=FALSE---------------------------------------------------------
#  browseVignettes('TreeExp')

## ----echo=FALSE, results='asis',message=FALSE----------------------------
library('TreeExp')
data('example_data')
knitr::kable(example_data)

## ---- eval=FALSE---------------------------------------------------------
#  taxa.objects = TEconstruct(ExpValueFP = system.file('extdata/primate_brain_expvalues.txt',
#  package = 'TreeExp'), taxa = "all", subtaxa = 'all')

## ----message=FALSE,results='hide'----------------------------------------
taxa.objects = TEconstruct(ExpValueFP = system.file('extdata/primate_brain_expvalues.txt',
package = 'TreeExp'), taxa = "all", subtaxa = c("ACC","CB"))

## ------------------------------------------------------------------------
print(taxa.objects, details = TRUE)

## ------------------------------------------------------------------------
print(taxa.objects[[1]], printlen = 6)

## ----eval=FALSE----------------------------------------------------------
#  taxa.objects[[6]]$exp_value[1:5,]

## ---- results='asis', echo=FALSE-----------------------------------------
knitr::kable(taxa.objects[[6]]$exp_value[1:5,])

## ---- echo=FALSE, fig.cap="Figure1. The schematic of transcriptome evolution along a phylogeny", out.width = '90%'----
knitr::include_graphics("Figure1.png")

## ---- warning=FALSE------------------------------------------------------
library(TreeExp)

## ------------------------------------------------------------------------
data("tetraExp")

## ---- message=FALSE------------------------------------------------------
dismat_pea <- expdist(tetraExp, taxa = "all",
                 subtaxa = "Brain",
                 method = "pea")
as.dist(dismat_pea)

## ---- message=FALSE------------------------------------------------------
dismat_sou <- expdist(tetraExp, taxa = "all",
                 subtaxa = "Brain",
                 method = "sou")
as.dist(dismat_sou)

## ---- message=FALSE------------------------------------------------------
dismat_sou_v <- expdist(tetraExp, taxa = "all",
                 subtaxa = "Brain",
                 method = "sou_v")
dismat_sou_v$pi
as.dist(dismat_sou_v$distance)

## ---- message=FALSE------------------------------------------------------
dismat_sou_v$pi

## ---- eval=FALSE---------------------------------------------------------
#  
#  expression_table <- exptabTE(tetraExp, taxa = "all",
#                              subtaxa = "Brain")
#  
#  dismat <- dist.pea(expression_table)
#  dismat <- dist.sou(expression_table)
#  dismat <- dist.sou_v(expression_table)

## ---- eval=FALSE---------------------------------------------------------
#  
#  dismat <- dist.pea(your_own_dataframe)
#  colnames(dismat) <- colnames(your_own_dataframe)
#  rownames(dismat) <- colnames(dismat)
#  

## ---- message=FALSE, warning=FALSE, results='hide', fig.height=4, fig.width=6----
library('ape')
 ##### build NJ tree
tr <- NJ(dismat_sou)
 ### estract expression table of Brain tissue from 'tetraExp'
expression_table <-exptabTE(objects = tetraExp,taxa = 'all',
                            subtaxa = 'Brain')
### root the tree using 'Chicken_Brain' as outgroup
tr <- root(tr, "Chicken_Brain", resolve.root = TRUE) 
###### generating bootstrap values of the tree
bs<-boot.exphy(phy=tr, x = expression_table, method = 'sou',
               outgroup = 'Chicken_Brain', B = 100) 
### assign bootstrap values to the "node.label" of the tree
tr$node.label = bs 
### plot the tree
plot(tr, show.node.label = TRUE)

## ---- warning=FALSE, message=FALSE---------------------------------------
#### load the tetraExp data firstly
data(tetraExp)
### extract the gene expression values of the selected genes from 'tetraExp' object.
exp_table <-exptabTE(tetraExp, taxa = 'all', subtaxa = 'Brain',rowindex = 200:800)

## ---- warning=FALSE, message=FALSE---------------------------------------
ztest <- RelaRate.test(expTable = exp_table, x = 'human', y = 'chimpanzee',
                       outgroup = 'macaque', alternative = 'greater')
ztest

## ----pressure, echo=FALSE, fig.cap="Figure2. The evolutionary scenario for comparative transcriptome analysis", out.width = '50%'----
knitr::include_graphics("Figure2.png")

## ---- eval=FALSE---------------------------------------------------------
#  library('TreeExp')

## ---- warning = FALSE, message = FALSE-----------------------------------
data('tetraExp')

## ---- warning = FALSE, message = FALSE-----------------------------------

species.group <- c("Human", "Chimpanzee", "Bonobo", "Gorilla",
"Macaque", "Mouse", "Opossum", "Platypus")
### all mammalian species

inv.corr.mat <- corrMatInv(tetraExp, taxa = species.group, subtaxa = "Brain")
inv.corr.mat

## ---- warning = FALSE, message = FALSE-----------------------------------
brain.exptable <- exptabTE(tetraExp, taxa = species.group, subtaxa = "Brain" ,logrithm = TRUE)
head(brain.exptable)

## ---- warning = FALSE, message = FALSE-----------------------------------
gamma.paras <- estParaGamma(exptable =brain.exptable, corrmatinv =inv.corr.mat)
## print the elements of gamma.paras
gamma.paras

## ---- warning = FALSE, message = FALSE-----------------------------------
brain.Q <- estParaQ(brain.exptable, corrmatinv = inv.corr.mat)
# with prior expression values and inversed correlation matrix
    
brain.post<- estParaWBayesian(brain.Q, gamma.paras)
brain.W <- brain.post$w # posterior selection pressures
brain.CI <- brain.post$ci95 # posterior expression 95% confidence interval


## ------------------------------------------------------------------------
names(brain.W) <- rownames(brain.exptable)
head(sort(brain.W, decreasing = TRUE)) #check a few genes with highest seletion pressure

## ----fig.height=4, fig.width=6-------------------------------------------
plot(density(brain.W))

## ---- echo=FALSE, fig.cap="Figure3. The schematic of ancestral expression inference along a phylogeny", out.width = '90%'----
knitr::include_graphics("Figure3.png")

## ---- warning=FALSE, message=FALSE---------------------------------------
data('tetraExp')
primate_group <-c('human','chimpanzee','gibbon','bonobo','gorilla','macaque')

dismat <- expdist(tetraExp, taxa = primate_group, 
                  subtaxa = "brain", method = "sou")
dismat

## ----warning=FALSE, message=FALSE----------------------------------------
primate_expT <-exptabTE(objects = tetraExp,taxa = primate_group,
                            subtaxa = 'Brain')
dismat <- dist.sou(primate_expT)
dismat

## ----warning=FALSE, message=FALSE----------------------------------------
primate_tree <- NJ(dismat)

## ----warning=FALSE, message=FALSE,results='hide', fig.height=4, fig.width=6----
primate_tree_root <- root(primate_tree, outgroup = "Macaque_Brain", resolve.root = TRUE)
bs<-boot.exphy(phy=primate_tree_root, x = primate_expT, method = 'sou', outgroup = 'Macaque_Brain',B = 100) 
primate_tree_root$node.label = bs 
plot(primate_tree_root, show.node.label = TRUE)

## ---- warning=FALSE, message=FALSE---------------------------------------
var_mat <- varMatInv(objects = tetraExp,phy = primate_tree,
                     taxa = primate_group, subtaxa = "Brain")
var_mat

## ---- warning=FALSE, message=FALSE---------------------------------------
    
primate_expT <-exptabTE(objects = tetraExp,taxa = primate_group,
                            subtaxa = 'Brain')

EMP1_expression <- primate_expT[which(rownames(primate_expT) == "ENSG00000134531"),]

## ---- warning=FALSE, message=FALSE---------------------------------------
EMP1_anc <- aee(x = EMP1_expression, phy = primate_tree, mat = var_mat)

## ---- warning=FALSE, message=FALSE, fig.height=4, fig.width=6------------
primate_tree$node.label <- sprintf("%.4f",EMP1_anc$est)
primate_tree$tip.label <- paste0(primate_tree$tip.label, "  ", 
                                 sprintf("%.4f", EMP1_expression))

plot(primate_tree, edge.color = "grey80", edge.width = 4, 
     show.node.label = TRUE, align.tip.label = TRUE)

