---
title: "README"
author: "Jingwen Yang"
date: "2019/4/27"
output: md_document
---

# TreeExp

*TreeExp* is an *R* package developed to provides useful phylogenetic tools applicable to RNA-seq data. The package can be applied to comparative expression evolution analysis, which includes but not liminited to:

* pairwise expression distance estimation;
* relative rate test for transcriptome evolution;
* the strength of expression conservation estimation;
* ancestral transcriptome inference.

Statistical methods implemented in the package was based on Ornstein-Uhlenbeck (OU) model of transcriptome evolution which claims that expression changes are constrained by stabilizing selection.

*TreeExp* package is under active developing, current developing version 2.0 is available at <https://github.com/jingwyang/TreeExp>.

A convenient way to install package from github is through *devtools* package:

```{r, eval=FALSE}
install.packages('devtools')
devtools::install_github("jingwyang/TreeExp")
```

Please refer to [Tutorial](https://jingwyang.github.io/) for startup guide of *TreeExp*.
