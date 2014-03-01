FDRreg
======

R package for false discovery rate regression (FDRR)

Many approaches for multiple testing begin with the assumption that all tests in a given study should be combined into a global false-discovery-rate analysis. But this may be inappropriate for many of today's large-scale screening problems, where auxiliary information about each test is often available, and where a combined analysis can lead to poorly calibrated error rates within different subsets of the experiment.

This package implements false-discovery-rate regression (FDRR), in which auxiliary covariate information is used to improve power while maintaining control over the global error rate. The method can be motivated by a hierarchical Bayesian model in which covariates are allowed to influence the local false discovery rate (or equivalently, the posterior probability that a given observation is a signal) via a logistic regression. 

To install the package in R, first install the devtools package, and then use the commands
`````````
library(devtools)
install_github('jgscott/FDRreg')
`````````

The method is described in the paper 'False discovery rate regression: an application to neural synchrony detection in primary visual cortex', available as [arXiv:1307.3495 (stat.ME)](http://arxiv.org/abs/1307.3495).
