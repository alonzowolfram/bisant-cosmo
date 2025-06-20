#!/usr/bin/env Rscript
# https://training.nextflow.io/advanced/structure/#bin

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Presto
devtools::install_github("immunogenomics/presto") # Faster marker identification.

# InSituType
devtools::install_github("https://github.com/Nanostring-Biostats/InSituType") # NanoString's CosMx cell-typing package.

# InSituCor
devtools::install_github("https://github.com/Nanostring-Biostats/InSituCor")
# ERROR: dependency ‘Rfast’ is not available for package ‘InSituCor’
# * removing ‘/rsrch6/home/genomic_med/lwfong/miniconda3/envs/cosmx/lib/R/library/InSituCor’
# Warning messages:
#   1: In i.p(...) :
#   installation of package ‘RcppGSL’ had non-zero exit status
# 2: In i.p(...) :
#   installation of package ‘RcppZiggurat’ had non-zero exit status
# 3: In i.p(...) : installation of package ‘Rfast’ had non-zero exit status
# 4: In i.p(...) :
#   installation of package ‘/tmp/RtmpdRkegl/file348731ba295d0/InSituCor_0.0.0.9000.tar.gz’ had non-zero exit status

# smiDE and dependencies
install.packages("spaMM") # Required by smiDE. Must be installed manually because not avaiable from conda.
library(devtools)
install_github("lhe17/nebula") # Required by smiDE.
install.packages("/rsrch6/home/genomic_med/lwfong/PRIME-TR/common/R-packages/INLA/INLA_23.05.30-1.tar.gz",
                 repos = NULL,
                 type = "source",
                 INSTALL_opts = c("--no-multiarch", "--clean", "--no-byte-compile")) # See note in "Version Compatibility" for which version to use.
install.packages("~/PRIME-TR/common/R-packages/smiDE/smiDE_0.0.2.00.tar.gz", repos = NULL, type = "source")

# Additional packages for misgdbr
install.packages('msigdbdf', repos = c('https://igordot.r-universe.dev', 'https://cloud.r-project.org'))

# spatialTIME
BiocManager::install("spatialTIME") # Analysis of spatial single-cell data.

# CellChat
devtools::install_github("jinworks/CellChat")

# DoRothEA
BiocManager::install("dorothea")
BiocManager::install("decoupleR")
BiocManager::install("OmnipathR")

# NicheNet
devtools::install_github("saeyslab/nichenetr")

# AUCell
BiocManager::install("AUCell")
# Other packages recommended for AUCell functionality:
# To support parallel execution:
BiocManager::install(c("doMC", "doSNOW"))
# For the examples in the follow-up section of the tutorial:
BiocManager::install(c("NMF", "dynamicTreeCut", "R2HTML")) # "d3heatmap" not available.

# BiocManager::install("SingleR")
# BiocManager::install("HGNChelper") # For ScType.
# BiocManager::install("scCustomize")
# install.packages("SCINA")
# BiocManager::install("progeny") # Pathway analysis of transcriptional data.
# install.packages("OpenRepGrid") # Vector recycling.
# BiocManager::install("sigmoid")
# BiocManager::install("randomcoloR")