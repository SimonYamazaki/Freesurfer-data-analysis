
#load packages used in all scripts
library(data.table)
library(ggplot2)
library(gridExtra)
library(lsmeans)
library(grid)
library(writexl)
library(car)
library(lsr)
library(ggseg3d)
library(rlist)
library(ggdist)
library(tidyr)
library(ggnewscale)
library(NCmisc)
library(readxl)
library(plotly)
library(htmlwidgets)
library(stringr)
library(filesstrings)
library(BayesFactor)

packinfo <- installed.packages(fields = c("Package", "Version"))
packinfo <- packinfo[,  c("Package","Version")]
packinfo <- packinfo[,  "Version", drop=F]

sink("/mrhome/simonyj/Freesurfer-data-analysis/R_package_versions.txt",append=TRUE)


"######### R SESSION INFO ##########"


sessionInfo()


"######### ALL PACKAGE VERSIONS INSTALLED ##########"


packinfo


sink()

