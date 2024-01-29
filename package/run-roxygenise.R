#!/usr/bin/env Rscript
setwd("PhyloSim")
system("rm man/*")
roxygen2::roxygenise()
system("find -type f \\( -not -name \"MD5\" \\) -exec md5sum '{}' \\; > MD5")
setwd("..")
