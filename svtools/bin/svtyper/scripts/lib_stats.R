#!/usr/bin/env Rscript

# print usage
usage <- function() {
    cat(
'usage: lib_stats.R <file>
lib_stats.R
author: Colby Chiang (colbychiang@wustl.edu)
description: Plot a library information for BAM file
positional arguments:
  <input JSON>               JSON file of library info, created by SVTyper
  <output PDF>               path to output pdf
')
}

# Draw a histogram from a text file
args <- commandArgs(trailingOnly=TRUE)
file <- args[1]
filename <- basename(args[1])
output <- args[2]

# Check input args
if (is.na(file) || is.na(output)) {
    usage()
    quit(save='no', status=1)
}

# install R packages if needed
options(repos=structure(c(CRAN="http://cran.wustl.edu")))
list.of.packages <- c("jsonlite")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(jsonlite)

bam <- fromJSON(file)

pdf(output, h=6, w=8)
for (sample in bam) {
    ins <- as.data.frame(t(sample$libraryArray$histogram))
    
    for (i in 1:ncol(ins)) {
        plot(as.numeric(row.names(ins)), ins[,i], type='h', xlab='Insert size (bp)', ylab='Number of read-pairs', main=paste0('LB: ', sample$libraryArray$library_name[i], '\nprevalence: ', round(100*sample$libraryArray$prevalence[i],1), "%"), col='steelblue', bty='l')
        abline(v=sample$libraryArray$mean[i], col='red', lty=1)
        legend('topright', c(paste0('mean: ', round(sample$libraryArray$mean[i],1)), paste0('sd: ', round(sample$libraryArray$sd[i],1))), lty=c(1, 0), col=c('red'))
    }
}
dev.off()

