library(memes)
suppressPackageStartupMessages(library(GenomicRanges))
library(magrittr)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)



data("example_peaks", package = "memes")

genome <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6

genes <- genes(txdb)




