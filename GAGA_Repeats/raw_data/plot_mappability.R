library(BSgenome)


setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/GAGA_Repeats/raw_data/")

import.single_repeats = function(file, qualthresh = 10, exclude.duplicates = FALSE)
{
  # this function will import paired-end data that has been mapped and duplicate-marked. It takes as input a path to a .bamfile. It gives the option to set the map-quality score cutoff as well as whether or not to exclude duplicates. This is a more generalized and versatile version of the importer for ATAC seq reads. It does not output a graph of fragment size distributions or coverage per chromosome.
  
  require(GenomicAlignments)
  require(BSgenome.Dmelanogaster.UCSC.dm6)
  dm6 = Dmelanogaster
  
  shortname = substring(file,tail(gregexpr('[/]',file)[[1]],1)+1,nchar(file))
  parent.dir = dirname(dirname(file))
  #	if(!dir.exists(paste0(parent.dir,'/sample_properties'))){dir.create(paste0(parent.dir,'/sample_properties')); cat('made a new directory (sample_properties) in the parental directory.\n')}
  
  cat(paste('Running import.single() on',file,'\n\n'))
  cat('Only importing from the canonical chromosomes, Drosophila genome version dm6. Reads only with above threshold map quality, non secondary mappings.\n\n')
  
  good.uns = c('chrX','chr2L','chr2R','chr3L','chr3R','chr4')
  
  cat('importing all reads and bam file metadata\n')
  
  if(exclude.duplicates){
    cat("Duplicates will be excluded. \n")
    allreads = readGAlignments(file, 
                               param = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE), 
                                                    what = c('flag','mrnm','mpos','mapq','isize')),
                               use.names = FALSE)
  }else{cat("Duplicates will not be excluded.\n")
    allreads = readGAlignments(file, 
                               param = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE), 
                                                    what = c('flag','mrnm','mpos','mapq','isize')),
                               use.names = FALSE)
  }
  
  cat('filtering out reads on non-canonical chromosomes.\n')
  init = length(allreads)
  allreads = allreads[seqnames(allreads) %in% good.uns]
  after = length(allreads)
  cat(paste('filtered', init-after, 'reads on non-canonical chromosomes.\n\n'))
  
  cat('filtering all reads by map quality score.\n')
  init = length(allreads)
  allreads = allreads[mcols(allreads)$mapq >= qualthresh]
  after = length(allreads)
  cat(paste('filtered', init-after, 'reads with low map quality.\n\n'))
  
  #cat('reconstructing paired end reads from filtered reads.\n')
  #pairs = makeGAlignmentPairs(allreads,use.names = TRUE,use.mcols = c('flag','mapq'))
  
  cat('coercing the data into a granges object.\n')
  single.interval = granges(allreads)
  
  cat('Assigning <genome> variable as <dm6>.\n')
  genome(single.interval) = 'dm6'
  
  seqlevels(single.interval) = seqlevels(dm6)
  seqinfo(single.interval) = seqinfo(dm6)
  single.interval = trim(single.interval)
  
  cat('Removing seqlevels not in use.\n')
  seqlevels(single.interval) = seqlevelsInUse(single.interval)
  
  names(single.interval) = NULL
  return(single.interval)
}



gaga_dm6 <- import.single_repeats("gaga_verysensitive_local.mapped.bam", qualthresh = 0, exclude.duplicates = FALSE)


#Load Drosophila Genome
bins = tileGenome(seqlengths(Dmelanogaster), tilewidth = 10000, cut.last.tile.in.chrom = TRUE)
bins = bins[seqnames(bins) %in% paste0('chr', c('2L','2R','3L','3R','4','X'))] # this can be commented out if you want to keep all chromosomes/contigs.
bins = bins[width(bins) == 10000] 


#Create a Matrix of the reads by genomic position
matrix <- countOverlaps(bins, gaga_dm6)

overlaps <- mergeByOverlaps(bins, gaga_dm6)


#Make the matrix Long

library(tidyverse)
library(reshape2)
library(ggplot2)

longdata <- melt(matrix) %>%  
  rownames_to_column() %>%  
  mutate(value = log(value, 10)) #%>%  
  #mutate(bin = as.numeric(rowname))


pdf("test_plot.pdf")
ggplot(longdata, aes(x =  rowname)) + 
  geom_path()
dev.off()


ggplot(longdata, aes(x =  rowname, y = value)) + 
  geom_point()

gaga_bam <- readGAlignments("gaga_verysensitive_local.mapped.bam")

















## Stuff With Heterochromatin Mappers ##



#Tile the Heterochromatin Reference
chrlens <- Biostrings::fasta.seqlengths("/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/chang2019/dmel_scaffold2_V5.fasta")
genomeobj <- Biostrings::readDNAStringSet("/Users/ryan/Documents/GitHub/sm_RNA_seq/sat_dna/chang2019/dmel_scaffold2_V5.fasta")
tileWidth = 10
gr <- GenomicRanges::tileGenome(seqlengths = chrlens, tilewidth = tileWidth,
                            cut.last.tile.in.chrom = TRUE)
gr <- gr[GenomicRanges::width(gr) == tileWidth]



#Load up Gaga 50mer mapped to heterochromatin reference
gaga_hetchrom<- readGAlignments("gaga_verysensitive_local_heterochrom.mapped.sam")



#Calculate Coverage
cover = coverage(gaga_hetchrom)

ba = binnedAverage(bins = gr, numvar = cover, varname = 'score') 
score(ba) = signif(score(ba) / length(gr) * 1e+6, 4) # this assumes that the length of the input genomic ranges object is the relevant denominator to do your normalization with.

mapped_regions <- as.data.table(ba) %>%  
  filter(score > 0 )
  #filter(seqnames == "3_scaffold2")

#What are the annotations for those regions












