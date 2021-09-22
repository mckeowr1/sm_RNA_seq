library(GenomicAlignments)
library(BSgenome.Dmelanogaster.UCSC.dm6)
#library(tidyverse) 

#Read in the BAM to a gRanges Object
obj <- import.single("SRR1187947/Mapped_Reads/SRR1187947_trimmed.fq.gz.mapped.md.bam", 
            qualthresh = 1, 
            exclude.duplicates = FALSE)

countBam("SRR1187947/Mapped_Reads/SRR1187947_trimmed.fq.gz.mapped.md.bam") %>% 
  .$records

#Look at Read Length Distribution 

hist(width(obj))


## Mean Coverage # BP per Tile ### - Function this later

#Tile the Genome - Binning the Reference into a GRanges Obj

bins = tileGenome(seqlengths(Dmelanogaster), tilewidth = 10, cut.last.tile.in.chrom = TRUE)
bins = bins[seqnames(bins) %in% paste0('chr', c('2L','2R','3L','3R','4','X'))] # this can be commented out if you want to keep all chromosomes/contigs.
bins = bins[width(bins) == 10] 

seqlevels(bins) = seqlevelsInUse(bins)

genome(bins) = 'dm6'



#Calculate Coverage 

cover = coverage(obj)



#Calculate the Binned Average of Coverage - Genomic Ranges Object
ba = binnedAverage(bins = bins, numvar = cover, varname = 'score')

#Lets Plot Some Coverage


#Export a Wiggle File 


wiggler = function(coverage.data, file, name = "short name", description = "description", color = 'forestgreen', smoothing = 2, y.axis.range = c(0,20), height = 50)
{
  # this function will take coverage data, like that produced by `basecount()` and turn it into a viewable track for display on the UCSC genome browser. Necessary input is the coverage data (make sure that all widths are equal), and the output filepath. Tracklines will be generated with default values unless they are overwritten. 
  require(rtracklayer)
  
  trackline = new('GraphTrackLine',
                  name = name,
                  description = description,
                  alwaysZero = TRUE,
                  type = 'wig',
                  color = col2rgb(color)[,1],
                  smoothingWindow = smoothing,
                  windowingFunction = 'mean',
                  viewLimits = y.axis.range,
                  graphType = 'bar',
                  autoScale = FALSE,
                  visibility = 'full',
                  maxHeightPixels = as(c(10, height, 100), 'integer')	
  )
  
  write(x = as(trackline, "character"), file = file)
  
  export(object = coverage.data, con = file, format = 'wig', append = TRUE)
}

wiggler(ba, "/Users/ryan/Documents/GitHub/sm_RNA_seq/SRR1187947.wig")



# values(ba) %>%#Acess the Score Column of the Coverage 
#   as.data.frame() %>% 
#   rownames_to_column() %>%  
#   #sample_n(10000) %>% 
#   ggplot() + 
#   geom_line(aes(x = rowname, y =score))









#How to display reads that map to multiple locations (non-unqie mappers randomly assigned)

#How does wiggle look if map quality is 1 or map quality is 0  when Importing 
  # 0 should not change 
  # How are reads flowing through bowtie


#What is the mapability of any given sliding window in the drosohpial genom,e 







