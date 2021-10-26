library(GenomicAlignments)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(bedr)
dm6 = Dmelanogaster


setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster")



generate_binned_readnames <- function(bin_size = 100, roi_readnames_bed , roi_chr, roi_start, roi_stop, proj_dir){
  ## Set up the Reference Genome ## 
  bins = tileGenome(seqlengths(Dmelanogaster), tilewidth = bin_size, cut.last.tile.in.chrom = TRUE)
  bins = bins[seqnames(bins) %in% paste0('chr', c('2L','2R','3L','3R','4','X'))] # this can be commented out if you want to keep all chromosomes/contigs.
  bins = bins[width(bins) == bin_size] 
  
  roi=GRanges(seqnames =c(roi_chr), 
                ranges=IRanges(start =c(roi_start), end=c(roi_stop))) 
  
  ## Subset the Region of Interest from the tiled Genome 
  
  roi_t <- subsetByOverlaps(bins, roi)
            
  #Give the Bins an id 
  roi_t$binID = 1:length(roi_t)
  
  ## Get the readnames in each bin
  
  #Read in Region BED 
  roi_bed <- bed_to_granges(roi_readnames_bed)
  
  #If a read is mapped to >1 bin, assign it to the first bin it occurs in 
  
  names(roi_bed) <- mcols(roi_bed)[,"id"] #set the names feature of the Granges obj to the readname 
  
  roi_bed_clean = roi_bed[unique(mcols(roi_bed)[,1])] #Run unique on the Readnames data of the gRanges obj
  
  overlaps <- mergeByOverlaps(roi_t, roi_bed_clean) #Join the clean readnames gRanges obj to the tiled ROI gRanges obj
   
  library(tidyverse)
  
  df <- overlaps[c(2, 4)] %>% as.data.frame() #Pull the bin ID and the readnames columns from the gRanges obj and convert to a DF
  
  #This list to character conversion is probably not needed but a later function expects the input to be in this format
  #Summarize the dataframe by binID
  reads_in_bin <- df %>% group_by(binID) %>% #Unlist the readnames 
    summarize(reads = paste0(as.character(id), collapse = ",")) #Put a summary column of read
  
  setwd(proj_dir)
  dir.create("binned.readnames")
  setwd("binned.readnames")
  
  write_readnames_binned <- function(reads_in_bin) { 
    bin_name <- reads_in_bin[1]
    reads <- as.character(reads_in_bin[,2]) %>%  
      str_split(pattern = ",") %>% 
      .[[1]] %>%
      sort()
    
    writeLines(reads, glue::glue("bin{bin_name}_readnames.txt"))
    
  } 
  
  for(read in 1:nrow(reads_in_bin)){ 
    # print(reads_in_bin[read, 1:2])
    df <- reads_in_bin[read,]
    
    write_readnames_binned(df)
    
    
  }
  
}



generate_binned_readnames(roi_readnames_bed = "/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/histonecluster_reads.bed",
                          roi_chr = "chr2L", 
                          roi_start = 21415940, 
                          roi_stop = 21543673, 
                          proj_dir = "/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster")






#### Setup the Reference Genome & The ROI Genome ####

bins = tileGenome(seqlengths(Dmelanogaster), tilewidth = 100, cut.last.tile.in.chrom = TRUE)
bins = bins[seqnames(bins) %in% paste0('chr', c('2L','2R','3L','3R','4','X'))] # this can be commented out if you want to keep all chromosomes/contigs.
bins = bins[width(bins) == 100] 


# Subset the Region of Interest From Tiled Genome # 

forty=GRanges(seqnames =c("chr2R"), 
              ranges=IRanges(start =c(6255432), end=c(6499291)))

#Get the tiled region of the genome
roi_t <- subsetByOverlaps(bins, forty)

#Name the bins 
roi_t$binID = 1:length(roi_t)

#### Assign Reads to BIN ####

#Read in Region BED 
roi_bed <- bed_to_granges("bed_roi/42AB.bed")

#sm_bed <- head(roi_bed, n= 1000)



## Deal with reads mapped to Region >1 ##

#Set the names to be Seq ID
names(roi_bed) <- mcols(roi_bed)[,"id"]

#Remove Duplicate Region Mappings
roi_bed_clean = roi_bed[unique(mcols(roi_bed)[,1])]


## Intersect the ROI_read_bed with the Tiled Genome Get a tiled ROI## 
overlaps <- mergeByOverlaps(roi_t, roi_bed_clean)

#### Get a Read Names File for each Bin #### 

library(tidyverse)


#
df <- overlaps[c(2, 4)] %>% as.data.frame()

reads_in_bin <- df %>% group_by(binID) %>% #Unlist the readnames 
  summarize(reads = paste0(as.character(id), collapse = ",")) 


#Write this out to a file (name = bin_ID, Text = Read names)

setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/binned_readnames/")
write_readnames_binned <- function(reads_in_bin) { 
  bin_name <- reads_in_bin[1]
  reads <- as.character(reads_in_bin[,2]) %>%  
    str_split(pattern = ",") %>% 
    .[[1]] %>%
    sort()
  
  writeLines(reads, glue::glue("bin{bin_name}_readnames.txt"))
  
} 

for(read in 1:nrow(reads_in_bin)){ 
  # print(reads_in_bin[read, 1:2])
  df <- reads_in_bin[read,]
  
  write_readnames_binned(df)

  
  }


# test_file <- reads_in_bin[10,]
# write_readnames_binned(test_file)







