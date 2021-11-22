library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(GenomicAlignments)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(tidyverse)

#Load Reference gene_set
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

genes <- genes(txdb) #Get the genes as a gRanges obj


histone <- data.table::fread("histone_cluster_meltedmatrix.tsv")






#Load ROI Beds into a list 

setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/histone_cluster_beds/")
bin_beds <- list.files()
beds <- list()
#Need to incorproate the strandedness!!
for(bed in bin_beds){
  df <- data.table::fread(bed)
  
  gr <- makeGRangesFromDataFrame(df, 
                                 seqnames.field = "V1", 
                                 start.field = "V2", 
                                 end.field = "V3")
  
  #Assign the list ID number to the bed name
  beds[[bed]] <- gr
  
}


#Load up the Genomic Bin DF

setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/analysis")

#Load Binned Genome DF
bins_df <- data.table::fread("binned_dm6_100000bins.tsv")


genome <- makeGRangesFromDataFrame(bins_df, 
                                   seqnames.field = "seqnames", 
                                   start.field = "start", 
                                   end.field = "end", 
                                   strand.field = "strand")


#Load Long Histone Matrix
hist_matrix <- data.table::fread("histone_cluster_meltedmatrix.tsv") 


#Test Gene List




## Start for gene loop here ##

#make a granges for a high mapping histone cluster region - bin260
high=GRanges(seqnames =c("chr2L"), 
             ranges=IRanges(start =c(10900001), end=c(11000000)))


high_genes<- subsetByOverlaps(genes, high)

#Set up vector for outputs
no_overlapping <- c()
gene_count <- 0
for(gene in 1:length(genes)){ 
  #print(head(high_genes[gene]))

gene_count <- gene_count + 1
print(gene_count )
##Bin the Genes##
tiled_gene <- tile(genes[gene], width = 100)

t_df <- as.data.frame(tiled_gene) %>% dplyr::select(seqnames, start, end, strand) #Get rid of gRanges group obj

t_gr <- makeGRangesFromDataFrame(t_df, 
                               seqnames.field = "seqnames", 
                               start.field = "start", 
                               end.field = "end", 
                               strand.field = "strand")


#Assign the gene bins a Genomic bin (dm6 with 100 kb bins) this gets us the bed files that we need to pull reads from
counts <- countOverlaps(genome, genes[gene])

genome_bins <- which(counts >0)


gene_bins <- hist_matrix %>%  dplyr::filter(GenomeBin %in% genome_bins & Log10_ReadCount > 0) %>%  
  .$HistoneClusterBin

if(length(gene_bins) == 0){ 
  append(no_overlapping, names(tiled_gene))
} else{

#Filter the list of beds to the gene specific beds
gene_beds <- beds[names(beds) %in% gene_bins] 


#Intersect the histone cluster bed with the gene bins
intersect_gene <- function(gr){
  c <- countOverlaps(t_gr, gr)
  
  return(c)
  
}



my_files = lapply(gene_beds, intersect_gene)


matrix = do.call("cbind", my_files)


long <- melt(matrix) %>% 
  rename("Gene_Bin" = "Var1", "ROI_bed" = "Var2", "num_reads" = "value")


#Now check to see if any portion of the gene is covered by a ROI mapping Read

check<- max(long$num_reads)

setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/analysis/test _gene")
if(check == 0) { 
  append(no_overlapping, names(tiled_gene))
} else {
  out_df <- long %>%
    group_by(Gene_Bin) %>%  
    summarise(reads = num_reads) %>%  
    filter(reads > 0)
  
  average <- mean(out_df$reads) %>% round()
  
  write_tsv(long, glue::glue("{names(tiled_gene)}__{average}.tsv")) 
  
}



}
} 

test_out <- long %>% 
  group_by(Gene_Bin) %>%  
  summarise(reads = num_reads)



head(high_genes[1])

length(high_genes)


test_out <- data.table::fread(input = "FBgn0261871.tsv")

sum <- test_out %>% group_by(Gene_Bin) %>%  
  summarise(reads = num_reads) %>%  
  filter(reads > 0)

avg <- mean(sum$reads) %>% round()

#Get Some Test Genes that are in a high mapping regions






#Overlap Gene with histone mapping reads





# Get the intergenic regions