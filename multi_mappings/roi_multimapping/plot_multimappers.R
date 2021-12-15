library(GenomicAlignments)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(data.table)
library(cowplot)

bins = tileGenome(seqlengths(Dmelanogaster), tilewidth = 100000, cut.last.tile.in.chrom = TRUE)
bins = bins[seqnames(bins) %in% paste0('chr', c('2L','2R','3L','3R','4','X'))] # this can be commented out if you want to keep all chromosomes/contigs.
bins = bins[width(bins) == 100000] 

make_vector<- function(file){ 
  print(file)
  
  df <- data.table::fread(file)
  
  
  gr <- makeGRangesFromDataFrame(df, 
                                 seqnames.field = "V1", 
                                 start.field = "V2", 
                                 end.field = "V3")
  c <- countOverlaps(bins, gr)
  
  return(c)                             

}


setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/histone_cluster_beds/")
bin_beds <- list.files()

library(stringr)
ordered_beds = str_sort(bin_beds, numeric = T)


myfiles = lapply(ordered_beds, make_vector)

names(myfiles) <- ordered_beds

matrix = do.call("cbind", myfiles)




heatmap(matrix, Colv = NA, Rowv = NA, na.rm = FALSE)



#GGplot a Matrix

library(reshape2)
library(ggplot2)
library(dplyr)

longdata <- melt(matrix)

longdata %>%
  ggplot(aes(x = value)) +
  geom_histogram()


#Read in the file to not have to load up and run everythin

setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/analysis")
longdata <- data.table::fread("histone_cluster_meltedmatrix.tsv") %>%  
  mutate(binid = str_extract(.$HistoneClusterBin, "[:digit:]"))


longdata2 <- mutate(longdata,value = log(value + 1, 10) ) 

ggplot(longdata, aes(x = GenomeBin, y = HistoneClusterBin, fill = Log10_ReadCount)) +
  geom_raster() +
  scale_fill_distiller(palette = "RdPu") +
  #scale_fill_gradient(low = "grey", high = "blue")+
  geom_vline(xintercept = 215, alpha = 0.05 )+
  #geom_vline(xintercept = no_extlow, alpha = 0.15, color = "#219ebc")+
  geom_vline(xintercept = 345:347, alpha = 0.15, color = "#219ebc")+
  geom_hline(yintercept = "bin1139.bed", alpha = 0.15, color = "#219ebc")+
  theme(
    axis.text = element_text(size = 10),
    axis.text.y = element_blank()) + 
  labs(x = "Genomic Bin ID (Bin size = 100000 Bp)", 
       y = "Histone Cluster Bin ID (Bin Size = 100 Bp")



cluster_loc<- ggplot(longdata, aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  scale_fill_distiller(palette = "RdPu") +
  #scale_fill_gradient(low = "grey", high = "blue")+
  geom_vline(xintercept = 215, alpha = 0.05 )
  #geom_vline(xintercept = no_extlow, alpha = 0.15, color = "#219ebc")
  #theme(
    #axis.text = element_text(size = 2))
 
highexp <- ggplot(longdata, aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  scale_fill_distiller(palette = "RdPu") +
  #scale_fill_gradient(low = "grey", high = "blue")+
  geom_vline(xintercept = 215, alpha = 0.05, color = "grey" )+
  geom_vline(xintercept = vhi_ovary, alpha = 0.15, color = "#219ebc")+
  theme(
    axis.text = element_text(size = 2),
    legend.position = "none")

pdf("histone_cluster.pdf")
plot_grid(highexp, lowexp)
dev.off()

pdf("histone_cluster_loc.pdf")
plot_grid(cluster_loc)
dev.off()
#Make annotation data into a list of Vlines
  
  
ggplot(longdata, aes(x = Var1, y = Var2, fill = value)) +
  geom_raster() +
  scale_fill_distiller(palette = "RdPu") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 1336))+
  #scale_fill_gradient(low = "grey", high = "blue")+
  geom_vline(xintercept = 215, alpha = 0.05, color = "grey" )+
  geom_vline(xintercept = vhi_ovary, alpha = 0.15, color = "#219ebc")+
  theme(
    axis.text = element_text(size = 2),
    legend.position = "none", 
    )
  
  
#Overlap the files with Binned References


### Make a rug plot with the chromosomes ###
library(tidyverse)
bins_df <- as.data.frame(bins) %>% rownames_to_column(var = "bin_id")
write_tsv(bins_df, "/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/analysis/binned_dm6_100000bins.tsv")

chromosomes  

for(f in unique(bins_df$seqnames)){
  print(f)
  
  #Make a vector for id
  v <- dplyr::filter(bins_df , seqnames == f) %>% 
    .$bin_id
  
  assign(f, v)  
  
}


longdata3 <- longdata2 %>%
  rename(GenomeBin = Var1, HistoneClusterBin = Var2, Log10_ReadCount = value) %>%  
  mutate(chr = ifelse(GenomeBin %in% chr2L, "chr2L", 
                      ifelse(GenomeBin %in% chr2R, "chr2R", 
                             ifelse(GenomeBin %in% chr3L, "chr3L", 
                                    ifelse(GenomeBin %in% chr3R, "chr3L", 
                                           ifelse(GenomeBin %in% chr4, "chr4", 
                                                 ifelse(GenomeBin %in% chrX, "chrX", NA)))))))

write_tsv(longdata3, "/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/analysis/histone_cluster_meltedmatrix.tsv")


ggplot(longdata3, aes(x = GenomeBin, y = HistoneClusterBin, fill = Log10_ReadCount)) +
  geom_raster() +
  scale_fill_distiller(palette = "RdPu") +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 1336))+
  #scale_fill_gradient(low = "grey", high = "blue")+
  #geom_vline(xintercept = 215, alpha = 0.05, color = "grey" )+
  #geom_rug(col="steelblue",alpha=0.1, size=1.5)+
  #geom_vline(xintercept = vhi_ovary, alpha = 0.15, color = "#219ebc")+
  facet_wrap( ~ chr, scales = "free_x")+
  theme(
    axis.text = element_text(size = 5)#,
    #legend.position = "none", 
  )



#Get a list of bins with 5 
high_mapping <- longdata3 %>%  
  dplyr::filter(Log10_ReadCount > 4)






#Annotation Data

forty=GRanges(seqnames =c("chr2L"), 
              ranges=IRanges(start =c(21415940), end=c(21543673)))
forty_matrix = countOverlaps(bins, forty)


setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/annotation_data/")

library(tidyverse)
ext_high_ov <- data.table::fread("extremly_highovary.txt") %>%  
  mutate(LOCATION_ARM = str_c("chr", .$LOCATION_ARM)) #Add chr to the location arm so that it can be intersected
  
hi_ov_gr <- makeGRangesFromDataFrame(high_ov, 
                         seqnames.field = "LOCATION_ARM", 
                         start.field = "LOCATION_MIN", 
                         end.field = "LOCATION_MAX")






v_high_ov <- data.table::fread("Veryhigh_ovary_expression.tsv") %>% 
  mutate(LOCATION_ARM = str_c("chr", .$LOCATION_ARM))


v_hi_ov_gr <- makeGRangesFromDataFrame(high_ov, 
                                       seqnames.field = "LOCATION_ARM", 
                                       start.field = "LOCATION_MIN", 
                                       end.field = "LOCATION_MAX")


get_expression_vec <- function(file){
  df <- data.table::fread(file) %>%  
    mutate(LOCATION_ARM = str_c("chr", .$LOCATION_ARM))
  
  gr <- makeGRangesFromDataFrame(df, 
                                 seqnames.field = "LOCATION_ARM", 
                                 start.field = "LOCATION_MIN", 
                                 end.field = "LOCATION_MAX")
  
  bin_IDs<- countOverlaps(bins, gr) %>%  
    as.data.frame() %>% 
    rownames_to_column %>%  
    filter(. == 1) %>% 
    pull(rowname)
  
  return(bin_IDs)
}


vhi_ovary <- as.numeric(get_expression_vec("Veryhigh_ovary_expression.tsv"))
hi_ovary <- get_expression_vec("high_ovary.tsv")
vlow_ovary <- as.numeric(get_expression_vec("verylow_ovary.tsv"))

no_extlow <- as.numeric(get_expression_vec("no_extmlow_ovaryexpression.tsv"))

#Make a character vector of bins with 42AB
forty_bins<- forty_matrix %>% as.data.frame() %>% rownames_to_column() %>% filter(. == 1) %>%  
  pull(rowname)




df<- countOverlaps(bins, hi_ov_gr) %>%  as.data.frame() %>%  
  rownames_to_column() %>%  
  mutate(color = ifelse( . >= 1, "#6a040f", NA)) %>%  
  mutate_cond(rowname %in% forty_bins, color = "black") %>%  
  #mutate_cond(rowname %in% vhi_ovary, color = "#9d0208") %>%  
  #mutate_cond(rowname %in% high_ovary, color = "#d00000" ) 
  mutate_cond(rowname %in% vlow_ovary, color = "#abc4ff")
  

colors <- df$color  


setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/plots")

pdf("42AB_multimapping_profile.pdf")

heatmap(matrix, Colv = NA, Rowv = NA, RowSideColors = colors , xlab = "42AB 100 bp Bin", ylab = "Genomic Pos 100000 bp Bin")
dev.off()




data




