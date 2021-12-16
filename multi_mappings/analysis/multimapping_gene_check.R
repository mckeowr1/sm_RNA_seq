library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(tidyverse)
library(plyranges)


histone_cluster <- data.table::fread("../histone_cluster/histone_cluster_meltedmatrix.tsv")

binned_genome <- data.table::fread("../histone_cluster/binned_dm6_100000bins.tsv") %>% 
  makeGRangesFromDataFrame(. ,seqnames.field = "seqnames", 
                           start.field = "start", 
                           end.field = "end", 
                           strand.field = "strand")

txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

#How to incorproate the reference genome ???

check_gene_overlaps <- function(gene_list, feature_type = "gene", roi_matrix_file, roi_beds, pad_gene = F, pad_size = 2000){

  df <- data.frame() # blank df to hold the output of the loop
  
  #Load and bin the reference genome 
  binned_genome = tileGenome(seqlengths(Dmelanogaster), tilewidth = 100000, cut.last.tile.in.chrom = TRUE)
  binned_genome = binned_genome[seqnames(binned_genome) %in% paste0('chr', c('2L','2R','3L','3R','4','X'))] # this can be commented out if you want to keep all chromosomes/contigs.
  binned_genome = binned_genome[width(binned_genome) == 100000] 
  
  roi_matrix <- data.table::fread(roi_matrix_file)
  
  for (gene in gene_list){
    name <- gene
    
    #### Find what genomic bin the gene overlaps with #### 
    #Filter genes GR to the specific gene or transcript
    
    if(pad_gene == T){
      if(feature_type == "gene"){
        genes <- genes(txdb)
        g1 <- plyranges::filter(genes, gene_id == gene)
        g <- stretch(anchor_center(g1), pad_size)
          
      }
      if(feature_type == "transcript"){
        transcripts <- transcripts(txdb)
        g1 <- plyranges::filter(transcripts, tx_name == gene)
        g <- stretch(anchor_center(g1), pad_size)
      }
    }
    else{
    
      if(feature_type == "gene"){
        genes <- genes(txdb)
        g <- plyranges::filter(genes, gene_id == gene)
     }
      if(feature_type == "transcript"){
        transcripts <- transcripts(txdb)
        g <- plyranges::filter(transcripts, tx_name == gene)
      }
    }
    
    #Count Overlaps of the gene obj with the binned genome  
    t <- suppressWarnings(countOverlaps(binned_genome, g))
    
    #Return the Genomic Bins Each gene overlaps with
    g_bins <- which(t > 0)
    
    #Check that the gene dosen't cover two genomic bins - if it does kill for this gene it (this is an extreme edge case since the genomic bins are so larg
    if(length(g_bins) > 1){ 
      print(glue::glue("The gene {name} spans two genomic bins. This is an edge case that this function cannot handle at the moment"))
      next }

    #### Zoom in on that genomic bin #### 
    genomic_region <- binned_genome[g_bins] #Filter the binned genome to just the bin the gene overlaps with
    zoom_gene <- tile_ranges(genomic_region, 20) #More finely tile that bin of the genome which is used later in the function.
    
    #### Find what ROI bins map to that genomic bin and find where in the Zoomed genomic bin they overlap ####
    
    #Make a list of all the ROI beds 
    all_beds <- list.files(glue::glue("{roi_beds}"))
    
    #Produce a Vector of ROI bins that map to that genomic loci
    roi_bin <- roi_matrix %>% 
      dplyr::filter(GenomeBin == g_bins) %>%  
      dplyr::filter(Log10_ReadCount > 0) %>%  
      .$HistoneClusterBin
    
    #If there are no histone cluster bins that map to skip the gene and return a message, and append a O to the df
    if(length(roi_bin) == 0){
      print(name)
      print("There are no histone cluster reads in this gene's genomic bin")
      num_overlaps <-c(0)
      #Create a df for the gene 
      t <- data.frame(name, num_overlaps)
      df <- rbind(df, t)
      
      next
    }
    #Filter the list of all beds to only include beds within the ROI
    bed <- all_beds[all_beds %in% roi_bin]
    
    #Load up those ROI beds and intersect them with the zoomed in genome
    
    make_vector<- function(file){ 
      
      
      df <- data.table::fread(glue::glue("{roi_beds}/{bed[1]}"))
      
      
      gr <- makeGRangesFromDataFrame(df, 
                                     seqnames.field = "V1", 
                                     start.field = "V2", 
                                     end.field = "V3", 
                                     strand.field = "V6")
      c <- suppressWarnings(countOverlaps(zoom_gene, gr))
      
      return(c)                             
      
    }
    
    myfiles = lapply(bed, make_vector)
  
  
  #Change the name of the columns in the matrix to the name of the ROI bed they came from
  for(i in 1:length(bed)) {
    
    names(myfiles)[[i]] = bed[i]
  }
  
  
  matrix = do.call("cbind", myfiles)
  
  
  out <- as.data.frame(matrix)

  #### Now intersect the zoomed in genomic region with the Gene GR obj ####
  
  #This variable comes from above when we filter the genes GR obj
  
  print(glue::glue("Calculating Zoomed bins for {name}:")) #A quick status update
  
  gene_overlaps <- suppressWarnings(countOverlaps(zoom_gene, g))
  
  
  #### Compare where the Gene is to where we have multimapping reads
  
  #Get a list of bins where ROI beds overlapped with Zoomed Genomic Bin 
  roi_bed_overlaps <- c()
  for(i in 1:length(bed)){
    #Get the zoomed gene bin that it overlaps with 
    overlaps <- (which(out[bed[i]] > 0))
    #Add that to a growing vector
    roi_bed_overlaps <- unique(append(roi_bed_overlaps, overlaps))
    #Make a Plot of the results in the future - gene zoom bin X gene_overlaps & ROI beds overlap
  }
  
  #Intersect the list of bins the gene maps to and the list of bins that have ROI mapped reads
  gene_zoom_bins <- which(gene_overlaps >0) # Filter to just bins where the gene is 
  num_overlaps <- length(intersect(gene_zoom_bins, roi_bed_overlaps)) #Get the number of bins that are the same in both lists 
  
  #Convert num_overlaps to numeric vector so we can append negative results to the dataframe 
  if(length(num_overlaps) == 0){ 
    num_overlaps <-c(0)
  }
  #Create a df for the gene 
  t <- data.frame(name, num_overlaps)
  
  df <- rbind(df, t)
  
  }
  return(df)
}


test_genes <- c("FBtr0072639")

check_gene_overlaps(gene_list = test_genes, 
                    feature_type = "transcript", 
                    roi_matrix_file = "/projects/b1059/projects/Ryan/mappings/histone_cluster/histone_cluster_meltedmatrix.tsv",
                    roi_beds = "/projects/b1059/projects/Ryan/mappings/histone_cluster/histone_cluster_beds", 
                    pad_gene = T, 
                    pad_size = 2000)