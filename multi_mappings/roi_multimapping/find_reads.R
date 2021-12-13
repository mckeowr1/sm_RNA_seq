library(data.table)
library(tidyverse)



find_reads <- function(search_bed_dir, proj_dir, bam_dir, index_file, sets){ 
  ## Load the required files
  
  bin_readnames <- list.files(glue::glue("{proj_dir}/binned.readnames")) #Get the binned readnames files
  
  index <- data.table::fread(glue::glue("{index_file}")) #Load the index file
  

  
  
  ## Split the files into sets 
  set_ranges <- split(bin_readnames,    
                      cut(seq_along(bin_readnames),
                          sets,
                          labels = FALSE))
  
  ## Search for the lines and write out lines file
  
  #Make a binned lines directory for the output of the function 
  dir.create(glue::glue("{proj_dir}/binned_lines"))
  
  get_bed_location <- function(readnamesfile, index) {
    
    #Name for the bin files
    bin_name <- str_split(readnamesfile, pattern = "_") %>% .[[1]] %>% .[1]
    
    dir.create(glue::glue("{proj_dir}/binned_lines/{set_name}/{bin_name}_lines"))
    
    
    #Set Working Directory to the binned lines in the project 
    setwd(glue::glue("{proj_dir}/binned.readnames"))
    #setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/projects/histone_cluster/binned.readnames") #For testing
    
    #Create a Vector of readnames
    read_names<- readLines(readnamesfile)
    
    
    #Filter the index to reads
    fil_index<- index %>% dplyr::filter(V4 %in% read_names)
    
    indexes <- dplyr::pull(fil_index, "index")
    
    #Create a vector of lines where reads are located for bin 
    lines <- unlist(str_split(indexes, pattern = ",")) %>%  
      as.numeric() %>%  
      base::sort()
    
    setwd(search_bed_dir)
    num_search_beds <- length(list.files()) #Find out how many search beds have been made 
    
    #Now subset the lines to which search file they need to be pulled from 
    for(i in 1:num_search_beds){
      letter <- letters[i]
      
      
      end <- i * 100000000
      start <- end - 100000000  + 1
      sub <- end - 100000000 
      
      
      if(i == 1) { 
        search_lines <- lines[lines <= 100000000]
      }
      if(i == max(1:num_search_beds)){ 
        search_lines <- lines[lines >= start] - sub
      } else
        
        search_lines <- lines[between(lines, start, end)] - sub
      
      
      
      
      setwd(glue::glue("{proj_dir}/binned_lines/{set_name}/{bin_name}_lines"))
      
      write(search_lines, glue::glue("{bin_name}_a{letter}.txt"), sep = "\n")
      
      setwd(glue::glue("{proj_dir}/binned.readnames"))
      #setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/projects/histone_cluster/binned.readnames") #For testing
      
    }
  }
  
  #Loop through all the sets
  for(set in 1:length(set_ranges)){ 
    set_name <- paste("set", names(set_ranges[set]), sep = "")
    print(set_name)
    dir.create(glue::glue("{proj_dir}/binned_lines/{set_name}"))
    
    for(bin in set_ranges[[set]]){
      get_bed_location(bin, index)
    }
  }  
  } 



find_reads( 
  search_bed_dir = "/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/",
  proj_dir = "/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/projects/testing", 
  bam_dir = "/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/BAMS",
  index_file = "/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/BAMS/SRR1187947_mapped_verysensitive_local.mapped_sorted.bedindex.tsv", 
  sets = 5 
  )



## For Troubleshootin ##


#Function Inputs

# search_bed_dir <- c("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/")
# 
# #proj_dir <- c("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/projects/testing")
# proj_dir <- c("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/projects/histone_cluster")
# 
# 
# bam_dir <- c("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/BAMS")
# 
# index_file <- c("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/BAMS/SRR1187947_mapped_verysensitive_local.mapped_sorted.bedindex.tsv")
# 
# sets <- 5




















