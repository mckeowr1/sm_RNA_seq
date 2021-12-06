library(data.table)
library(tidyverse)
library(feather)
library(ggplot2)


#Sets is the number of times you want to run the find lines c program. 

find_readlines <- function(proj_dir , sets) { 
  
  #Load up out Index File
  setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/") #Set the Directory to the loaction of the Index
  index <- data.table::fread("SRR1187947_mapped_verysensitive_local.mapped_sorted.bedindex.tsv")
  
  bin_readnames <- list.files(glue::glue("{proj_dir}/binned.readnames"))
  
  #Break the readnames into chunks to be run individually to not blow out computer 
  
  num_bins <- lenth(bin_readnames)
  
  set_size <- num_bins / sets
  
  set_ranges <- split(bin_readnames,          
                       cut(seq_along(bin_readnames),
                           set_size,
                           labels = FALSE))
  counter <- 0 
  for(i in set_ranges) { 
    counter <- counter + 1 
    set_id <- paste("set", counter, sep = "")
    assign(set_id, i)
      
  }  
      
  get_bed_location <- function(readnamesfile, index) {
    
    #Name for the bin files
    bin_name <- str_split(readnamesfile, pattern = "_") %>% .[[1]] %>% .[1]
    
    dir.create(glue::glue("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/binned_lines4/{bin_name}_lines"))
    
    #Create a Vector of readnames
    read_names<- readLines(readnamesfile)
    
    
    #Filter the index to reads
    fil_index<- index %>% dplyr::filter(V4 %in% read_names)
    
    indexes <- dplyr::pull(fil_index, "index")
    
    #Create a vector of lines where reads are loacted for bin 
    lines <- unlist(str_split(indexes, pattern = ",")) %>%  
      as.numeric() %>%  
      base::sort()
    
    #Now subset the lines to which file they need to be pulled from 
    
    aa <- lines[lines <= 100000000]
    ab <- lines[between(lines, 100000001 , 200000000)] - 100000000
    ac <- lines[between(lines, 200000001 , 300000000)] - 200000000
    ad <- lines[between(lines, 300000001 , 400000000)] - 300000000
    ae <- lines[between(lines, 400000001 , 500000000)] - 400000000
    af <- lines[between(lines, 500000001 , 600000000)] - 500000000
    ag <- lines[between(lines, 600000001 , 700000000)] - 600000000
    ah <- lines[lines >= 700000001] - 700000000
    
    setwd(glue::glue("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/binned_lines4/{bin_name}_lines"))
    
    write(aa, glue::glue("{bin_name}_aa.txt"), sep = "\n")
    write(ab, glue::glue("{bin_name}_ab.txt"), sep = "\n")
    write(ac, glue::glue("{bin_name}_ac.txt"), sep = "\n")
    write(ad, glue::glue("{bin_name}_ad.txt"), sep = "\n")
    write(ae, glue::glue("{bin_name}_ae.txt"), sep = "\n")
    write(af, glue::glue("{bin_name}_af.txt"), sep = "\n")
    write(ag, glue::glue("{bin_name}_ag.txt"), sep = "\n")
    write(ah, glue::glue("{bin_name}_ah.txt"), sep = "\n")
    
    setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/binned_readnames")
    
    
    
  }
  
  
  
  for(f in bin_readnames){ 
    get_bed_location(f , index)
  }
  
  }

test <- c(1:10) 
cut_interval(test, n =2)

#Stupid function
test_ranges <- split(test,             # Applying split() function
      cut(seq_along(test),
          2,
          labels = FALSE))
vec <- test_ranges[1][[1]]
vec <- test_ranges[2][[1]]

counter <- 0 
for(set in test_ranges){ 
  print(set)
  counter <- counter + 1
  print(counter)
  
  }

#Set the Directory to the loaction of the Index
setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/")

index <- data.table::fread("SRR1187947_mapped_verysensitive_local.mapped_sorted.bedindex.tsv")


#Set the Directory to the Projects Binned Readnames File 
setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/binned.readnames")


#Get a list of the binned readnames files

bin_readnames <- list.files("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/binned.readnames")

set1 <- bin_readnames[1:600]
set2 <- bin_readnames[601:1200]
set3 <- bin_readnames[1201:1800]
set4 <- bin_readnames[1801:2400]
set5 <- bin_readnames[2401:2439]

setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_test/sets")

writeLines(set5, "set5.txt")

setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/binned.readnames")

#Load up function to pull out index position
get_bed_location <- function(readnamesfile, index) {
  
  #Name for the bin files
  bin_name <- str_split(readnamesfile, pattern = "_") %>% .[[1]] %>% .[1]
  
  dir.create(glue::glue("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/binned_lines/{bin_name}_lines"))
  
  #Create a Vector of readnames
  read_names<- readLines(readnamesfile)
  
  
  #Filter the index to reads
  fil_index<- index %>% dplyr::filter(V4 %in% read_names)
  
  indexes <- dplyr::pull(fil_index, "index")
  
  #Create a vector of lines where reads are loacted for bin 
  lines <- unlist(str_split(indexes, pattern = ",")) %>%  
    as.numeric() %>%  
    base::sort()
  
  #Now subset the lines to which file they need to be pulled from 
  
  aa <- lines[lines <= 100000000]
  ab <- lines[between(lines, 100000001 , 200000000)] - 100000000
  ac <- lines[between(lines, 200000001 , 300000000)] - 200000000
  ad <- lines[between(lines, 300000001 , 400000000)] - 300000000
  ae <- lines[between(lines, 400000001 , 500000000)] - 400000000
  af <- lines[between(lines, 500000001 , 600000000)] - 500000000
  ag <- lines[between(lines, 600000001 , 700000000)] - 600000000
  ah <- lines[lines >= 700000001] - 700000000
  
  setwd(glue::glue("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/binned_lines/{bin_name}_lines"))
  
  write(aa, glue::glue("{bin_name}_aa.txt"), sep = "\n")
  write(ab, glue::glue("{bin_name}_ab.txt"), sep = "\n")
  write(ac, glue::glue("{bin_name}_ac.txt"), sep = "\n")
  write(ad, glue::glue("{bin_name}_ad.txt"), sep = "\n")
  write(ae, glue::glue("{bin_name}_ae.txt"), sep = "\n")
  write(af, glue::glue("{bin_name}_af.txt"), sep = "\n")
  write(ag, glue::glue("{bin_name}_ag.txt"), sep = "\n")
  write(ah, glue::glue("{bin_name}_ah.txt"), sep = "\n")
  
  setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster/binned.readnames")
  
  
  
}


for(f in bin_readnames){ 
  get_bed_location(f , index)
}
