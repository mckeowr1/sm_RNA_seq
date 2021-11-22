library(data.table)
library(tidyverse)

setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB/")


index <- data.table::fread("SRR1187947_mapped_verysensitive_local.mapped_sorted.bedindex.tsv")



#Filter the index based on Read Names Input

#Region


read_names <- readLines("files/42AB_multimapping_readnames.txt")


setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB_Lines/")
for(read in read_names){ 
  get_bed_location(index, read)
  
  } 


#Function to Search through a Bed Index and Pull out line loactions 
get_bed_location <- function(bedindex, 
                             readname)
{
  ## At some point want to figure out how to run this for multiple reads## 
  # #Read in read names to character vector 
  # read_names <- readLines(read_name )
  
  read <- bedindex %>% dplyr::filter(V4 == readname) 


  #Never Use this! It onlly take the 1st part of the list
  # lines <- as.character(read[,2]) %>%
  #   str_split(pattern = ",") %>%  #For some reason this converts to a list
  #   .[[1]] %>% #Make it not a list 
  #   sort()

    writeLines(lines, glue::glue("{readname}_lines.txt"))
  
  
} 


setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/42AB/")


#Function to look for multiple reads listed in a text file 


##get_multiread_bed_location 

#Make this a loop through a list of files

#For a file 

get_bined_bed <- function(readnamesfile, index){ 

read_names<- readLines(readnamesfile)

#Has to be a better way to bane the bin
bin_name <- str_split(readnamesfile, pattern = "_") %>% .[[1]] %>% .[1]

#Filter the index to reads
fil_index<- index %>% dplyr::filter(V4 %in% read_names)


#Make search files 
aa <- read %>% filter(fil_index)
ab <- read %>% filter(fil_index)
ac 
ad 
ae 
af 
ag 
ah 



#Pull Bed line numbers from the index
indexes <- dplyr::pull(fil_index, "index")

#Break apart multi position reads
lines <- unlist(str_split(indexes, pattern = ",")) %>%  
  as.numeric() %>%  
  base::sort()
  
#Write out the lines
write(lines, glue::glue("{bin_name}_lines.txt"), sep = "\n")

} 

#Mutate the index ? if >1000

system.time(get_bined_bed("bin10_readnames.txt", index))


read_names<- readLines("bin10_readnames.txt")

fil_index<- index %>% dplyr::filter(V4 %in% read_names)

lines <- unlist(str_split(indexes, pattern = ",")) %>%  
  as.numeric() %>%  
  base::sort()

indexes <- dplyr::pull(fil_index, "index")

aa <- lines[lines <= 100000000]
ab <- lines[between(lines, 100000001 , 200000000)] - 100000000
ac <- lines[between(lines, 200000001 , 300000000)] - 200000000


write(aa, "bin10_lines_aa.txt", sep = "\n")



#as.numeric(lines) %>% sort()












