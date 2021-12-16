import.single = function(file, qualthresh = 10, exclude.duplicates = FALSE)
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
                               param = ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = FALSE, isUnmappedQuery = FALSE, isDuplicate = FALSE), 
                                                    what = c('flag','mrnm','mpos','mapq','isize')),
                               use.names = FALSE)
  }else{cat("Duplicates will not be excluded.\n")
    allreads = readGAlignments(file, 
                               param = ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = FALSE, isUnmappedQuery = FALSE), 
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
