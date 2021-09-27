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

allreads = readGAlignments("SRR1187947.57.bam", 
                           param = ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = TRUE, isUnmappedQuery = FALSE), 
                                                what = c('flag','mpos','mapq','isize')),
                           use.names = TRUE)


cat('filtering out reads on non-canonical chromosomes.\n')
init = length(allreads)
allreads = allreads[seqnames(allreads) %in% good.uns]
after = length(allreads)


cat('coercing the data into a granges object.\n')
multi_reads = granges(allreads)


read_names <- unique(names(multi_reads)) 
ranges <- ranges(multi_reads)


bins = tileGenome(seqlengths(Dmelanogaster), tilewidth = 100, cut.last.tile.in.chrom = TRUE)
bins = bins[seqnames(bins) %in% paste0('chr', c('2L','2R','3L','3R','4','X'))] # this can be commented out if you want to keep all chromosomes/contigs.
bins = bins[width(bins) == 100] 



# 
# test <- union(bins, multi_reads)
# test2<- findOverlaps(multi_reads, bins)
# 
# ranges(test)


#Make a matrix for the reads

z <- matrix(0, length(bins), 1) #Define the matrix length length of binned genomes
z[subjectHits(findOverlaps(multi_reads, bins)),1] <- 1 #set the matrix value based
colnames(z) <- paste0("gr.", 1)
mcols(bins) <- z
show(bins)

#Check the work 

bins[849578]

granges(multi_reads)


cat('Assigning <genome> variable as <dm6>.\n')
genome(multi_reads) = 'dm6'


seqlevels(multi_reads) = seqlevels(dm6)
seqinfo(multi_reads) = seqinfo(dm6)
multi_reads = trim(multi_reads)

seqlevels(multi_reads) = seqlevelsInUse(multi_reads)

names(multi_reads) = NULL

seqnames(multi_reads)
