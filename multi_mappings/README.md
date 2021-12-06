# Scripts

# process_multimapping

# roi_multimapping
This section contains information about the scripts used to analyze multimappings in a particular ROI (ex Histone Cluster). They are listed in the order they are used.

`region_multi_mappers.sh`
Takes an input multimapped BAM file and A ROI Bed (ex chr2R 10 1000) and returns a 
bam file for just that region of interest. Also, can return a list of reads that map to ROI. Uses samtools.

`bamtobed.sh`
Takes an ROI bam and returns a bed file for that region

`bin_roi.R`
Takes the ROI bed file, and bins the region of interest. Then Returns a list of readnames for each bin. These readnames are stored in the file binned.readnames in the project directory. Ex) bin1_readnames.txt

`find_reads.R`
Pulls the line location of a list of reads 

Inputs: 
Sorted Bed Index 
Directory with readnames text files

Output:
Lines Directory


### parallel_linesearch.sh
Runs Filterline program across all search beds https://github.com/miku/filterline (must be complied out of normal conda env)

### Cleanup 

#### check_linesearch.sh 
This checks that all the bed folders have the right number of beds. This may not actually be that useful to keep around but it was useful during development.

#### organize_beds.sh 
The Filterline program runs inside the lines directory for each ROI bin (ex bin1) and outputs a bed file specific to each search bed. We beed to concatentate them and move them to thier own directory for downstream analyses. 

#### clean_beds.sh 
We also have to delete the beds that are output from the FilterLine Program this script does that. At some point they need to be combined 