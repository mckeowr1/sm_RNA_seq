# Scripts

#### process_multimapping

### roi_multimapping
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


`parallel_linesearch.sh`
Runs Filterline program across all search beds https://github.com/miku/filterline (must be complied out of normal conda env)

### Cleanup 

`check_linesearch.sh`
This checks that all the bed folders have the right number of beds. This may not actually be that useful to keep around but it was useful during development.

`organize_beds.sh` 
The Filterline program runs inside the lines directory for each ROI bin (ex bin1) and outputs a bed file specific to each search bed. We beed to concatentate them and move them to thier own directory for downstream analyses. 

`clean_beds.sh` 
We also have to delete the beds that are output from the FilterLine Program this script does that. At some point they need to be combined 

### Analysis 

`plot_multimappers.R` 
This scripts takes the ROI binned beds with thier other reads locations and plots them as a matrix.


# Workflow 
1) Clone the github repo 

2) Move your multimapped BAM file, it's indexs `.bai` and index created in process_multimapping steps to the BAMS file

3) Within the multimappings folder there is a projects directory. Within projects create a file for the ROI you wish to analyze (ex. Histone Cluster). 

4) Create a bed file using the template.bed file in the projects directory. Edit the genomic coordinates to your ROI. Save the file with the name of your ROI. ex) `histonecluster.bed` Move the bed file into your project dir. 

5) Run `regions_multi_mappers.sh`. The first command line arument for this script is the directory that contains your multimapped BAM files. The second argument is the path to your ROI.bed file. 

    ````
    regions_multi_mappers.sh path/to/multimapped/bams ~/multi_mappings/projects/your_roi/roi.bed    

    ````
6) The above step should sucessfuly create a `ROI.bam` file along with a `ROI_reads.bed` file. This is a reduced verison of the BAM file that was also created. Now that we have this bed file, we want to break it into bins using the `bin_roi.R` script. Load the function `generate_binned_readnames()` and set the working directory to your project directory. This will create the `binned.readnames` directory in your ROI folder. 
    ````
    setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster")
    
    
    generate_binned_readnames(roi_readnames_bed = "~/multi_mappings/projects/ROI/ROI.bed",
                          roi_chr = "chr2L", 
                          roi_start = 21415940, 
                          roi_stop = 21543673, 
                          proj_dir = "~multi_mappings/projects/ROI")
    ````
    
7) Now that we have our `binned.readnames` file we want to get all of the locations of those reads. To do that we use the `bedindex.tsv` that we created when we processed our multimapping BAM and the `find_reads.R` script. Sets is how many times you plan on running the `parallel_linesearch.sh` program that proceeds this script. Breaking the program into multuiple sets can be good since I'm not sure what the limit of the computer is. I have done it all in one set though. 

    ````
    find_reads( 
    search_bed_dir = "/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/",
    proj_dir = "/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/projects/testing", 
    bam_dir = "/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/BAMS",
    index_file = "/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/BAMS/SRR1187947_mapped_verysensitive_local.mapped_sorted.bedindex.tsv", 
    sets = 5 
    )

    ````

8) 



