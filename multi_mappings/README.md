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
The Filterline program runs inside the lines directory for each ROI bin (ex bin1) and outputs a bed file specific to each search bed. We beed to concatentate them and move them to their own directory for downstream analyses. 

`clean_beds.sh` 
We also have to delete the beds that are output from the FilterLine Program this script does that. At some point they need to be combined 

### Analysis 

`plot_multimappers.R` 
This scripts takes the ROI binned beds with thier other reads locations and plots them as a matrix.

# Set up local environments 
Things that need to go into the environment:
- bedtools

To get bedtools program we need to create a conda environment
```
#Put Channels In the Right order for Bioconda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

#Create the environment
conda create -n myenv

#Activate the environment
conda activate myenv

#Install Bedtools
conda install -c bioconda bedtools

```
For all subsequent steps we need to have bedtools. Make sure the environment is active when


# Process Multimapping SAM
This section describes how to create a database from a multimapped BAM file. This is a relatively simple process however, its computationally intensive so at the moment these steps have to be performed on Quest or another compute cluster. All scripts for this section are found in the `process_multimapping` directory.

1) The first step is to convert your BAM file to a BED file using the bedtools function `bamtobed` and then to sort the file using the bedtools function `sort` which sorts by chromosome and start position. The file `process_multimapping.sh` contains a slurm script to perform these two steps on Quest.

2) Next you want to build a read key from the bed file. This is just a TSV that contains the line number for each read in the bed file. This is built using the `build_bed_index.R` script. To use the script update the path to your bed file. Then call the `build_bed_index()` function. This should generate a large index file for your beds. 

3) The final step is to create search beds to optimize the search for reads in later steps. This is done by using the UNIX command `split`. It is crucial to note the number you use for the lines parameter since this is will affect the `find_reads.R` script later on. `split_bed.sh` is a slurm script to split your bed file into 100000000 files. 

# Analyze an ROI
1) Clone the github repo 

2) Move your multimapped BAM file, it's indexes `.bai` and index created in process_multimapping steps to the BAMS file

3) Within the multimappings folder there is a projects directory. Within projects create a file for the ROI you wish to analyze (ex. Histone Cluster). 

4) Create a bed file using the template.bed file in the projects directory. Edit the genomic coordinates to your ROI. Save the file with the name of your ROI. ex) `histonecluster.bed` Move the bed file into your project dir. 

5) Run `regions_multi_mappers.sh`. The first command line argument for this script is the directory that contains your multimapped BAM files. The second argument is the path to your ROI.bed file. 

    ````
    regions_multi_mappers.sh path/to/multimapped/bams ~/multi_mappings/projects/your_roi/roi.bed    

    ````
6) The above step should successfully create a `ROI.bam` file along with a `ROI_reads.bed` file. This is a reduced version of the BAM file that was also created. Now that we have this bed file, we want to break it into bins using the `bin_roi.R` script. Load the function `generate_binned_readnames()` and set the working directory to your project directory. This will create the `binned.readnames` directory in your ROI folder. 
    ````
    setwd("/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/histone_cluster")
    
    
    generate_binned_readnames(roi_readnames_bed = "~/multi_mappings/projects/ROI/ROI.bed",
                          roi_chr = "chr2L", 
                          roi_start = 21415940, 
                          roi_stop = 21543673, 
                          proj_dir = "~multi_mappings/projects/ROI")
    ````
    
7) Now that we have our `binned.readnames` file we want to get all of the locations of those reads. To do that we use the `bedindex.tsv` that we created when we processed our multimapping BAM and the `find_reads.R` script. The Sets parameter is how many times you plan on running the `parallel_linesearch.sh` program in step 8. Breaking the program into multuiple sets can be good since I'm not sure what the limit of the computer is. Though I did end up doing the histone cluster in one set. Set the following parameters in in R script and run the function. 

    ````
    find_reads( 
    search_bed_dir = "/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/search_beds/",
    proj_dir = "/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/projects/testing", 
    bam_dir = "/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/BAMS",
    index_file = "/Users/ryan/Documents/GitHub/sm_RNA_seq/multi_mappings/BAMS/SRR1187947_mapped_verysensitive_local.mapped_sorted.bedindex.tsv", 
    sets = 5 
    )

    ````
    This script will spit out a new directory in your project dir called `binned_lines`. Within `binned_lines` there will be a number of set directories ex `set1` based on however many sets you specified in the `find_reads.R` script. Each set contains a directory for each bin within the set. ex `bin1_lines`. Within the directory for that bin are text files that contain the lines for a particular search file. 

    ````
   
    ├── binned_lines
        ├── set1
            ├── bin1_lines
                ├── bin1_aa.txt 
                ├── bin1_ab.txt
                ├── bin1_ac.txt 
                ├── ... 
            ├── bin2_lines

    ````
    

8) We now can search for our reads! To do this we'll use the `parallel_linesearch.sh`. This script uses a C program called `filterline` to filter our search beds by line number extremely quickly (much faster than awk). To install filterline switch to your desired installation directory and run the following commands:
    ````
    #switch to install directory 

    cd path/to/install/dir

    #Clone and compile the repo 

    git clone https://github.com/miku/filterline.git
    cd filterline
    make
    
    ````
    The filterline program takes two inputs. 1) A file with the line numbers you wish to return and 2) a file you wish to pull those lines from.
    ````
    path/to/install/dir/filterline lines_file search_file

    ````
    We want to search through lots of files at the same time so we spin up a bunch of instances of this search program using the `parallel_linesearch.sh` script. At the moment this script has to be manually updated for each search we want to perform (working on a better function). There are 2 things that need to be edited. 
    
    1st edit the paths to the search_beds and the lines_directory
    ````
    search_path=~/multimappings/search_beds
    
    lines_dir=~/multimappings/projects/your_roi/binned_lines/set1
    
    ````
    

    2nd edit the `filterline` path to your install directory

    ````

    /Users/ryan/filterline/filterline

    ````

    If you have more than 8 search beds make sure to define an additional search_bed, lines file, and search program.

    Once all the changes are made and saved then you should be ready to run the program. It's best to run overnight while you won't be using your computer. To run the program execute the following command.

    ````
    $ bash parallel_linesearch.sh

    ````
    When we come back a few hours later, the file search should be complete. I like to verify that no more `filterline` programs are running by opening activity monitor on mac OS and searching for filterline. If any processes are still running, give it a few more hours. If they are still going after that, there was probably an error and the program has spun out (this has never happened). 

    This program will add a new folder `beds` to each `bin*_lines` folder.
    ````
    ├── binned_lines
        ├── set1
            ├── bin1_lines
                ├── beds 
                    ├── aa.bed
                    ├── ab.bed
                    ├── ac.bed
                    ├── ... 
                ├── bin1_aa.txt 
                ├── bin1_ab.txt
                ├── bin1_ac.txt 
                ├── ... 
            ├── bin2_lines
    ````

9) Now we should have a bunch of small bed files for each bin that contain the read ID and and the locations that its mapped across the genome. We want to clean up this directory structure a bit to make it easier to work with. To do this we'll use the `organize_beds.sh` script. This program goes into each line directory and concatenates all the beds. Then it move them to an output directory of your choosing. Best practice would be within your project directory ex) `~/multimappings/projects/your_projects/roi_beds`. Similar to the `parallel_linesearch.sh` script we have to update two variables manually. 
    ````
    lines_dir=~/multimappings/projects/your_project/binned_lines/set1
    out_dir=~/multimappings/projects/your_project/roi_beds
    ````
    Now our directory structure should look like this

    ````
    ├── binned_lines
        ├── set1
            ├── bin1_lines
                ├── beds 
                ├── bin1_aa.txt 
                ├── bin1_ab.txt
                ├── bin1_ac.txt 
                ├── ... 
    ├── roi_beds
        ├── binned_lines

    ````

10) We now have our data and it's time to analyze it! One way to look at the data is by generating a multi_mapped matrix which shows which parts of the genome our histone cluster bins map to. This matrix is generated by the script `multi_mapped_matrix.R`. It takes a our bed directory created in the step above and intersects each bed with a dm6 genome that has been binned (bin size set by the genome_bin_size parameter) to generate a vector or reads per genomic bin. The vectors are combined into a matrix and the matrix is melted into a long `.tsv` file that can be plotted or analyzed in different scripts. All of the following analysis code
    ````
    multi_map_matrix(roi_name = my_roi, 
                    genome_bin_size = 100000, 
                    bed_dir = "~/multimappings/projects/your_project/roi_beds", 
                    out_dir = "~multimappings/projects/your_project" 

                     )
    ````


# Analysis 

## Look for multimapping in a list of genes or transcripts 
To look for multimapping in a list of genes or transcripts we use the script `multimapping_gene_check.R` which uses the `check_gene_overlaps()` function. 

1) This script first finds the genomic bin a gene is within. 

2) It then more finely bins that genomic bin (20 bp windows)

3) We then filter for ROI bins that map to that genomic bin. We take the bed files for those bins and `countOverlaps` with the more finely binned genomic bin. This gives us better resolution of where in the genomic bin reads are multimapped to. 

4) We can now also ask where our gene is in this genomic region.

5) Then we compare the bins that the gene overlaps to the bins that histone cluster reads map to.

### Usage
````
check_gene_overlaps(gene_list = test_genes, 
                    feature_type = "transcript", 
                    roi_matrix_file = "/projects/b1059/projects/Ryan/mappings/histone_cluster/histone_cluster_meltedmatrix.tsv",
                    roi_beds = "/projects/b1059/projects/Ryan/mappings/histone_cluster/histone_cluster_beds", 
                    pad_gene = T, 
                    pad_size = 2000) # 1000 bp is added to start and end 
````
The function will return a dataframe with the name of each gene in the list and the number of 20bp gene bins that also have multimapped reads from the ROI.