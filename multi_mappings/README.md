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

10) We now have our data and it's time to analyze it! One way to visualize the data is by generating a multi_mapped matrix which shows which parts of the genome our histone cluster bins map to. This matrix is generated by the script `multi_mapped_matrix.R`. It takes a our bed directory created in the step above and intersects each bed with a dm6 genome that has been binned (bin size set by the genome_bin_size parameter) to generate a vector or reads per genomic bin. The vectors are combined into a matrix and the matrix is melted into a long `.tsv` file that can be plotted or analyzed in different scripts. 
    ````
    multi_map_matrix(roi_name = my_roi, 
                    genome_bin_size = 100000, 
                    bed_dir = "~/multimappings/projects/your_project/roi_beds", 
                    out_dir = "~multimappings/projects/your_project" 

                     )
    ````




