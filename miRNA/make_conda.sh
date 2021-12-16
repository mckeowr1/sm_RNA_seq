
#Put Channels In the Right order for Bioconda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

#Create the environment
conda create -n mirdeep

#Activate the environment
conda activate mirdeep

#Install mirdeep2
conda install -c bioconda mirdeep2