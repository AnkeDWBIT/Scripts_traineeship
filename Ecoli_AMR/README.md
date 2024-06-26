# Ecoli_AMR
Inside the Git repository "Scripts_traineeship", the folder "Ecoli_AMR" contains all scripts I made during my traineeship for the ba-na-ba bioinformatics program (04/2024 - 06/2024) for the E. coli study that compares bioinformatics tools for AMR detection.
This README file describes the order in which the scripts can be used, what the scripts do and how they can be used with examples of command line options.

## Input files
E. coli genomic FASTA-files.

## Bioinformatics tools for AMR detection
### AMRFinderPlus
#### Installing the tool
The tool can be installed in a conda-environment as following:
```
$ conda create -y -c conda-forge -c bioconda -n amrfinder ncbi-amrfinderplus
$ conda activate amrfinder
$ amrfinder --database_version
$ amrfinder -u
```
The version I had installed: Software version: 3.12.8 & Database version: 2024-05-02.2

#### AMRFinderPlus.sh
Script that loops through an input directory and runs AMRFinderPlus serially for all FASTA-files. The output directory "AMR_output" where the resulting tsv-files are saved, is created inside the input directory.
```
$ bash AMRFinderPlus.sh [1]
[1] = Full path to input directory with FASTA files 
```

#### tsv_to_excel.py
Script that makes an output directory "AMR_output_excel" inside the input directory and copies the tsv-files into Excel-files.
```
$ python tsv_to_excel.py [1]
[1] = Full path to input directory (=folder containing only tsv-files)
```

#### AMRFinder_results.py
! DRAFT VERSION  - script still has some errors !
Script to write results of AMRFinderPlus to an Excel-file
```
$ python AMRFinder_results.py [1] [2]
[1] = Full path to directory with AMRFinderPlus results
[2] = Full path to Excel file with AMR tools summary
```
#### RGI_CARD.sh
Script that loops through an input directory and executes RGI CARD on all fasta-files.
```
$ bash RGI_CARD.sh[1]
[1] = Full path to input directory with fasta-files.
```