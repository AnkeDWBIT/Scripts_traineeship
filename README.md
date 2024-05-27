# Scripts_traineeship
The Git repository named "Scripts_traineeship" contains all scripts I adapted or made during my traineeship for the ba-na-ba bioinformatics program (04/2024 - 06/2024). This README file describes the order in which the scripts can be used, what the script do and how they can be used with examples of command line options.

## Input files
Genomic fasta-files & gbff-files, can be downloaded using [NCBI Command Line Tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).

The tool can be installed in a conda-environment as following:
```
$ conda create -n ncbi_dataset 
$ conda activate ncbi_dataset 
$ conda install -c conda-forge ncbi-datasets-cli 

```
Sequences can be downloaded using following commands.\
For example, downloading whole-genome sequences (WGS) of *P.aeruginosa* that are annotated, assembled by RefSeq and exclude atypical genomes:
```
$ datasets download genome taxon 'Pseudomonas aeruginosa' --annotated --assembly-level complete --exclude-atypical --assembly-source ‘RefSeq’ --include genome,gbff

$ unzip ncbi_dataset.zip -d /home/guest/BIT11_Traineeship/Paeruginosa_genome_gbff

```
## STEP 1 - Quality analysis of WGS
### 1.16S rRNA sequence analysis

#### 16SrRNAExtractor.py
  
- Looks in subdirectories of the input directory for ".gbff" files.\
Renames files "genomic.gbff" to [RefSeq ID].gbff.\
Extracts 16S rRNA sequences (with length 1400-1600 bp).\
Pastes sequences in "rRNAs.fa".
```
$ python 16SrRNAExtractor.py [1]
[1] = Full path to directory with .gbff files to extract 16S rRNA sequences from
```
#### ClustalO
- *This section is still under development*
- Perform multiple sequence alignment of sequences in "rRNAs.fa"

### 2. Busco
#### genomes_fasta_dir.py
- In order to run Busco in batch mode, all fasta-files need to be in one directory. \
The script looks in each subdirectory of the specified input directory for fast-files and copies them to the output directory (will be made if it does not exist yet).\
```
$ python genomes_fasta_dir.py [1] [2]
[1] = Full path to input directory (= folder containing all subfolders with downloaded fasta-files)
[2] = Full path to output directory (will make new directory if it doesn't exist yet)
```
#### Busco analysis
- Busco is a tool that assesses assembly completeness using universal single copy orthologs. \
! Attention, the tool starts looking for input and output directories relative to the current directory ! \
When the tool is installed, it can be run with following command:
```
$ busco -i [1] -o [2] -m genome -l [3] --force -c [4]
[1] = Input sequence file in fasta format
[2] = Name of output folder
[3] = Specify the name of the BUSCO lineage to be used e.g. pseudomonadales_odb10
[4] = Specify the number of threads/cores to use
```
- More information on the tool, it's installation and usage via this [link](https://busco.ezlab.org/).
  
#### busco_summary_txtfile.py 
- Makes a new text-file with custom name. \
Iterates through input directory (containing Busco output) and looks for "short_summary.json" files. \
Writes RefSeq ID and one_line_summary to the new text-file.
```
$ python busco_summary_txtfile.py [1] [2] [3]
[1] = Full path to directory with BUSCO output files
[2] = Choose an identifier to add to the file-name (e.g. Paeurginosa_1 -> busco_summary_Paeurginosa_1.txt)
[3] = Full path to output directory
```
### 3. CheckM
- *This section is still under development*

## STEP 2 - fastANI analysis
#### make_list.sh
- Looks for fast-files in the input directory and writes their full paths to the output file. \
Ouput file will be overwritten if it exists already.
```
$ bash make_list.sh [inputdirectory] [outputfile]
```
#### fastANI analysis
- fastANI computes whole-genome Average Nucleotide Identity (ANI) between genomes. \
When the tool in installed, it can be run for all genomes with one command:
```
$ fastANI --ql [1] --rl [2] -o [3] --matrix -t [4]
[1] = File containing list of query genome files, one genome per line
[2] = File containing list of reference genome files, one genome per line
[3] = Output file name
[4] = Thread count for parallel execution
```
- In addition to the standard output file, a lower triangluar matrix will be constructed from these ANI values.
- More information on the tool, it's installation and usage via this [link](https://github.com/ParBLiSS/FastANI).

## STEP 3 - MLST analysis
#### MLST_batch_run
- Runs classical MLST for multiple genomes at once.
```
$ bash MLST_batch_run.sh [1] [2] [3] [4]
[1] Full path to folder with input files (fasta) to run MLST on (files can't be in subdirectories)
[2] Specify the species (e.g. paeruginosa).\n
[3] Full path to output folder (will create if it doens't exist yet).\n
[4] Create extended output (y/n)?";
```
- Extended output includes:
"results.txt", "results_tab.tsv", "Hit_in_genome_seq.fsa", "MLST_allele_seq.fsa".
- More information on the tool, it's installation and usage via this [link](https://github.com/genomicepidemiology/mlst).

## STEP 4 - Data processing
#### ANI_matrix_excel.py
- Uses ".matrix" file from fastANI to construct a full matrix in Excel (worksheet name = target species specified as input argument). \
Uses ".json" files in the MLST output-folder to make a new worksheet named "MLST" containing genome names and accompanying STs.
```
$ python ANI_matrix_excel.py [1] [2] [3]
[1] = Full path to FastANI matrix output file
[2] = Target species (_ as delimiter, e.g. Pseudomonas_aeruginosa)
[3] = Full path to directory with JSON files from MLST
```
#### QC_ANI_matrix.py
- Opens the Exce file and looks at worksheets with ANI matrix & MLST results (specified via input-arguments).
- Looks for genomes with average %identity <95% (indicative that the genome is a different species than rest of the genomes in the dataset).\
Detected genomes and their average ANI value are added to "genomes_to_remove.txt".\
ANI & MLST worksheets are copied to new worksheets named "ANI_matrix_reviewed" & "MLST_reviewed". \
Removes genomes with doubtful species identity from these new worksheets.
```
$ python QC_ANI_matrix.py [1] [2] [3]
[1] = Full path to Excel file with ANI matrix and MLST results
[2] = Worksheet with ANI matrix in Excel file [1] e.g. Pseudomonas_aeruginosa
[3] = Worksheet with MLST results in Excel file [1] e.g. MLST
```
*! It is advised to use the "ANI_matrix_reviewed" & "MLST_reviewed" worksheets from now on !*

#### subset_matrix.py
- In the Excel file, creates new worksheets "Subset_ANI_matrix" & "Subset_MLST_results". \
Copies data from ANI & MLST worksheets (specified via input-arguments), but only includes data of genomes having certain STs (also specified via input-argument).
```
$ python subset_matrix.py [1] [2] [3] [4] ...
[1] = Full path to Excel file with ANI matrix and MLST results
[2] = Worksheet with ANI matrix in Excel file [1] e.g. Pseudomonas_aeruginosa or ANI_matrix_reviewed
[3] = Worksheet with MLST results in Excel file [1] e.g. MLST or MLST_reviewed
[4] ... = Specify one or more strain types (ST) (e.g. 262 1076 Unknown), all genomes of those STs will be selected from the Excel file
```
#### subset_min3STs.py
- In the Excel file, creates new worksheets "Subset_ANI_matrix_3STs" & "Subset_MLST_3STs". \
Copies data from ANI & MLST worksheets (specified via input-arguments), but only includes data of genomes having an STs that occurs at least 3 times in the dataset.
```
$ python subset_min3STs.py [1] [2] [3]
    [1] = Full path to Excel file with ANI matrix and MLST results
    [2] = Worksheet with ANI matrix in Excel file [1] e.g. Pseudomonas_aeruginosa or ANI_matrix_reviewed
    [3] = Worksheet with MLST results in Excel file [1] e.g. MLST or MLST_reviewed
```

## STEP 5 - Clustering ANI values
#### ANI_clustering.py
- Converts ANI %identity matrix to a distance matrix. \
Performs Single-Linkage Agglomerative Clustering of ANI (distance) values  with different values for distance_threshold (specified via input-argument).\
Makes a new worksheet "MLST_ordered" where MLST results are ordered to match the order of the genomes in the ANI matrix.\
Writes clustering results to this new worksheet (one column per tested distance_threshold value).
```
$ python ANI_clustering.py [1] [2] [3] [4] ...
[1] = Full path to Excel file with ANI matrix and MLST results
[2] = Worksheet with ANI matrix in Excel file [1] e.g. Pseudomonas_aeruginosa or ANI_matrix_reviewed
[3] = Worksheet with MLST results in Excel file [1] e.g. MLST or MLST_reviewed
[4] ... = Values to test for the distance threshold (e.g. 0.5 0.05)
```
## STEP 6 - Clustering processing and visualisation
### Comparing partitioning methods
- The [Comparing Partitions Website](http://www.comparingpartitions.info/index.php?link=Tool) can be used to line up the results of MLST analysis and Agglomerative clustering of the data. \
It calculates values such as: \
Simpson's Index of Diversity, to measure the discriminatory ability of typing systems &\
Rand Index, to give the global congruence between two typing methods.
### Visualization of ANI clusters
- The [Correlation Analysis website](https://bioit.shinyapps.io/ManiR/) can be used to make heatmaps of the ANI matrix and use MLST results as metadata. \
Since the Excel-file can be quite large, it is advised to only upload worksheets with a subset of the data (these can be constructed using the "subset_matrix.py" script). \
The heatmap can be used to compare visual clusters with clusters made by Agglomerative clustering.