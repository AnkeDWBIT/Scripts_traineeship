#!/usr/bin/env python
#This script will copy the fasta files downloaded with NCBI Datasets command line tool to a directory without subfolders
#When all the genomic-fasta files are in this new folder, it can be used as input to run BUSCO in batch mode
import os
import shutil
# Specify the input directory (=folder containing all subfolders with downloaded genomic-fasta files)
#input_directory = input("Specify the input directory (folder containing all genomic-fasta files): ")
input_dir = "/home/guest/BIT11_Traineeship/01_Paeruginosa_refseq_genomes/ncbi_dataset/data/"

#Make a new folder as output-folder to store the genomic fasta-files together
#new_output_dir = input("Full path and name of the new directory: ")
new_output_dir = "/home/guest/BIT11_Traineeship/02_QC_genomes/BUSCO/Paeruginosa_genomes_fasta/"
if not os.path.exists(new_output_dir):
    os.makedirs(new_output_dir)

#Look for the genomic fasta-files and store them in the new output directory
# Iterate over subdirectories
for subdir in os.listdir(input_dir):
	subdir_path = os.path.join(input_dir, subdir)
	if os.path.isdir(subdir_path):
		# Iterate over files in each subdirectory
		for filename in os.listdir(subdir_path):
			if ".fna" in filename:
				# Construct the new filename with the subdirectory name
				input_dir_path = os.path.join(subdir_path, filename)
				shutil.copy(input_dir_path, new_output_dir)