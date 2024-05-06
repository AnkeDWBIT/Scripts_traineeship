#!/usr/bin/env python
# Script to copy the fasta files downloaded with NCBI Datasets command line tool to a directory without subfolders
# Usage example : need all genomic-fasta files directly in one folder to run BUSCO (tool) in batch mode

import sys
import os
import shutil

# STEP 1 : Checking the command line arguments
############################################################################################################
# Check if any command-line arguments have been provided
if len(sys.argv) < 2:
	print("Error: Need to provide command-line arguments.")
	print("Usage: python scriptname.py [1] [2]")
	print("\t[1] = Full path to input directory -> folder containing all subfolders with downloaded genomic-fasta files e.g. .../ncbi_dataset/data/)")
	print("\t[2] = Full path to output directory (will make new directory if it doesn't exist et)")
	sys.exit(1)

# Store the command-line argument(s) in an object
input_dir = sys.argv[1]
new_output_dir = sys.argv[2]

# Specify the input directory (=folder containing all subfolders with downloaded genomic-fasta files)
#input_directory = input("Specify the input directory (folder containing all genomic-fasta files): ")
#input_dir = "/home/guest/BIT11_Traineeship/01_Paeruginosa_refseq_genomes/ncbi_dataset/data/"


# STEP 2 : Making a new directory if specified output directory doesn't exist yet
############################################################################################################
#new_output_dir = input("Full path and name of the new directory: ")
#new_output_dir = "/home/guest/BIT11_Traineeship/02_QC_genomes/BUSCO/Paeruginosa_genomes_fasta/"
if not os.path.exists(new_output_dir):
    os.makedirs(new_output_dir)

# STEP 1 : Look for the genomic fasta-files and store them in the new output directory
############################################################################################################
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