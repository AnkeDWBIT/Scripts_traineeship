#!/usr/bin/env python
# Script to rename .gbff files to RefSeq identifier & extract 16S rRNA sequences from .gbff files

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os


output_handle=open("rRNAs.fa","w")

# STEP 1 : Checking the command line arguments
############################################################################################################
# Check if any command-line arguments have been provided
if len(sys.argv) < 1:
    print("Error: Need to provide a command-line argument.")
    print("Usage: python scriptname.py [1]")
    print("\t[1] = Full path to directory with .gbff files to extract 16S rRNA sequences from")
    sys.exit(1)

# Store the command-line argument(s) in an object
input_directory = sys.argv[1]

# Specify the input directory
#input_directory = input("Specify the input directory: ")
#input_directory = "/home/guest/BIT11_Traineeship/01_Paeruginosa_refseq_genomes/ncbi_dataset/data/"


# STEP 2 : Renaming the .gbff files by adding the refseq accession
############################################################################################################

"""
# Iterate over subdirectories
for subdir in os.listdir(input_directory):
	subdir_path = os.path.join(input_directory, subdir)
	print(subdir_path)
	if os.path.isdir(subdir_path):
		# Iterate over files in each subdirectory
		for filename in os.listdir(subdir_path):
			if ".gbff" in filename:
				# Construct the new filename with the subdirectory name
				new_filename = os.path.join(subdir_path, subdir + ".gbff")
				# Rename the file
				os.rename(os.path.join(subdir_path, filename), new_filename)
"""

# STEP 3 : Create a file containing all 16S rRNA sequences
############################################################################################################

sixteen_s_records=[]

for root, dirs, files in os.walk(input_directory):
	for file in files:
		if file.endswith(".gbff"):
			name = os.path.splitext(os.path.basename(file))[0]
			gbff_file = os.path.join(input_directory,name,file)
			with open(gbff_file, 'r') as gbff:
				for record in SeqIO.parse(gbff, 'genbank'):
					for feature in record.features:
						if(feature.type == "rRNA"):
							if '16S ribosomal RNA' in feature.qualifiers['product'][0]:
								sequence = feature.extract(record.seq)
								if (len(sequence) <= 1600) and (len(sequence) > 1400):
									sixteen_s_records.append(SeqRecord(sequence, id=f"{name}_16S_{record.id}", description=""))

SeqIO.write(sixteen_s_records, output_handle, "fasta")
output_handle.close()