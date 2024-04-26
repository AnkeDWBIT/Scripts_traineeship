#!/usr/bin/env python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os


output_handle=open("rRNAs.fa","w")

# Specify the input directory
input_directory = "/home/guest/BIT11_Traineeship/Paeruginosa_refseq_genomes/ncbi_dataset/data"
sixteen_s_records=[]
for gbff_file in os.listdir(input_directory):
	if os.path.isfile(os.path.join(input_directory, gbff_file)):
		name = os.path.splitext(os.path.basename(gbff_file))[0]
		gbff_file = os.path.join(input_directory, gbff_file)
		with open(gbff_file, 'r') as gbff:
			for record in SeqIO.parse(gbff, 'genbank'):
				for feature in record.features:
					if(feature.type == "rRNA"):
						if '16S ribosomal RNA' in feature.qualifiers['product'][0]:
							sequence = feature.extract(record.seq)
							if (len(sequence) <= 1604) and (len(sequence) > 1451):
								sixteen_s_records.append(SeqRecord(sequence, id=f"{name}_16S_{record.id}", description=""))

SeqIO.write(sixteen_s_records, output_handle, "fasta")
output_handle.close()
