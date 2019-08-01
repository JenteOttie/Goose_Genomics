# Script to extract genes that overlap with list of windows
# Written by Jente Ottenburghs on 2018-08-16

import sys

# Input files
# 1 = shared windows
# 2 = GFF-file with gene coordinates

# Open input and output-file
GFF_file = open(sys.argv[2], "r")
output = open("Genes_Shared.txt", "a")

# Put scaffold names into a list. This will later be used to compare with the gene file.
Window_file = open(sys.argv[1], "r")
scaffold_names = []
for line in Window_file:
	data = line.split()
	scaffold_names.append(data[0])
Window_file.close()

# Check which genes are located in a window
for line in GFF_file:
	if "#" not in line:
		contents = line.split()
		
		# Only use lines of the class "gene"
		if contents[2] == "gene":
			scaffold = contents[0]
			start = contents[3]

			if scaffold in scaffold_names:
				# If the scaffold is in the list of scaffold, open the file and see if coordinates match
				Window_file = open(sys.argv[1], "r")

				# Check if start of gene is within the window
				for window in Window_file:
					window_start = window.split()[1]
					window_end = window.split()[2]

					# Write to output-file
					if start > window_start and start < window_end and scaffold == window.split()[0]:
						output.write(line)

output.close()
GFF_file.close()

# Get gene names for GO-analysis
Genes = open("Genes_Shared.txt", "r")
Gene_names = open("Genes_Shared_Names.txt", "a")

for line in Genes:
	name = line.replace(";","\t").replace("=","\t").split()[13]
	if "LOC" in name:
		next
	else:
		Gene_names.write(name + '\n')

Genes.close()
Gene_names.close()
