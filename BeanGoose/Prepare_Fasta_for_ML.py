# Script to prepare data for the Machine Learning analyses
# Written by Jente Ottenburghs on 2019-05-16

# Step1 - Choose scaffold and make window file
# Step2 - Extract data from Fasta-files
        # Reference genome -> anc.fa
        # Anser fabalis samples -> fabalis.fa
        # Anser serrirostris samples -> rossicus.fa

import re

############################################################
### STEP 1 - Choose scaffold and window size (automate later) ###
Scaffold_Names_File = open("Scaffolds.txt", "r")

for scaffold  in Scaffold_Names_File:
        scaffold_name = scaffold.strip() # Get rid of any whitespaces in the name

        # Extract windows for particular scaffold and save in temporary file
        Windows_file = open("BeanGoose_Simulations/dataToClassify/Window_Coordinates.bed", "r") # This file contains all the windows
        Scaffold_windows = open("BeanGoose_Simulations/dataToClassify/" + scaffold_name + ".bed", "a") # This output-file will contain only windows for particular scaffold

        for line in Windows_file:
                if scaffold_name in line:
                        Scaffold_windows.write(line)

        Windows_file.close()
        Scaffold_windows.close()

### STEP 2 - Extract data from Fasta-files ###
# Reference genome
        Reference = open("/proj/sllstore2017033/nobackup/work/jente/Reference_Genome/ansCyg.fa", 'r')
        Anc_file = open("BeanGoose_Simulations/dataToClassify/" + scaffold_name + "_anc.fa", "a")


        for line in Reference:
                if ">" in line:
                        header = line
                        if scaffold_name in header:
                                name = ">" + scaffold_name + "\n"
                                Anc_file.write(name)
                else:
                        sequence = line

                        if scaffold_name in header:
                                sequence_with_Ns = re.sub('[a-z]', 'N', sequence) # Replace lower case nucleotides with Ns
                                Anc_file.write(sequence_with_Ns)

        Reference.close()
        Anc_file.close()

# Fabalis samples
        List_of_Fabalis_samples = ["AnFa01U01","AnFa01U02","AnFa01U03","AnFa01U04","AnFa01U05","AnFa01U06","AnFa01U07","AnFa01U08","AnFa01U09"]

        for sample in List_of_Fabalis_samples:
                File_name = "/proj/sllstore2017033/nobackup/work/jente/Bean_Goose_Project/Modelling/" + sample + ".fa"

                Fabalis_in = open(File_name, "r")
                Fabalis_out = open("BeanGoose_Simulations/dataToClassify/" + scaffold_name + "_fabalis.fa", "a")

                for line in Fabalis_in:
                        if ">" in line:
                                header = line
                                if scaffold_name in header:
                                        name = ">" + sample + "\n"
                                        Fabalis_out.write(name)
                        else:
                                sequence = line

                                if scaffold_name in header:
                                        sequence_with_Ns = re.sub('[a-z]', 'N', sequence) # Replace lower case nucleotides with Ns
                                        Fabalis_out.write(sequence_with_Ns)

                Fabalis_in.close()
                Fabalis_out.close()

# Rossicus samples
        List_of_Rossicus_samples = ["AnRo01U01","AnRo01U02","AnRo01U03","AnRo01U04","AnRo01U05","AnRo01U06","AnRo01U07","AnRo01F01","AnRo01M01"]

        for sample in List_of_Rossicus_samples:
                File_name = "/proj/sllstore2017033/nobackup/work/jente/Bean_Goose_Project/Modelling/" + sample + ".fa"

                Rossicus_in = open(File_name, "r")
                Rossicus_out = open("BeanGoose_Simulations/dataToClassify/" + scaffold_name + "_rossicus.fa", "a")

                for line in Rossicus_in:
                        if ">" in line:
                                header = line
                                if scaffold_name in header:
                                        name = ">" + sample + "\n"
                                        Rossicus_out.write(name)
                        else:
                                sequence = line

                                if scaffold_name in header:
                                        sequence_with_Ns = re.sub('[a-z]', 'N', sequence) # Replace lower case nucleotides with Ns
                                        Rossicus_out.write(sequence_with_Ns)

                Rossicus_in.close()
                Rossicus_out.close()

Scaffold_Names_File.close(
