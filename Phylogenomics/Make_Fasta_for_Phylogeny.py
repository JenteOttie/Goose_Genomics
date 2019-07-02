# Script to make fasta files for phylogeny estimation.
# Written by Jente Ottenburghs on 2019-07-02
# Workflow:
# 1. Determine window to extract (will be randomized)
# 2. Locate scaffold and save sequence
# 3. Extract window from sequence
# 4. Check number of Ns in extracted windows - if too high (threshold), discard
# 5. Save window to fasta-file

##############################################

# Determine window to extract
scaffold = "NW_013185654.1"
start = 1000
end = 2000

# Set threshold for number of Ns (in percentage)
N_threshold = 25

# Loop through list of samples and extract the window
List_of_samples = ["AnFa01U01","AnFa01U02","AnFa01U03","AnFa01U04","AnFa01U05","AnFa01U06","AnFa01U07","AnFa01U08","AnFa01U09",\
                "AnRo01U01","AnRo01U02","AnRo01U03","AnRo01U04","AnRo01U05","AnRo01U06","AnRo01U07","AnRo01F01","AnRo01M01"]

for sample in List_of_samples:
        # Get file with genome sequence
        File_name = "/proj/sllstore2017033/nobackup/work/jente/Bean_Goose_Project/Modelling/" + sample + ".fa"
        print("Analyzing file: " + File_name)

        # Set input-file and output-file
        Input_file = open(File_name, "r")
        Output_file = open("Raw_Sequences.fa", "a")

        sequence = ""
        for line in Input_file:
                if ">" in line:
                        header = line
                        if scaffold in header:
                                name = ">" + sample + "\n"
                                Output_file.write(name)
                else:
                        if scaffold in header:
                                sequence = sequence + line.strip() # Save DNA-sequence in string while header equals scaffold
        print("Saved scaffold " + scaffold + " sequence. Now extracting window.\n")

        # Extract specified window from DNA-sequence
        window = sequence[start:end] + "\n"

        # Count number of Ns and discard window if bigger than threshold
        Number_of_Ns = window.count("n")
        Percentage_of_Ns = Number_of_Ns / float(len(window)) *100
        print("The percentage of Ns in this string is " + str(round(Percentage_of_Ns,2)) + "%")

        if Percentage_of_Ns > N_threshold:
                print("Too many Ns, discarding sequence")
        else:
                print("Saving sequence\n***************\n")
                Output_file.write(window)

        Input_file.close()
        Output_file.close()
