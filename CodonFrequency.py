# Mapping amino acids to their respective codons
codon_mapping = {
    'Phe': ['UUU', 'UUC'], 'Leu': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'Ser': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], 'Tyr': ['UAU', 'UAC'],
    'Cys': ['UGU', 'UGC'], 'Trp': ['UGG'], 'Pro': ['CCU', 'CCC', 'CCA', 'CCG'],
    'His': ['CAU', 'CAC'], 'Gln': ['CAA', 'CAG'], 'Arg': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'Ile': ['AUU', 'AUC', 'AUA'], 'Met': ['AUG'], 'Thr': ['ACU', 'ACC', 'ACA', 'ACG'],
    'Asn': ['AAU', 'AAC'], 'Lys': ['AAA', 'AAG'], 'Val': ['GUU', 'GUC', 'GUA', 'GUG'],
    'Ala': ['GCU', 'GCC', 'GCA', 'GCG'], 'Asp': ['GAU', 'GAC'], 'Glu': ['GAA', 'GAG'],
    'Gly': ['GGU', 'GGC', 'GGA', 'GGG']
}

# Dictionary for single-letter amino acid codes
amino_acid_names = {
    'F': 'Phe', 'L': 'Leu', 'S': 'Ser', 'Y': 'Tyr', 'C': 'Cys', 'W': 'Trp', 'P': 'Pro',
    'H': 'His', 'Q': 'Gln', 'R': 'Arg', 'I': 'Ile', 'M': 'Met', 'T': 'Thr', 'N': 'Asn',
    'K': 'Lys', 'V': 'Val', 'A': 'Ala', 'D': 'Asp', 'E': 'Glu', 'G': 'Gly'
}

# Function to get a valid amino acid input from the user
def get_amino_acid_input():
    while True:
        user_input = input("Input Amino Acid: ").upper()
        if all(aa in amino_acid_names for aa in user_input):
            return user_input
        else:
            print("Detected invalid amino acid code(s). Please enter valid codes.")

# Collect possible codons for the given amino acid sequence
def collect_possible_codons(amino_acid_seq):
    codons_list = []
    for aa in amino_acid_seq:
        full_name = amino_acid_names.get(aa)
        if full_name:
            codons = codon_mapping.get(full_name, [])
            codons_list.append(codons)
        else:
            print(f"Invalid amino acid code: {aa}")
    return codons_list

# Function to display possible codons for each amino acid
def display_codons(amino_acid_seq, possible_codons):
    print("\nmRNA Codon Possibilities:")
    for index, codons in enumerate(possible_codons):
        amino_acid = amino_acid_seq[index]
        full_name = amino_acid_names[amino_acid]
        print(f"Amino Acid {amino_acid} (Full Name: {full_name}):")
        print(f"Possible Codons: {', '.join(codons)}")
        print()  # Add a line break for formatting

# Function to get a valid mRNA sequence from the user
def get_mRNA_sequence():
    while True:
        mRNA_input = input("Enter the mRNA sequence: ").upper()
        if all(nucleotide in 'AUCG' for nucleotide in mRNA_input):
            return mRNA_input
        else:
            print("Invalid mRNA sequence")

# Main logic
amino_acid_sequence = get_amino_acid_input()
possible_codons = collect_possible_codons(amino_acid_sequence)
display_codons(amino_acid_sequence, possible_codons)

mRNA_sequence = get_mRNA_sequence()

# Calculate and display codon frequencies in the mRNA sequence
codon_frequencies = {}
for codons in possible_codons:
    for codon in codons:
        codon_frequencies[codon] = mRNA_sequence.count(codon)

# Display codon frequencies
print("\nCodon Frequencies in the mRNA sequence:")
for codon, count in codon_frequencies.items():
    print(f"{codon} = {count}")
