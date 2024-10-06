# Mapping from DNA nucleotides to mRNA nucleotides
nucleotide_mapping = {
    'A': 'U',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}

# Codon table mapping mRNA codons to amino acids
codon_dictionary = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UGU': 'Cys', 'UGC': 'Cys',
    'UGG': 'Trp', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu',
    'CUG': 'Leu', 'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro',
    'CCG': 'Pro', 'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln',
    'CAG': 'Gln', 'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg',
    'CGG': 'Arg', 'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
    'AUG': 'Met', 'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr',
    'ACG': 'Thr', 'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys',
    'AAG': 'Lys', 'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg',
    'AGG': 'Arg', 'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val',
    'GUG': 'Val', 'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala',
    'GCG': 'Ala', 'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu',
    'GAG': 'Glu', 'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly',
    'GGG': 'Gly'
}

# Function to find the mRNA complement from a DNA sequence
def get_mRNA_complement(dna_seq):
    mRNA_seq = ""
    for base in dna_seq:
        if base in nucleotide_mapping:
            mRNA_seq += nucleotide_mapping[base]
        else:
            return None
    return mRNA_seq

# Function to convert mRNA sequence to a protein sequence
def translate_mRNA_to_protein(mRNA_seq):
    protein_chain = []
    for idx in range(0, len(mRNA_seq), 3):
        codon = mRNA_seq[idx:idx + 3]
        if codon in codon_dictionary:
            protein_chain.append(codon_dictionary[codon])
        else:
            protein_chain.append("Unknown")
    return protein_chain

# Main loop to get user input and process the DNA sequence
while True:
    user_input = input("Please provide the DNA sequence: ").upper()

    # Validate the DNA sequence
    if all(nucleotide in 'ATCG' for nucleotide in user_input):
        mRNA_result = get_mRNA_complement(user_input)
        if mRNA_result is not None:
            print("Resulting mRNA Sequence:", mRNA_result)

            protein_result = translate_mRNA_to_protein(mRNA_result)
            print("Corresponding Amino Acid Sequence:", " - ".join(protein_result))
            break
        else:
            print("The provided DNA sequence is invalid.")
    else:
        print("Invalid DNA Sequence")
