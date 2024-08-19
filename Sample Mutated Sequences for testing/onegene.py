from Bio import SeqIO
from Bio.Seq import Seq

def introduce_stop_codon(sequence, position):
    """
    Introduce a premature stop codon (TAA) at the specified position.
    """
    return sequence[:position-1] + 'TAA' + sequence[position+2:]

def process_fasta(input_fasta, output_fasta, gene_id, start, position_relative_to_gene_start):
    # Read the original sequences from the fasta file
    sequences = SeqIO.parse(input_fasta, "fasta")

    # Open the output file to write the mutated sequence
    with open(output_fasta, "w") as output_handle:
        for seq_record in sequences:
            if seq_record.id == gene_id:
                # Calculate the exact position for the stop codon
                absolute_position = start + position_relative_to_gene_start - 1
                # Introduce the stop codon in the sequence
                mutated_sequence = introduce_stop_codon(str(seq_record.seq), absolute_position)
                # Create a new Seq object with the mutated sequence
                seq_record.seq = Seq(mutated_sequence)
                # Write the mutated sequence to the output file
                SeqIO.write(seq_record, output_handle, "fasta")
            else:
                # Write the original sequence to the output file
                SeqIO.write(seq_record, output_handle, "fasta")

# File paths
input_fasta = "/home/javeria/Downloads/onegenemutation/17428_844.contigs_velvet.fa"
output_fasta = "/home/javeria/Downloads/onegenemutation/17428_844.mutated.fa"
gene_id = ".17428_8_44.1"  # ID of the sequence containing the gene
start = 41402 # Start position of the gene
position_relative_to_gene_start = 52  # Position to introduce the premature stop codon, relative to the start of the gene

# Process the fasta file to introduce the mutation
process_fasta(input_fasta, output_fasta, gene_id, start, position_relative_to_gene_start)
