from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

def extract_gene_sequence(gff_file, fasta_file, gene_name, output_file):
    """
    Extract the sequence of a specific gene from a FASTA file using a GFF file.
    
    Args:
    gff_file (str): Path to the GFF file.
    fasta_file (str): Path to the FASTA file.
    gene_name (str): The name of the gene to extract.
    output_file (str): Path to the output file where the extracted sequence will be saved.
    """
    # Parse the GFF file
    gff_data = pd.read_csv(gff_file, sep='\t', comment='#', header=None,
                           names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    
    # Drop rows with missing 'attributes' to avoid errors
    gff_data = gff_data.dropna(subset=['attributes'])
    
    # Extract the gene information
    gene_info = gff_data[gff_data['attributes'].str.contains(f"Name={gene_name}", na=False)]
    
    if gene_info.empty:
        print(f"Gene {gene_name} not found in the GFF file.")
        return
    
    # Assuming there is only one occurrence of the gene
    gene_info = gene_info.iloc[0]
    seqid = gene_info['seqid']
    start = int(gene_info['start']) - 1  # Convert to 0-based index for Python
    end = int(gene_info['end'])
    strand = gene_info['strand']
    
    # Read the FASTA file and extract the sequence
    fasta_sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    
    if seqid not in fasta_sequences:
        print(f"Sequence ID {seqid} not found in the FASTA file.")
        return
    
    sequence = fasta_sequences[seqid].seq[start:end]
    
    # If the gene is on the negative strand, reverse complement the sequence
    if strand == '-':
        sequence = sequence.reverse_complement()
    
    # Write the extracted sequence to the output file
    with open(output_file, 'w') as output_handle:
        output_handle.write(f">{gene_name} extracted from {seqid}:{start+1}-{end} ({strand} strand)\n")
        output_handle.write(str(sequence) + "\n")
    
    print(f"Sequence for gene {gene_name} has been extracted and written to {output_file}.")

# Example usage
gff_file = "/home/javeria/Downloads/onegenemutation/17428_844.gff"  # Path to your GFF file
fasta_file = "/home/javeria/Downloads/onegenemutation/17428_844.mutated.fa"  # Path to your FASTA file
gene_name = "gmk"  # Name of the gene to extract
output_file = "/home/javeria/Downloads/onegenemutation/17428_844_output_sequence.fa"  # Path for your output file

extract_gene_sequence(gff_file, fasta_file, gene_name, output_file)

