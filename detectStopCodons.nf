process detectStopCodons {

    input:
    path snps_consensus_fasta

    output:
    path("${snps_consensus_fasta.baseName}_premature_stop_codons.txt")

    script:
    """
    python3 << 'EOF'
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable

def detect_premature_stop_codons(consensus_fasta, output_file, input_format, stop_codon_symbol):
    alignment = SeqIO.parse(consensus_fasta, input_format)
    
    premature_stop_codons = []
    
    for record in alignment:
        sequence = str(record.seq)
        gene_name = record.id
        
        # Translate the sequence
        protein_sequence = str(Seq(sequence).translate())
        
        stop_positions = []
        
        for i, aa in enumerate(protein_sequence):
            if aa == stop_codon_symbol and i < len(protein_sequence) - 1:  
                stop_codon_position = (i + 1) * 3
                stop_positions.append(stop_codon_position)
        
        if stop_positions:
            for pos in stop_positions:
                premature_stop_codons.append((gene_name, pos, sequence))
    
    # Write all detected premature stop codons to the output file
    with open(output_file, 'w') as out_file:
        out_file.write("Gene_Name\\tStop_Codon_Position\\tOriginal_Sequence\\n")
        for gene, pos, seq in premature_stop_codons:
            out_file.write(f"{gene}\\t{pos}\\t{seq}\\n")
    
    if premature_stop_codons:
        print(f"Premature stop codons detected in {len(premature_stop_codons)} genes. Details written to {output_file}.")
    else:
        print("No premature stop codons detected.")

consensus_fasta = "${snps_consensus_fasta}"
output_file = "${snps_consensus_fasta.baseName}_premature_stop_codons.txt"
input_format = "${params.detectStopCodons.input_format}"
stop_codon_symbol = "${params.detectStopCodons.stop_codon_symbol}"
detect_premature_stop_codons(consensus_fasta, output_file, input_format, stop_codon_symbol)
EOF
    """
}

