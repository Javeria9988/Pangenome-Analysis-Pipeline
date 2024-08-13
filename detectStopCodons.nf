process detectStopCodons {

    input:
    path alignment_file

    output:
    path("${alignment_file.baseName}_premature_stop_codons.txt")

    publishDir "detectStopCodons_results", mode: 'copy'

    script:
    """
    python3 << 'EOF'
from Bio import SeqIO
from Bio.Seq import Seq

def detect_premature_stop_codons(alignment_file, output_file, input_format, stop_codon_symbol):
    alignment = SeqIO.parse(alignment_file, input_format)
    
    premature_stop_codons = []
    
    for record in alignment:
        sequence = str(record.seq).replace("-", "").replace("N", "") 
        gene_name = record.id  
        
        protein_sequence = Seq(sequence).translate(to_stop=False)
        
        stop_codon_index = protein_sequence.find(stop_codon_symbol)
        
        if 0 <= stop_codon_index < len(protein_sequence) - 1:
            stop_codon_position = (stop_codon_index + 1) * 3
            premature_stop_codons.append((gene_name, stop_codon_position, sequence))
    
    with open(output_file, 'w') as out_file:
        out_file.write("Gene_Name\\tStop_Codon_Position\\tOriginal_Sequence\\n")
        for gene, pos, seq in premature_stop_codons:
            out_file.write(f"{gene}\\t{pos}\\t{seq}\\n")
    
    if premature_stop_codons:
        print(f"Premature stop codons detected in {len(premature_stop_codons)} genes. Details written to {output_file}.")
    else:
        print("No premature stop codons detected.")

alignment_file = "${alignment_file}"
output_file = "${alignment_file.baseName}_premature_stop_codons.txt"
input_format = "${params.detectStopCodons.input_format}"
stop_codon_symbol = "${params.detectStopCodons.stop_codon_symbol}"
detect_premature_stop_codons(alignment_file, output_file, input_format, stop_codon_symbol)
EOF
    """
}

