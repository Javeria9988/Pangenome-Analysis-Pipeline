process genedescription {
    input:
    path gene_presence_absence_csv
    path detect_stop_codon_txt

    output:
    path "${detect_stop_codon_txt.baseName}_gene_description.txt"

    script:
    """
    python3 << EOF
import pandas as pd

# Load gene annotations from the gene_presence_absence.csv
gene_annotations = pd.read_csv("${gene_presence_absence_csv}")

# Replace empty annotations with 'No annotation'
gene_annotations['Annotation'] = gene_annotations['Annotation'].fillna('No annotation')

# Load stop codon file
stop_codons = pd.read_csv("${detect_stop_codon_txt}", sep='\\t', header=None, names=['Gene_Name', 'Stop_Codon_Position', 'Original_Sequence'])

# Merge the two dataframes on the gene name (assuming the first column in stop_codons is the gene name)
merged_data = pd.merge(stop_codons, gene_annotations[['Gene', 'Annotation']], left_on='Gene_Name', right_on='Gene', how='left')

# Replace NaN values in annotations with 'No annotation'
merged_data['Annotation'] = merged_data['Annotation'].fillna('No annotation')

# Create a new column with the desired output format
merged_data['Output_Line'] = merged_data['Gene_Name'] + "   " + merged_data['Stop_Codon_Position'].astype(str) + "  " + merged_data['Annotation'] + " " + merged_data['Original_Sequence']

# Select only the Output_Line column to save the final result
final_output = merged_data['Output_Line']

# Write the result to a new file
final_output.to_csv("${detect_stop_codon_txt.baseName}_gene_description.txt", sep='\\t', index=False, header=False)

# Now read the file and replace 'No annotation' with 'Product' in the first line if present
with open("${detect_stop_codon_txt.baseName}_gene_description.txt", 'r') as file:
    lines = file.readlines()

# Replace 'No annotation' with 'Product' in the first line
if lines:
    lines[0] = lines[0].replace('No annotation', 'Product')

# Write the modified lines back to the file
with open("${detect_stop_codon_txt.baseName}_gene_description.txt", 'w') as file:
    file.writelines(lines)

EOF
    """
}

