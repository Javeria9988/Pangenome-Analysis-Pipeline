# Pangenome-Analysis-Pipeline

## Overview
This pipeline performs a comprehensive pangenome analysis, starting from genome annotation, through pangenome identification, fastq reads generation, SNP calling, and finally, detecting premature stop codons and addition of gene description. The pipeline leverages various bioinformatics tools, integrated within a Nextflow workflow to streamline the analysis of multiple sequences.

## Pipeline Structure
The pipeline consists of the following steps:

1. **Annotate:** Annotates the input sequences using prokka.
2. **Panaroo:** Identifies core and accessory genes across the annotated genomes and generates relevant outputs.
3. **DWGSIM:** Generates error-free fastq reads from the input sequences.
4. **Snippy:** Aligns sequences and calls SNPs based on the simulated reads.
5. **Detect Stop Codons:** Identifies premature stop codons from the SNP data using python code.
6. **GeneDescription:** Adds gene product details using python code.

## Workflow
![rb2 workflow](https://github.com/user-attachments/assets/8b6bc41c-3dd7-4b7f-9860-4ee32742dc1b)

## Installation

### Prerequisites
- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/) for containerized execution of processes.

### Setup
Clone the repository and navigate to the directory:
```bash
git clone https://github.com/Javeria9988/pangenome-pipeline.git
cd pangenome-pipeline
```
### Directory structure for running pipeline
![Screenshot from 2024-08-20 21-21-11](https://github.com/user-attachments/assets/b1d4c0b9-df32-4e07-902d-e387013f464b)

### Usage
### Running the Pipeline
To run the pipeline, use the following command:
nextflow run main.nf -c nextflow.config -resume

### Parameters
params.sequences: Path to the input sequences in .fa format. The default is set to 'sequences/*.fa'.

### Input Data
Sequences: The pipeline expects input sequences in .fa format, located in the directory specified by the params.sequences parameter.

### Output Data
Annotated Genomes: GFF files generated by the annotation step.
Panaroo Outputs: Core and accessory gene information in CSV and FASTA formats.
fastq Reads generation: Error-free reads generated by DWGSIM.
SNPs: SNP calling results from Snippy.
Stop Codon Detection: A report detailing any detected premature stop codons in each sample.
GeneDescription: Adds protein details to output files generated in above step.

## Docker Integration
This pipeline utilizes Docker containers to ensure reproducibility and consistency across different environments. Each step of the pipeline is associated with a specific Docker container, as outlined below:

### Containers Used:

1. **Annotate:**
   - **Container:** `staphb/prokka:latest`
   - **Description:** This container includes Prokka, a tool used for rapid annotation of prokaryotic genomes.

2. **Panaroo:**
   - **Container:** `staphb/panaroo:latest`
   - **Description:** This container includes Panaroo, a tool for pangenome analysis, particularly focused on bacterial genomes.

3. **DWGSIM:**
   - **Container:** `dwgsim-with-ps:latest`
   - **Description:** This customized container includes DWGSIM, a tool for getting error-free reads from genomic FASTA sequences and ps, which is required for genration of nextflow report.

4. **Snippy:**
   - **Container:** `staphb/snippy:4.6.0`
   - **Description:** This container includes Snippy, a tool for rapid bacterial SNP calling and variant detection.

5. **Detect Stop Codons:**
   - **Container:** `python-biopython-pandas-ps:latest`
   - **Description:** This customized container includes Python along with the ps, Biopython and Pandas libraries, used for detecting premature stop codons and processing genomic data.

6. **Gene Description:**
   - **Container:** `python-biopython-pandas-ps:latest`
   - **Description:** This customized container includes Python along with the ps, Biopython and Pandas libraries, used for adding gene product details from panaroo generated gene_presence_absence.csv file.

### Usage
To ensure that Docker is enabled, the pipeline is configured with Docker integration:

```groovy
docker {
    enabled = true
}
