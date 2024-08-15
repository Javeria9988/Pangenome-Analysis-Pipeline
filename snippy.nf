process snippy {

    input:
    path reference_file
    tuple val(sample_name), path(fq1), path(fq2)

    output:
    path "snippy_output/${sample_name}/${sample_name}_snps.consensus.fa"

    script:
    """
    snippy --outdir snippy_output/${sample_name} \
           --ref ${reference_file} \
           --R1 ${fq1} \
           --R2 ${fq2} \
           --mincov ${params.snippy.mincov} \
           --minfrac ${params.snippy.minfrac} \
           --maxsoft ${params.snippy.maxsoft} \
           --minqual ${params.snippy.minqual} \
           --cpus ${params.snippy.cpus}
    
    # Rename the snps.consensus.fa file to include the sample name
    mv snippy_output/${sample_name}/snps.consensus.fa snippy_output/${sample_name}/${sample_name}_snps.consensus.fa
    """
}

