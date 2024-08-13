process panaroo {

    input:
    path gff_files

    output:
    tuple val(gff_files), path("panaroo_output/gene_presence_absence.csv"), path("panaroo_output/pan_genome_reference.fa")

    script:
    """
    panaroo -i ${gff_files.join(' ')} \
            -o panaroo_output \
            --clean-mode ${params.panaroo.clean_mode} \
            -a ${params.panaroo.analysis_mode} \
            --aligner ${params.panaroo.aligner} \
            -t ${params.panaroo.threads}
    """
}

