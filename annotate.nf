process annotate {

    input:
    path contigfile

    output:
    path "${contigfile.simpleName}.gff", emit: gff_files

    script:
    """
    prokka --kingdom ${params.annotate.kingdom} --genus ${params.annotate.genus} --species ${params.annotate.species} --prefix ${contigfile.simpleName} --outdir ${contigfile.simpleName} ${contigfile}
    cp ${contigfile.simpleName}/${contigfile.simpleName}.gff ${contigfile.simpleName}.gff
    """
}

