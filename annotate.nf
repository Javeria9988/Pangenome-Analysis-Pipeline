process annotate {

    input:
    path contigfile

    output:
    path "${contigfile.simpleName}.gff", emit: gff_files

    script:
    """
    prokka --prefix ${contigfile.simpleName} --outdir ${contigfile.simpleName} ${contigfile}
    cp ${contigfile.simpleName}/${contigfile.simpleName}.gff ${contigfile.simpleName}.gff
    """
}

