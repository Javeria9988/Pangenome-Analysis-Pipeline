process annotate {

    input:
    path contigfile

    output:
    path "${contigfile.simpleName}.gff", emit: gff_files
    path "${contigfile.simpleName}.txt", emit: txt_files

    script:
    """
    prokka --prefix ${contigfile.simpleName} --outdir ${contigfile.simpleName} ${contigfile}

    # Copy the generated GFF file to the desired location
    cp ${contigfile.simpleName}/${contigfile.simpleName}.gff ${contigfile.simpleName}.gff
    
    # Copy the generated TXT file to the desired location
    cp ${contigfile.simpleName}/${contigfile.simpleName}.txt ${contigfile.simpleName}.txt
    """
}

