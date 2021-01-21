rule get_fourfold_zerofold:
    input:
        ref = REFERENCE_GENOME,
        gff = GFF_FILE
    output:
        four = '{0}/4fold_0fold/Trepens_4fold.bed.gz'.format(PROGRAM_RESOURCE_DIR),
        zero = '{0}/4fold_0fold/Trepens_0fold.bed.gz'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/4fold_0fold/get_fourfold_zerofold.log'
    conda: '../envs/ref.yaml'
    params:
        outpath = '{0}/4fold_0fold/'.format(PROGRAM_RESOURCE_DIR)
    resources:
        mem_mb = 4000,
        time = '04:00:00'
    shell:
        """
        ( git clone https://github.com/James-S-Santangelo/Degeneracy.git &&
        cd Degeneracy &&
        get_4fold_sites.sh {input.gff} {input.ref} {params.outpath} &&
        gzip {output.zero} && {output.four} &&
        rm -rf Degeneracy ) 2> {log}
        """

