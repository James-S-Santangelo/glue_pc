rule get_fourfold_zerofold:
    input:
        ref = REFERENCE_GENOME,
        gff = GFF_FILE
    output:
        '{0}/4fold_0fold/Trepens_{{site}}.bed'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/4fold_0fold/{{site}}_get_fourfold_zerofold.log'
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
        rm -rf Degeneracy ) 2> {log}
        """

