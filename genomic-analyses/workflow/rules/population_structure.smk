rule pcangsd:
    input:
        rules.concat_angsd_gl.output
    output:
        '{0}/pcangsd/allSamples_{{site}}_maf{{maf}}_pcangsd.cov'.format(POP_STRUC_DIR),
    log: 'logs/pcangsd/{site}_maf{maf}_pcangsd.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:pcangsd_v0.99'
    threads: 10
    params:
        out = '{0}/pcangsd/allSamples_{{site}}_maf{{maf}}_pcangsd'.format(POP_STRUC_DIR)
    wildcard_constraints:
        site = '4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        python3 /opt/pcangsd/pcangsd.py \
            -beagle {input} \
            -o {params.out} \
            -threads {threads} \
            > {log}
        """

rule pop_structure_done:
    input:
        expand(rules.pcangsd.output, site = '4fold', maf = '0.05')
    output:
        '{0}/population_structure.done'.format(POP_STRUC_DIR)
    shell:
        """
        touch {output}
        """
