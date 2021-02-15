rule angsd_single_sample_saf:
    input:
        bam = rules.samtools_markdup.output.bam,
        index = rules.index_bam.output,
        ref = REFERENCE_GENOME
    output:
        saf = '{0}/{{sample}}/{{sample}}_CM019101.1.saf.gz'.format(SS_SFS_DIR),
        saf_pos = '{0}/{{sample}}/{{sample}}_CM019101.1.saf.pos.gz'.format(SS_SFS_DIR),
        saf_idx = '{0}/{{sample}}/{{sample}}_CM019101.1.saf.idx'.format(SS_SFS_DIR)
    log: 'logs/angsd_single_sample_saf/{sample}_ss_sfs.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    params:
        out = '{0}/{{sample}}/{{sample}}_CM019101.1'.format(SS_SFS_DIR)
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doSaf 1 \
            -anc {input.ref} \
            -r CM019101.1 \
            -i {input.bam} 2> {log}
        """

rule angsd_estimate_sfs_single_sample:
    input:
        rules.angsd_single_sample_saf.output.saf_idx
    output:
        '{0}/{{sample}}/{{sample}}_CM019101.1.sfs'.format(SS_SFS_DIR)
    log: 'logs/angsd_estimate_sfs_single_sample/{sample}_sfs.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        realSFS {input} -P {threads} > {output} 2> {log}
        """

rule single_sample_sfs_done:
    input:
        expand(rules.angsd_estimate_sfs_single_sample.output, sample=SAMPLES)
    output:
        '{0}/single_sample_sfs.done'.format(SS_SFS_DIR)
    shell:
        """
        touch {output}
        """




