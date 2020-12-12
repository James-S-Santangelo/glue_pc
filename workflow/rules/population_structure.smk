rule angsd_gl:
    input:
        bams = rules.create_bam_list.output
    output:
        '{0}/angsd_gl/{{chrom}}/{{chrom}}_genolike_allSamples_maf0.05.beagle.gz'.format(POP_STRUC_DIR),
        '{0}/angsd_gl/{{chrom}}/{{chrom}}_genolike_allSamples_maf0.05.mafs.gz'.format(POP_STRUC_DIR)
    log: 'logs/angsd_gl/{chrom}_angsd_gl.log'
    conda: '../envs/population_structure.yaml'
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 5000,
        time = '06:00:00'
    shell:
        """
        angsd -GL 1 \
            -out {0}/angsd_gl/{{wildcards.chrom}}/{{wildcards.chrom}}_genolike_allSamples_maf0.05 \
            -nThreads {{threads}} \
            -doGlf 2 \
            -doMajorMinor 1 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -minQ 20 \
            -minMapQ 30 \
            -minMaf 0.05 \
            -r {{wildcards.chrom}} \
            -bam {{input.bams}} 2> {{log}}
        """.format(POP_STRUC_DIR)
