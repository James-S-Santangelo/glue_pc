rule angsd_gl:
    input:
        bams = rules.create_bam_list.output
    output:
        gls = temp('{0}/angsd_gl/full/{{chrom}}/{{chrom}}_genolike_allSamples.beagle.gz'.format(POP_STRUC_DIR)),
        mafs = temp('{0}/angsd_gl/full/{{chrom}}/{{chrom}}_genolike_allSamples.mafs.gz'.format(POP_STRUC_DIR))
    log: 'logs/angsd_gl/{chrom}_angsd_gl.log'
    conda: '../envs/population_structure.yaml'
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '06:00:00'
    shell:
        """
        angsd -GL 1 \
            -out {0}/angsd_gl/{{wildcards.chrom}}/{{wildcards.chrom}}_genolike_allSamples \
            -nThreads {{threads}} \
            -doGlf 2 \
            -doMajorMinor 1 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -minQ 20 \
            -minMapQ 30 \
            -r {{wildcards.chrom}} \
            -bam {{input.bams}} 2> {{log}}
        """.format(POP_STRUC_DIR)

rule concat_angsd_gl:
    input:
        expand(rules.angsd_gl.output.gls, chrom=CHROMOSOMES)
    output:
        '{0}/angsd_gl/full/genolike_allSamples.beagle.gz'.format(POP_STRUC_DIR)
    log: 'logs/concat_angsd_gl/concat_angsd_gl.log'
    conda: '../envs/population_structure.yaml'
    shell:
        """
        first=1
        for f in {input}; do
            if [ "$first"  ]; then
                zcat "$f"
                first=
            else
                zcat "$f"| tail -n +2
            fi
        done | bgzip -c > {output} 2> {log}
        """

rule concat_angsd_mafs:
    input:
        expand(rules.angsd_gl.output.mafs, chrom=CHROMOSOMES)
    output:
        '{0}/angsd_gl/full/genolike_allSamples.mafs.gz'.format(POP_STRUC_DIR)
    log: 'logs/concat_angsd_mafs/concat_angsd_mafs.log'
    conda: '../envs/population_structure.yaml'
    shell:
        """
        first=1
        for f in {input}; do
            if [ "$first"  ]; then
                zcat "$f"
                first=
            else
                zcat "$f"| tail -n +2
            fi
        done | bgzip -c > {output} 2> {log}
        """

rule ld_prune_angsd_gl:
    input:
        rules.concat_angsd_gl.output
    output:
        '{0}/angsd_gl/pruned/{{chrom}}/{{chrom}}_genolike_allSamples_prunned.beagle.gz'.format(POP_STRUC_DIR)
    log: 'logs/ld_prune_angsd_gl/{chrom}_ld_prune.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:ngsld_v1.1.1'
    threads: 8
    shell:
        """
        NUM_SITES=$(( $(zcat {input} | wc -l) -1 )) &&
        ngsLD --geno {input} \
            --n_ind 120 \
            --n_sites $NUM_SITES
            --n_threads {threads} \
            --out {output} 2> {log}
        """
