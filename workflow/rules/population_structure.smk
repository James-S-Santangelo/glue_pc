rule angsd_gl:
    input:
        bams = rules.create_bam_list.output
    output:
        gls = temp('{0}/angsd_gl/full/{{chrom}}/{{chrom}}_genolike_allSamples.beagle.gz'.format(POP_STRUC_DIR)),
        mafs = temp('{0}/angsd_gl/full/{{chrom}}/{{chrom}}_genolike_allSamples.mafs.gz'.format(POP_STRUC_DIR)),
        saf = temp('{0}/angsd_gl/full/{{chrom}}/{{chrom}}_genolike_allSamples.saf.gz'.format(POP_STRUC_DIR)),
        saf_idx = temp('{0}/angsd_gl/full/{{chrom}}/{{chrom}}_genolike_allSamples.saf.idx'.format(POP_STRUC_DIR)),
        saf_pos = temp('{0}/angsd_gl/full/{{chrom}}/{{chrom}}_genolike_allSamples.saf.pos.gz'.format(POP_STRUC_DIR))
    log: 'logs/angsd_gl/{chrom}_angsd_gl.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933' 
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '06:00:00'
    shell:
        """
        angsd -GL 1 \
            -out {0}/angsd_gl/full/{{wildcards.chrom}}/{{wildcards.chrom}}_genolike_allSamples \
            -nThreads {{threads}} \
            -doGlf 2 \
            -doMajorMinor 1 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -doCounts 1 \
            -setMinDepthInd 3 \
            -setMaxDepth 4500 \
            -baq 2 \
            -ref {1} \
            -minInd 96 \
            -minQ 20 \
            -minMapQ 30 \
            -doSaf 1 \
            -anc {1} \
            -r {{wildcards.chrom}} \
            -bam {{input.bams}} 2> {{log}}
        """.format(POP_STRUC_DIR, REFERENCE_GENOME)

rule merge_safs:
    input:
        expand(rules.angsd_gl.output.saf_idx, chrom=CHROMOSOMES)
    output:
        saf = '{0}/angsd_gl/full/genolike_allSamples.saf.gz'.format(POP_STRUC_DIR),
        saf_idx = '{0}/angsd_gl/full/genolike_allSamples.saf.idx'.format(POP_STRUC_DIR),
        saf_pos = '{0}/angsd_gl/full/genolike_allSamples.saf.pos.gz'.format(POP_STRUC_DIR)
    log: 'logs/merge_safs/merge_safs.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933' 
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '02:00:00'
    shell:
        """
        realSFS cat {{input}} \
            -outnames {0}/angsd_gl/full/genolike_allSamples 2> {{log}}
        """.format(POP_STRUC_DIR)

rule concat_angsd_gl:
    input:
        expand(rules.angsd_gl.output.gls, chrom=CHROMOSOMES)
    output:
        '{0}/angsd_gl/full/genolike_allSamples.beagle.gz'.format(POP_STRUC_DIR)
    log: 'logs/concat_angsd_gl/concat_angsd_gl.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933' 
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
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933' 
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

rule create_pos_file_for_ngsLD:
    input:
        rules.angsd_gl.output.mafs
    output:
        '{0}/angsd_gl/full/{{chrom}}/{{chrom}}_angsdGL.pos'.format(POP_STRUC_DIR)
    log: 'logs/create_pos_file_for_ngsLD/{chrom}_pos.log'
    shell:
        """
        zcat {input} | cut -f 1,2 | tail -n +2 > {output} 2> {log}
        """

rule calc_ld_angsd_gl:
    input:
        pos = rules.create_pos_file_for_ngsLD.output,
        gls = rules.angsd_gl.output.gls
    output:
        '{0}/angsd_gl/ld/{{chrom}}/{{chrom}}_genolike_allSamples.ld.gz'.format(POP_STRUC_DIR)
    log: 'logs/calc_ld_angsd_gl/{chrom}_calc_ld.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:ngsld_v1.1.1'
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '48:00:00'
    shell:
        """
        ( NUM_SITES=$(cat {input.pos} | wc -l) &&
          ngsLD --geno {input.gls} \
            --pos {input.pos} \
            --n_ind 120 \
            --n_sites $NUM_SITES \
            --probs \
            --n_threads {threads} \
            --max_kb_dist 25 | gzip --best > {output} ) 2> {log}
        """
