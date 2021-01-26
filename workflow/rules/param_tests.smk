rule test_angsd_baq_GL:
    input:
        bams = rules.create_bam_list.output
    output:
        saf = '{0}/test_params/GL{{GL}}_baq{{baq}}/CM019112.1_allSamples_allSites_GL{{GL}}_baq{{baq}}.saf.gz'.format(ANGSD_DIR),
        saf_idx = '{0}/test_params/GL{{GL}}_baq{{baq}}/CM019112.1_allSamples_allSites_GL{{GL}}_baq{{baq}}.saf.idx'.format(ANGSD_DIR),
        saf_pos = '{0}/test_params/GL{{GL}}_baq{{baq}}/CM019112.1_allSamples_allSites_GL{{GL}}_baq{{baq}}.saf.pos.gz'.format(ANGSD_DIR)
    log: 'logs/test_angsd_baq_GL/CM019112.1_GL{GL}_baq{baq}.log'
    conda: '../envs/angsd.yaml'
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE,
        time = '12:00:00'
    shell:
        """
        angsd -GL {{wildcards.GL}} \
            -out {0}/test_params/GL{{wildcards.GL}}_baq{{wildcards.baq}}/CM019112.1_allSamples_allSites_GL{{wildcards.GL}}_baq{{wildcards.baq}} \
            -nThreads {{resources.ntasks}} \
            -doCounts 1 \
            -dumpCounts 2 \
            -setMinDepthInd 1 \
            -setMaxDepth 4500 \
            -baq {{wildcards.baq}} \
            -ref {1} \
            -minInd 60 \
            -minQ 20 \
            -minMapQ 30 \
            -doSaf 1 \
            -anc {1} \
            -r CM019112.1 \
            -bam {{input.bams}} 2> {{log}}
        """.format(ANGSD_DIR, REFERENCE_GENOME)

rule test_angsd_sfs_baq_GL:
    input:
        rules.test_angsd_baq_GL.output.saf_idx 
    output:
        '{0}/test_params/GL{{GL}}_baq{{baq}}/CM019112.1_allSamples_allSites_GL{{GL}}_baq{{baq}}.sfs'.format(ANGSD_DIR)
    log: 'logs/sfs_allSites/CM019112.1_GL{GL}_baq{baq}.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 25000,
        time = '06:00:00'
    shell:
        """
        realSFS {input} -P {threads} -fold 1 > {output} 2> {log}
        """
