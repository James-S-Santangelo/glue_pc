rule angsd_withInvar:
    input:
        bams = rules.create_bam_list.output
    output:
        saf = temp('{0}/withInvar/{{chrom}}/{{chrom}}_genolike_allSamples_withInvar.saf.gz'.format(ANGSD_DIR)),
        saf_idx = temp('{0}/withInvar/{{chrom}}/{{chrom}}_genolike_allSamples_withInvar.saf.idx'.format(ANGSD_DIR)),
        saf_pos = temp('{0}/withInvar/{{chrom}}/{{chrom}}_genolike_allSamples_withInvar.saf.pos.gz'.format(ANGSD_DIR))
    log: 'logs/angsd_withInvar/{chrom}_angsd_withInvar.log'
    conda: '../envs/angsd.yaml'
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE * 2,
        time = '12:00:00'
    shell:
        """
        angsd -GL 1 \
            -out {0}/withInvar/{{wildcards.chrom}}/{{wildcards.chrom}}_genolike_allSamples_withInvar \
            -nThreads {{resources.ntasks}} \
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
        """.format(ANGSD_DIR, REFERENCE_GENOME)

rule angsd_withMaf:
    input:
        bams = rules.create_bam_list.output
    output:
        gls = temp('{0}/withMaf/{{chrom}}/{{chrom}}_genolike_allSamples_withMaf{{maf}}.beagle.gz'.format(ANGSD_DIR)),
        mafs = temp('{0}/withMaf/{{chrom}}/{{chrom}}_genolike_allSamples_withMaf{{maf}}.mafs.gz'.format(ANGSD_DIR)),
    log: 'logs/angsd_withMaf/{chrom}_angsd_maf{maf}.log'
    conda: '../envs/angsd.yaml'
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE * 2,
        time = '12:00:00'
    shell:
        """
        angsd -GL 1 \
            -out {0}/withMaf/{{wildcards.chrom}}/{{wildcards.chrom}}_genolike_allSamples_withMaf{{wildcards.maf}} \
            -nThreads {{resources.ntasks}} \
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
            -minMaf {{wildcards.maf}} \
            -r {{wildcards.chrom}} \
            -bam {{input.bams}} 2> {{log}}
        """.format(ANGSD_DIR, REFERENCE_GENOME)
 
rule sfs_allSites:
    input:
        rules.angsd_withInvar.output.saf_idx 
    output:
        '{0}/sfs/{{chrom}}/{{chrom}}_allSamples_allSites.sfs'.format(ANGSD_DIR)
    log: 'logs/sfs_allSites/{chrom}_sfs_allSites.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 25000,
        time = '06:00:00'
    shell:
        """
        realSFS {input} -P {threads} -fold 1 > {output} 2> {log}
        """

rule thetas_per_site:
    input:
        saf_idx = rules.angsd_withInvar.output.saf_idx,
        sfs = rules.sfs_allSites.output
    output:
        idx = '{0}/summary_stats/thetas/{{chrom}}/{{chrom}}_allSamples_perSite.thetas.idx'.format(ANGSD_DIR),
        thet = '{0}/summary_stats/thetas/{{chrom}}/{{chrom}}_allSamples_perSite.thetas.gz'.format(ANGSD_DIR)
    log: 'logs/thetas_per_site/{chrom}_thetas.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '06:00:00'
    shell:
        """
        realSFS saf2theta {{input.saf_idx}} \
            -P {{threads}} \
            -sfs {{input.sfs}} \
            -outname {0}/summary_stats/thetas/{{wildcards.chrom}}/{{wildcards.chrom}}_allSamples_perSite 2> {{log}}
        """.format(ANGSD_DIR)

rule diversity_neutrality_byChrom:
    input:
        rules.thetas_per_site.output.idx
    output:
        '{0}/summary_stats/thetas/{{chrom}}/{{chrom}}_allSamples_allSites_diversityNeutrality.thetas.idx.pestPG'.format(ANGSD_DIR)
    log: 'logs/diversity_neutrality_byChrom/{{chrom}}_allSites_diversity_neutrality.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '06:00:00'
    shell:
        """
        thetaStat do_stat {input} 2> {log}
        """

# rule concat_angsd_gl:
#     input:
#         get_angsd_gl_toConcat
#     output:
#         '{0}/genolike_allSamples_withMaf{{minMaf}}.beagle.gz'.format(ANGSD_DIR)
#     log: 'logs/concat_angsd_gl/concat_angsd_gl_withMaf{{minMaf}}.log'
#     container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933' 
#     shell:
#         """
#         first=1
#         for f in {input}; do
#             if [ "$first"  ]; then
#                 zcat "$f"
#                 first=
#             else
#                 zcat "$f"| tail -n +2
#             fi
#         done | bgzip -c > {output} 2> {log}
#         """
# 
# rule concat_angsd_mafs:
#     input:
#         get_angsd_maf_toConcat
#     output:
#         '{0}/genolike_allSamples_withMaf{{minMaf}}.mafs.gz'.format(ANGSD_DIR)
#     log: 'logs/concat_angsd_mafs/concat_angsd_mafs_withMaf{{minMaf}}.log'
#     container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933' 
#     shell:
#         """
#         first=1
#         for f in {input}; do
#             if [ "$first"  ]; then
#                 zcat "$f"
#                 first=
#             else
#                 zcat "$f"| tail -n +2
#             fi
#         done | bgzip -c > {output} 2> {log}
#         """

