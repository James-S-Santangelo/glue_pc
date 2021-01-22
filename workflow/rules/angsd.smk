rule angsd_depth:
    input:
        bams = rules.create_bam_list.output
    output:
        sam = '{0}/depth/{{chrom}}/{{chrom}}_allSamples_allSites.depthSample'.format(ANGSD_DIR),
        glo = '{0}/depth/{{chrom}}/{{chrom}}_allSamples_allSites.depthGlobal'.format(ANGSD_DIR)
    log: 'logs/angsd_depth/{chrom}_angsd_depth.log'
    conda: '../envs/angsd.yaml'
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE,
        time = '12:00:00'
    shell:
        """
        angsd -bam {{input.bams}} \
            -nThreads {{resources.ntasks}} \
            -doDepth 1 \
            -doCounts 1 \
            -r {{wildcards.chrom}} \
            -minMapQ 30 \
            -minQ 20 \
            -maxDepth 2000 \
            -out {0}/depth/{{wildcards.chrom}}/{{wildcards.chrom}}_allSamples_allSites
        """.format(ANGSD_DIR)

rule convert_sites_for_angsd:
    input:
        rules.get_fourfold_zerofold.output
    output:
        '{0}/angsd_sites/Trepens_{{site}}.sites'.format(PROGRAM_RESOURCE_DIR), 
    log: 'logs/convert_sites_for_angsd/convert_{site}_for_angsd.log'
    conda: '../envs/angsd.yaml'
    shell:
        """
        touch {output}
        ( awk '{{print $1"\t"$2+1"\t"$3}}' {input} > {output} && \
            angsd sites index {output} ) 2> {log}
        """

rule angsd_saf_likelihood:
    input:
        bams = rules.create_bam_list.output,
        sites = rules.convert_sites_for_angsd.output
    output:
        saf = temp('{0}/sfs/{{site}}/{{chrom}}/{{chrom}}_allSamples_{{site}}.saf.gz'.format(ANGSD_DIR)),
        saf_idx = temp('{0}/sfs/{{site}}/{{chrom}}/{{chrom}}_allSamples_{{site}}.saf.idx'.format(ANGSD_DIR)),
        saf_pos = temp('{0}/sfs/{{site}}/{{chrom}}/{{chrom}}_allSamples_{{site}}.saf.pos.gz'.format(ANGSD_DIR))
    log: 'logs/angsd_{site}/{chrom}_angsd_{site}.log'
    conda: '../envs/angsd.yaml'
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE,
        time = '12:00:00'
    shell:
        """
        if [ wildcards.site = 'allSites' ]
        then
            angsd -GL 1 \
                -out {0}/sfs/{{wildcards.site}}/{{wildcards.chrom}}/{{wildcards.chrom}}_allSamples_{{wildcards.site}} \
                -nThreads {{resources.ntasks}} \
                -doCounts 1 \
                -dumpCounts 2 \
                -setMinDepthInd 1 \
                -setMaxDepth 4500 \
                -baq 2 \
                -ref {1} \
                -minInd 60 \
                -minQ 20 \
                -minMapQ 30 \
                -doSaf 1 \
                -anc {1} \
                -r {{wildcards.chrom}} \
                -bam {{input.bams}} 2> {{log}}
        else
            angsd -GL 1 \
                -out {0}/sfs/{{wildcards.site}}/{{wildcards.chrom}}/{{wildcards.chrom}}_allSamples_{{wildcards.site}} \
                -nThreads {{resources.ntasks}} \
                -sites {{input.sites}} \
                -doCounts 1 \
                -dumpCounts 2 \
                -setMinDepthInd 1 \
                -setMaxDepth 4500 \
                -baq 2 \
                -ref {1} \
                -minInd 60 \
                -minQ 20 \
                -minMapQ 30 \
                -doSaf 1 \
                -anc {1} \
                -r {{wildcards.chrom}} \
                -bam {{input.bams}} 2> {{log}}
        """.format(ANGSD_DIR, REFERENCE_GENOME)

rule angsd_gl_withMaf:
    input:
        bams = rules.create_bam_list.output
    output:
        gls = temp('{0}/gl/withMaf/{{chrom}}/{{chrom}}_genolike_allSamples_withMaf{{maf}}.beagle.gz'.format(ANGSD_DIR)),
        mafs = temp('{0}/gl/withMaf/{{chrom}}/{{chrom}}_genolike_allSamples_withMaf{{maf}}.mafs.gz'.format(ANGSD_DIR)),
    log: 'logs/angsd_gl_withMaf/{chrom}_angsd_maf{maf}.log'
    conda: '../envs/angsd.yaml'
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE,
        time = '12:00:00'
    shell:
        """
        angsd -GL 1 \
            -out {0}/gl/withMaf/{{wildcards.chrom}}/{{wildcards.chrom}}_genolike_allSamples_withMaf{{wildcards.maf}} \
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
 
rule angsd_estimate_sfs:
    input:
        rules.angsd_saf_likelihood.output.saf_idx 
    output:
        '{0}/sfs/{{site}}/{{chrom}}/{{chrom}}_allSamples_{{site}}.sfs'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_sfs/{chrom}_{site}_sfs.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 25000,
        time = '06:00:00'
    shell:
        """
        realSFS {input} -P {threads} -fold 1 > {output} 2> {log}
        """

# rule angsd_estimate_thetas:
#     input:
#         saf_idx = rules.angsd_saf_likelihood.output.saf_idx,
#         sfs = rules.angsd_estimate_sfs.output
#     output:
#         idx = '{0}/summary_stats/thetas/{{chrom}}/{{chrom}}_allSamples_{{site}}.thetas.idx'.format(ANGSD_DIR),
#         thet = '{0}/summary_stats/thetas/{{chrom}}/{{chrom}}_allSamples_{{site}}.thetas.gz'.format(ANGSD_DIR)
#     log: 'logs/angsd_estimate_thetas/{chrom}_{site}_thetas.log'
#     container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
#     threads: 10
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 10000,
#         time = '06:00:00'
#     shell:
#         """
#         realSFS saf2theta {{input.saf_idx}} \
#             -P {{threads}} \
#             -sfs {{input.sfs}} \
#             -outname {0}/summary_stats/thetas/{{wildcards.chrom}}/{{wildcards.chrom}}_allSamples_{{wildcards.site}} 2> {{log}}
#         """.format(ANGSD_DIR)
# 
# rule angsd_diversity_neutrality_stats:
#     input:
#         rules.angsd_estimate_thetas.output.idx
#     output:
#         '{0}/summary_stats/thetas/{{chrom}}/{{chrom}}_allSamples_{{site}}.thetas.idx.pestPG'.format(ANGSD_DIR)
#     log: 'logs/angsd_diversity_neutrality_stats/{{chrom}}_{site}_diversity_neutrality.log'
#     container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 4000,
#         time = '06:00:00'
#     shell:
#         """
#         thetaStat do_stat {input} 2> {log}
#         """
# 
# rule concat_angsd_stats:
#     input:
#         expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES)
#     output:
#         '{0}/summary_stats/thetas/allSamples_{{site}}_diversityNeutrality.thetas.idx.pestPG'.format(ANGSD_DIR)
#     log: 'logs/concat_angsd_stats/concat_angsd_stats.log'
#     shell:
#         """
#         first=1
#         for f in {input}; do
#             if [ "$first"  ]; then
#                 cat "$f"
#                 first=
#             else
#                 cat "$f"| tail -n +2
#             fi
#         done > {output} 2> {log}
#         """
# 
# rule concat_sfs_allSites:
#     input:
#         expand(rules.sfs_allSites.output, chrom=CHROMOSOMES)
#     output:
#         '{0}/sfs/allSites/allSamples_allSites_allChroms.sfs'.format(ANGSD_DIR)
#     log: 'logs/concat_sfs_allSites/concat_sfs_allSites.log'
#     run:
#         shell('cat {input} > temp.txt 2> {log}')
#         import pandas as pd
#         sfs_allChroms = pd.read_table('temp.txt', delimiter = '\t')
#         sfs_sum = sfs_allChroms.sum(axis=0) 
#         sfs_sum.to_csv(output[0], sep = '\t', header = None)
#         shell('rm temp.txt')

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

