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
        awk '{{print $1"\t"$2+1"\t"$3}}' {input} > {output} 2> {log}
        """

rule split_angsd_sites_byChrom:
    input:
        rules.convert_sites_for_angsd.output
    output:
        '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}.sites'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/split_angsd_sites_byChrom/{{chrom}}/{{chrom}}_{{site}}_split_angsd_sites.log'
    shell:
        """
        grep {wildcards.chrom} {input} > {output} 2> {log}
        """

rule angsd_index_sites:
    input:
        rules.split_angsd_sites_byChrom.output
    output:
        binary = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}.sites.bin'.format(PROGRAM_RESOURCE_DIR),
        idx = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}.sites.idx'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/angsd_index_sites/{chrom}_{site}_index.log'
    conda: '../envs/angsd.yaml'
    shell:
        """
        angsd sites index {input} 2> {log}
        """

rule angsd_saf_likelihood_allSites:
    input:
        bams = rules.create_bam_list.output
    output:
        saf = temp('{0}/sfs/allSites/{{chrom}}/{{chrom}}_allSamples_allSites.saf.gz'.format(ANGSD_DIR)),
        saf_idx = temp('{0}/sfs/allSites/{{chrom}}/{{chrom}}_allSamples_allSites.saf.idx'.format(ANGSD_DIR)),
        saf_pos = temp('{0}/sfs/allSites/{{chrom}}/{{chrom}}_allSamples_allSites.saf.pos.gz'.format(ANGSD_DIR)),
        pos = temp('{0}/sfs/allSites/{{chrom}}/{{chrom}}_allSamples_allSites.pos.gz'.format(ANGSD_DIR)),
        counts = temp('{0}/sfs/allSites/{{chrom}}/{{chrom}}_allSamples_allSites.counts.gz'.format(ANGSD_DIR))
    log: 'logs/angsd_saf_likelihood_allSites/{chrom}_allSites_angsd_saf.log'
    conda: '../envs/angsd.yaml'
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE,
        time = '12:00:00'
    wildcard_constraints:
        site='allSites'
    shell:
        """
        angsd -GL 1 \
            -out {0}/sfs/allSites/{{wildcards.chrom}}/{{wildcards.chrom}}_allSamples_allSites \
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
        """.format(ANGSD_DIR, REFERENCE_GENOME)

rule angsd_saf_likelihood_specificSites:
    input:
        bams = rules.create_bam_list.output,
        sites = rules.split_angsd_sites_byChrom.output,
        idx = rules.angsd_index_sites.output
    output:
        saf = temp('{0}/sfs/{{site}}/{{chrom}}/{{chrom}}_allSamples_{{site}}.saf.gz'.format(ANGSD_DIR)),
        saf_idx = temp('{0}/sfs/{{site}}/{{chrom}}/{{chrom}}_allSamples_{{site}}.saf.idx'.format(ANGSD_DIR)),
        saf_pos = temp('{0}/sfs/{{site}}/{{chrom}}/{{chrom}}_allSamples_{{site}}.saf.pos.gz'.format(ANGSD_DIR)),
        pos = temp('{0}/sfs/{{site}}/{{chrom}}/{{chrom}}_allSamples_{{site}}.pos.gz'.format(ANGSD_DIR)),
        counts = temp('{0}/sfs/{{site}}/{{chrom}}/{{chrom}}_allSamples_{{site}}.counts.gz'.format(ANGSD_DIR))
    log: 'logs/angsd_saf_likelihood_specificSites/{chrom}_{site}_angsd_saf.log'
    conda: '../envs/angsd.yaml'
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE,
        time = '12:00:00'
    wildcard_constraints:
        site='0fold|4fold'
    shell:
        """
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

rule angsd_gl_allSites:
    input:
        bams = rules.create_bam_list.output
    output:
        gls = temp('{0}/gl/allSites/{{chrom}}/{{chrom}}_allSamples_allSites_maf{{maf}}.beagle.gz'.format(ANGSD_DIR)),
        mafs = temp('{0}/gl/allSites/{{chrom}}/{{chrom}}_allSamples_allSites_maf{{maf}}.mafs.gz'.format(ANGSD_DIR)),
    log: 'logs/angsd_gl_allSites/{chrom}_allSites_maf{maf}_angsd_gl.log'
    conda: '../envs/angsd.yaml'
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE,
        time = '12:00:00'
    wildcard_constraints:
        site='allSites'
    shell:
        """
        angsd -GL 1 \
            -out {0}/gl/allSites/{{wildcards.chrom}}/{{wildcards.chrom}}_allSamples_allSites_maf{{wildcards.maf}} \
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
 
rule angsd_gl_specificSites:
    input:
        bams = rules.create_bam_list.output,
        sites = rules.split_angsd_sites_byChrom.output,
        idx = rules.angsd_index_sites.output
    output:
        gls = temp('{0}/gl/{{site}}/{{chrom}}/{{chrom}}_allSamples_{{site}}_maf{{maf}}.beagle.gz'.format(ANGSD_DIR)),
        mafs = temp('{0}/gl/{{site}}/{{chrom}}/{{chrom}}_allSamples_{{site}}_maf{{maf}}.mafs.gz'.format(ANGSD_DIR)),
    log: 'logs/angsd_gl_specificSites/{chrom}_{site}_maf{maf}_angsd_gl.log'
    conda: '../envs/angsd.yaml'
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE,
        time = '12:00:00'
    wildcard_constraints:
        site='0fold|4fold'
    shell:
        """
        angsd -GL 1 \
            -out {0}/gl/{{wildcards.site}}/{{wildcards.chrom}}/{{wildcards.chrom}}_allSamples_{{wildcards.site}}_maf{{wildcards.maf}} \
            -nThreads {{resources.ntasks}} \
            -sites {{input.sites}} \
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

rule angsd_estimate_sfs_allSites:
    input:
        rules.angsd_saf_likelihood_allSites.output.saf_idx 
    output:
        '{0}/sfs/allSites/{{chrom}}/{{chrom}}_allSamples_allSites.sfs'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_sfs_allSites/{chrom}_allSites_sfs.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 10
    wildcard_constraints:
        site='allSites'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 30000,
        time = '06:00:00'
    shell:
        """
        realSFS {input} -P {threads} -fold 1 > {output} 2> {log}
        """

rule angsd_estimate_sfs_specificSites:
    input:
        rules.angsd_saf_likelihood_specificSites.output.saf_idx 
    output:
        '{0}/sfs/{{site}}/{{chrom}}/{{chrom}}_allSamples_{{site}}.sfs'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_sfs_specificSites/{chrom}_{site}_sfs.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 10
    wildcard_constraints:
        site='0fold|4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 30000,
        time = '06:00:00'
    shell:
        """
        realSFS {input} -P {threads} -fold 1 > {output} 2> {log}
        """

rule angsd_estimate_thetas_allSites:
    input:
        saf_idx = rules.angsd_saf_likelihood_allSites.output.saf_idx,
        sfs = rules.angsd_estimate_sfs_allSites.output
    output:
        idx = '{0}/summary_stats/thetas/allSites/{{chrom}}/{{chrom}}_allSamples_allSites.thetas.idx'.format(ANGSD_DIR),
        thet = '{0}/summary_stats/thetas/allSites/{{chrom}}/{{chrom}}_allSamples_allSites.thetas.gz'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_thetas_allSites/{chrom}_allSites_thetas.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 10
    wildcard_constraints:
        site='allSites'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '06:00:00'
    shell:
        """
        realSFS saf2theta {{input.saf_idx}} \
            -P {{threads}} \
            -sfs {{input.sfs}} \
            -outname {0}/summary_stats/thetas/allSites/{{wildcards.chrom}}/{{wildcards.chrom}}_allSamples_allSites 2> {{log}}
        """.format(ANGSD_DIR)

rule angsd_estimate_thetas_specificSites:
    input:
        saf_idx = rules.angsd_saf_likelihood_specificSites.output.saf_idx,
        sfs = rules.angsd_estimate_sfs_specificSites.output
    output:
        idx = '{0}/summary_stats/thetas/{{site}}/{{chrom}}/{{chrom}}_allSamples_{{site}}.thetas.idx'.format(ANGSD_DIR),
        thet = '{0}/summary_stats/thetas/{{site}}/{{chrom}}/{{chrom}}_allSamples_{{site}}.thetas.gz'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_thetas_specificSites/{chrom}_{site}_thetas.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 10
    wildcard_constraints:
        site='0fold|4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '06:00:00'
    shell:
        """
        realSFS saf2theta {{input.saf_idx}} \
            -P {{threads}} \
            -sfs {{input.sfs}} \
            -outname {0}/summary_stats/thetas/{{wildcards.site}}/{{wildcards.chrom}}/{{wildcards.chrom}}_allSamples_{{wildcards.site}} 2> {{log}}
        """.format(ANGSD_DIR)

rule angsd_diversity_neutrality_stats_allSites:
    input:
        rules.angsd_estimate_thetas_allSites.output.idx
    output:
        '{0}/summary_stats/thetas/allSites/{{chrom}}/{{chrom}}_allSamples_allSites.thetas.idx.pestPG'.format(ANGSD_DIR)
    log: 'logs/angsd_diversity_neutrality_stats_allSites/{chrom}_allSites_diversity_neutrality.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '06:00:00'
    wildcard_constraints:
        site='allSites'
    shell:
        """
        thetaStat do_stat {input} 2> {log}
        """
 
rule angsd_diversity_neutrality_stats_specificSites:
    input:
        rules.angsd_estimate_thetas_specificSites.output.idx
    output:
        '{0}/summary_stats/thetas/{{site}}/{{chrom}}/{{chrom}}_allSamples_{{site}}.thetas.idx.pestPG'.format(ANGSD_DIR)
    log: 'logs/angsd_diversity_neutrality_stats_specificSites/{chrom}}_{site}_diversity_neutrality.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '06:00:00'
    wildcard_constraints:
        site='0fold|4fold'
    shell:
        """
        thetaStat do_stat {input} 2> {log}
        """

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

