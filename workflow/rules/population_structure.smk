# rule angsd_full:
#     input:
#         bams = rules.create_bam_list.output
#     output:
#         gls = temp('{0}/full/{{chrom}}/{{chrom}}_genolike_allSamples_full.beagle.gz'.format(ANGSD_DIR)),
#         mafs = temp('{0}/full/{{chrom}}/{{chrom}}_genolike_allSamples_full.mafs.gz'.format(ANGSD_DIR)),
#         saf = temp('{0}/full/{{chrom}}/{{chrom}}_genolike_allSamples_full.saf.gz'.format(ANGSD_DIR)),
#         saf_idx = temp('{0}/full/{{chrom}}/{{chrom}}_genolike_allSamples_full.saf.idx'.format(ANGSD_DIR)),
#         saf_pos = temp('{0}/full/{{chrom}}/{{chrom}}_genolike_allSamples_full.saf.pos.gz'.format(ANGSD_DIR))
#     log: 'logs/angsd_full/{chrom}_angsd_full.log'
#     conda: '../envs/angsd.yaml'
#     resources:
#         nodes = 1,
#         ntasks = CORES_PER_NODE * 2,
#         time = '12:00:00'
#     shell:
#         """
#         angsd -GL 1 \
#             -out {0}/full/{{wildcards.chrom}}/{{wildcards.chrom}}_genolike_allSamples_full \
#             -nThreads {{resources.ntasks}} \
#             -doGlf 2 \
#             -doMajorMinor 1 \
#             -SNP_pval 1e-6 \
#             -doMaf 1 \
#             -doCounts 1 \
#             -setMinDepthInd 3 \
#             -setMaxDepth 4500 \
#             -baq 2 \
#             -ref {1} \
#             -minInd 96 \
#             -minQ 20 \
#             -minMapQ 30 \
#             -doSaf 1 \
#             -anc {1} \
#             -r {{wildcards.chrom}} \
#             -bam {{input.bams}} 2> {{log}}
#         """.format(ANGSD_DIR, REFERENCE_GENOME)

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

# rule angsd_withMaf:
#     input:
#         bams = rules.create_bam_list.output
#     output:
#         gls = temp('{0}/withMaf/{{chrom}}/{{chrom}}_genolike_allSamples_withMaf{{maf}}.beagle.gz'.format(ANGSD_DIR)),
#         mafs = temp('{0}/withMaf/{{chrom}}/{{chrom}}_genolike_allSamples_withMaf{{maf}}.mafs.gz'.format(ANGSD_DIR)),
#     log: 'logs/angsd_withMaf/{chrom}_angsd_maf{maf}.log'
#     conda: '../envs/angsd.yaml'
#     resources:
#         nodes = 1,
#         ntasks = CORES_PER_NODE * 2,
#         time = '12:00:00'
#     shell:
#         """
#         angsd -GL 1 \
#             -out {0}/withMaf/{{wildcards.chrom}}/{{wildcards.chrom}}_genolike_allSamples_withMaf{{wildcards.maf}} \
#             -nThreads {{resources.ntasks}} \
#             -doGlf 2 \
#             -doMajorMinor 1 \
#             -SNP_pval 1e-6 \
#             -doMaf 1 \
#             -doCounts 1 \
#             -setMinDepthInd 3 \
#             -setMaxDepth 4500 \
#             -baq 2 \
#             -ref {1} \
#             -minInd 96 \
#             -minQ 20 \
#             -minMapQ 30 \
#             -minMaf {{wildcards.maf}} \
#             -r {{wildcards.chrom}} \
#             -bam {{input.bams}} 2> {{log}}
#         """.format(ANGSD_DIR, REFERENCE_GENOME)
# 
# rule merge_safs:
#     input:
#         expand(rules.angsd_full.output.saf_idx, chrom=CHROMOSOMES)
#     output:
#         saf = '{0}/sfs/genolike_allSamples_full.saf.gz'.format(ANGSD_DIR),
#         saf_idx = '{0}/sfs/genolike_allSamples_full.saf.idx'.format(ANGSD_DIR),
#         saf_pos = '{0}/sfs/genolike_allSamples_full.saf.pos.gz'.format(ANGSD_DIR)
#     log: 'logs/merge_safs/merge_safs.log'
#     container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933' 
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 4000,
#         time = '02:00:00'
#     shell:
#         """
#         realSFS cat {{input}} \
#             -outnames {0}/sfs/genolike_allSamples_full 2> {{log}}
#         """.format(ANGSD_DIR)
# 
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
# 
# rule create_pos_file_for_ngsLD:
#     input:
#         rules.angsd_withMaf.output.mafs
#     output:
#         '{0}/ngsld_pos/{{chrom}}_angsdGL_withMaf{{maf}}.pos'.format(PROGRAM_RESOURCE_DIR)
#     log: 'logs/create_pos_file_for_ngsLD/{chrom}_withMaf{maf}_pos.log'
#     shell:
#         """
#         zcat {input} | cut -f 1,2 | tail -n +2 > {output} 2> {log}
#         """
# 
# rule calc_ld_angsd_gl:
#     input:
#         pos = rules.create_pos_file_for_ngsLD.output,
#         gls = rules.angsd_withMaf.output.gls
#     output:
#         '{0}/{{chrom}}/{{chrom}}_genolike_allSamples_withMaf{{maf}}.ld.gz'.format(NGSLD_DIR)
#     log: 'logs/calc_ld_angsd_gl/{chrom}_withMaf{maf}_calc_ld.log'
#     container: 'shub://James-S-Santangelo/singularity-recipes:ngsld_v1.1.1'
#     threads: 16
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 8000,
#         time = '16:00:00'
#     shell:
#         """
#         ( NUM_SITES=$(cat {{input.pos}} | wc -l) &&
#           ngsLD --geno {{input.gls}} \
#             --pos {{input.pos}} \
#             --n_ind {0} \
#             --n_sites $NUM_SITES \
#             --probs \
#             --min_maf {{wildcards.maf}} \
#             --n_threads {{threads}} \
#             --max_kb_dist 100 | gzip --best > {{output}} ) 2> {{log}}
#         """.format(len(SAMPLES))
# 
# rule global_sfs:
#     input:
#         rules.merge_safs.output.saf_idx 
#     output:
#         '{0}/sfs/allSamples_global.sfs'.format(ANGSD_DIR)
#     log: 'logs/global_sfs/global_sfs.log'
#     container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
#     threads: 16
#     resources:
#         mem_mb = 30000,
#         time = '01:00:00'
#     shell:
#         """
#         realSFS {input} -P {threads} -fold 1 > {output} 2> {log}
#         """
# 
# rule thetas_per_site:
#     input:
#         saf_idx = rules.merge_safs.output.saf_idx,
#         sfs = rules.global_sfs.output
#     output:
#         idx = '{0}/summary_stats/thetas/allSamples_perSite.thetas.idx'.format(ANGSD_DIR),
#         thet = '{0}/summary_stats/thetas/allSamples_perSite.thetas.gz'.format(ANGSD_DIR)
#     log: 'logs/thetas_per_site/thetas.log'
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
#             -outname {0}/summary_stats/thetas/allSamples_perSite 2> {{log}}
#         """.format(ANGSD_DIR)
# 
# rule theta_stat_wholeGenome:
#     input:
#         rules.thetas_per_site.output.idx
#     output:
#         '{0}/summary_stats/thetas/genome-wide/allSamples_perSite.thetas.idx.pestPG'.format(ANGSD_DIR)
#     log: 'logs/theta_stat_wholeGenome.log'
#     container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 4000,
#         time = '06:00:00'
#     shell:
#         """
#         thetaStat do_stat {input} 2> {output}
#         """
#         
# 
# rule prune_ld:
#     input:
#         rules.calc_ld_angsd_gl.output
#     output:
#         '{0}/pruned/{{chrom}}/{{chrom}}_withMaf{{maf}}_pruned.id'.format(NGSLD_DIR)
#     log: 'logs/prune_ld/{chrom}_withMaf{maf}_prune_ld.log'
#     container: 'shub://James-S-Santangelo/singularity-recipes:ngsld_v1.1.1'
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 40000,
#         time = '24:00:00'
#     shell:
#         """
#         ( zcat {input} | perl /opt/bin/prune_graph.pl \
#             --max_kb_dist 20 \
#             --min_weight 0.2 | sort -V > {output} ) 2> {log}
#         """
