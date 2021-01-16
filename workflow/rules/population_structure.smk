rule angsd_full:
    input:
        bams = rules.create_bam_list.output
    output:
        gls = temp('{0}/full/{{chrom}}/{{chrom}}_genolike_allSamples_full.beagle.gz'.format(ANGSD_DIR)),
        mafs = temp('{0}/full/{{chrom}}/{{chrom}}_genolike_allSamples_full.mafs.gz'.format(ANGSD_DIR)),
        saf = temp('{0}/full/{{chrom}}/{{chrom}}_genolike_allSamples_full.saf.gz'.format(ANGSD_DIR)),
        saf_idx = temp('{0}/full/{{chrom}}/{{chrom}}_genolike_allSamples_full.saf.idx'.format(ANGSD_DIR)),
        saf_pos = temp('{0}/full/{{chrom}}/{{chrom}}_genolike_allSamples_full.saf.pos.gz'.format(ANGSD_DIR))
    log: 'logs/angsd_full/{chrom}_angsd_full.log'
    conda: '../envs/angsd.yaml'
    #container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE * 2,
        time = '12:00:00'
    shell:
        """
        angsd -GL 1 \
            -out {0}/full/{{wildcards.chrom}}/{{wildcards.chrom}}_genolike_allSamples_full \
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
            -doSaf 1 \
            -anc {1} \
            -r {{wildcards.chrom}} \
            -bam {{input.bams}} 2> {{log}}
        """.format(ANGSD_DIR, REFERENCE_GENOME)

rule angsd_withMaf:
    input:
        bams = rules.create_bam_list.output
    output:
        gls = temp('{0}/withMaf/{{chrom}}/{{chrom}}_genolike_allSamples_maf{{maf}}.beagle.gz'.format(ANGSD_DIR)),
        mafs = temp('{0}/withMaf/{{chrom}}/{{chrom}}_genolike_allSamples_maf{{maf}}.mafs.gz'.format(ANGSD_DIR)),
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
# 
# rule merge_safs:
#     input:
#         expand(rules.angsd_full.output.saf_idx, chrom=CHROMOSOMES)
#     output:
#         saf = '{0}/angsd/full/genolike_allSamples_full.saf.gz'.format(ANGSD_DIR),
#         saf_idx = '{0}/angsd/full/genolike_allSamples_full.saf.idx'.format(ANGSD_DIR),
#         saf_pos = '{0}/angsd/full/genolike_allSamples_full.saf.pos.gz'.format(ANGSD_DIR)
#     log: 'logs/merge_safs/merge_safs.log'
#     container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933' 
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 4000,
#         time = '02:00:00'
#     shell:
#         """
#         realSFS cat {{input}} \
#             -outnames {0}/angsd/genolike_allSamples_full 2> {{log}}
#         """.format(ANGSD_DIR)

# rule concat_angsd_gl:
#     input:
#         get_angsd_gl_toConcat
#     output:
#         '{0}/angsd/genolike_allSamples_{{minMaf}}.beagle.gz'.format(ANGSD_DIR)
#     log: 'logs/concat_angsd_gl/concat_angsd_gl_{{minMaf}}.log'
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

# rule concat_angsd_mafs:
#     input:
#         get_angsd_maf_toConcat
#     output:
#         '{0}/angsd/genolike_allSamples_{{minMaf}}.mafs.gz'.format(ANGSD_DIR)
#     log: 'logs/concat_angsd_mafs/concat_angsd_mafs_{{minMaf}}.log'
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

# rule create_pos_file_for_ngsLD:
#     input:
#         rules.angsd_withMaf.output.mafs
#     output:
#         '{0}/withMaf/{{chrom}}/{{chrom}}_angsdGL_maf{{maf}}.pos'.format(ANGSD_DIR)
#     log: 'logs/create_pos_file_for_ngsLD/{chrom}_maf{maf}_pos.log'
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
#         '{0}/{{chrom}}/{{chrom}}_genolike_allSamples_maf{{maf}}.ld.gz'.format(NGSLD_DIR)
#     log: 'logs/calc_ld_angsd_gl/{chrom}_maf{maf}_calc_ld.log'
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

# rule prune_ld:
#     input:
#         rules.calc_ld_angsd_gl.output
#     output:
#         '{0}/ld/pruned/{{chrom}}/{{chrom}}_pruned.id'.format(NGSLD_DIR)
#     log: 'logs/prune_ld/{chrom}_prune_ld.log'
#     container: 'shub://James-S-Santangelo/singularity-recipes:ngsld_v1.1.1'
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 16000,
#         time = '24:00:00'
#     shell:
#         """
#         ( zcat {input} | perl /opt/bin/prune_graph.pl \
#             --max_kb_dist 25 \
#             --min_weight 0.5 | sort -V > {output} ) 2> {log}
#         """
