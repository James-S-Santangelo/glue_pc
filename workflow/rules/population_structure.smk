rule create_pos_file_for_ngsLD:
    input:
        rules.angsd_gl_withMaf.output.mafs
    output:
        '{0}/ngsld_pos/{{chrom}}_angsdGL_withMaf{{maf}}.pos'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/create_pos_file_for_ngsLD/{chrom}_withMaf{maf}_pos.log'
    shell:
        """
        zcat {input} | cut -f 1,2 | tail -n +2 > {output} 2> {log}
        """

rule calc_ld_angsd_gl:
    input:
        pos = rules.create_pos_file_for_ngsLD.output,
        gls = rules.angsd_gl_withMaf.output.gls
    output:
        '{0}/{{chrom}}/{{chrom}}_genolike_allSamples_withMaf{{maf}}.ld.gz'.format(NGSLD_DIR)
    log: 'logs/calc_ld_angsd_gl/{chrom}_withMaf{maf}_calc_ld.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:ngsld_v1.1.1'
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '16:00:00'
    shell:
        """
        ( NUM_SITES=$(cat {{input.pos}} | wc -l) &&
          ngsLD --geno {{input.gls}} \
            --pos {{input.pos}} \
            --n_ind {0} \
            --n_sites $NUM_SITES \
            --probs \
            --min_maf {{wildcards.maf}} \
            --n_threads {{threads}} \
            --max_kb_dist 100 | gzip --best > {{output}} ) 2> {{log}}
        """.format(len(SAMPLES))

# rule prune_ld:
#     input:
#         rules.calc_ld_angsd_gl.output
#     output:
#         '{0}/pruned/{{chrom}}/{{chrom}}_withMaf{{maf}}_pruned.id'.format(NGSLD_DIR)
#     log: 'logs/prune_ld/{chrom}_withMaf{maf}_prune_ld.log'
#     container: 'shub://James-S-Santangelo/singularity-recipes:ngsld_v1.1.1'
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 50000,
#         time = '24:00:00'
#     shell:
#         """
#         ( zcat {input} | perl /opt/bin/prune_graph.pl \
#             --max_kb_dist 20 \
#             --min_weight 0.2 | sort -V > {output} ) 2> {log}
#         """
