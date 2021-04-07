rule urban_rural_toronto_bam_lists:
    input:
        TOR_BAMS
    output:
        '{0}/bam_lists/toronto_pi_fst_test/{{group}}_{{habitat}}_{{ss}}ind_bams.list'.format(PROGRAM_RESOURCE_DIR)
    wildcard_constraints:
        habitat = 'urban|rural',
        group = 'highCov',
        ss = str(max(SS_PI_FST_TEST))
    run:
        import os 
        import glob
        import random
        random.seed(42)
        if wildcards.habitat == 'urban':
            pop = '41'
        elif wildcards.habitat == 'rural':
            pop = '83'
        bams = glob.glob(input[0] + '/s_{0}_*.bam'.format(pop))
        bams_relate_remove = [x for x in bams if 's_83_11' not in x]
        subset_bams = random.sample(bams_relate_remove, int(wildcards.ss))
        with open(output[0], 'w') as fout:
            for bam in subset_bams:
                fout.write('{0}\n'.format(bam))

rule downsample_urban_rural_toronto_pi_fst_test_bams:
    input:
        get_toronto_bam_pi_fst_test
    output:
        '{0}/toronto_urban_rural_pi_fst/{{tor_test_sample}}_CM019101.1_merged_sorted_dupsMarked_downsampled.bam'.format(BAM_DIR)
    log: 'logs/downsample_urban_rural_toronto_bams/{tor_test_sample}_downsample.log'
    conda: '../envs/mapping.yaml'
    threads: 6
    params:
        region = 'CM019101.1'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '03:00:00'
    shell:
        """
        samtools view -hb -s 0.1 -@ {threads} {input} {params.region} > {output} 2> {log}
        """

rule index_urban_rural_toronto_pi_fst_test_bam:
    input:
        rules.downsample_urban_rural_toronto_pi_fst_test_bams.output
    output:
        '{0}/toronto_urban_rural_pi_fst/{{tor_test_sample}}_CM019101.1_merged_sorted_dupsMarked_downsampled.bam.bai'.format(BAM_DIR)
    conda: '../envs/mapping.yaml'
    log: 'logs/index_urban_rural_toronto_pi_fst_bams/{tor_test_sample}_index.log'
    resources:
        mem_mb = 1000,
        time = '01:00:00'
    shell:
        """
        samtools index {input} 2> {log}
        """

rule urban_rural_toronto_downsampled_bam_lists:
    input:
        expand(rules.downsample_urban_rural_toronto_pi_fst_test_bams.output, tor_test_sample=TOR_SAMPLES_PI_FST_TEST)
    output:
        '{0}/bam_lists/toronto_pi_fst_test/{{group}}_{{habitat}}_{{ss}}ind_bams.list'.format(PROGRAM_RESOURCE_DIR)
    wildcard_constraints:
        habitat = 'urban|rural',
        group = 'lowCov'
    run:
        import os
        import glob
        import random
        random.seed(42)
        if wildcards.habitat == 'urban':
            pop = '41'
        elif wildcards.habitat == 'rural':
            pop = '83'
        possible_bams = [x for x in input if os.path.basename(x).split('_')[1] == pop]
        subset_bams = random.sample(possible_bams, int(wildcards.ss))
        with open(output[0], 'w') as fout:
            for bam in subset_bams:
                fout.write('{0}\n'.format(bam))


# rule angsd_saf_likelihood_toronto_pi_fst_test:
#     input:
#         bams = rules.st_byCity_byHabitat.output,
#         sites_idx = rules.angsd_index_degenerate_allChroms.output,
#         sites = rules.convert_sites_for_angsd.output, 
#         ref = REFERENCE_GENOME,
#         chroms = config['chromosomes']
#     output:
#         saf = temp('{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.saf.gz'.format(ANGSD_DIR)),
#         saf_idx = temp('{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.saf.idx'.format(ANGSD_DIR)),
#         saf_pos = temp('{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.saf.pos.gz'.format(ANGSD_DIR))
#     log: 'logs/angsd_saf_likelihood_byCity_byHabitat/{city}_{habitat}_{site}_saf.log'
#     container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
#     params:
#         out = '{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}'.format(ANGSD_DIR)
#     threads: 10
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 8000,
#         time = '06:00:00'
#     wildcard_constraints:
#         site='4fold'
#     shell:
#         """
#         NUM_IND=$( wc -l < {input.bams} );
#         MIN_IND=$(( NUM_IND*50/100 ));
#         if [[ $MIN_IND -eq 0 ]]; then MIN_IND=1; fi;
#         angsd -GL 1 \
#             -out {params.out} \
#             -nThreads {threads} \
#             -doMajorMinor 4 \
#             -baq 2 \
#             -ref {input.ref} \
#             -minInd $MIN_IND \
#             -sites {input.sites} \
#             -minQ 20 \
#             -minMapQ 30 \
#             -doSaf 1 \
#             -anc {input.ref} \
#             -rf {input.chroms} \
#             -bam {input.bams} 2> {log}
#         """



rule pi_fst_sample_size_test_done:
    input:
        expand(rules.urban_rural_toronto_bam_lists.output, habitat=['urban', 'rural'], group = 'highCov', ss='7'),
        expand(rules.index_urban_rural_toronto_pi_fst_test_bam.output, tor_test_sample=TOR_SAMPLES_PI_FST_TEST),
        expand(rules.urban_rural_toronto_downsampled_bam_lists.output, habitat=['urban', 'rural'], ss=SS_PI_FST_TEST, group = 'lowCov')
    output:
        '{0}/pi_fst_sample_size_test/pi_fst_sample_size_test.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """
