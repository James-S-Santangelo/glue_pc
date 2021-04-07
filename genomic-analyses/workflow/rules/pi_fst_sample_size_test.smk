rule urban_rural_toronto_bam_lists:
    input:
        TOR_BAMS
    output:
        '{0}/bam_lists/toronto_pi_fst_test/{{habitat}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    run:
        import os
        import glob
        if wildcards.habitat == 'urban':
            pop = '41'
        elif wildcards.habitat == 'rural':
            pop = '83'
        bams = glob.glob(input[0] + '/*.bam')
        with open(output[0], 'w') as fout:
            for bam in bams:
                base = os.path.basename(bam)
                population = base.split('_')[1]
                if population == pop:
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

rule pi_fst_sample_size_test_done:
    input:
        expand(rules.urban_rural_toronto_bam_lists.output, habitat=['urban', 'rural']),
        expand(rules.index_urban_rural_toronto_pi_fst_test_bam.output, tor_test_sample=TOR_SAMPLES_PI_FST_TEST)
    output:
        '{0}/pi_fst_sample_size_test/pi_fst_sample_size_test.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """
