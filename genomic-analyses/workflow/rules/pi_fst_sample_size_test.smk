rule urban_rural_toronto_bam_lists:
    input:
        all_bams = TOR_BAMS
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
            

rule pi_fst_sample_size_test_done:
    input:
        expand(rules.urban_rural_toronto_bam_lists.output, habitat=['urban', 'rural'])
    output:
        '{0}/pi_fst_sample_size_test/pi_fst_sample_size_test.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """
