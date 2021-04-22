# Rules to perform small-scale test of the effects of sample size (i.e., # individuals)
# on estimating Pi and Fst. Compares Pi and Fst estimated using combinations of 1, 3, 5, or 7
# individuals in each habitat (sampled without replacement). "Control" is pi and Fst estimated
# from 7 urban and 7 rural individuals from high-coverage data (~12X).

rule urban_rural_toronto_bam_lists:
    """
    Generates bam lists with differing numbers of urban and rural individuals from Toronto populations.
    Samples either 1, 3, 5, or 7 individuals without replacement. Writes paths to raw BAMs (i.e., not downsampled)
    """
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
    """
    Downsamples all Toronto BAMs from select urban (population = 41) and rural (population = 83)
    populations. Samples 10% of reads.
    """
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
    """
    Indexes downsampled Toronto BAMs.
    """
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
    """
    Generates bam lists with differing numbers of urban and rural individuals from Toronto populations.
    Samples either 1, 3, 5, or 7 individuals without replacement. Writes paths to downsampled_bams.
    """
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


rule angsd_saf_likelihood_toronto_pi_fst_test:
    """
    Estimate Site Allele Frequency (SAF) likelihood file for urban and rural sample size combinations. 
    Includes "control" (i.e., max urban and rural sample size and high coverage data)
    """
    input:
        unpack(toronto_pi_fst_test_saf_input),
        ref = REFERENCE_GENOME
    output:
        saf = temp('{0}/pi_fst_sample_size_test/{{group}}/{{group}}_{{habitat}}_{{ss}}ind.saf.gz'.format(ANGSD_DIR)),
        saf_idx = temp('{0}/pi_fst_sample_size_test/{{group}}/{{group}}_{{habitat}}_{{ss}}ind.saf.idx'.format(ANGSD_DIR)),
        saf_pos = temp('{0}/pi_fst_sample_size_test/{{group}}/{{group}}_{{habitat}}_{{ss}}ind.saf.pos.gz'.format(ANGSD_DIR)),
    log: 'logs/angsd_saf_likelihood_toronto_pi_fst_test/{group}_{habitat}_{ss}ind_saf.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    params:
        out = '{0}/pi_fst_sample_size_test/{{group}}/{{group}}_{{habitat}}_{{ss}}ind'.format(ANGSD_DIR)
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '03:00:00'
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND*50/100 ));
        if [[ $MIN_IND -eq 0 ]]; then MIN_IND=1; fi;
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 4 \
            -baq 2 \
            -ref {input.ref} \
            -minInd $MIN_IND \
            -sites {input.sites} \
            -minQ 20 \
            -minMapQ 30 \
            -doSaf 1 \
            -anc {input.ref} \
            -r CM019101.1 \
            -bam {input.bams} 2> {log}
        """


rule angsd_estimate_joint_sfs_toronto_fst_test:
    """
    Estimate 2D urban-rural folded SFS from SAF files for sample size test.
    """
    input:
        get_saf_joint_sfs_pi_fst_test
    output:
        '{0}/pi_fst_sample_size_test/{{group}}/sfs/{{group}}_{{joint_sfs_ss}}.2dsfs'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_joint_sfs_toronto_fst_test/{group}_{joint_sfs_ss}.2dsfs.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '01:00:00'
    shell:
        """
        realSFS {input} -maxIter 2000 -seed 42 -fold 1 -P {threads} > {output} 2> {log}
        """

rule angsd_fst_index_toronto_fst_test:
    """
    Estimate per-site alphas (numerator) and betas (denominator) for Hudson's Fst
    """
    input: 
        saf_idx = get_saf_joint_sfs_pi_fst_test,
        joint_sfs = rules.angsd_estimate_joint_sfs_toronto_fst_test.output
    output:
        fst = '{0}/pi_fst_sample_size_test/stats/fst/{{group}}/{{group}}_{{joint_sfs_ss}}.fst.gz'.format(ANGSD_DIR),
        idx = '{0}/pi_fst_sample_size_test/stats/fst/{{group}}/{{group}}_{{joint_sfs_ss}}.fst.idx'.format(ANGSD_DIR)
    log: 'logs/angsd_fst_index_toronto_fst_test/{group}_{joint_sfs_ss}_index.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 4
    resources:
        mem_mb = 4000,
        time = '01:00:00'
    params:
        fstout = '{0}/pi_fst_sample_size_test/stats/fst/{{group}}/{{group}}_{{joint_sfs_ss}}'.format(ANGSD_DIR)
    shell:
        """
        realSFS fst index {input.saf_idx} -sfs {input.joint_sfs} -fold 1 -P {threads} -whichFst 1 -fstout {params.fstout} 2> {log}
        """

rule angsd_fst_readable_toronto_fst_test:
    """
    Generate readable version of per-site Fst. Required due to file format of raw per-site Fst
    written by realSFS fst index.
    """
    input:
        rules.angsd_fst_index_toronto_fst_test.output.idx
    output:
        '{0}/pi_fst_sample_size_test/stats/fst/{{group}}/{{group}}_{{joint_sfs_ss}}_readable.fst'.format(ANGSD_DIR)
    log: 'logs/angsd_fst_readable_toronto_fst_test/{group}_{joint_sfs_ss}_readable.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    shell:
        """
        realSFS fst print {input} > {output} 2> {log}
        """

rule angsd_estimate_sfs_toronto_pi_test:
    """
    Estimate 1D folded SFS from SAF files for sample size test. 
    """
    input:
        rules.angsd_saf_likelihood_toronto_pi_fst_test.output.saf_idx
    output:
        '{0}/pi_fst_sample_size_test/{{group}}/sfs/{{group}}_{{habitat}}_{{ss}}ind.sfs'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_sfs_toronto_pi_test/{group}_{habitat}_{ss}ind_sfs.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '01:00:00'
    shell:
        """
        realSFS {input} -P {threads} -fold 1 -maxIter 2000 -seed 42 > {output} 2> {log}
        """

rule angsd_estimate_thetas_torontp_pi_test:
    """
    Estimate per-site thetas for urban/rural habitats under varying sample size combinations.
    """
    input:
        saf_idx = rules.angsd_saf_likelihood_toronto_pi_fst_test.output.saf_idx,
        sfs = rules.angsd_estimate_sfs_toronto_pi_test.output
    output:
        idx = '{0}/pi_fst_sample_size_test/stats/thetas/{{group}}/{{group}}_{{habitat}}_{{ss}}ind.thetas.idx'.format(ANGSD_DIR),
        thet = '{0}/pi_fst_sample_size_test/stats/thetas/{{group}}/{{group}}_{{habitat}}_{{ss}}ind.thetas.gz'.format(ANGSD_DIR),
    log: 'logs/angsd_estimate_thetas_torontp_pi_test/{group}_{habitat}_{ss}_thetas.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 4
    params:
        out = '{0}/pi_fst_sample_size_test/stats/thetas/{{group}}/{{group}}_{{habitat}}_{{ss}}ind'.format(ANGSD_DIR)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        realSFS saf2theta {input.saf_idx} \
            -P {threads} \
            -fold 1 \
            -sfs {input.sfs} \
            -outname {params.out} 2> {log}
        """

rule angsd_diversity_neutrality_stats_toronto_pi_test:
    """
    Estimate pi, Waterson's theta, Tajima's D, etc., for urban/rural habitat under varying sample size combinations. 
    """
    input:
        rules.angsd_estimate_thetas_torontp_pi_test.output.idx
    output:
       '{0}/pi_fst_sample_size_test/stats/thetas/{{group}}/{{group}}_{{habitat}}_{{ss}}ind.thetas.idx.pestPG'.format(ANGSD_DIR)
    log: 'logs/angsd_diversity_neutrality_stats_toronto_pi_test/{group}_{habitat}_{ss}ind_diversity_neutrality.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        thetaStat do_stat {input} 2> {log}
        """

rule pi_fst_sample_size_test_done:
    """
    Generates empty flag file signaling successful completion of test of the effects of sample sizes on
    estimation of pi and Fst
    """
    input:
        expand(rules.urban_rural_toronto_bam_lists.output, habitat=['urban', 'rural'], group = 'highCov', ss='7'),
        expand(rules.index_urban_rural_toronto_pi_fst_test_bam.output, tor_test_sample=TOR_SAMPLES_PI_FST_TEST),
        expand(rules.urban_rural_toronto_downsampled_bam_lists.output, habitat=['urban', 'rural'], ss=SS_PI_FST_TEST, group = 'lowCov'),
        expand(rules.angsd_fst_readable_toronto_fst_test.output, group = 'highCov', joint_sfs_ss = 'u7_r7'),
        expand(rules.angsd_fst_readable_toronto_fst_test.output, group = 'lowCov', joint_sfs_ss = JOINT_SFS_WILDCARDS_TOR_TEST),
        expand(rules.angsd_diversity_neutrality_stats_toronto_pi_test.output, group = 'highCov', habitat = ['urban', 'rural'], ss = '7'),
        expand(rules.angsd_diversity_neutrality_stats_toronto_pi_test.output, group = 'lowCov', habitat = ['urban', 'rural'], ss = SS_PI_FST_TEST)
    output:
        '{0}/pi_fst_sample_size_test/pi_fst_sample_size_test.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """

rule toronto_pi_fst_test_notebook:
    """
    Interactive exploration of the effects of varying sample sizes on estimation of Pi and Fst. 
    """
    input:
        rules.pi_fst_sample_size_test_done.output
    output:
        '{0}/pi_fst_sample_size_test/notebook.done'.format(ANGSD_DIR)
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/pi_fst_sample_size_test.r.ipynb"


