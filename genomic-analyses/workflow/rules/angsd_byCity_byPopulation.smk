# Rules for estimating SFS and summary stats (e.g. pairwise Fst) for all populations within a city

###############################
#### SFS AND SUMMARY STATS ####
###############################

checkpoint populations_byCity_byHabitat:
    """
    Checkpoint to create file with paths to BAMs for each population in a city that still contains individuals 
    post-filtering. Checkpoint required since some populations will not have any individuals, so there may be a 
    variable number of files created. 
    """
    input:
        rules.concat_habitat_bamLists_withinCities.output
    output:
        directory('{0}/bam_lists/by_city/{{city}}/by_pop/'.format(PROGRAM_RESOURCE_DIR))
    run:
        import os
        import re
        import pandas as pd
        shell('mkdir -p {output}')
        df = pd.read_table(config['samples'], sep = '\t')
        bams = open(input[0], 'r').readlines()
        pops = []
        for line in bams:
            pop = str(re.findall('(\d+)(?=_\d+)', line)[0])
            df_sub = df[(df['city'] == wildcards.city) & (df['pop'].astype(str) == pop)]
            samples_city_habitat = df_sub['sample'].tolist()
            bam_list_out = output[0] + '/{0}_{1}_bams.list'.format(wildcards.city, pop)
            with open(bam_list_out, 'w') as f:
                for bam in bams:
                    sample = os.path.basename(bam).split('_merged')[0]
                    if sample in samples_city_habitat:
                        f.write('{0}'.format(bam))
            
rule angsd_saf_likelihood_byCity_byPopulation:
    """
    Generate Site Allele Frequency (SAF) likelihood file for each population in each city using ANGSD. 
    Uses only 4fold sites.
    """
    input:
        unpack(get_files_saf_estimation_byPopulation)
    output:
        saf = temp('{0}/sfs/by_city/{{city}}/by_pop/{{city}}_{{popu}}_{{site}}.saf.gz'.format(ANGSD_DIR)),
        saf_idx = temp('{0}/sfs/by_city/{{city}}/by_pop/{{city}}_{{popu}}_{{site}}.saf.idx'.format(ANGSD_DIR)),
        saf_pos = temp('{0}/sfs/by_city/{{city}}/by_pop/{{city}}_{{popu}}_{{site}}.saf.pos.gz'.format(ANGSD_DIR))
    log: 'logs/angsd_saf_likelihood_byCity_byPopulation/{city}_{popu}_{site}_saf.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    params:
        out = '{0}/sfs/by_city/{{city}}/by_pop/{{city}}_{{popu}}_{{site}}'.format(ANGSD_DIR)
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '02:00:00'
    wildcard_constraints:
        site='4fold',
        chrom='CM019101.1'
    shell:
        """
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 4 \
            -baq 2 \
            -ref {input.ref} \
            -sites {input.sites} \
            -minQ 20 \
            -minMapQ 30 \
            -doSaf 1 \
            -anc {input.ref} \
            -r CM019101.1 \
            -bam {input.bams} 2> {log}
        """

rule angsd_estimate_sfs_byCity_byPopulation:
    """
    Estimate folded SFS separately for each population in each city (i.e., 1D SFS) using realSFS. 
    """
    input:
        rules.angsd_saf_likelihood_byCity_byPopulation.output.saf_idx
    output:
        '{0}/sfs/by_city/{{city}}/{{city}}_{{popu}}_{{site}}.sfs'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_sfs_byCity_byPopulation/{city}_{popu}_{site}_sfs.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '01:00:00'
    shell:
        """
        realSFS {input} -P {threads} -fold 1 -maxIter 2000 -seed 42 > {output} 2> {log}
        """

rule angsd_estimate_thetas_byCity_byPopulation:
    """
    Generate per-site thetas in each population for each city from 1DSFS
    """
    input:
        saf_idx = rules.angsd_saf_likelihood_byCity_byPopulation.output.saf_idx,
        sfs = rules.angsd_estimate_sfs_byCity_byPopulation.output
    output:
        idx = '{0}/summary_stats/thetas/by_city/{{city}}/by_pop/{{city}}_{{popu}}_{{site}}.thetas.idx'.format(ANGSD_DIR),
        thet = '{0}/summary_stats/thetas/by_city/{{city}}/by_pop/{{city}}_{{popu}}_{{site}}.thetas.gz'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_thetas_byCity_byPopulation/{city}_{popu}_{site}_thetas.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 4
    params:
        out = '{0}/summary_stats/thetas/by_city/{{city}}/by_pop/{{city}}_{{popu}}_{{site}}'.format(ANGSD_DIR)
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

rule angsd_diversity_neutrality_stats_byCity_byPopulation:
    """
    Estimate pi, Waterson's theta, Tajima's D, etc. in each population in each city.
    """
    input:
        rules.angsd_estimate_thetas_byCity_byPopulation.output.idx
    output:
       '{0}/summary_stats/thetas/by_city/{{city}}/by_pop/{{city}}_{{popu}}_{{site}}.thetas.idx.pestPG'.format(ANGSD_DIR)
    log: 'logs/angsd_diversity_neutrality_stats_byCity_byPopulation/{city}_{popu}_{site}_diversity_neutrality.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        thetaStat do_stat {input} 2> {log}
        """

rule aggregate_theta:
    """
    Aggregate outputs from checkpoint
    """
    input: 
        aggregate_input_theta
    output:
        '{0}/sfs/by_city/{{city}}/by_pop/{{city}}_{{site}}_theta.done'.format(ANGSD_DIR)
    shell:
        """
        echo {input} > {output}
        """

rule angsd_byCity_byPopulation_done:
    """
    Generate empty flag file signalling successful completion on Pi and pairwise Fst estimation among populations
    within a city. 
    """
    input:
        expand('{0}/sfs/by_city/{{city}}/by_pop/{{city}}_{{site}}_theta.done'.format(ANGSD_DIR), city='Albuquerque', site='4fold')
    output:
        '{0}/angsd_byCity_byPopulation.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """

