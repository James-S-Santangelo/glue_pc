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

rule aggregate:
    """
    Aggregate outputs from checkpoint
    """
    input: 
        aggregate_input
    output:
        '{0}/sfs/by_city/{{city}}/by_pop/{{city}}_{{site}}.done'.format(ANGSD_DIR)
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
        expand('{0}/sfs/by_city/{{city}}/by_pop/{{city}}_{{site}}.done'.format(ANGSD_DIR), city='Albuquerque', site='4fold')
    output:
        '{0}/angsd_byCity_byPopulation.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """

