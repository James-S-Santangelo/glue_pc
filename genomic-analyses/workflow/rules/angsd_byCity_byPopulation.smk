# Rules for estimating SFS and summary stats (e.g. pairwise Fst) for all populations within a city

###############################
#### SFS AND SUMMARY STATS ####
###############################

checkpoint populations_byCity_byHabitat:
    input:
        rules.create_bam_list_byCity_byHabitat.output
    output:
        touch('{0}/bam_lists/by_city/{{city}}/by_pop/{{habitat}}_flag'.format(PROGRAM_RESOURCE_DIR))
    params:
        '{0}/bam_lists/by_city/{{city}}/by_pop/{{habitat}}'.format(PROGRAM_RESOURCE_DIR)
    run:
        import os
        import re
        import pandas as pd
        shell('mkdir -p {params}')
        df = pd.read_table(config['samples'], sep = '\t')
        bams = open(input[0], 'r').readlines()
        pops = []
        for line in bams:
            pop = str(re.findall('(\d+)(?=_\d+)', line)[0])
            df_sub = df[(df['city'] == wildcards.city) & (df['site'] == wildcards.habitat) & (df['pop'].astype(str) == pop)]
            print(df_sub.head())
            samples_city_habitat = df_sub['sample'].tolist()
            bam_list_out = params[0] + '/{0}_{1}_{2}_bams.list'.format(wildcards.city, wildcards.habitat, pop)
            with open(bam_list_out, 'w') as f:
                for bam in bams:
                    sample = os.path.basename(bam).split('_merged')[0]
                    if sample in samples_city_habitat:
                        f.write('{0}'.format(bam))
            
rule append_hello:
    input: 
        aggregate_input
    output:
        '{0}/bam_lists/by_city/{{city}}/by_pop/{{habitat}}/processed.txt'.format(PROGRAM_RESOURCE_DIR)
    shell:
        """
        echo {input} 'hello' > {output}
        """

rule angsd_byCity_byPopulation_done:
    input:
        expand('{0}/bam_lists/by_city/{{city}}/by_pop/{{habitat}}/processed.txt'.format(PROGRAM_RESOURCE_DIR), city='Albuquerque', habitat='u')
    output:
        '{0}/angsd_byCity_byPopulation.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """




# rule create_bam_list_byCity_byPopulation:
#     """
#     Create text file with paths to BAM files in each population by city. Uses only "finalSample_lowCovRemoved"
#     sample set.
#     """
#     input:
#         rules.create_bam_list_finalSamples_lowCovRemoved.output
#     output:
#         '{0}/bam_lists/by_city/{{city}}/by_pop/{{city}}_{{habitat}}_{{pop}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
#     run:
#         import os
#         import pandas as pd
#         df = pd.read_table(config['samples'], sep = '\t')
#         df_sub = df[(df['city'] == wildcards.city) & (df['site'] == wildcards.habitat) & (df['pop'] == wildcards.pop)]
#         samples_city_habitat = df_sub['sample'].tolist()
#         bams = open(input[0], 'r').readlines()
#         with open(output[0], 'w') as f:
#             for bam in bams:
#                 sample = os.path.basename(bam).split('_merged')[0]
#                 if sample in samples_city_habitat:
#                     f.write('{0}'.format(bam))

