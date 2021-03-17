rule create_bam_list_byCity_byHabitat:
    input:
        rules.create_bam_list_finalSamples_lowCovRemoved.output
    output:
        '{0}/bam_lists/by_city/{{city}}/{{city}}_{{habitat}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/create_bam_list/{city}_{habitat}.log'
    run:
        import os
        import pandas as pd
        df = pd.read_table(config['samples'], sep = '\t')
        df_sub = df[(df['city'] == wildcards.city) & (df['site'] == wildcards.habitat)]
        samples_city_habitat = df_sub['sample'].tolist()
        bams = open(input[0], 'r').readlines()
        with open(output[0], 'w') as f:
            for bam in bams:
                sample = os.path.basename(bam).split('_merged')[0]
                if sample in samples_city_habitat:
                    f.write('{0}'.format(bam))



rule angsd_pairwise_done:
    input:
        expand(rules.create_bam_list_byCity_byHabitat.output, city=CITIES, habitat=HABITATS)
    output:
        '{0}/angsd_pairwise.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """
