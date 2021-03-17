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

rule angsd_index_degenerate_allChroms:
    input:
        rules.convert_sites_for_angsd.output
    output:
        binary = '{0}/angsd_sites/Trepens_{{site}}.sites.bin'.format(PROGRAM_RESOURCE_DIR),
        idx = '{0}/angsd_sites/Trepens_{{site}}.sites.idx'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/angsd_index/allChroms_{site}_index.log'
    conda: '../envs/angsd.yaml'
    shell:
        """
        angsd sites index {input} 2> {log}
        """

rule angsd_saf_likelihood_byCity_byHabitat:
    input:
        bams = rules.create_bam_list_byCity_byHabitat.output,
        sites = rules.angsd_index_degenerate_allChroms.output,
        ref = REFERENCE_GENOME
    output:
        saf = temp('{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.saf.gz'.format(ANGSD_DIR)),
        saf_idx = temp('{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.saf.idx'.format(ANGSD_DIR)),
        saf_pos = temp('{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.saf.pos.gz'.format(ANGSD_DIR))
    log: 'logs/angsd_saf_likelihood_byCity_byHabitat/{city}_{habitat}_{site}_saf.log'
    conda: '../envs/angsd.yaml'
    params:
        out = '{0}/sfs/by_city/{{city}}/{{city}}/_{{habitat}}_{{site}}'.format(ANGSD_DIR)
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '06:00:00'
    wildcard_constraints:
        site='4fold'
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
            -bam {input.bams} 2> {log}
        """

rule angsd_pairwise_done:
    input:
        expand(rules.angsd_saf_likelihood_byCity_byHabitat.output, city=CITIES, habitat=HABITATS, site=['4fold'])
    output:
        '{0}/angsd_pairwise.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """
