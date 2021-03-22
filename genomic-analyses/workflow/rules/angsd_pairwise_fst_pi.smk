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
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    shell:
        """
        angsd sites index {input} 2> {log}
        """

rule angsd_saf_likelihood_byCity_byHabitat:
    input:
        bams = rules.create_bam_list_byCity_byHabitat.output,
        sites_idx = rules.angsd_index_degenerate_allChroms.output,
        sites = rules.convert_sites_for_angsd.output, 
        ref = REFERENCE_GENOME,
        chroms = config['chromosomes']
    output:
        saf = temp('{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.saf.gz'.format(ANGSD_DIR)),
        saf_idx = temp('{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.saf.idx'.format(ANGSD_DIR)),
        saf_pos = temp('{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.saf.pos.gz'.format(ANGSD_DIR))
    log: 'logs/angsd_saf_likelihood_byCity_byHabitat/{city}_{habitat}_{site}_saf.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    params:
        out = '{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}'.format(ANGSD_DIR)
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
            -rf {input.chroms} \
            -bam {input.bams} 2> {log}
        """

rule angsd_estimate_joint_sfs_byCity:
    input:
        get_habitat_saf_files_byCity
    output:
        '{0}/sfs/by_city/{{city}}/{{city}}_{{site}}_r_u.2dsfs'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_2dsfs_byCity/{city}_{site}.2dsfs.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '01:00:00'
    shell:
        """
        realSFS {input} -maxIter 2000 -seed 42 -fold 1 -P {threads} > {output} 2> {log}
        """

rule angsd_fst_index:
    input: 
        saf_idx = get_habitat_saf_files_byCity,
        joint_sfs = rules.angsd_estimate_joint_sfs_byCity.output
    output:
        fst = '{0}/summary_stats/fst/fst{{fst}}/{{city}}/{{city}}_{{site}}_r_u_fst{{fst}}.fst.gz'.format(ANGSD_DIR),
        idx = '{0}/summary_stats/fst/fst{{fst}}/{{city}}/{{city}}_{{site}}_r_u_fst{{fst}}.fst.idx'.format(ANGSD_DIR)
    log: 'logs/angsd_fst_index/{city}_{site}_fst{fst}_index.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 4
    resources:
        mem_mb = 4000,
        time = '01:00:00'
    params:
        fstout = '{0}/summary_stats/fst/fst{{fst}}/{{city}}/{{city}}_{{site}}_r_u_fst{{fst}}'.format(ANGSD_DIR)
    shell:
        """
        realSFS fst index {input.saf_idx} -sfs {input.joint_sfs} -fold 1 -P {threads} -whichFst {wildcards.fst} -fstout {params.fstout} 2> {log}
        """

rule angsd_fst_readable:
    input:
        rules.angsd_fst_index.output.idx
    output:
        '{0}/summary_stats/fst/fst{{fst}}/{{city}}/{{city}}_{{site}}_r_u_fst{{fst}}_readable.fst'.format(ANGSD_DIR)
    log: 'logs/angsd_fst_readable/{city}_{site}_fst{fst}_readable.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    shell:
        """
        realSFS fst print {input} > {output} 2> {log}
        """

rule angsd_estimate_sfs_byCity_byHabitat:
    input:
        rules.angsd_saf_likelihood_byCity_byHabitat.output.saf_idx
    output:
        '{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.sfs'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_sfs_byCity_byHabitat/{city}_{habitat}_{site}_sfs.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '01:00:00'
    shell:
        """
        realSFS {input} -P {threads} -fold 1 -maxIter 2000 -seed 42 > {output} 2> {log}
        """

rule angsd_estimate_thetas_byCity_byHabitat:
    input:
        saf_idx = rules.angsd_saf_likelihood_byCity_byHabitat.output.saf_idx,
        sfs = rules.angsd_estimate_sfs_byCity_byHabitat.output
    output:
        idx = '{0}/summary_stats/thetas/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.thetas.idx'.format(ANGSD_DIR),
        thet = '{0}/summary_stats/thetas/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.thetas.gz'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_thetas_byCity_byHabitat/{city}_{habitat}_{site}_thetas.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 4
    params:
        out = '{0}/summary_stats/thetas/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}'.format(ANGSD_DIR)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        realSFS saf2theta {input.saf_idx} \
            -P {threads} \
            -sfs {input.sfs} \
            -outname {params.out} 2> {log}
        """

rule angsd_diversity_neutrality_stats_byCity_byHabitat:
    input:
        rules.angsd_estimate_thetas_byCity_byHabitat.output.idx
    output:
       '{0}/summary_stats/thetas/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.thetas.idx.pestPG'.format(ANGSD_DIR)
    log: 'logs/angsd_diversity_neutrality_stats_byCity_byHabitat/{city}_{habitat}_{site}_diversity_neutrality.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        thetaStat do_stat {input} 2> {log}
        """

rule angsd_pairwise_done:
    input:
        expand(rules.angsd_fst_readable.output, city=CITIES, site=['4fold'], fst=['0', '1']),
        expand(rules.angsd_diversity_neutrality_stats_byCity_byHabitat.output, city=CITIES, habitat=HABITATS, site=['4fold'])
    output:
        '{0}/angsd_pairwise.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """

rule pairwise_pi_fst_notebook:
    input:
        rules.angsd_pairwise_done.output
    output:
        '{0}/supplemental/fst/wc_hudson_fst_corr.pdf'.format(FIGURES_DIR)
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/pairwise_pi_fst.r.ipynb"
