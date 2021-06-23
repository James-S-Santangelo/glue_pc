# Rules to permute urban and rural samples by city to estimate null distribution of pi and Fst
# Used to test whether urban-rural difference in Pi and Fst are different than expected by chance.

rule create_random_bam_list_byCity_byHabitat:
    """
    Create text files with paths to BAM files for each city. 'urban' and 'rural' samples are randomly
    selected from all possible bams for a city. Sample sizes are the same as those in the observed samples.
    """
    input:
        unpack(get_urban_rural_bam_lists)
    output:
        urban = '{0}/bam_lists/by_city/{{city}}/randomized/{{city}}_randU_seed{{seed}}_bams.list'.format(PROGRAM_RESOURCE_DIR),
        rural = '{0}/bam_lists/by_city/{{city}}/randomized/{{city}}_randR_seed{{seed}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    run:
        import random
        urban_bams = list(open(input.urban_bams, 'r'))
        rural_bams = list(open(input.rural_bams, 'r'))
        urban_n = len(urban_bams)
        rural_n = len(rural_bams)
        all_bams = urban_bams + rural_bams
        random.seed(int(wildcards.seed))
        randU = random.sample(all_bams, urban_n)
        randR = [bam for bam in all_bams if not bam in randU] 
        with open(output.urban, 'w') as uout:
            for bam in randU:
                uout.write(bam)
        with open(output.rural, 'w') as rout:
            for bam in randR:
                rout.write(bam)

rule angsd_permuted_saf_likelihood_byCity_byHabitat:
    """
    Generate Site Allele Frequency (SAF) likelihood file for each habitat in each city using ANGSD. 
    Uses only 4fold sites. Uses permuted urban and rural samples 
    """
    input:
        unpack(get_files_for_permuted_saf_estimation)
    output:
        saf = temp('{0}/sfs/by_city/{{city}}/randomized/{{city}}_{{habitat}}_{{site}}_seed{{seed}}.saf.gz'.format(ANGSD_DIR)),
        saf_idx = temp('{0}/sfs/by_city/{{city}}/randomized/{{city}}_{{habitat}}_{{site}}_seed{{seed}}.saf.idx'.format(ANGSD_DIR)),
        saf_pos = temp('{0}/sfs/by_city/{{city}}/randomized/{{city}}_{{habitat}}_{{site}}_seed{{seed}}.saf.pos.gz'.format(ANGSD_DIR))
    log: 'logs/angsd_permuted_saf_likelihood_byCity_byHabitat/{city}_{habitat}_{site}_seed{seed}_saf.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    params:
        out = '{0}/sfs/by_city/{{city}}/randomized/{{city}}_{{habitat}}_{{site}}_seed{{seed}}'.format(ANGSD_DIR)
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    wildcard_constraints:
        site='4fold'
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

rule angsd_estimate_permuted_joint_sfs_byCity:
    """
    Estimated folded, two-dimensional urban-rural SFS for each city using realSFS. Uses 4fold sites.
    Uses permuted urban and rural samples.
    """
    input:
        get_habitat_saf_files_byCity_permuted
    output:
        '{0}/sfs/by_city/{{city}}/randomized/{{city}}_{{site}}_seed{{seed}}_r_u.2dsfs'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_permuted_2dsfs_byCity/{city}_{site}_seed{{seed}}.2dsfs.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        realSFS {input} -maxIter 2000 -seed 42 -fold 1 -P {threads} > {output} 2> {log}
        """

rule angsd_permuted_fst_index:
    """
    Estimate per-site alphas (numerator) and betas (denominator) for Hudson's Fst estimator.
    Uses permuted urban and rural samples.
    """
    input: 
        saf_idx = get_habitat_saf_files_byCity_permuted,
        joint_sfs = rules.angsd_estimate_permuted_joint_sfs_byCity.output
    output:
        fst = '{0}/summary_stats/fst/fst1/{{city}}/randomized/{{city}}_{{site}}_seed{{seed}}_r_u.fst.gz'.format(ANGSD_DIR),
        idx = '{0}/summary_stats/fst/fst1/{{city}}/randomized/{{city}}_{{site}}_seed{{seed}}_r_u.fst.idx'.format(ANGSD_DIR)
    log: 'logs/angsd_permuted_fst_index/{city}_{site}_seed{seed}_index.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 4
    resources:
        mem_mb = 4000,
        time = '01:00:00'
    params:
        fstout = '{0}/summary_stats/fst/fst1/{{city}}/randomized/{{city}}_{{site}}_seed{{seed}}_r_u'.format(ANGSD_DIR)
    shell:
        """
        realSFS fst index {input.saf_idx} -sfs {input.joint_sfs} -fold 1 -P {threads} -whichFst 1 -fstout {params.fstout} 2> {log}
        """

rule angsd_permuted_fst_readable:
    """
    Create readable Fst files. Required due to format of realSFS fst index output files. Uses permuted 
    urban and rural samples.
    """
    input:
        rules.angsd_permuted_fst_index.output.idx
    output:
        '{0}/summary_stats/fst/fst1/{{city}}/randomized/{{city}}_{{site}}_seed{{seed}}_r_u_readable.fst'.format(ANGSD_DIR)
    log: 'logs/angsd_permuted_fst_readable/{city}_{site}_seed{{seed}}_readable.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    shell:
        """
        realSFS fst print {input} > {output} 2> {log}
        """

rule angsd_estimate_permuted_sfs_byCity_byHabitat:
    """
    Estimate folded SFS separately for each habitat in each city (i.e., 1D SFS) using realSFS. Uses permuted 
    urban and rural samples.
    """
    input:
        rules.angsd_permuted_saf_likelihood_byCity_byHabitat.output.saf_idx
    output:
        '{0}/sfs/by_city/{{city}}/randomized/{{city}}_{{habitat}}_{{site}}_seed{{seed}}.sfs'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_permuted_sfs_byCity_byHabitat/{city}_{habitat}_{site}_seed{{seed}}_sfs.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        realSFS {input} -P {threads} -fold 1 -maxIter 2000 -seed 42 > {output} 2> {log}
        """

rule angsd_estimate_permuted_thetas_byCity_byHabitat:
    """
    Generate per-site thetas in each habitat for each city from 1DSFS. Uses permuted urban and rural samples.
    """
    input:
        saf_idx = rules.angsd_permuted_saf_likelihood_byCity_byHabitat.output.saf_idx,
        sfs = rules.angsd_estimate_permuted_sfs_byCity_byHabitat.output
    output:
        idx = '{0}/summary_stats/thetas/by_city/{{city}}/randomized/{{city}}_{{habitat}}_{{site}}_seed{{seed}}.thetas.idx'.format(ANGSD_DIR),
        thet = '{0}/summary_stats/thetas/by_city/{{city}}/randomized/{{city}}_{{habitat}}_{{site}}_seed{{seed}}.thetas.gz'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_permuted_thetas_byCity_byHabitat/{city}_{habitat}_{site}_seed{{seed}}_thetas.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 4
    params:
        out = '{0}/summary_stats/thetas/by_city/{{city}}/randomized/{{city}}_{{habitat}}_{{site}}_seed{{seed}}'.format(ANGSD_DIR)
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

rule angsd_permuted_diversity_neutrality_stats_byCity_byHabitat:
    """
    Estimate pi, Waterson's theta, Tajima's D, etc. in each habitat in each city. Uses permuted urban and rural samples
    """
    input:
        rules.angsd_estimate_permuted_thetas_byCity_byHabitat.output.idx
    output:
       '{0}/summary_stats/thetas/by_city/{{city}}/randomized/{{city}}_{{habitat}}_{{site}}_seed{{seed}}.thetas.idx.pestPG'.format(ANGSD_DIR)
    log: 'logs/angsd_permuted_diversity_neutrality_stats_byCity_byHabitat/{city}_{habitat}_{site}_seed{{seed}}_diversity_neutrality.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        thetaStat do_stat {input} 2> {log}
        """

rule angsd_byCity_byHabitat_permuted_done:
    """
    Generate empty flag file signalling successful completion of pairwise pi and Fst analysis
    """
    input:
        expand(rules.angsd_permuted_fst_readable.output, city=CITIES, site=['4fold'], seed=BOOT_SEEDS),
        expand(rules.angsd_permuted_diversity_neutrality_stats_byCity_byHabitat.output, city=CITIES, habitat=HABITATS, site=['4fold'], seed=BOOT_SEEDS)
    output:
        '{0}/angsd_byCity_byHabitat_permuted.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """

rule pi_fst_byCity_byHabitat_notebook:
    """
    Interactive exploration of pairwise urban-rural pi and Fst analysis.
    """
    input:
        rules.angsd_byCity_byHabitat_done.output,
        rules.angsd_byCity_byHabitat_permuted_done.output
    output:
        '{0}/pi_fst_byCity_byHabitat_notebook.done'.format(FIGURES_DIR)
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/pi_fst_byCity_byHabitat.r.ipynb"
