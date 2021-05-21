# Perform Principal Components Analysis from genotype likelihoods using PCAngsd

rule pcangsd:
    """
    Perform PCA using genome-wide 4fold dengenerate sites. PCA is performed o both sample sets for comparison.  
    """
    input:
        rules.concat_angsd_gl.output
    output:
        '{0}/pcangsd/{{sample_set}}_{{site}}_maf{{maf}}_pcangsd.cov'.format(POP_STRUC_DIR),
    log: 'logs/pcangsd/{sample_set}_{site}_maf{maf}_pcangsd.log'
    container: 'library://james-s-santangelo/pcangsd/pcangsd:0.99'
    threads: 10
    params:
        out = '{0}/pcangsd/{{sample_set}}_{{site}}_maf{{maf}}_pcangsd'.format(POP_STRUC_DIR)
    wildcard_constraints:
        site = '4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '02:00:00'
    shell:
        """
        python3 /opt/pcangsd-v.0.99/pcangsd.py \
            -beagle {input} \
            -o {params.out} \
            -minMaf {wildcards.maf} \
            -threads {threads} \
            > {log}
        """

rule ngsrelate:
    """
    Estimate pairwise relatedness among samples within cities from genotype likelihoods
    """
    input:
        bams = rules.concat_habitat_bamLists_withinCities.output,
        gls = rules.angsd_gl_byCity_binary.output.gls,
        freq = rules.convert_freq_forNGSrelate.output
    output:
        '{0}/ngsrelate/{{city}}/{{city}}_{{site}}_ngsrelate_maf{{maf}}.out'.format(POP_STRUC_DIR)
    log: 'logs/ngsrelate/{city}_{site}_ngsRelate_maf{maf}.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:ngsrelate_vlatest'
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    wildcard_constraints:
        site = '4fold'
    shell:
        """
        N=$( wc -l < {input.bams} );
        ngsRelate -f {input.freq} \
            -O {output} \
            -g {input.gls} \
            -p {threads} \
            -n $N 2> {log}
        """

rule ngsadmix:
    input:
        rules.angsd_gl_byCity_beagle.output.gls
    output:
        '{0}/ngsadmix/{{city}}/K{{k}}/{{city}}_ngsadmix_K{{k}}_{{site}}_maf{{maf}}_seed{{seed}}.fopt.gz'.format(POP_STRUC_DIR),
        '{0}/ngsadmix/{{city}}/K{{k}}/{{city}}_ngsadmix_K{{k}}_{{site}}_maf{{maf}}_seed{{seed}}.qopt'.format(POP_STRUC_DIR) 
    log: 'logs/ngsadmix/{city}/K{k}/{city}_{site}_maf{maf}_K{k}_seed{seed}_ngsadmix.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 10
    params:
        out = '{0}/ngsadmix/{{city}}/K{{k}}/{{city}}_ngsadmix_K{{k}}_{{site}}_maf{{maf}}_seed{{seed}}'.format(POP_STRUC_DIR)
    wildcard_constraints:
        site = '4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '02:00:00'
    shell:
        """
        NGSadmix -likes {input} \
            -K {wildcards.k} \
            -seed {wildcards.seed} \
            -P {threads} \
            -outfiles {params.out} 2> {log}
        """

rule pop_structure_done:
    """
    Generate empty flag file signaling successful completion of PCAngsd
    """
    input:
        expand(rules.pcangsd.output, site = '4fold', maf = ['0.005', '0.01', '0.05'], sample_set=['highErrorRemoved','finalSamples_lowCovRemoved']),
        expand(rules.ngsrelate.output, site = '4fold', maf = '0.05', city = CITIES),
        expand(rules.ngsadmix.output, city = CITIES, k = NGSADMIX_K, seed=NGSADMIX_SEEDS, site = '4fold', maf = '0.05')
    output:
        '{0}/population_structure.done'.format(POP_STRUC_DIR)
    shell:
        """
        touch {output}
        """

rule global_depth_pi_sfs_theta_notebook:
    """
    Interactive exploration of global ANGSD analysis and population structure. Include analysis of depth
    along chromosome 1, and patterns of diversity and structure across all samples. 
    """
    input:
        rules.angsd_done.output,
        rules.pop_structure_done.output
    output:
        plot1 = '{0}/supplemental/angsd/depth/depthGlobal_histogram.pdf'.format(FIGURES_DIR),
        plot2 = '{0}/supplemental/angsd/sfs/sfs_allChroms_bySite_bySampleSet.pdf'.format(FIGURES_DIR),
        plot3 = '{0}/supplemental/angsd/sfs/sfs_allChroms_bySite_finalSamples_lowCovRemoved.pdf'.format(FIGURES_DIR),
        table = '{0}/tables/allChroms_diversity_bySite_bySampleSet.txt'.format(FIGURES_DIR),
        plot4 = '{0}/supplemental/population_structure/highErrorRemoved_PCA_byHabitat_bySampleSet_byMAF.pdf'.format(FIGURES_DIR),
        plot5 = '{0}/supplemental/population_structure/highErrorRemoved_PCA_byCity_byMAF.pdf'.format(FIGURES_DIR),
        plot6 = '{0}/supplemental/population_structure/highErrorRemoved_PCA_byContinent_byMAF.pdf'.format(FIGURES_DIR),
        plot7 = '{0}/supplemental/population_structure/highErrorRemoved_PCA_byRange_byMAF.pdf'.format(FIGURES_DIR),
        plot8 = '{0}/supplemental/population_structure/highErrorRemoved_PCA_byCity_maf0.05.pdf'.format(FIGURES_DIR),
        plot9 = '{0}/supplemental/population_structure/highErrorRemoved_PCA_byContinent_maf0.05.pdf'.format(FIGURES_DIR),
        plot10 = '{0}/supplemental/population_structure/highErrorRemoved_PCA_byRange_maf0.05.pdf'.format(FIGURES_DIR)
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/global_depth_pi_sfs_theta.r.ipynb"

rule relatedness_notebook:
    input:
        rules.pop_structure_done.output
    output:
        '{0}/ngsrelate/relatedness_analysis.done'.format(POP_STRUC_DIR)
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/relatedness_analysis.r.ipynb"

