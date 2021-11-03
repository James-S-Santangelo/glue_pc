# Population structure analyses. 

################
#### GLOBAL ####
################

rule pcangsd_allSamples:
    """
    Perform PCA using genome-wide 4fold dengenerate sites using all samples from all cities.
    """
    input:
        rules.concat_angsd_gl.output
    output:
        '{0}/pcangsd/allSamples/allSamples_{{site}}_maf{{maf}}_pcangsd.cov'.format(POP_STRUC_DIR),
    log: 'logs/pcangsd_allSamples/allSamples_{site}_maf{maf}_pcangsd.log'
    container: 'library://james-s-santangelo/pcangsd/pcangsd:0.99'
    threads: 10
    params:
        out = '{0}/pcangsd/allSamples/allSamples_{{site}}_maf{{maf}}_pcangsd'.format(POP_STRUC_DIR)
    wildcard_constraints:
        site = '4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 12000,
        time = '03:00:00'
    shell:
        """
        python3 /opt/pcangsd-v.0.99/pcangsd.py \
            -beagle {input} \
            -o {params.out} \
            -minMaf {wildcards.maf} \
            -threads {threads} \
            -iter 5000 \
            > {log}
        """

###################
#### ADMIXTURE ####
###################

rule pcangsd_byCity:
    """
    Perform PCA by city using genome-wide 4fold dengenerate sites and estimate admixture proportions
    """
    input:
        rules.angsd_gl_byCity.output.gls
    output:
        cov = '{0}/pcangsd/by_city/{{city}}/{{city}}_{{site}}_maf{{maf}}_pcangsd.cov'.format(POP_STRUC_DIR),
        pi = '{0}/pcangsd/by_city/{{city}}/{{city}}_{{site}}_maf{{maf}}_pcangsd.pi.npy'.format(POP_STRUC_DIR)
    log: 'logs/pcangsd_byCity/{city}_{site}_maf{maf}_pcangsd.log'
    container: 'library://james-s-santangelo/pcangsd/pcangsd:0.99'
    threads: 6
    params:
        out = '{0}/pcangsd/by_city/{{city}}/{{city}}_{{site}}_maf{{maf}}_pcangsd'.format(POP_STRUC_DIR)
    wildcard_constraints:
        site = '4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '02:00:00'
    shell:
        """
        python3 /opt/pcangsd/pcangsd.py \
            -beagle {input} \
            -o {params.out} \
            -minMaf {wildcards.maf} \
            -threads {threads} \
            -admix \
            -admix_seed 42 \
            -iter 5000 \
            > {log}
        """
                
rule pop_structure_done:
    """
    Generate empty flag file signaling successful completion of PCAngsd
    """
    input:
        expand(rules.pcangsd_allSamples.output, site = '4fold', maf = ['0.05']),
        expand(rules.pcangsd_byCity.output, site = '4fold', maf = '0.05', city = CITIES)
    output:
        '{0}/population_structure.done'.format(POP_STRUC_DIR)
    shell:
        """
        touch {output}
        """

rule admixture_notebook:
    input:
        rules.pop_structure_done.output
    output:
        '{0}/admixture_notebook.done'.format(POP_STRUC_DIR)
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/admixture.r.ipynb"

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
