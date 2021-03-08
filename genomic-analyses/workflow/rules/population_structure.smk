rule pcangsd:
    input:
        rules.concat_angsd_gl.output
    output:
        '{0}/pcangsd/{{sample_set}}_{{site}}_maf{{maf}}_pcangsd.cov'.format(POP_STRUC_DIR),
    log: 'logs/pcangsd/{sample_set}_{site}_maf{maf}_pcangsd.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:pcangsd_v0.99'
    threads: 10
    params:
        out = '{0}/pcangsd/{{sample_set}}_{{site}}_maf{{maf}}_pcangsd'.format(POP_STRUC_DIR)
    wildcard_constraints:
        site = '4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        python3 /opt/pcangsd/pcangsd.py \
            -beagle {input} \
            -o {params.out} \
            -threads {threads} \
            > {log}
        """

rule pop_structure_done:
    input:
        expand(rules.pcangsd.output, site = '4fold', maf = '0.05', sample_set=['highErrorRemoved','finalSamples_lowCovRemoved'])
    output:
        '{0}/population_structure.done'.format(POP_STRUC_DIR)
    shell:
        """
        touch {output}
        """

rule global_depth_pi_sfs_theta_notebook:
    input:
        rules.angsd_done.output,
        rules.pop_structure_done.output
    output:
        plot1 = '{0}/supplemental/angsd/depth/depthGlobal_histogram.pdf'.format(FIGURES_DIR),
        plot2 = '{0}/supplemental/angsd/sfs/allSamples_allChroms_sfs_bySite.pdf'.format(FIGURES_DIR),
        plot3 = '{0}/supplemental/population_structure/allSamples_PCA_byHabitat.pdf'.format(FIGURES_DIR),
        plot4 = '{0}/supplemental/population_structure/allSamples_PCA_byCity.pdf'.format(FIGURES_DIR),
        plot5 = '{0}/supplemental/population_structure/allSamples_PCA_byContinent.pdf'.format(FIGURES_DIR),
        plot6 = '{0}/supplemental/population_structure/allSamples_PCA_byRange.pdf'.format(FIGURES_DIR)
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/global_depth_pi_sfs_theta.r.ipynb"
