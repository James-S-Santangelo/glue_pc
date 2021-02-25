rule global_depth_pi_sfs_theta_notebook:
    input:
        rules.angsd_done.output
    output:
        '{0}/supplemental/angsd/depth/depthGlobal_histogram.png'.format(FIGURES_DIR),
        '{0}/supplemental/angsd/sfs/allSamples_allChroms_sfs_bySite.png'.format(FIGURES_DIR),
        '{0}/supplemental/population_structure/allSamples_PCA_byHabitat.png'.format(FIGURES_DIR),
        '{0}/supplemental/population_structure/allSamples_PCA_byCity.png'.format(FIGURES_DIR),
        '{0}/supplemental/population_structure/allSamples_PCA_byContinent.png'.format(FIGURES_DIR)
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/global_depth_pi_sfs_theta.r.ipynb"
