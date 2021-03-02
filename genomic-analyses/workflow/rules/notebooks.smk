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

rule qc_analysis_notebook:
    input:
        rules.multiqc.output,
        rules.single_sample_sfs_done.output
    output:
        plot1 = '{0}/supplemental/qc/mean_coverage_histogram.pdf'.format(FIGURES_DIR),
        plot2 = '{0}/supplemental/qc/general_error_rate_histogram.pdf'.format(FIGURES_DIR),
        plot3 = '{0}/supplemental/qc/aligned_vs_coverage_by_propVar_highErrorRemoved.pdf'.format(FIGURES_DIR),
        plot4 = '{0}/supplemental/qc/aligned_vs_coverage_by_propVar_highQualOnly.pdf'.format(FIGURES_DIR),
        table1 = '{0}/tables/lowQualitySamples_by_city_and_habitat.txt'.format(FIGURES_DIR),
        error_df = '{0}/highErrorRate_toRemove.txt'.format(PROGRAM_RESOURCE_DIR),
        low_cov_df = '{0}/lowCoverageSamples_toRemove.txt'.format(PROGRAM_RESOURCE_DIR)
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/qc_analysis.r.ipynb"
