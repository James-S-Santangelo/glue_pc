# Rules to perform a small-scale test of the effects of number of sites on SFS estimation

rule angsd_saf_nSites_test:
    """
    Estimate Site Allele Frequency (SAF) likelihood file using 2 Mb region of chromosome 1. 
    Manipulate sites by changing the proportion of individuals required to have reads
    (lower proportion = less missing data = more sites)
    """
    input:
        bams = get_bams_for_angsd,
        ref = REFERENCE_GENOME
    output:
        saf = '{0}/nSites_test/{{sample_set}}/{{prop}}/{{sample_set}}_{{prop}}_allSites.saf.gz'.format(ANGSD_DIR),
        saf_idx = '{0}/nSites_test/{{sample_set}}/{{prop}}/{{sample_set}}_{{prop}}_allSites.saf.idx'.format(ANGSD_DIR),
        saf_pos = '{0}/nSites_test/{{sample_set}}/{{prop}}/{{sample_set}}_{{prop}}_allSites.saf.pos.gz'.format(ANGSD_DIR)
    log: 'logs/angsd_saf_nSites_test/{sample_set}_{prop}_allSites_angsd_saf.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    params:
        out = '{0}/nSites_test/{{sample_set}}/{{prop}}/{{sample_set}}_{{prop}}_allSites'.format(ANGSD_DIR),
        max_dp = ANGSD_MAX_DP
    threads: 6
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND*{wildcards.prop}/100 ));
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doCounts 1 \
            -setMaxDepth {params.max_dp} \
            -baq 2 \
            -ref {input.ref} \
            -minInd $MIN_IND \
            -minQ 20 \
            -minMapQ 30 \
            -doSaf 1 \
            -anc {input.ref} \
            -r CM019101.1:2000000-4000000 \
            -bam {input.bams} 2> {log}
        """

rule angsd_sfs_nSites_test:
    """
    Estimate folded SFS from SAF files using realSFS.
    """
    input:
        rules.angsd_saf_nSites_test.output.saf_idx
    output:
        '{0}/nSites_test/{{sample_set}}/{{prop}}/{{sample_set}}_{{prop}}_allSites.sfs'.format(ANGSD_DIR)
    log: 'logs/angsd_sfs_nSitest_test/{sample_set}_{prop}_sfs.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 6
    shell:
        """
        realSFS {input} -P {threads} -fold 1 -seed 42 -maxIter 2000 > {output} 2> {log}
        """

rule angsd_nSites_test_done:
    """
    Generate empty flag file signalling successful completion of small-scale test of the number of sites
    effects on the SFS estimation.
    """
    input:
        expand(rules.angsd_sfs_nSites_test.output, sample_set=['finalSamples_lowCovRemoved', 'highErrorRemoved'], prop=['30','40','50','60','70'])
    output:
        '{0}/nSites_test/angsd_nSites_test.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """

rule angsd_nSites_test_notebook:
    """
    Interactive exploration of the effects of the number of sites on SFS estimation in ANGSD
    """
    input:
        rules.angsd_nSites_test_done.output
    output:
        '{0}/supplemental/angsd_nSites_test/sfs_by_sampleSet_and_propInd.pdf'.format(FIGURES_DIR)
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/angsd_nSites_test.r.ipynb"
