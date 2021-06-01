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
    """
    Estimate admixture proportions directly from genotype likelihoods.

    """
    input:
        rules.angsd_gl_byCity_beagle.output.gls
    output:
        fopt = '{0}/ngsadmix/{{city}}/K{{k}}/{{city}}_ngsadmix_K{{k}}_{{site}}_maf{{maf}}_seed{{seed}}.fopt.gz'.format(POP_STRUC_DIR),
        qopt = '{0}/ngsadmix/{{city}}/K{{k}}/{{city}}_ngsadmix_K{{k}}_{{site}}_maf{{maf}}_seed{{seed}}.qopt'.format(POP_STRUC_DIR),
        lf = '{0}/ngsadmix/{{city}}/K{{k}}/{{city}}_ngsadmix_K{{k}}_{{site}}_maf{{maf}}_seed{{seed}}.log'.format(POP_STRUC_DIR)
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

rule pcangsd_byCity:
    """
    Perform PCA by city using genome-wide 4fold dengenerate sites and estimate admixture proportions
    """
    input:
        rules.angsd_gl_byCity_beagle.output.gls
    output:
        cov = '{0}/pcangsd/by_city/{{city}}/{{city}}_{{site}}_maf{{maf}}_pcangsd.cov'.format(POP_STRUC_DIR),
        adm = '{0}/pcangsd/by_city/{{city}}/{{city}}_{{site}}_maf{{maf}}_pcangsd.admix.Q.npy'.format(POP_STRUC_DIR)
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
        python3 /opt/pcangsd-v.0.99/pcangsd.py \
            -beagle {input} \
            -o {params.out} \
            -minMaf {wildcards.maf} \
            -threads {threads} \
            -admix \
            -admix_seed 42 \
            > {log}
        """

rule logfile_for_clumpak:
    """
    Create Inputfile for CLUMPAK containing Log likelihood values of NGSadmix runs for each K
    """
    input:
        get_ngsadmix_logfiles_byCity
    output:
        '{0}/for_clumpak/{{city}}_ngsadmix_logfile_for_clumpak.txt'.format(PROGRAM_RESOURCE_DIR)
    run:
        import re
        with open(output[0], 'w') as fout:
            for lf in input:
                # Get K
                m1 = re.search('(?<=_K)(\d+)', lf)
                k = m1.group(1)
                # Get likelihood
                line = open(lf, 'r').readlines()[-1]  # Likelihood always on last line
                m2 = re.search('(?<=like=)(-?\d+.\d+)', line)
                like = m2.group(1)
                fout.write('{0}\t{1}\n'.format(k, like))

rule clumpak_best_k_by_evanno:
    """
    Find optimal K value by city using Evanno method, as implemented in CLUMPAK
    """
    input:
        rules.logfile_for_clumpak.output
    output:
        directory('{0}/bestKbyEvanno/{{city}}'.format(POP_STRUC_DIR))
    log: 'logs/clumpak_best_k_by_evanno/{city}_evanno.log'
    container: 'library://james-s-santangelo/clumpak/clumpak:1.1'
    params:
        outdir = lambda wildcards: '{0}/bestKbyEvanno/{1}'.format(POP_STRUC_DIR, wildcards.city)
    wildcard_constraints:
        city='|'.join(CITIES)
    resources:
        mem_mb = 1000,
        time = '01:00:00'
    shell:
        """
        perl /opt/bin/BestKByEvanno.pl --id {wildcards.city}_out \
            --d {params.outdir} \
            --f {input} \
            --inputtype lnprobbyk 2>&1 > {log}
        """

rule bestK_byCity:
    input:
        rules.clumpak_best_k_by_evanno.output
    output:
        '{0}/bestKbyEvanno/bestK_files/{{city}}_evanno_bestK.txt'.format(POP_STRUC_DIR)
    run:
        import re
        with open (output[0], 'w') as fout:
            lines = open(input[0] + '/output.log', 'r').readlines()
            for l in lines:
                if 'Optimal K by Evanno' in l:
                    m = re.search('(?<=: )(\d+)', l)
                    k = m.group(1)
                    fout.write('city\tbest_k\n'.format(wildcards.city))  # Header
                    fout.write('{0}\t{1}\n'.format(wildcards.city, k))
                
rule pop_structure_done:
    """
    Generate empty flag file signaling successful completion of PCAngsd
    """
    input:
        expand(rules.pcangsd.output, site = '4fold', maf = ['0.005', '0.01', '0.05'], sample_set=['highErrorRemoved','finalSamples_lowCovRemoved']),
        expand(rules.ngsrelate.output, site = '4fold', maf = '0.05', city = CITIES),
        expand(rules.pcangsd_byCity.output, site = '4fold', maf = '0.05', city = CITIES),
        expand(rules.bestK_byCity.output, city=CITIES)
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

rule relatedness_notebook:
    input:
        rules.pop_structure_done.output
    output:
        '{0}/ngsrelate/relatedness_analysis.done'.format(POP_STRUC_DIR)
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/relatedness_analysis.r.ipynb"

