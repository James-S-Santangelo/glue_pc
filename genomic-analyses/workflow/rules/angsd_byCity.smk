# Rules to get genotype likelihoods across all individuals per city. 
# Used for population structure analyses. 

##############################
#### GENOTYPE LIKELIHOODS ####
##############################

rule concat_habitat_bamLists_withinCities:
    """
    Concatenate urban and rural sample BAM lists within cities. Generates a single file with
    the paths to all of the BAM files for samples within a city
    """
    input:
        get_bamLists_toConcat
    output:
        '{0}/bam_lists/by_city/{{city}}/{{city}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/concat_habitat_bamLists_withinCities/{city}_concat.log'
    shell:
        """
        cat {input} > {output} 2> {log}
        """

rule angsd_index_degenerate_allChroms:
    """
    Indexes ANGSD sites files containing genome-wide positions. Since only 4fold sites are going to be
    used here, we can generate the SAF file across all chromosomes simultaneously.
    """
    input:
        rules.convert_sites_for_angsd.output
    output:
        binary = '{0}/angsd_sites/Trepens_{{site}}.sites.bin'.format(PROGRAM_RESOURCE_DIR),
        idx = '{0}/angsd_sites/Trepens_{{site}}.sites.idx'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/angsd_index/allChroms_{site}_index.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    wildcard_constraints:
        site='4fold'
    shell:
        """
        angsd sites index {input} 2> {log}
        """

rule angsd_gl_byCity_binary:
    """
    Estimate genotype likelihoods jointly for all samples within a city. Output is binary format for use
    with NGSrelate
    """
    input:
        bams = rules.concat_habitat_bamLists_withinCities.output,
        sites = rules.convert_sites_for_angsd.output,
        sites_idx = rules.angsd_index_degenerate_allChroms.output,
        ref = REFERENCE_GENOME,
        chroms = config['chromosomes']
    output:
        gls = '{0}/gls/by_city/{{city}}/{{city}}_{{site}}_binaryGLs_maf{{maf}}.glf.gz'.format(ANGSD_DIR),
        mafs = '{0}/gls/by_city/{{city}}/{{city}}_{{site}}_binaryGLs_maf{{maf}}.mafs.gz'.format(ANGSD_DIR),
        pos = '{0}/gls/by_city/{{city}}/{{city}}_{{site}}_binaryGLs_maf{{maf}}.glf.pos.gz'.format(ANGSD_DIR)
    log: 'logs/angsd_gl_byCity_binary/{city}_{site}_maf{maf}_binaryGL.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    params:
        out = '{0}/gls/by_city/{{city}}/{{city}}_{{site}}_binaryGLs_maf{{maf}}'.format(ANGSD_DIR)
    threads: 6
    wildcard_constraints:
        site = '4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 5000,
        time = '06:00:00'
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND * 60/100 ));
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doGlf 3 \
            -doMajorMinor 1 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -doCounts 1 \
            -baq 2 \
            -ref {input.ref} \
            -minInd $MIN_IND \
            -minQ 20 \
            -minMapQ 30 \
            -minMaf {wildcards.maf} \
            -sites {input.sites} \
            -rf {input.chroms} \
            -bam {input.bams} 2> {log}
        """

rule angsd_gl_byCity_beagle:
    """
    Estimate genotype likelihoods jointly for all samples within a city. Output is beagle format for use
    with NGSadmix
    """
    input:
        bams = rules.concat_habitat_bamLists_withinCities.output,
        sites = rules.convert_sites_for_angsd.output,
        sites_idx = rules.angsd_index_degenerate_allChroms.output,
        ref = REFERENCE_GENOME,
        chroms = config['chromosomes']
    output:
        gls = '{0}/gls/by_city/{{city}}/{{city}}_{{site}}_beagleGLs_maf{{maf}}.beagle.gz'.format(ANGSD_DIR),
        mafs = '{0}/gls/by_city/{{city}}/{{city}}_{{site}}_beagleGLs_maf{{maf}}.mafs.gz'.format(ANGSD_DIR)
    log: 'logs/angsd_gl_byCity_beagle/{city}_{site}_maf{maf}_beagleGL.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    params:
        out = '{0}/gls/by_city/{{city}}/{{city}}_{{site}}_beagleGLs_maf{{maf}}'.format(ANGSD_DIR)
    threads: 6
    wildcard_constraints:
        site = '4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 5000,
        time = '06:00:00' 
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND * 60/100 ));
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doGlf 2 \
            -doMajorMinor 1 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -doCounts 1 \
            -baq 2 \
            -ref {input.ref} \
            -minInd $MIN_IND \
            -minQ 20 \
            -minMapQ 30 \
            -minMaf {wildcards.maf} \
            -sites {input.sites} \
            -rf {input.chroms} \
            -bam {input.bams} 2> {log}
        """

rule convert_freq_forNGSrelate:
    input:
        rules.angsd_gl_byCity_binary.output.mafs
    output:
        '{0}/gls/by_city/{{city}}/{{city}}_{{site}}_ngsRelate_maf{{maf}}.freqs'.format(ANGSD_DIR)
    log: 'logs/convert_freq_forNGSrelate/{city}_{site}_maf{maf}_convert_freqs.log'
    shell:
        """
        zcat {input} | cut -f6 | sed 1d > {output} 2> {log}
        """

##############
#### POST ####
##############

rule angsd_byCity_done:
    input:
        expand(rules.angsd_gl_byCity_binary.output, city=CITIES, site='4fold', maf='0.05'),
        expand(rules.convert_freq_forNGSrelate.output, city=CITIES, site='4fold', maf='0.05')
    output:
        '{0}/angsd_byCity.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """
        
