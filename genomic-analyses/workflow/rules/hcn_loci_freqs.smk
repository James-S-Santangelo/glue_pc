# Rules to get allele frequencies at Ac and Li loci from read count data
# Also uses ANGSD to estimate Fst of SNPs along chromosomes containing Ac and Li
# Estimates 4fold allele frequencie sin urban and rural habitats for HCN differentiation test

###############
#### SETUP ####
###############

rule split_angsd_sites_byChrom:
    """
    Split ANGSD sites file into separate sites files by chromosome
    """
    input:
        rules.convert_sites_for_angsd.output
    output:
        sites = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}.sites'.format(PROGRAM_RESOURCE_DIR),
    log: 'logs/split_random_angsd_sites_byChrom/{chrom}_{site}_split_angsd_sites_random.log'
    shell:
        """
        grep {wildcards.chrom} {input} > {output.sites} 2> {log}
        """

rule index_chromosomal_angsd_sites:
    """
    Index chromosomal ANGSD sites files for use with ANGSD
    """
    input:
        rules.split_angsd_sites_byChrom.output
    output:
        binary = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}.sites.bin'.format(PROGRAM_RESOURCE_DIR),
        idx = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}.sites.idx'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/index_random_chromosomal_angsd_sites/{chrom}_{site}_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    shell:
        """
        angsd sites index {input} 2> {log}
        """

####################################
#### HCN LOCI ALLELE FREQUENIES ####
####################################

rule read_count_data:
    """
    Calculate number of reads overlapping target region in Ac or Li locus for each sample
    """
    input:
        ancient(get_all_bams(BAM_DIR))
    output:
        '{0}/{{gene}}_read_counts.txt'.format(HCN_LOCI_DIR)
    log: 'logs/read_count_data/{gene}_counts.log'
    conda: '../envs/hcn_loci_freqs.yaml'
    params: 
        region = lambda w: 'CM019108.1:30218214-30229250 CM019108.1:30230911-30247247' if w.gene == 'li' else 'CM019103.1:19559221-19573344'
    resources:
        mem_mb = 8000,
        time = '03:00:00'
    shell:
        """
        ( for bam in {input}; do
            count=( $( samtools view $bam {params.region} | wc -l ) )
            printf "%s\t%s\n" $bam $count >> {output}
        done ) 2> {log}
        """

rule calculate_hcn_loci_frequencies:
    """
    Estimate per-sample genotype likelihoods and per-city deletion frequencies from read counts
    at Ac or Li locus. Writes two dataframes to disk per locus.
    """
    input:
        multiqc = rules.glue_dnaSeqQC_multiqc.output,
        qual_datafile = '{0}/multiqc/multiqc_data/multiqc_bamtools_stats_bamtools_stats.txt'.format(QC_DIR),
        counts = rules.read_count_data.output
    output:
        freqs = '{0}/{{gene}}_freqs.txt'.format(HCN_LOCI_DIR),
        likes = '{0}/{{gene}}_GLs.txt'.format(HCN_LOCI_DIR)
    log: 'logs/calculate_hcn_loci_frequencies/{gene}_freqs.log'
    conda: '../envs/hcn_loci_freqs.yaml'
    resources:
        mem_mb = 4000,
        time = '03:00:00'
    script:
        "../scripts/python/hcn_loci_GLs_freqs.py"

########################################
#### FST ALONG HCN LOCI CHROMOSOMES ####
########################################

# Used as null distribution against which to compare Ac and Li frequencies estimated above

rule angsd_saf_likelihood_snps_hcn_chroms:
    """
    Estimate Site Allele Frequency Likelihood (SAF) by city and habitat for SNPs along same chromosomes
    as Ac and Li loci
    """
    input:
        unpack(get_files_for_saf_estimation_snps_hcn_chroms)
    output:
        saf = temp('{0}/sfs/by_city/{{city}}/{{gene}}/{{city}}_{{gene}}_{{habitat}}_{{site}}.saf.gz'.format(ANGSD_DIR)),
        saf_idx = temp('{0}/sfs/by_city/{{city}}/{{gene}}/{{city}}_{{gene}}_{{habitat}}_{{site}}.saf.idx'.format(ANGSD_DIR)),
        saf_pos = temp('{0}/sfs/by_city/{{city}}/{{gene}}/{{city}}_{{gene}}_{{habitat}}_{{site}}.saf.pos.gz'.format(ANGSD_DIR))
    log: 'logs/angsd_saf_likelihood_snps_hcn_chroms/{city}_{gene}_{habitat}_{site}_saf.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933' 
    params:
        out = '{0}/sfs/by_city/{{city}}/{{gene}}/{{city}}_{{gene}}_{{habitat}}_{{site}}'.format(ANGSD_DIR),
        region = lambda w: 'CM019108.1' if w.gene == 'li' else 'CM019103.1' 
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 4 \
            -doMaf 1 \
            -SNP_pval 1e-6 \
            -baq 2 \
            -ref {input.ref} \
            -sites {input.sites} \
            -minQ 20 \
            -minMapQ 30 \
            -doSaf 1 \
            -anc {input.ref} \
            -r {params.region} \
            -bam {input.bams} 2> {log}
        """

rule angsd_estimate_joint_sfs_snps_hcn_chroms:
    """
    Estimated folded, two-dimensional urban-rural SFS for each city using SNPs along HCN chroms.
    """
    input:
        saf = get_habitat_saf_files_byCity_hcn_chroms,
        sites = get_sites_for_hcn_chrom_analysis
    output:
        '{0}/sfs/by_city/{{city}}/{{gene}}/{{city}}_{{gene}}_{{site}}_r_u.2dsfs'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_2dsfs_snps_hcn_chroms/{city}_{gene}_{site}.2dsfs.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933' 
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '03:00:00'
    shell:
        """
        realSFS {input.saf} \
            -sites {input.sites}\
            -maxIter 2000 \
            -seed 42 \
            -fold 1 \
            -P {threads} > {output} 2> {log}
        """

rule angsd_fst_index_snps_hcn_chroms:
    """
    Estimate per-site alphas (numerator) and betas (denominator) for Fst estimation. Done separately using 
    both Weir and Cockeram and Hudson's Fst
    """
    input: 
        saf_idx = get_habitat_saf_files_byCity_hcn_chroms,
        sites = get_sites_for_hcn_chrom_analysis,
        joint_sfs = rules.angsd_estimate_joint_sfs_snps_hcn_chroms.output
    output:
        fst = '{0}/summary_stats/fst/fst1/{{city}}/{{gene}}/{{city}}_{{gene}}_{{site}}_r_u_fst1.fst.gz'.format(ANGSD_DIR),
        idx = '{0}/summary_stats/fst/fst1/{{city}}/{{gene}}/{{city}}_{{gene}}_{{site}}_r_u_fst1.fst.idx'.format(ANGSD_DIR)
    log: 'logs/angsd_fst_index_snps_hcn_chroms/{city}_{gene}_{site}_fst1_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933' 
    threads: 4
    resources:
        mem_mb = 4000,
        time = '01:00:00'
    params:
        fstout = '{0}/summary_stats/fst/fst1/{{city}}/{{gene}}/{{city}}_{{gene}}_{{site}}_r_u_fst1'.format(ANGSD_DIR)
    shell:
        """
        realSFS fst index {input.saf_idx} \
            -sites {input.sites} \
            -sfs {input.joint_sfs} \
            -fold 1 \
            -P {threads} \
            -whichFst 1 \
            -fstout {params.fstout} 2> {log}
        """

rule angsd_fst_readable_snps_hcn_chroms:
    """
    Create readable Fst files. Required due to format of realSFS fst index output files. 
    """
    input:
        rules.angsd_fst_index_snps_hcn_chroms.output.idx
    output:
        '{0}/summary_stats/fst/fst1/{{city}}/{{gene}}/{{city}}_{{gene}}_{{site}}_r_u_fst1_readable.fst'.format(ANGSD_DIR)
    log: 'logs/angsd_fst_readable_snps_hcn_chroms/{city}_{gene}_{site}_fst1_readable.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933' 
    shell:
        """
        realSFS fst print {input} > {output} 2> {log}
        """

######################################
#### 4FOLD SNP ALLELE FREQUENCIES ####
######################################

# Used to estimate null distribution against which to compare observed differentiation in HCN

rule angsd_alleleFreqs_byCity:
    """
    Estimate minor allele frequency of 4fold SNPs across all samples in a city.
    Used to get polymophic SNPs across all samples in a city.
    """
    input:
        sites = rules.convert_sites_for_angsd.output,
        sites_idx = rules.angsd_index_allDegenerateSites.output.idx,
        ref = rules.glue_dnaSeqQC_unzip_reference.output,
        bams = rules.concat_habitat_bamLists_withinCities.output,
        chroms = config['chromosomes']
    output:
        afs = '{0}/afs/by_city/{{city}}/{{city}}_{{site}}.mafs.gz'.format(ANGSD_DIR)
    log: 'logs/angsd_alleleFreq_byCity/{city}_{site}_maf.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    params:
        out = '{0}/afs/by_city/{{city}}/{{city}}_{{site}}'.format(ANGSD_DIR)
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '08:00:00'
    shell:
        """
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 4 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -baq 2 \
            -ref {input.ref} \
            -sites {input.sites} \
            -minQ 20 \
            -minMapQ 30 \
            -anc {input.ref} \
            -rf {input.chroms} \
            -bam {input.bams} 2> {log}
        """

rule snps_forAlleleFreqs_byCity_byHabitat:
    """
    Extract positions of 4fold SNPs across all samples in a city
    """
    input:
        rules.angsd_alleleFreqs_byCity.output
    output:
        '{0}/angsd_sites/alleleFreqs/{{city}}_{{site}}_snps.sites'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/snps_forAlleleFreqs_byCity_byHabitat/{city}_{site}.log'
    shell:
        """
        zcat {input} | tail -n +2 | cut -f1,2 > {output} 2> {log}
        """

rule angsd_index_city_snps:
    """
    Index within-city 4fold SNP positions for use with ANGSD.
    """
    input:
        rules.snps_forAlleleFreqs_byCity_byHabitat.output
    output:
        idx = '{0}/angsd_sites/alleleFreqs/{{city}}_{{site}}_snps.sites.idx'.format(PROGRAM_RESOURCE_DIR),
        binary = '{0}/angsd_sites/alleleFreqs/{{city}}_{{site}}_snps.sites.bin'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/angsd_index_city_snps/{city}_{site}_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    shell:
        """
        angsd sites index {input} 2> {log}
        """

rule angsd_alleleFreqs_byCity_byHabitat:
    """
    Estimate frequency of 4fold SNPs separately in urban and rural habits for each city.
    Uses city-wide 4fold SNPs as -sites input to include SNPs that are fixed in one habitat or the other.
    """
    input:
        unpack(get_files_for_alleleFreq_estimation_byCity_byHabitat)
    output:
        afs = '{0}/afs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.mafs.gz'.format(ANGSD_DIR)
    log: 'logs/angsd_alleleFreq_byCity_byHabitat/{city}_{habitat}_{site}_saf.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    params:
        out = '{0}/afs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}'.format(ANGSD_DIR)
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '06:00:00'
    shell:
        """
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 4 \
            -doMaf 1 \
            -baq 2 \
            -ref {input.ref} \
            -sites {input.sites} \
            -minQ 20 \
            -minMapQ 30 \
            -bam {input.bams} 2> {log}
        """

rule hcn_loci_freq_done:
    """
    Generate empty flag file to signal successful completion of HCN and Ac/Li differentiation test rules.
    """
    input:
        expand(rules.calculate_hcn_loci_frequencies.output, gene=['li', 'ac']),
        expand(rules.angsd_fst_readable_snps_hcn_chroms.output, city=CITIES, gene=['li','ac'], site='4fold'),
        expand(rules.angsd_alleleFreqs_byCity_byHabitat.output, city=CITIES, habitat=HABITATS, site=['4fold'])
    output:
        '{0}/hcn_loci_freq.done'.format(HCN_LOCI_DIR)
    shell:
        """
        touch {output}
        """
