# Rules to get allele frequencies at Ac and Li loci from read count data
# Also uses ANGSD to estimate Fst of SNPs along chromosomes containing Ac and Li

####################################
#### HCN LOCI ALLELE FREQUENIES ####
####################################

rule read_count_data:
    """
    Calculate number of reads overlapping target region in Ac or Li locus for each sample
    """
    input:
        expand(rules.samtools_markdup.output.bam, sample = SAMPLES)
    output:
        '{0}/{{gene}}_read_counts.txt'.format(HCN_LOCI_DIR)
    log: 'logs/read_count_data/{gene}_counts.log'
    conda: '../envs/hcn_loci_freqs.yaml'
    params: 
        region = lambda w: 'CM019108.1:30218214-30229250 CM019108.1:30230911-30247247' if w.gene == 'li' else 'CM019103.1:19559221-19573344' 
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
        multiqc = rules.multiqc.output,
        counts = rules.read_count_data.output
    output:
        freqs = '{0}/{{gene}}_freqs.txt'.format(HCN_LOCI_DIR),
        likes = '{0}/{{gene}}_GLs.txt'.format(HCN_LOCI_DIR)
    log: 'logs/calculate_hcn_loci_frequencies/{gene}_freqs.log'
    conda: '../envs/hcn_loci_freqs.yaml'
    script:
        "../scripts/python/hcn_loci_GLs_freqs.py"

########################################
#### FST ALONG HCN LOCI CHROMOSOMES ####
########################################

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
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    params:
        out = '{0}/sfs/by_city/{{city}}/{{gene}}/{{city}}_{{gene}}_{{habitat}}_{{site}}'.format(ANGSD_DIR),
        region = lambda w: 'CM019108.1' if w.gene == 'li' else 'CM019103.1' 
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
        get_habitat_saf_files_byCity_hcn_chroms
    output:
        '{0}/sfs/by_city/{{city}}/{{gene}}/{{city}}_{{gene}}_{{site}}_r_u.2dsfs'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_2dsfs_snps_hcn_chroms/{city}_{gene}_{site}.2dsfs.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '01:00:00'
    shell:
        """
        realSFS {input} -maxIter 2000 -seed 42 -fold 1 -P {threads} > {output} 2> {log}
        """

rule angsd_fst_index_snps_hcn_chroms:
    """
    Estimate per-site alphas (numerator) and betas (denominator) for Fst estimation. Done separately using 
    both Weir and Cockeram and Hudson's Fst
    """
    input: 
        saf_idx = get_habitat_saf_files_byCity_hcn_chroms, 
        joint_sfs = rules.angsd_estimate_joint_sfs_snps_hcn_chroms.output
    output:
        fst = '{0}/summary_stats/fst/fst1/{{city}}/{{gene}}/{{city}}_{{gene}}_{{site}}_r_u_fst1.fst.gz'.format(ANGSD_DIR),
        idx = '{0}/summary_stats/fst/fst1/{{city}}/{{gene}}/{{city}}_{{gene}}_{{site}}_r_u_fst1.fst.idx'.format(ANGSD_DIR)
    log: 'logs/angsd_fst_index_snps_hcn_chroms/{city}_{gene}_{site}_fst1_index.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 4
    resources:
        mem_mb = 4000,
        time = '01:00:00'
    params:
        fstout = '{0}/summary_stats/fst/fst1/{{city}}/{{gene}}/{{city}}_{{gene}}_{{site}}_r_u_fst1'.format(ANGSD_DIR)
    shell:
        """
        realSFS fst index {input.saf_idx} -sfs {input.joint_sfs} -fold 1 -P {threads} -whichFst 1 -fstout {params.fstout} 2> {log}
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
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    shell:
        """
        realSFS fst print {input} > {output} 2> {log}
        """

rule hcn_loci_freq_done:
    """
    Generate empty flag file to signal successful completion of GL and deletion frequency estimation
    """
    input:
        expand(rules.calculate_hcn_loci_frequencies.output, gene=['li', 'ac']),
        expand(rules.angsd_fst_readable_snps_hcn_chroms.output, city=CITIES, gene=['li','ac'], site='4fold')
    output:
        '{0}/hcn_loci_freq_done'.format(HCN_LOCI_DIR)
    shell:
        """
        touch {output}
        """

rule hcn_loci_notebook:
    input:
        rules.hcn_loci_freq_done.output
    output:
        '{0}/hcn_loci_notebook.done'.format(HCN_LOCI_DIR)
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/hcn_loci.r.ipynb"
