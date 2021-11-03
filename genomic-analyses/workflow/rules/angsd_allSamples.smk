# Uses ANGSD to estimate SFS and genotype likelihoods across the genome using all samples
# as input (i.e., a global dataset). Only really used to perform PCA of global samples. 

###############
#### SETUP ####
###############

rule subset_bams_degeneracy:
    """
    Subset BAMs for all samples around 4fold sites to speed up ANGSD computations.
    """
    input:
        unpack(get_subset_bams_degeneracy_input)
    output:
        '{0}/{{site}}/{{sample}}_{{site}}.bam'.format(BAM_DIR)
    log: 'logs/subset_bams_degenerate/{sample}_{site}_subset.log'
    conda: '../envs/degeneracy.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time = '01:00:00'
    shell:
        """
        samtools view -bh -L {input.regions} {input.bam} > {output} 2> {log}
        """

rule index_degenerate_bam:
    input:
        rules.subset_bams_degeneracy.output
    output:
        '{0}/{{site}}/{{sample}}_{{site}}.bam.bai'.format(BAM_DIR)
    log: 'logs/index_degenerate_bam/{sample}_{site}_index.log'
    conda: '../envs/degeneracy.yaml'
    shell:
        """
        samtools index {input} 2> {log}
        """

rule create_samples_to_remove:
    """
    Writes file with sample names for those with high alignment error rates and another 
    with samples with low coverage. Thresholds were assessed through exploratory analysis
    of QC data. 
    """
    input:
        flag = rules.glue_dnaSeqQC_multiqc.output,
        qc_data = '{0}/multiqc/multiqc_data/multiqc_qualimap_bamqc_genome_results_qualimap_bamqc.txt'.format(QC_DIR)
    output: 
        error_df = '{0}/highErrorRate_toRemove.txt'.format(PROGRAM_RESOURCE_DIR)
    run:
        import pandas as pd
        qc_data = pd.read_table(input.qc_data, sep = '\t')
        qc_data['sample'] = qc_data['Sample'].str.extract('(\w+_\d+_\d+)')
        cols = ['sample', 'mean_coverage', 'general_error_rate']
        qc_data = qc_data[cols]
        # Samples with high ealignment errors have error rates > 0.04
        highError_samples = qc_data[qc_data['general_error_rate'] >= 0.03]
        highError_samples['sample'].to_csv(output.error_df, header = None, index = None)

rule create_bam_list_finalSamples:
    """
    Create text file with paths to BAMs, excluding samples with high alignment error rates
    """
    input:
        bams = get_all_bams,
        highErr = rules.create_samples_to_remove.output.error_df
    output:
        '{0}/bam_lists/finalSamples_{{site}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/create_bam_list/finalSamples_{site}_bam_list.log'
    run:
        import os
        import re
        import pandas as pd
        from contextlib import redirect_stderr
        with open(log[0], 'w') as stderr, redirect_stderr(stderr):
            bad_samples = pd.read_table(input.highErr, header=None).iloc[:,0].tolist()
            with open(output[0], 'w') as f:
                for bam in input.bams:
                    search = re.search('^(.+)(?=_\w)', os.path.basename(bam))
                    sample = search.group(1)
                    if sample not in bad_samples:
                        f.write('{0}\n'.format(bam))

rule convert_sites_for_angsd:
    """
    Converts 0-based BED files for 0fold and 4fold sites output from Degeneracy to 1-based 
    site coordinates required by ANGSD
    """
    input:
        get_bed
    output:
        '{0}/angsd_sites/Trepens_{{site}}.sites'.format(PROGRAM_RESOURCE_DIR) 
    log: 'logs/convert_sites_for_angsd/convert_{site}_for_angsd.log'
    wildcard_constraints:
        site='0fold|4fold'
    shell:
        """
        awk '{{print $1"\t"$2+1}}' {input} > {output} 2> {log}
        """

rule select_random_degenerate_sites:
    """
    Random select `params.nSites` degenerate (i.e., 4fold or 0fold) sites from across the genome
    """
    input:
        rules.convert_sites_for_angsd.output
    output:
        '{0}/angsd_sites/Trepens_{{site}}_random.sites'.format(PROGRAM_RESOURCE_DIR)
    params:
        nSites = 250000
    run:
        import random
        sites = open(input[0], 'r').readlines()
        random.seed(42)
        rand_sites = sorted(random.sample(sites, k = int(params.nSites)))
        with open(output[0], 'w') as fout:
            for site in rand_sites:
                fout.write(site)

rule angsd_index_random_degen_sites:
    """
    Index randomly selected genome-wide degenerate sites for use with ANGSD
    """
    input:
        rules.select_random_degenerate_sites.output
    output:
        binary = '{0}/angsd_sites/Trepens_{{site}}_random.sites.bin'.format(PROGRAM_RESOURCE_DIR),
        idx = '{0}/angsd_sites/Trepens_{{site}}_random.sites.idx'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/angsd_index_random_degen_sites/random_{site}_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    shell:
        """
        angsd sites index {input} 2> {log}
        """

rule split_random_angsd_sites_byChrom:
    """
    Split ANGSD sites file into separate sites files by chromosome
    """
    input:
        rules.select_random_degenerate_sites.output
    output:
        sites = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}_random.sites'.format(PROGRAM_RESOURCE_DIR),
    log: 'logs/split_random_angsd_sites_byChrom/{chrom}_{site}_split_angsd_sites_random.log'
    shell:
        """
        grep {wildcards.chrom} {input} > {output.sites} 2> {log}
        """

rule index_random_chromosomal_angsd_sites:
    """
    Index chromosomal ANGSD sites files for use with ANGSD
    """
    input:
        rules..output
    output:
        binary = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}_random.sites.bin'.format(PROGRAM_RESOURCE_DIR),
        idx = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}_random.sites.idx'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/index_random_chromosomal_angsd_sites/{chrom}_{site}_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933' 
    shell:
        """
        angsd sites index {input} 2> {log}
        """

##############################
#### GENOTYPE LIKELIHOODS ####
##############################
        
rule angsd_gl_allSamples:
    """
    Estimate genotype likelihoods for all samples separately for each of 16 chromosomes using ANGSD.
    """
    input:
        bams = rules.create_bam_list_finalSamples.output,
        ref = rules.glue_dnaSeqQC_unzip_reference.output,
        sites = rules.split_random_angsd_sites_byChrom.output,
        idx = rules.index_random_chromosomal_angsd_sites.output
    output:
        gls = temp('{0}/gls/allSamples/{{site}}/{{chrom}}/{{chrom}}_{{site}}.beagle.gz'.format(ANGSD_DIR)),
        mafs = temp('{0}/gls/allSamples/{{site}}/{{chrom}}/{{chrom}}_{{site}}.mafs.gz'.format(ANGSD_DIR)),
    log: 'logs/angsd_gl_allSamples/{chrom}_{site}_angsd_gl.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933' 
    params:
        out = '{0}/gls/allSamples/{{site}}/{{chrom}}/{{chrom}}_{{site}}'.format(ANGSD_DIR),
        max_dp = ANGSD_MAX_DP
    threads: 12
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 30000,
        time = '12:00:00'
    wildcard_constraints:
        site='4fold'
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND*50/100 ))
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doGlf 2 \
            -doMajorMinor 1 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -doCounts 1 \
            -setMaxDepth {params.max_dp} \
            -baq 2 \
            -ref {input.ref} \
            -minInd $MIN_IND \
            -minQ 20 \
            -minMapQ 30 \
            -sites {input.sites} \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """
 
rule concat_angsd_gl:
    """
    Concatenated GLs from all 16 chromosomes into single file. Done separately for each site type.
    """
    input:
        get_angsd_gl_toConcat
    output:
        '{0}/gls/allSamples/{{site}}/allChroms_{{site}}.beagle.gz'.format(ANGSD_DIR)
    log: 'logs/concat_angsd_gl/allSamples_{site}_concat.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933' 
    shell:
        """
        first=1
        for f in {input}; do
            if [ "$first"  ]; then
                zcat "$f"
                first=
            else
                zcat "$f"| tail -n +2
            fi
        done | bgzip -c > {output} 2> {log}
        """

rule concat_angsd_mafs:
    """
    Concatenate MAF files for each of 16 chromosomes into single file. Done separately for each site type.
    """
    input:
        get_angsd_maf_toConcat
    output:
        '{0}/gls/allSamples/{{site}}/allChroms_{{site}}.mafs.gz'.format(ANGSD_DIR)
    log: 'logs/concat_angsd_mafs/allSamples_{site}_concat.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933' 
    shell:
        """
        first=1
        for f in {input}; do
            if [ "$first"  ]; then
                zcat "$f"
                first=
            else
                zcat "$f"| tail -n +2
            fi
        done | bgzip -c > {output} 2> {log}
        """

rule angsd_allSamples_done:
    """
    Generate empty flag file signalling successful completion of global SFS and GL estimation across all samples. 
    """
    input:
        expand(rules.concat_angsd_gl.output, site=['4fold']),
        expand(rules.concat_angsd_mafs.output, site=['4fold'])
    output:
        '{0}/angsd_allSamples.done'.format(ANGSD_DIR)
    shell:
        "touch {output}"

