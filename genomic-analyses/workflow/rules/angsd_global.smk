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
        error_df = '{0}/highErrorRate_toRemove.txt'.format(PROGRAM_RESOURCE_DIR),
        lowCov_df = '{0}/lowCoverageSamples_toRemove.txt'.format(PROGRAM_RESOURCE_DIR)
    run:
        import pandas as pd
        qc_data = pd.read_table(input.qc_data, sep = '\t')
        qc_data['sample'] = qc_data['Sample'].str.extract('(\w+_\d+_\d+)')
        cols = ['sample', 'mean_coverage', 'general_error_rate']
        qc_data = qc_data[cols]
        # Samples with high ealignment errors have error rates > 0.04
        highError_samples = qc_data[qc_data['general_error_rate'] >= 0.03]
        highError_samples['sample'].to_csv(output.error_df, header = None, index = None)
        # Samples with low coverage are those with mean coverage < 0.31X
        lowCov_samples = qc_data[qc_data['mean_coverage'] < 0.31]
        lowCov_samples['sample'].to_csv(output.lowCov_df, header = None, index = None)

rule create_bam_list_highErrorRemoved:
    """
    Create text file with paths to BAMs, excluding samples with high alignment error rates
    """
    input:
        bams = get_all_bams,
        highErr = rules.create_samples_to_remove.output.error_df
    output:
        '{0}/bam_lists/highErrorSamplesRemoved_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/create_bam_list/highErrorSamples_bam_list.log'
    run:
        import os
        import pandas as pd
        bad_samples = pd.read_table(input.highErr, header=None).iloc[:,0].tolist()
        with open(output[0], 'w') as f:
            for bam in input.bams:
                if wildcards.sample not in bad_samples:
                    f.write('{0}\n'.format(bam))

rule create_bam_list_finalSamples_lowCovRemoved:
    """
    Create text file with paths to BAMs, excluding an additional 83 samples with low coverage
    """
    input:
        allSamples = rules.create_bam_list_highErrorRemoved.output,
        lowCov = rules.create_samples_to_remove.output.lowCov_df
    output:
        '{0}/bam_lists/finalSamples_lowCovRemoved_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/create_bam_list/finalSamples_lowCovRemoved_bam_list.log'
    run:
        import os
        import pandas as pd
        lowCov_samples = pd.read_table(input.lowCov, header=None).iloc[:,0].tolist()
        bams = open(input.allSamples[0], 'r').readlines()
        with open(output[0], 'w') as f:
            for bam in bams:
                if wildcards.sample not in lowCov_samples:
                    f.write('{0}'.format(bam))

rule convert_sites_for_angsd:
    """
    Converts 0-based BED files for 0fold and 4fold sites output from Degeneracy to 1-based 
    site coordinates required by ANGSD
    """
    input:
        get_bed_to_subset
    output:
        '{0}/angsd_sites/Trepens_{{site}}.sites'.format(PROGRAM_RESOURCE_DIR) 
    log: 'logs/convert_sites_for_angsd/convert_{site}_for_angsd.log'
    wildcard_constraints:
        site='0fold|4fold'
    shell:
        """
        awk '{{print $1"\t"$2+1}}' {input} > {output} 2> {log}
        """

rule split_angsd_sites_byChrom:
    """
    Split ANGSD sites file by chromosome for SFS estimation. 
    """
    input:
        rules.convert_sites_for_angsd.output
    output:
        sites = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}.sites'.format(PROGRAM_RESOURCE_DIR),
    log: 'logs/split_angsd_sites_byChrom/{chrom}_{site}_split_angsd_sites.log'
    shell:
        """
        grep {wildcards.chrom} {input} > {output.sites} 2> {log}
        """

rule index_angsd_sites:
    """
    Index ANGSD sites files for 0fold and 4fold sites. Using sites files requires the binary representation 
    of sites files generated by this command. 
    """
    input:
        rules.split_angsd_sites_byChrom.output
    output:
        binary = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}.sites.bin'.format(PROGRAM_RESOURCE_DIR),
        idx = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}.sites.idx'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/angsd_index/{chrom}_{site}_index.log'
    conda: '../envs/angsd.yaml'
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
        bams = rules.subset_bams_degeneracy.output,
        ref = rules.glue_dnaSeqQC_unzip_reference.output,
        sites = rules.split_angsd_sites_byChrom.output,
        idx = rules.index_angsd_sites.output
    output:
        gls = temp('{0}/gls/allSamples/{{site}}/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}.beagle.gz'.format(ANGSD_DIR)),
        mafs = temp('{0}/gls/allSamples/{{site}}/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}.mafs.gz'.format(ANGSD_DIR)),
    log: 'logs/angsd_gl_allSamples/{chrom}_{site}_maf{maf}_angsd_gl.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933' 
    params:
        out = '{0}/gls/allSamples/{{site}}/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}'.format(ANGSD_DIR),
        max_dp = ANGSD_MAX_DP
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
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
            -minMaf {wildcards.maf} \
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
        '{0}/gls/{{sample_set}}/{{site}}/allChroms_{{sample_set}}_{{site}}_maf{{maf}}.beagle.gz'.format(ANGSD_DIR)
    log: 'logs/concat_angsd_gl/concat_angsd_gl_{sample_set}_{site}_maf{maf}.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
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
        '{0}/gls/{{sample_set}}/{{site}}/allChroms_{{sample_set}}_{{site}}_maf{{maf}}.mafs.gz'.format(ANGSD_DIR)
    log: 'logs/concat_angsd_mafs/concat_angsd_mafs_{sample_set}_{site}_maf{maf}.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
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

# rule extract_sample_angsd:
#     """
#     Create text file with order of samples in ANGSD output files. Order is the same as the sample
#     order in the BAM list used as input for ANGSD.
#     """
#     input:
#         get_bams_for_angsd
#     output:
#         '{0}/angsd_{{sample_set}}_order.txt'.format(PROGRAM_RESOURCE_DIR)
#     run:
#         with open(output[0], 'w') as fout:
#             with open(input[0], 'r') as fin:
#                 for line in fin:
#                     sline = line.strip().split('/')
#                     bam = sline[-1]
#                     sample = bam.split('_merged')[0]
#                     fout.write('{0}\n'.format(sample))
# 
rule angsd_done:
    """
    Generate empty flag file signalling successful completion of global SFS and GL estimation across all samples. 
    """
    input:
        expand('{0}/{{site}}/{{sample}}_{{site}}.bam.bai'.format(BAM_DIR), sample=SAMPLES, site=['4fold'])
#         expand(rules.angsd_depth.output, chrom='CM019101.1', sample_set=['highErrorRemoved','finalSamples_lowCovRemoved']),
#         expand(rules.concat_angsd_stats.output, site=['allSites','0fold','4fold'], sample_set=['highErrorRemoved','finalSamples_lowCovRemoved']),
#         expand(rules.sum_sfs.output, site=['allSites','0fold','4fold'], sample_set=['highErrorRemoved','finalSamples_lowCovRemoved']),
#         expand(rules.concat_angsd_gl.output, sample_set=['highErrorRemoved','finalSamples_lowCovRemoved'], site=['allSites','0fold','4fold'], maf=['0.005','0.01','0.05']),
#         expand(rules.concat_angsd_mafs.output, sample_set=['highErrorRemoved','finalSamples_lowCovRemoved'], site=['allSites','0fold','4fold'], maf=['0.005','0.01','0.05']),
#         expand(rules.extract_sample_angsd.output, sample_set=['highErrorRemoved','finalSamples_lowCovRemoved'])
    output:
        '{0}/angsd.done'.format(ANGSD_DIR)
    shell:
        "touch {output}"

