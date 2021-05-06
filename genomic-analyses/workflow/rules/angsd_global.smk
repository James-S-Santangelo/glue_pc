# Uses ANGSD to estimate SFS and genotype likelihoods across the genome using all samples
# as input (i.e., a global dataset). Only really used to perform PCA of global samples. 

rule create_samples_to_remove:
    """
    Writes file with sample names for those with high alignment error rates and another 
    with samples with low coverage. Thresholds were assessed through exploratory analysis
    of QC data. 
    """
    input:
        flag = rules.multiqc.output,
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
        highError_samples = qc_data[qc_data['general_error_rate'] > 0.04]
        highError_samples['sample'].to_csv(output.error_df, header = None, index = None)
        # Samples with low coverage are those with mean coverage < 0.31X
        lowCov_samples = qc_data[qc_data['mean_coverage'] < 0.31]
        lowCov_samples['sample'].to_csv(output.lowCov_df, header = None, index = None)

rule create_bam_list_highErrorRemoved:
    """
    Create text file with paths to BAMs, excluding 5 samples with high alignment error rates
    (see qc_analysis_notebook file). These BAMs make up the "highErrorRemoved" sample set. 
    TODO: Should rename this to "full" sample set to match manuscript. 
    """
    input:
        bams = get_all_bams,
        highErr = rules.create_samples_to_remove.output.error_df,
        flag = rules.downsample_toronto_done.output
    output:
        '{0}/bam_lists/highErrorSamplesRemoved_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/create_bam_list/highErrorSamples_bam_list.log'
    run:
        import os
        import pandas as pd
        bad_samples = pd.read_table(input.highErr, header=None).iloc[:,0].tolist()
        with open(output[0], 'w') as f:
            for bam in input.bams:
                sample = os.path.basename(bam).split('_merged')[0]
                if sample not in bad_samples:
                    f.write('{0}\n'.format(bam))

rule create_bam_list_finalSamples_lowCovRemoved:
    """
    Create text file with paths to BAMs, excluding an additional 83 samples with low coverage
    (see qc_analysis_notebook file). These BAMs make up the "finalSamples_lowCovRemoved" sample set.
    TODO: Should rename this to "Reduced" sample set to match manuscript.
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
                sample = os.path.basename(bam).split('_merged')[0]
                if sample not in lowCov_samples:
                    f.write('{0}'.format(bam))

rule angsd_depth:
    """
    Calculate total depth at all sites that pass filters for each sample (.depthSample) and across all samples
    (.depthGlobal). Used to assess whether max depth cutoff used during variant calling is reasonable.
    Done separately for each of 16 chromosomes, though we only analyse chromosome 1 (see qc_analysis_notebook file).
    """
    input:
        bams = get_bams_for_angsd
    output:
        sam = '{0}/depth/{{sample_set}}/{{chrom}}/{{chrom}}_{{sample_set}}_allSites.depthSample'.format(ANGSD_DIR),
        glo = '{0}/depth/{{sample_set}}/{{chrom}}/{{chrom}}_{{sample_set}}_allSites.depthGlobal'.format(ANGSD_DIR)
    log: 'logs/angsd_depth/{sample_set}_{chrom}_angsd_depth.log'
    conda: '../envs/angsd.yaml'
    #container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    params:
        out = '{0}/depth/{{sample_set}}/{{chrom}}/{{chrom}}_{{sample_set}}_allSites'.format(ANGSD_DIR)
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE,
        time = '12:00:00'
    wildcard_constraints:
        chrom = 'CM019101.1'
    shell:
        """
        angsd -bam {input.bams} \
            -nThreads {resources.ntasks} \
            -doDepth 1 \
            -doCounts 1 \
            -r {wildcards.chrom} \
            -minMapQ 30 \
            -minQ 20 \
            -maxDepth 10000 \
            -out {params.out} 2> {log}
        """

rule angsd_saf_likelihood_allSites:
    """
    Estimate Site Allele Frequency (SAF) likelihood file separately for each of 16 chromosomes using ANGSD. 
    """
    input:
        bams = get_bams_for_angsd,
        ref = REFERENCE_GENOME
    output:
        saf = temp('{0}/sfs/{{sample_set}}/allSites/{{chrom}}/{{chrom}}_{{sample_set}}_allSites.saf.gz'.format(ANGSD_DIR)),
        saf_idx = temp('{0}/sfs/{{sample_set}}/allSites/{{chrom}}/{{chrom}}_{{sample_set}}_allSites.saf.idx'.format(ANGSD_DIR)),
        saf_pos = temp('{0}/sfs/{{sample_set}}/allSites/{{chrom}}/{{chrom}}_{{sample_set}}_allSites.saf.pos.gz'.format(ANGSD_DIR)),
        pos = '{0}/sfs/{{sample_set}}/allSites/{{chrom}}/{{chrom}}_{{sample_set}}_allSites.pos.gz'.format(ANGSD_DIR),
        counts = '{0}/sfs/{{sample_set}}/allSites/{{chrom}}/{{chrom}}_{{sample_set}}_allSites.counts.gz'.format(ANGSD_DIR)
    log: 'logs/angsd_saf_likelihood_allSites/{sample_set}_{chrom}_allSites_angsd_saf.log'
    conda: '../envs/angsd.yaml'
    #container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    params:
        out = '{0}/sfs/{{sample_set}}/allSites/{{chrom}}/{{chrom}}_{{sample_set}}_allSites'.format(ANGSD_DIR),
        max_dp = ANGSD_MAX_DP
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE,
        time = '12:00:00'
    wildcard_constraints:
        site='allSites'
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND*50/100 ));
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {resources.ntasks} \
            -doCounts 1 \
            -dumpCounts 2 \
            -setMaxDepth {params.max_dp} \
            -baq 2 \
            -ref {input.ref} \
            -minInd $MIN_IND \
            -minQ 20 \
            -minMapQ 30 \
            -doSaf 1 \
            -anc {input.ref} \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """

rule angsd_gl_allSites:
    """
    Estimate genotype likelihoods for all samples separately for each of 16 chromosomes using ANGSD.
    """
    input:
        bams = get_bams_for_angsd,
        ref = REFERENCE_GENOME
    output:
        gls = temp('{0}/gls/{{sample_set}}/allSites/{{chrom}}/{{chrom}}_{{sample_set}}_allSites_maf{{maf}}.beagle.gz'.format(ANGSD_DIR)),
        mafs = temp('{0}/gls/{{sample_set}}/allSites/{{chrom}}/{{chrom}}_{{sample_set}}_allSites_maf{{maf}}.mafs.gz'.format(ANGSD_DIR)),
    log: 'logs/angsd_gl_allSites/{chrom}_{sample_set}_allSites_maf{maf}_angsd_gl.log'
    conda: '../envs/angsd.yaml'
    #container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    params:
        out = '{0}/gls/{{sample_set}}/allSites/{{chrom}}/{{chrom}}_{{sample_set}}_allSites_maf{{maf}}'.format(ANGSD_DIR),
        max_dp = ANGSD_MAX_DP
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE,
        time = '12:00:00'
    wildcard_constraints:
        site='allSites'
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND*50/100 ))
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {resources.ntasks} \
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
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """
 
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

rule extract_angsd_allSites:
    """
    Generates ANGSD sites file for all sites used in estimation of SAF file. This avoids
    having to write separate rules or an if statement for SFS estimation of degenerate sites
    or all sites.
    """
    input:
        rules.angsd_saf_likelihood_allSites.output.pos
    output:
        sites = '{0}/angsd_sites/{{chrom}}/{{sample_set}}_{{chrom}}_Trepens_allSites.sites'.format(PROGRAM_RESOURCE_DIR),
    log: 'logs/extract_angsd_allSites/{sample_set}_{chrom}_extract_allSites.log'
    shell:
        """
        zcat {input} | tail -n +2 | cut -f1,2 > {output.sites} 2> {log} 
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

rule angsd_index_allSites:
    """
    Index ANGSD sites files for "allSites". Using sites files requires the binary representation 
    of sites files generated by this command. 
    """
    input:
        rules.extract_angsd_allSites.output
    output:
        binary = '{0}/angsd_sites/{{chrom}}/{{sample_set}}_{{chrom}}_Trepens_allSites.sites.bin'.format(PROGRAM_RESOURCE_DIR),
        idx = '{0}/angsd_sites/{{chrom}}/{{sample_set}}_{{chrom}}_Trepens_allSites.sites.idx'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/angsd_index/{sample_set}_{chrom}_allSites_index.log'
    conda: '../envs/angsd.yaml'
    shell:
        """
        angsd sites index {input} 2> {log}
        """

rule angsd_index_degenerate:
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

rule angsd_estimate_sfs:
    """
    Estimate folded SFS using realSFS separately for each chromosome.
    Done for "allSites", "0fold", and "4fold" sites.
    """
    input:
        unpack(angsd_sfs_input) 
    output:
        temp('{0}/sfs/{{sample_set}}/{{site}}/{{chrom}}/{{chrom}}_{{sample_set}}_{{site}}.sfs'.format(ANGSD_DIR))
    log: 'logs/angsd_estimate_sfs/{chrom}_{sample_set}_{site}_sfs.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 60000,
        time = '03:00:00'
    shell:
        """
        realSFS {input.saf_idx} -P {threads} -sites {input.sites} -fold 1 -maxIter 2000 -seed 42 > {output} 2> {log}
        """

rule angsd_estimate_thetas:
    """
    Estimate per-site thetas for all site types and seaprately for all 16 chromosomes.
    """
    input:
        unpack(angsd_estimate_thetas_input)
    output:
        idx = '{0}/summary_stats/thetas/{{sample_set}}/{{site}}/{{chrom}}/{{chrom}}_{{sample_set}}_{{site}}.thetas.idx'.format(ANGSD_DIR),
        thet = '{0}/summary_stats/thetas/{{sample_set}}/{{site}}/{{chrom}}/{{chrom}}_{{sample_set}}_{{site}}.thetas.gz'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_thetas/{chrom}_{sample_set}_{site}_thetas.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 10
    params:
        out = '{0}/summary_stats/thetas/{{sample_set}}/{{site}}/{{chrom}}/{{chrom}}_{{sample_set}}_{{site}}'.format(ANGSD_DIR)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '03:00:00'
    shell:
        """
        realSFS saf2theta {input.saf_idx} \
            -fold 1 \
            -P {threads} \
            -sfs {input.sfs} \
            -sites {input.sites} \
            -outname {params.out} 2> {log}
        """

rule angsd_diversity_neutrality_stats:
    """
    Estimate pi, Waterson's theta, Tajima's D, etc., for all site types and seaprately for each of 16 chromosomes.
    """
    input:
        rules.angsd_estimate_thetas.output.idx
    output:
       temp( '{0}/summary_stats/thetas/{{sample_set}}/{{site}}/{{chrom}}/{{chrom}}_{{sample_set}}_{{site}}.thetas.idx.pestPG'.format(ANGSD_DIR))
    log: 'logs/angsd_diversity_neutrality_stats/{chrom}_{sample_set}_{site}_diversity_neutrality.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '03:00:00'
    shell:
        """
        thetaStat do_stat {input} 2> {log}
        """
 
rule concat_angsd_stats:
    """
    Concatenate diversity and neutrality statistics for each of 16 chromosomes into single file.
    """
    input:
        get_angsd_stats_toConcat
    output:
        '{0}/summary_stats/thetas/{{sample_set}}/{{site}}/allChroms_{{sample_set}}_{{site}}_diversityNeutrality.thetas.idx.pestPG'.format(ANGSD_DIR)
    log: 'logs/concat_angsd_stats_specificSites/{sample_set}_{site}_concat_angsd_stats.log'
    shell:
        """
        first=1
        for f in {input}; do
            if [ "$first"  ]; then
                cat "$f"
                first=
            else
                cat "$f"| tail -n +2
            fi
        done > {output} 2> {log}
        """

rule concat_sfs:
    """
    Concatenate folded SFS files for each of 16 chromosomes into single file (one chromosome per line)
    """
    input:
        get_angsd_sfs_toConcat
    output:
        '{0}/sfs/{{sample_set}}/{{site}}/{{sample_set}}_{{site}}_allChroms.cat'.format(ANGSD_DIR)
    log: 'logs/concat_sfs/{sample_set}_{site}_concat_sfs.log'
    shell:
        """
        cat {input} > {output} 2> {log}
        """

rule sum_sfs:
    """
    Generate genome-wide folded SFS by summing allele frequency bins across each of 16 chromosomes.
    (e.g., sum all singletons, sum all doubletons, etc.)
    """
    input:
        rules.concat_sfs.output
    output:
        '{0}/sfs/{{sample_set}}/{{site}}/{{sample_set}}_{{site}}_allChroms.sfs'.format(ANGSD_DIR)
    log: 'logs/sum_sfs/{sample_set}_{site}.log'
    run:
        import pandas as pd
        import sys
        with open(log[0], 'w') as logfile:
            sys.stderr = logfile
            sfs_allChroms = pd.read_table(input[0], sep = ' ', header = None)
            sfs_sum = sfs_allChroms.sum(axis=0) 
            sfs_sum.to_csv(output[0], sep = '\t', header = None)
        

rule files_for_angsd_site_subset:
    """
    Generate files for later subsetting of MAF and GL files by site type (i.e., allSites, 0fold, 4fold). 
    One file to subset GLs and one for MAFs (positions are encoded separately for these two filetypes)
    """
    input:
        rules.split_angsd_sites_byChrom.output
    output:
        gl = '{0}/angsd_sites/{{chrom}}/{{chrom}}_{{site}}_gl.positions'.format(PROGRAM_RESOURCE_DIR),
        maf = '{0}/angsd_sites/{{chrom}}/{{chrom}}_{{site}}_maf.positions'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/files_for_angsd_site_subset/{chrom}_{site}_subsetFile.log'
    shell:
        """
        ( sed 's/\t/_/g' {input} > {output.gl};
        cut -f2 {input} > {output.maf} ) 2> {log}
        """

rule subset_angsd_gl:
    """
    Subset GLs for each site type (allSites, 0fold, 4fold). Done separately for each of 16 chromosomes.
    """
    input:
        sites = rules.split_angsd_sites_byChrom.output,
        subset = rules.files_for_angsd_site_subset.output.gl,
        gl = rules.angsd_gl_allSites.output.gls
    output:
        '{0}/gls/{{sample_set}}/{{site}}/{{chrom}}/{{chrom}}_{{sample_set}}_{{site}}_maf{{maf}}.beagle.gz'.format(ANGSD_DIR)
    log: 'logs/subset_angsd_gl/{chrom}_{sample_set}_{site}_maf{maf}_subset_gl.log'
    params:
        unzip_out='{0}/gls/{{sample_set}}/{{site}}/{{chrom}}/{{chrom}}_{{sample_set}}_{{site}}_maf{{maf}}.beagle'.format(ANGSD_DIR)
    wildcard_constraints:
        site='0fold|4fold'
    shell:
        """
        ( zcat {input.gl} | sed -n 1p > {params.unzip_out} &&
        zcat {input.gl} | awk 'NR==FNR{{A[$1]; next}} $1 in A' {input.subset} - >> {params.unzip_out} &&
        gzip {params.unzip_out} ) 2> {log}
        """

rule subset_angsd_maf:
    """
    Subset MAFs for each site type (allSites, 0fold, 4fold). Done separately for each of 16 chromosomes.
    """
    input:
        sites = rules.split_angsd_sites_byChrom.output,
        subset = rules.files_for_angsd_site_subset.output.maf,
        mafs = rules.angsd_gl_allSites.output.mafs
    output:
        '{0}/gls/{{sample_set}}/{{site}}/{{chrom}}/{{chrom}}_{{sample_set}}_{{site}}_maf{{maf}}.mafs.gz'.format(ANGSD_DIR)
    log: 'logs/subset_angsd_maf/{chrom}_{sample_set}_{site}_maf{maf}_subset_maf.log'
    params:
        unzip_out='{0}/gls/{{sample_set}}/{{site}}/{{chrom}}/{{chrom}}_{{sample_set}}_{{site}}_maf{{maf}}.mafs'.format(ANGSD_DIR)
    wildcard_constraints:
        site='0fold|4fold'
    shell:
        """
        ( zcat {input.mafs} | sed -n 1p  > {params.unzip_out} &&
        zcat {input.mafs} | awk 'NR==FNR{{A[$1]; next}} $2 in A' {input.subset} - >> {params.unzip_out} &&
        gzip {params.unzip_out} ) 2> {log}
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

rule extract_sample_angsd:
    """
    Create text file with order of samples in ANGSD output files. Order is the same as the sample
    order in the BAM list used as input for ANGSD.
    """
    input:
        get_bams_for_angsd
    output:
        '{0}/angsd_{{sample_set}}_order.txt'.format(PROGRAM_RESOURCE_DIR)
    run:
        with open(output[0], 'w') as fout:
            with open(input[0], 'r') as fin:
                for line in fin:
                    sline = line.strip().split('/')
                    bam = sline[-1]
                    sample = bam.split('_merged')[0]
                    fout.write('{0}\n'.format(sample))

rule angsd_done:
    """
    Generate empty flag file signalling successful completion of global SFS and GL estimation across all samples. 
    """
    input:
        expand(rules.angsd_depth.output, chrom='CM019101.1', sample_set=['highErrorRemoved','finalSamples_lowCovRemoved']),
        expand(rules.concat_angsd_stats.output, site=['allSites','0fold','4fold'], sample_set=['highErrorRemoved','finalSamples_lowCovRemoved']),
        expand(rules.sum_sfs.output, site=['allSites','0fold','4fold'], sample_set=['highErrorRemoved','finalSamples_lowCovRemoved']),
        expand(rules.concat_angsd_gl.output, sample_set=['highErrorRemoved','finalSamples_lowCovRemoved'], site=['allSites','0fold','4fold'], maf=['0.005','0.01','0.05']),
        expand(rules.concat_angsd_mafs.output, sample_set=['highErrorRemoved','finalSamples_lowCovRemoved'], site=['allSites','0fold','4fold'], maf=['0.005','0.01','0.05']),
        expand(rules.extract_sample_angsd.output, sample_set=['highErrorRemoved','finalSamples_lowCovRemoved'])
    output:
        '{0}/angsd.done'.format(ANGSD_DIR)
    shell:
        "touch {output}"

