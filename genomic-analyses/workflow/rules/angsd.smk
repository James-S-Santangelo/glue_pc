rule create_bam_list_highErrorRemoved:
    input:
        glue_bams = expand(rules.samtools_markdup.output.bam, sample=SAMPLES),
        tor_bams = expand(rules.downsample_toronto_bam.output, tor_sample=TOR_SAMPLES),
        highErr = rules.qc_analysis_notebook.output.error_df,
        flag = rules.downsample_toronto_done.output
    output:
        '{0}/bam_lists/highErrorSamplesRemoved_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/create_bam_list/highErrorSamples_bam_list.log'
    run:
        import os
        import pandas as pd
        bad_samples = pd.read_table(input.highErr, header=None).iloc[:,0].tolist()
        with open(output[0], 'w') as f:
            for bam in input.glue_bams:
                sample = os.path.basename(bam).split('_merged')[0]
                if sample not in bad_samples:
                    f.write('{0}\n'.format(bam))
            for bam in input.tor_bams:
                sample = os.path.basename(bam).split('_merged')[0]
                f.write('{0}\n'.format(bam))

rule create_bam_list_finalSamples_lowCovRemoved:
    input:
        allSamples = rules.create_bam_list_highErrorRemoved.output,
        lowCov = rules.qc_analysis_notebook.output.low_cov_df
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
    input:
        bams = get_bams_for_angsd
    output:
        sam = '{0}/depth/{{sample_set}}/{{chrom}}/{{chrom}}_{{sample_set}}_allSites.depthSample'.format(ANGSD_DIR),
        glo = '{0}/depth/{{sample_set}}/{{chrom}}/{{chrom}}_{{sample_set}}_allSites.depthGlobal'.format(ANGSD_DIR)
    log: 'logs/angsd_depth/{sample_set}_{chrom}_angsd_depth.log'
    conda: '../envs/angsd.yaml'
    params:
        out = '{0}/depth/{{sample_set}}/{{chrom}}/{{chrom}}_{{sample_set}}_allSites'.format(ANGSD_DIR)
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE,
        time = '12:00:00'
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
    input:
        bams = get_bams_for_angsd,
        ref = REFERENCE_GENOME
    output:
        gls = temp('{0}/gls/{{sample_set}}/allSites/{{chrom}}/{{chrom}}_{{sample_set}}_allSites_maf{{maf}}.beagle.gz'.format(ANGSD_DIR)),
        mafs = temp('{0}/gls/{{sample_set}}/allSites/{{chrom}}/{{chrom}}_{{sample_set}}_allSites_maf{{maf}}.mafs.gz'.format(ANGSD_DIR)),
    log: 'logs/angsd_gl_allSites/{chrom}_{sample_set}_allSites_maf{maf}_angsd_gl.log'
    conda: '../envs/angsd.yaml'
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
    input:
        rules.angsd_saf_likelihood_allSites.output.pos
    output:
        sites = '{0}/angsd_sites/{{chrom}}/{{sample_set}}_{{chrom}}_Trepens_allSites.sites'.format(PROGRAM_RESOURCE_DIR),
        binary = '{0}/angsd_sites/{{chrom}}/{{sample_set}}_{{chrom}}_Trepens_allSites.sites.bin'.format(PROGRAM_RESOURCE_DIR),
        idx = '{0}/angsd_sites/{{chrom}}/{{sample_set}}_{{chrom}}_Trepens_allSites.sites.idx'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/extract_angsd_allSites/{sample_set}_{chrom}_extract_allSites.log'
    conda: '../envs/angsd.yaml'
    wildcard_constraints:
        site = 'allSites'
    shell:
        """
        ( zcat {input} | tail -n +2 | cut -f1,2 > {output.sites} &&\
            sleep 2m &&\
            angsd sites index {output.sites} ) 2> {log}
        """

rule split_angsd_sites_byChrom:
    input:
        rules.convert_sites_for_angsd.output
    output:
        sites = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}.sites'.format(PROGRAM_RESOURCE_DIR),
        binary = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}.sites.bin'.format(PROGRAM_RESOURCE_DIR),
        idx = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}.sites.idx'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/split_angsd_sites_byChrom/{chrom}_{site}_split_angsd_sites.log'
    conda: '../envs/angsd.yaml'
    shell:
        """
        ( grep {wildcards.chrom} {input} > {output.sites} &&\
            sleep 2m &&\
            angsd sites index {output.sites} ) 2> {log}
        """

rule angsd_estimate_sfs:
    input:
        unpack(angsd_sfs_input) 
    output:
        temp('{0}/{{sample_set}}/sfs/{{site}}/{{chrom}}/{{chrom}}_{{sample_set}}_{{site}}.sfs'.format(ANGSD_DIR))
    log: 'logs/angsd_estimate_sfs/{chrom}_{sample_set}_{site}_sfs.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 60000,
        time = '03:00:00'
    shell:
        """
        realSFS {input.saf_idx} -P {threads} -sites {input.sites} -fold 1 > {output} 2> {log}
        """

rule angsd_estimate_thetas:
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
            -P {threads} \
            -sfs {input.sfs} \
            -sites {input.sites} \
            -outname {params.out} 2> {log}
        """

rule angsd_diversity_neutrality_stats:
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
    input:
        get_angsd_stats_toConcat
    output:
        '{0}/summary_stats/thetas/{{sample_set}}/{{site}}/{{sample_set}}_{{site}}_diversityNeutrality.thetas.idx.pestPG'.format(ANGSD_DIR)
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
    input:
        rules.concat_sfs.output
    output:
        '{0}/sfs/{{sample_set}}/{{site}}/{{sample_set}}_{{site}}_allChroms.sfs'.format(ANGSD_DIR)
    log: 'logs/sum_sfs/{sample_set}_{site}.log'
    run:
        import pandas as pd
        import sys
        with open(log, 'w') as logfile:
            sys.stderr = logfile
            sfs_allChroms = pd.read_table(input[0], sep = ' ', header = None)
            sfs_sum = sfs_allChroms.sum(axis=0) 
            sfs_sum.to_csv(output[0], sep = '\t', header = None)
        

rule files_for_angsd_site_subset:
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
    input:
        expand(rules.angsd_depth.output, chrom=CHROMOSOMES, sample_set=['highErrorRemoved','finalSamples_lowCovRemoved']),
        expand(rules.concat_angsd_stats.output, site=['allSites','0fold','4fold'], sample_set=['highErrorRemoved','finalSamples_lowCovRemoved']),
        expand(rules.sum_sfs.output, site=['allSites','0fold','4fold'], sample_set=['highErrorRemoved','finalSamples_lowCovRemoved']),
        expand(rules.concat_angsd_gl.output, sample_set=['highErrorRemoved','finalSamples_lowCovRemoved'], site=['allSites','0fold','4fold'], maf=['0.05']),
        expand(rules.concat_angsd_mafs.output, sample_set=['highErrorRemoved','finalSamples_lowCovRemoved'], site=['allSites','0fold','4fold'], maf=['0.05']),
        expand(rules.extract_sample_angsd.output, sample_set=['highErrorRemoved','finalSamples_lowCovRemoved'])
    output:
        '{0}/angsd.done'.format(ANGSD_DIR)
    shell:
        "touch {output}"

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
