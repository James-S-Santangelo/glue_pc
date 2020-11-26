rule create_regions_equal_coverage:
    input: 
        get_representative_bam
    output:
        '{0}/{{chrom}}_forFreebayes.regions'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/create_regions_equal_cov/{chrom}_create_regions_equal_cov.log'
    conda: '../envs/variant_calling.yaml'
    threads: 8
    resources:
        time = '06:00:00'
    shell:
        """
        ( samtools view --threads {{threads}} -b -s 0.20 {{input}} {{wildcards.chrom}} |\
            bamtools coverage |\
            coverage_to_regions.py {0}.fai {1} > {{output}} ) 2> {{log}} 
        """.format(REFERENCE_GENOME, NUM_REGIONS_PER_CHROM)

rule create_bam_list:
    input:
        expand('{0}/final/{{sample}}_merged_sorted_dupsMarked.bam'.format(BAM_DIR), sample=SAMPLES)
    output:
        '{0}/freebayes_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/create_bam_list/create_bam_list.log'
    run:
        with open(output[0], 'w') as f:
            for bam in input:
                f.write('{0}\n'.format(bam))

rule freebayes_call_variants:
    input:
        bams = rules.create_bam_list.output,
        #regions = rules.concat_regions_forFreebayes.output
        regions = '../resources/{chrom}_forFreebayes.regions'
    output:
        temp('{0}/vcf/{{chrom}}/{{chrom}}_allSamples.vcf'.format(VARIANT_DIR))
    log: 'logs/freebayes/{chrom}_freebayes.log'
    conda: '../envs/variant_calling.yaml'
    threads: 5
    resources:
        time = '24:00:00'
    shell:
        """
        ( freebayes-parallel {{input.regions}} {{threads}} \
            --fasta-reference {0} \
            --bam-list {{input.bams}} \
            --use-best-n-alleles 4 \
            --report-monomorphic \
            --max-complex-gap 1 \
            --haplotype-length 1 \
            --genotype-qualities > {{output}} ) 2> {{log}}
        """.format(REFERENCE_GENOME)

rule bgzip_vcf:
    input:
        rules.freebayes_call_variants.output
    output:
        temp('{0}/vcf/{{chrom}}/{{chrom}}_allSamples.vcf.gz'.format(VARIANT_DIR))
    log: 'logs/bgzip/{chrom}_bgzip.log'
    conda: '../envs/variant_calling.yaml',
    threads: 8
    resources:
        time = '01:00:00'
    shell:
        """
        bgzip -@ {threads} {input}
        """

rule bcftools_sort:
    input:
        rules.bgzip_vcf.output
    output:
        '{0}/vcf/{{chrom}}/{{chrom}}_allSamples_sorted.vcf.gz'.format(VARIANT_DIR)
    log: 'logs/bcftools_sort/{chrom}_bcftools_sort.log'
    conda: '../envs/variant_calling.yaml'
    resources:
        time = '12:00:00'
    shell:
        """
        mkdir {0}/{{wildcards.chrom}};
        bcftools sort -O z -o {{output}} -T {0}/{{wildcards.chrom}} {{input}} 2> {{log}}
        """.format(TMPDIR)

rule bcftools_split_variants:
    input:
        vcf = rules.bcftools_sort.output
    output:
        '{0}/vcf/{{chrom}}/{{chrom}}_allSamples_{{site_type}}_sorted.vcf.gz'.format(VARIANT_DIR)
    log: 'logs/bcftools_split_variants/{chrom}_bcftools_split_variants_{site_type}.log'
    conda: '../envs/variant_calling.yaml'
    wildcard_constraints:
        site_type='snps|indels|invariant|mnps|other'
    threads: 8
    resources:
        time = '12:00:00'
    shell:
        """
        if [ {{wildcards.site_type}} = 'invariant' ]; then
            bcftools view --threads {{threads}} -O z --include 'N_ALT = 0' {{input}} > {{output}} 2> {{log}}
        elif [ {{wildcards.site_type}} = 'snps' ]; then
            ( bcftools view --threads {{threads}} -O v --types {{wildcards.site_type}} {{input}} |\
            vcfallelicprimitives --keep-info --keep-geno |\
            bcftools view --threads {{threads}} --types {{wildcards.site_type}} --min-alleles 2 --max-alleles 2 |\
            bcftools sort -O z -T {0}/{{wildcards.chrom}} -o {{output}} ) 2> {{log}}
        else
            bcftools view --threads {{threads}} -O z --types {{wildcards.site_type}} {{input}} > {{output}} 2> {{log}}
        fi
        """.format(TMPDIR)

rule tabix_vcf:
    input: get_tabix_files
    output:
        '{0}/vcf/{{chrom}}/{{chrom}}_allSamples_{{site_type}}_sorted.vcf.gz.tbi'.format(VARIANT_DIR)
    log: 'logs/tabix/{chrom}_tabix_{site_type}.log'
    conda: '../envs/variant_calling.yaml'
    shell:
        """
        tabix {input}
        """

rule vcf_to_zarr:
    input:
        rules.bcftools_split_variants.output
    output:
        directory('{0}/zarr_db/{{chrom}}/{{chrom}}_allSamples_{{site_type}}_sorted.zarr'.format(VARIANT_DIR))
    log: 'logs/vcf_to_zarr/{chrom}_vcf_to_zarr_{site_type}.log'
    conda: '../envs/vcf_to_zarr.yaml'
    wildcard_constraints:
        site_type='snps|invariant'
    threads: 1
    resources:
        time = '01:00:00'
    script:
        "../scripts/python/vcf_to_zarr.py"
