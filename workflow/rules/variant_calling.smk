rule create_regions_equal_coverage:
    input: 
        get_representative_bam
    output:
        "../results/program_resources/{chrom}_forFreebayes.regions"
    log: "logs/create_regions_equal_cov/{chrom}_create_regions_equal_cov.log"
    conda: "../envs/variant_calling.yaml"
    shell:
        """
        ( samtools view -b -s 0.20 {{input}} {{wildcards.chrom}} |\
            bamtools coverage |\
            coverage_to_regions.py {0}.fai {1} > {{output}} ) 2> {{log}} 
        """.format(REFERENCE_GENOME, NUM_REGIONS_PER_CHROM)

rule concat_regions_forFreebayes:
    input:
        expand("../results/program_resources/{chrom}_forFreebayes.regions", chrom=CHROMOSOMES)
    output:
        "../results/program_resources/wholeGenome_forFreebayes.regions"
    log: "logs/concat_regions_forFreebayes/concat_regions_forFreebayes.log"
    shell:
        """
        cat {input} >> {output} 2> {log}
        """

rule create_bam_list:
    input:
        expand('../results/bams/final/{sample}_merged_sorted_dupsMarked.bam', sample=SAMPLES)
    output:
        '../results/program_resources/freebayes_bams.list'
    log: 'logs/create_bam_list/create_bam_list.log'
    run:
        with open(output[0], 'w') as f:
            for bam in input:
                f.write('{0}\n'.format(bam))

rule freebayes_call_variants:
    input:
        bams = rules.create_bam_list.output,
        #regions = rules.concat_regions_forFreebayes.output
        regions = '../resources/test.regions'
    output:
        '../results/vcf/wholeGenome_allSamples_allSites.vcf'
    log: 'logs/freebayes/freebayes.log'
    conda: '../envs/variant_calling.yaml'
    threads: 10 
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
        '../results/vcf/wholeGenome_allSamples_allSites.vcf.gz'
    log: 'logs/bgzip/bgzip.log'
    conda: '../envs/variant_calling.yaml'
    shell:
        """
        bgzip {input}
        """

rule bcftools_sort:
    input:
        rules.bgzip_vcf.output
    output:
        '../results/vcf/wholeGenome_allSamples_allSites_sorted.vcf.gz'
    log: 'logs/bcftools_sort/bcftools_sort.log'
    conda: '../envs/variant_calling.yaml'
    shell:
        """
        bcftools sort -O z -o {{output}} -T {0} {{input}} 2> {{log}}
        """.format(TMPDIR)

rule bcftools_split_variants:
    input:
        vcf = rules.bcftools_sort.output
    output:
        '../results/vcf/wholeGenome_allSamples_{site_type}_sorted.vcf.gz'
    log: 'logs/bcftools_split_variants/bcftools_split_variants_{site_type}.log'
    #conda: '../envs/variant_calling.yaml'
    wildcard_constraints:
        site_type='snps|indels|invariant|mnps|other'
    run:
        if wildcards.site_type == 'invariant':
            shell("bcftools view -O z --include 'N_ALT = 0' {input} > {output} 2> {log}")
        elif wildcards.site_type == 'snps':
            shell("( bcftools view -O v --types {{wildcards.site_type}} {{input}} |\
                vcfallelicprimitives --keep-info --keep-geno |\
                bcftools view --types {{wildcards.site_type}} --min-alleles 2 --max-alleles 2 |\
                bcftools sort -O z -T {0} -o {{output}} ) 2> {{log}}".format(TMPDIR))
        else:
            shell("bcftools view -O z --types {wildcards.site_type} {input} > {output} 2> {log}")

rule tabix_vcf:
    input: get_tabix_files
    output:
        '../results/vcf/wholeGenome_allSamples_{site_type}_sorted.vcf.gz.tbi'
    log: 'logs/tabix/tabix_{site_type}.log'
    conda: '../envs/variant_calling.yaml'
    shell:
        """
        tabix {input}
        """
