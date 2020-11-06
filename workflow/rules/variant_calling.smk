rule create_regions_equal_coverage:
    input: 
        get_representative_bam
    output:
        "../resources/{chrom}_forFreebayes.regions"
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
        expand("../resources/{chrom}_forFreebayes.regions", chrom=CHROMOSOMES)
    output:
        "../resources/wholeGenome_forFreebayes.regions"
    log: "logs/concat_regions_forFreebayes/concat_regions_forFreebayes.log"
    shell:
        """
        cat {input} >> {output} 2> {log}
        """

rule create_bam_list:
    input:
        expand('../results/bams/final/{sample}_merged_sorted_dupsMarked.bam', sample=SAMPLES)
    output:
        '../resources/freebayes_bams.list'
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
        '../results/vcf/wholeGenome_allSamples.vcf'
    log: 'logs/freebayes/freebayes.log'
    conda: '../envs/variant_calling.yaml'
    threads: 5
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
        '../results/vcf/wholeGenome_allSamples.vcf.gz'
    log: 'logs/bgzip/bgzip.log'
    conda: '../envs/variant_calling.yaml'
    shell:
        """
        bgzip {input}
        """

rule tabix_vcf:
    input:
        rules.bgzip_vcf.output
    output:
        '../results/vcf/wholeGenome_allSamples.vcf.gz.tbi'
    log: 'logs/tabix/tabix.log'
    conda: '../envs/variant_calling.yaml'
    shell:
        """
        tabix {input}
        """
