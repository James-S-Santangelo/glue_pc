rule create_chrom_file:
    input: get_representative_bam()
    output:
        "../resources/chromosome_file.txt"
    log: "logs/create_chrom_file/create_chrom_file.log"
    conda: "../envs/variant_calling.yaml"
    shell:
        """
        samtools idxstats {{input}} | cut -f1 | grep {0} > {{output}} 
        """.format(CHROM_PREFIX)

rule create_regions_equal_coverage:
    input: 
        chrom_file = '../resources/chromosome_file.txt',
        bam = get_representative_bam()
    output:
        "../resources/{chrom}_forFreebayes.regions"
    log: "logs/create_regions_equal_cov/{chrom}_create_regions_equal_cov.log"
    conda: "../envs/variant_calling.yaml"
    shell:
        """
        ( samtools view -b -s 0.20 {{input.bam}} {{wildcards.chrom}} |\
            bamtools coverage |\
            coverage_to_regions.py {0}.fai {1} > {{output}} ) 2> {{log}} 
        """.format(REFERENCE_GENOME, NUM_REGIONS_PER_CHROM)

rule concat_regions_forFreebayes:
    input:
        expand('../resources/{chrom}_forFreebayes.regions', chrom=get_chrom_list())
    output:
        "../resources/wholeGenome_forFreebayes.regions"
    log: "logs/concat_regions_forFreebayes/concat_regions_forFreebayes.log"
    conda: "../envs/variant_calling.yaml"
    shell:
        """
        cat {input} >> {output} 2> {log}
        """
