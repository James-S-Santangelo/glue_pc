rule create_chrom_file:
    input: glob.glob('../results/bams/final/{0}_*.bam'.format(REPRESENTATIVE_SAMPLE))[0]
    output:
        "../resources/chromosome_file.txt"
    log: "logs/create_chrom_file/create_chrom_file.log"
    conda: "../envs/variant_calling.yaml"
    shell:
        """
        samtools idxstats {{input}} | cut -f1 | grep {0} > {{output}} 
        """.format(CHROM_PREFIX)
