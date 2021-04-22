# Downsample 20 Toronto samples for inclusion in GLUE
# This is kind of hacky and not very reproducible since BAMs will not be distributed.
# Need to incorporate raw reads for these samples instead

rule downsample_toronto_bam:
    """
    Sample 25% of reads from Toronto sample using samtools. 
    """
    input:
        get_toronto_bam
    output:
        '{0}/toronto_bams/{{tor_sample}}_merged_sorted_dupsMarked_downsampled.bam'.format(BAM_DIR)
    conda: '../envs/mapping.yaml'
    log: 'logs/downsample_toronto_bam/{tor_sample}_downsample.log'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '03:00:00'
    shell:
        """
        samtools view -hb -s 0.25 {input} > {output} 2> {log}
        """

rule index_toronto_bam:
    """
    Index downsampled Toronto BAM using samtools. Required for variant calling using ANGSD.
    """
    input:
        rules.downsample_toronto_bam.output
    output:
        '{0}/toronto_bams/{{tor_sample}}_merged_sorted_dupsMarked_downsampled.bam.bai'.format(BAM_DIR)
    conda: '../envs/mapping.yaml'
    log: 'logs/index_toronto_bam/{tor_sample}_index.log'
    resources:
        mem_mb = 1000,
        time = '01:00:00'
    shell:
        """
        samtools index {input} 2> {log}
        """

rule downsample_toronto_done:
    """
    Generate empty flag file signalling successful completion of Toronto BAM downsampling.
    """
    input:
        expand(rules.index_toronto_bam.output, tor_sample=TOR_SAMPLES)
    output:
        '{0}/downsample_toronto.done'.format(BAM_DIR)
    shell:
        """
        touch {output}
        """
