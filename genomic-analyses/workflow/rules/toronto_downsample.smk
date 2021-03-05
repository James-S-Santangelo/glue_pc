rule downsample_toronto_bam:
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
    input:
        expand(rules.index_toronto_bam.output, tor_sample=TOR_SAMPLES)
    output:
        '{0}/downsample_toronto.done'.format(BAM_DIR)
    shell:
        """
        touch {output}
        """
