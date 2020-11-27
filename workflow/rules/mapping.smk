rule bwa_map_unpaired:
    input:
        fastqc_trimmed_done = rules.fastqc_trimmed_reads.output.html1,
        unp = rules.fastp_trim.output.unp
    output:
        temp('{0}/unpaired/{{sample}}_unpaired_sorted.bam'.format(BAM_DIR))
    params:
        r"-R '@RG\tID:${sample}\tCN:NOVOGENE\tPL:ILLUMINA\tPM:NOVASEQ.S4\tSM:${sample}'"
    conda: '../envs/bwa_mapping.yaml'
    log: 'logs/bwa_map_unpaired/{sample}_bwa_map.unpaired.log'
    threads: 8
    resources:
        mem_mb = lambda wildcards, input: int(input.size_mb),
        time = '03:00:00'
    shell:
        """
        ( bwa mem -t {{threads}} {0} {{input.unp}} {{params}} |\
            samtools view -hb -o {{output}} - ) 2> {{log}}
        """.format(REFERENCE_GENOME, TMPDIR)

rule bwa_map_paired:
    input:
        fastqc_trimmed_done = rules.fastqc_trimmed_reads.output.html1,
        r1 = rules.fastp_trim.output.r1_trim,
        r2 = rules.fastp_trim.output.r2_trim
    output:
        temp('{0}/paired/{{sample}}_paired_sorted.bam'.format(BAM_DIR))
    params:
        r"-R '@RG\tID:${sample}\tCN:NOVOGENE\tPL:ILLUMINA\tPM:NOVASEQ.S4\tSM:${sample}'"
    conda: '../envs/bwa_mapping.yaml'
    log: 'logs/bwa_map_paired/{sample}_bwa_map.paired.log'
    threads: 8
    resources:
        mem_mb = lambda wildcards, input: int(input.size_mb),
        time = '06:00:00'
    shell:
        """
        ( bwa mem -t {{threads}} {0} {{input.r1}} {{input.r2}} {{params}} |\
            samtools view -hb -o {{output}} - ) 2> {{log}}
        """.format(REFERENCE_GENOME, TMPDIR)

rule merge_bams:
    input:
        unp = rules.bwa_map_unpaired.output,
        pair = rules.bwa_map_paired.output
    output:
        temp('{0}/merged/{{sample}}_merged_sorted.bam'.format(BAM_DIR))
    conda: '../envs/bwa_mapping.yaml'
    log: 'logs/merge_bams/{sample}_merge_bams.log'
    threads: 8
    resources:
        time = '03:00:00'
    shell:
        """
        ( samtools cat --threads {{threads}} {{input.pair}} {{input.unp}} |\
            samtools collate --threads {{threads}} -o {{output}} - {0}/{{wildcards.sample}}_merged ) 2> {{log}}
        """.format(TMPDIR)

rule samtools_markdup:
    input:
        rules.merge_bams.output
    output:
        bam = '{0}/final/{{sample}}_merged_sorted_dupsMarked.bam'.format(BAM_DIR),
        stats = '{0}/duplication_stats/{{sample}}_dupStats.txt'.format(QC_DIR)
    conda: '../envs/bwa_mapping.yaml'
    log: 'logs/samtools_markdup/{sample}_samtools_markdup.log'
    threads: 8
    resources:
        mem_mb = lambda wildcards, input: int(input.size_mb),
        time = '12:00:00'
    shell:
        """
        ( samtools fixmate --threads {{threads}} -m {{input}} - |\
            samtools sort --threads {{threads}} -T {0}/{{wildcards.sample}} -o - |\
            samtools markdup --threads {{threads}} -T {0} -f {{output.stats}} - {{output.bam}} ) 2> {{log}}
        """.format(TMPDIR)

checkpoint index_bam:
    input:
        rules.samtools_markdup.output.bam
    output:
        '{0}/final/{{sample}}_merged_sorted_dupsMarked.bam.bai'.format(BAM_DIR),
    conda: '../envs/bwa_mapping.yaml'
    log: 'logs/index_bam/{sample}_index_bam.log'
    threads: 8
    resources:
        time = '01:00:00'
    shell:
        """
        samtools index -@ {threads} {input} 2> {log}
        """

