rule bwa_map_unpaired:
    input:
        trimmed_reads_done_file = 'fastqc_trimmed_reads.done',
        unp = rules.fastp_trim.output.unp
    output:
        "../results/bams/unpaired/{sample}_unpaired_sorted.bam"
    params:
        r"-R '@RG\tID:${sample}\tCN:NOVOGENE\tPL:ILLUMINA\tPM:NOVASEQ.S4\tSM:${sample}'"
    conda: "../envs/bwa_mapping.yaml"
    log: "logs/bwa_map_unpaired/{sample}_bwa_map.unpaired.log"
    shell:
        #"touch {output}"
        """
        ( bwa mem {0} {{input.unp}} {{params}} |\
            samtools view -hb -o {{output}} - ) 2> {{log}}
        """.format(REFERENCE_GENOME, TMPDIR)

rule bwa_map_paired:
    input:
        trimmed_reads_done_file = 'fastqc_trimmed_reads.done',
        r1 = rules.fastp_trim.output.r1_trim,
        r2 = rules.fastp_trim.output.r2_trim
    output:
        "../results/bams/paired/{sample}_paired_sorted.bam"
    params:
        r"-R '@RG\tID:${sample}\tCN:NOVOGENE\tPL:ILLUMINA\tPM:NOVASEQ.S4\tSM:${sample}'"
    conda: "../envs/bwa_mapping.yaml"
    log: "logs/bwa_map_paired/{sample}_bwa_map.paired.log"
    shell:
        #"touch {output}"
        """
        ( bwa mem {0} {{input.r1}} {{input.r2}} {{params}} |\
            samtools view -hb -o {{output}} - ) 2> {{log}}
        """.format(REFERENCE_GENOME, TMPDIR)

rule merge_bams:
    input:
        unp = rules.bwa_map_unpaired.output,
        pair = rules.bwa_map_paired.output
    output:
        "../results/bams/merged/{sample}_merged_sorted.bam"
    conda: "../envs/bwa_mapping.yaml"
    log: "logs/merge_bams/{sample}_merge_bams.log"
    shell:
        #"touch {output}"
        """
        ( samtools cat {{input.pair}} {{input.unp}} |\
            samtools collate -o {{output}} - {0}/{{wildcards.sample}}_merged ) 2> {{log}}
        """.format(TMPDIR)

rule samtools_markdup:
    input:
        rules.merge_bams.output
    output:
        bam = "../results/bams/final/{sample}_merged_sorted_dupsMarked.bam",
        stats = "../results/duplication_stats/{sample}_dupStats.txt"
    conda: "../envs/bwa_mapping.yaml"
    log: "logs/samtools_markdup/{sample}_samtools_markdup.log"
    shell:
        """
        ( samtools fixmate -m {{input}} - |\
            samtools sort -T {0}/{{wildcards.sample}} -o - |\
            samtools markdup -T {0} -f {{output.stats}} - {{output.bam}} ) 2> {{log}}
        """.format(TMPDIR)

checkpoint index_bam:
    input:
        rules.samtools_markdup.output.bam
    output:
        "../results/bams/final/{sample}_merged_sorted_dupsMarked.bam.bai",
    conda: "../envs/bwa_mapping.yaml"
    log: "logs/index_bam/{sample}_index_bam.log"
    shell:
        """
        samtools index {input} 2> {log}
        """

