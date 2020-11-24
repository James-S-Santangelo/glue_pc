rule bwa_map_unpaired:
    input:
        trimmed_reads_done_file = '../results/flag_files/fastqc_trimmed_reads.done',
        unp = rules.fastp_trim.output.unp
    output:
        temp('{0}/unpaired/{{sample}}_unpaired_sorted.bam'.format(BAM_DIR))
    params:
        r"-R '@RG\tID:${sample}\tCN:NOVOGENE\tPL:ILLUMINA\tPM:NOVASEQ.S4\tSM:${sample}'"
    conda: '../envs/bwa_mapping.yaml'
    log: 'logs/bwa_map_unpaired/{sample}_bwa_map.unpaired.log'
    resources:
        cpus = 4,
        mem_mb = lambda wildcards, input: int(input.size_mb)
    shell:
        """
        ( bwa mem -t {{resources.cpus}} {0} {{input.unp}} {{params}} |\
            samtools view -hb -o {{output}} - ) 2> {{log}}
        """.format(REFERENCE_GENOME, TMPDIR)

rule bwa_map_paired:
    input:
        trimmed_reads_done_file = '{0}/fastqc_trimmed_reads.done'.format(FLAG_FILES_DIR),
        r1 = rules.fastp_trim.output.r1_trim,
        r2 = rules.fastp_trim.output.r2_trim
    output:
        temp('{0}/paired/{{sample}}_paired_sorted.bam'.format(BAM_DIR))
    params:
        r"-R '@RG\tID:${sample}\tCN:NOVOGENE\tPL:ILLUMINA\tPM:NOVASEQ.S4\tSM:${sample}'"
    conda: '../envs/bwa_mapping.yaml'
    log: 'logs/bwa_map_paired/{sample}_bwa_map.paired.log'
    resources:
        cpus = 4,
        mem_mb = lambda wildcards, input: int(input.size_mb)
    shell:
        """
        ( bwa mem -t {{resources.cpus}} {0} {{input.r1}} {{input.r2}} {{params}} |\
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
    resources:
        cpus = 4
    shell:
        """
        ( samtools cat --threads {{resources.cpus}} {{input.pair}} {{input.unp}} |\
            samtools collate --threads {{resources.cpus}} -o {{output}} - {0}/{{wildcards.sample}}_merged ) 2> {{log}}
        """.format(TMPDIR)

rule samtools_markdup:
    input:
        rules.merge_bams.output
    output:
        bam = '{0}/final/{{sample}}_merged_sorted_dupsMarked.bam'.format(BAM_DIR),
        stats = '{0}/duplication_stats/{{sample}}_dupStats.txt'.format(QC_DIR)
    conda: '../envs/bwa_mapping.yaml'
    log: 'logs/samtools_markdup/{sample}_samtools_markdup.log'
    resources:
        cpus = 4,
        mem_mb = lambda wildcards, input: int(input.size_mb)
    shell:
        """
        ( samtools fixmate --threads {{resources.cpus}} -m {{input}} - |\
            samtools sort --threads {{resources.cpus}} -T {0}/{{wildcards.sample}} -o - |\
            samtools markdup --threads {{resources.cpus}} -T {0} -f {{output.stats}} - {{output.bam}} ) 2> {{log}}
        """.format(TMPDIR)

checkpoint index_bam:
    input:
        rules.samtools_markdup.output.bam
    output:
        '{0}/final/{{sample}}_merged_sorted_dupsMarked.bam.bai'.format(BAM_DIR),
    conda: '../envs/bwa_mapping.yaml'
    log: 'logs/index_bam/{sample}_index_bam.log'
    resources:
        cpus = 4
    shell:
        """
        samtools index -@ {resources.cpus} {input} 2> {log}
        """

