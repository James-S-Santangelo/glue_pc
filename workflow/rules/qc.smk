rule fastqc_raw_reads:
    input:
        tmp = rules.create_tmp_dir.output,
        read1 = lambda wildcards: raw_read_dict[wildcards.sample]['R1'],
        read2 = lambda wildcards: raw_read_dict[wildcards.sample]['R2']
    output:
        html1 = '{0}/fastqc_raw_reads/{{sample}}_1_fastqc.html'.format(QC_DIR),
        html2 = '{0}/fastqc_raw_reads/{{sample}}_2_fastqc.html'.format(QC_DIR),
        zip1 = '{0}/fastqc_raw_reads/{{sample}}_1_fastqc.zip'.format(QC_DIR),
        zip2 = '{0}/fastqc_raw_reads/{{sample}}_2_fastqc.zip'.format(QC_DIR)
    conda: '../envs/fastqc.yaml'
    log: 'logs/fastqc_raw_reads/{sample}_fastqc_raw_reads.log'
    threads: 2
    resources: 
        mem_mb = 1000, 
        time = '01:00:00'
    shell:
        """
        ( fastqc --threads {{threads}} --outdir {0}/fastqc_raw_reads --noextract --quiet --dir {{input.tmp}} {{input.read1}} {{input.read2}} &&

        mv {0}/fastqc_raw_reads/{{wildcards.sample}}_*_1_fastqc.html {0}/fastqc_raw_reads/{{wildcards.sample}}_1_fastqc.html &&
        mv {0}/fastqc_raw_reads/{{wildcards.sample}}_*_1_fastqc.zip {0}/fastqc_raw_reads/{{wildcards.sample}}_1_fastqc.zip &&
        mv {0}/fastqc_raw_reads/{{wildcards.sample}}_*_2_fastqc.html {0}/fastqc_raw_reads/{{wildcards.sample}}_2_fastqc.html &&
        mv {0}/fastqc_raw_reads/{{wildcards.sample}}_*_2_fastqc.zip {0}/fastqc_raw_reads/{{wildcards.sample}}_2_fastqc.zip ) 2> {{log}}
        """.format(QC_DIR)

rule fastqc_trimmed_reads:
    input:
        tmp = rules.create_tmp_dir.output,
        read1 = rules.fastp_trim.output.r1_trim,
        read2 = rules.fastp_trim.output.r2_trim
    output:
        html1 = '{0}/fastqc_trimmed_reads/{{sample}}_trimmed_1_fastqc.html'.format(QC_DIR),
        html2 = '{0}/fastqc_trimmed_reads/{{sample}}_trimmed_2_fastqc.html'.format(QC_DIR),
        zip1 = '{0}/fastqc_trimmed_reads/{{sample}}_trimmed_1_fastqc.zip'.format(QC_DIR),
        zip2 = '{0}/fastqc_trimmed_reads/{{sample}}_trimmed_2_fastqc.zip'.format(QC_DIR)
    conda: '../envs/fastqc.yaml'
    log: 'logs/fastqc_trimmed_reads/{sample}_fastqc_trimmed_reads.log'
    threads: 2
    resources:
        mem_mb = 1000,
        time = '01:00:00'
    shell:
        """
        fastqc --threads {{threads}} --outdir {0}/fastqc_trimmed_reads --noextract --quiet --dir {1} {{input.read1}} {{input.read2}} 2> {{log}}
        """.format(QC_DIR, TMPDIR)

rule qualimap_bam_qc:
    input:
        get_bam
    output:
        temp(directory('{0}/qualimap/{{sample}}_qualimap_bamqc'.format(QC_DIR)))
    log: 'logs/qualimap/{sample}_bamqc.log'
    conda: '../envs/qualimap.yaml'
    threads: 8
    resources:
        mem_mb = lambda wildcards, input: 4 * int(input.size_mb),
        time = '04:00:00'
    shell:
        """
        unset DISPLAY;
        qualimap bamqc -bam {{input}} \
            --paint-chromosome-limits \
            --collect-overlap-pairs \
            -nt {{threads}} \
            -outdir {0}/qualimap/{{wildcards.sample}}_qualimap_bamqc \
            -outformat html \
            --java-mem-size={{resources.mem_mb}}M >> {{log}} 2>&1
        """.format(QC_DIR)

rule multiqc:
    input:
       fastqc_raw = expand('{0}/fastqc_raw_reads/{{sample}}_{{read}}_fastqc.zip'.format(QC_DIR), sample=SAMPLES, read=['1', '2']),
       fastqc_trim = expand('{0}/fastqc_trimmed_reads/{{sample}}_trimmed_{{read}}_fastqc.zip'.format(QC_DIR), sample=SAMPLES, read=['1', '2']),
       fastp = expand('{0}/fastp_trim_reports/{{sample}}_fastp.json'.format(QC_DIR), sample=SAMPLES),
       qualimap = expand('{0}/qualimap/{{sample}}_qualimap_bamqc'.format(QC_DIR), sample=SAMPLES),
       bamstats = expand('{0}/bamtools_stats/{{sample}}_bamtools.stats'.format(QC_DIR), sample=SAMPLES)
    output:
        '{0}/multiqc/multiqc_report.html'.format(QC_DIR)
    conda: '../envs/multiqc.yaml'
    log: 'logs/multiqc/multiqc.log'
    shell:
        """
        multiqc --verbose \
            --dirs \
            --outdir {0}/multiqc \
            --config ../config/multiqc_config.yaml \
            {{input}} 2> {{log}}
        """.format(QC_DIR)

rule bamutil_validate:
    input:
        get_bam
    output:
        '{0}/bamutil_validate/{{sample}}_validation.txt'.format(QC_DIR)
    log: 'logs/bamutil_validate/{sample}_validation.log'
    conda: '../envs/bamutil.yaml'
    resources:
        mem = lambda wildcards, input: int(input.size_mb),
        time = '04:00:00'
    shell:
        """
        bam validate --in {input} \
            --so_coord \
            --verbose 2> {output}
        """

rule bamtools_stats:
    input:
        rules.samtools_markdup.output.bam
    output:
        '{0}/bamtools_stats/{{sample}}_bamtools.stats'.format(QC_DIR)
    conda: '../envs/bamtools.yaml'
    log: 'logs/bamtools_stats/{sample}_bamtools_stats.log'
    resources:
        time = '01:00:00'
    shell:
        """
        bamtools stats -in {input} > {output} 2> {log}
        """
