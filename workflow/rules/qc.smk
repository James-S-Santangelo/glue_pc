rule multiqc:
    input:
       #fastqc_raw = expand('{0}/fastqc_raw_reads/{{sample}}_{{read}}_fastqc.zip'.format(QC_DIR), sample=SAMPLES, read=['1', '2']),
       #fastqc_trim = expand('{0}/fastqc_trimmed_reads/{{sample}}_trimmed_{{read}}_fastqc.zip'.format(QC_DIR), sample=SAMPLES, read=['1', '2']),
       #fastp = expand('{0}/fastp_trim_reports/{{sample}}_fastp.json'.format(QC_DIR), sample=SAMPLES),
       #qualimap = expand('{0}/qualimap/{{sample}}_qualimap_bamqc'.format(QC_DIR), sample=SAMPLES),
       #flagstat = expand('{0}/samtools_flagstat/{{sample}}_flagstat.tsv'.format(QC_DIR), sample=SAMPLES)
       '{0}'.format(QC_DIR)
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
