rule fastp_trim:
    input:
        read1 = lambda wildcards: raw_read_dict[wildcards.sample]['R1'],
        read2 = lambda wildcards: raw_read_dict[wildcards.sample]['R2']
    output:
        r1_trim = '{0}/{{sample}}/{{sample}}_trimmed_1.fq.gz'.format(TRIMMED_READ_DIR),
        r2_trim = '{0}/{{sample}}/{{sample}}_trimmed_2.fq.gz'.format(TRIMMED_READ_DIR),
        unp = '{0}/{{sample}}/{{sample}}_trimmed_unpaired.fq.gz'.format(TRIMMED_READ_DIR),
        html = '{0}/fastp_trim_reports/{{sample}}_fastp.html'.format(QC_DIR)
    conda: '../envs/fastp.yaml'
    log: 'logs/fastp_trim/{sample}_fastp.log'
    threads: 4
    resources:
        mem_mb = lambda wildcards, input, attempt: attempt * int(input.size_mb / 2),
        time = '03:00:00'
    shell:
        """
        fastp --in1 {input.read1} \
            --in2 {input.read2} \
            --out1 {output.r1_trim} \
            --out2 {output.r2_trim} \
            --unpaired1 {output.unp} \
            --unpaired2 {output.unp} \
            --html {output.html} \
            --thread {threads} \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --overrepresentation_analysis 2> {log}
        """ 
