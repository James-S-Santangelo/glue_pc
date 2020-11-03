rule fastp_trim:
    input:
        raw_read_qc_done = 'fastqc_raw_reads.done'
        r1 = lambda wc: str(glob.glob('{0}/{1}/{1}_*_1.fq.gz'.format(RAW_READ_DIR, wc.sample))[0]),
        r2 = lambda wc: str(glob.glob('{0}/{1}/{1}_*_2.fq.gz'.format(RAW_READ_DIR, wc.sample))[0]),
    output:
        r1_trim = "{0}/{{sample}}/{{sample}}_trimmed_1.fq.gz".format(TRIMMED_READ_DIR),
        r2_trim = "{0}/{{sample}}/{{sample}}_trimmed_2.fq.gz".format(TRIMMED_READ_DIR),
        unp = "{0}/{{sample}}/{{sample}}_unpaired.fq.gz".format(TRIMMED_READ_DIR),
        html = "../results/fastp_trim_reports/{sample}_fastp.html"
    conda: "../envs/fastp.yaml"
    log: "logs/fastp_trim/{sample}_fastp.log"  
    shell:
        #"touch {output}"
        """
        fastp --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1_trim} \
            --out2 {output.r2_trim} \
            --unpaired1 {output.unp} \
            --unpaired2 {output.unp} \
            --html {output.html} \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --overrepresentation_analysis 2> {log}
        """ 
