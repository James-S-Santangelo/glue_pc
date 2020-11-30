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
