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

rule write_qualimap_multiqc_file:
    input:
        expand(rules.qualimap_bam_qc.output, sample=SAMPLES)
    output:
        '{0}/bamqcPaths_forQualimap_multiQC.txt'.format(PROGRAM_RESOURCE_DIR)
    run:
        with open(output[0], 'w') as f:
            for sample in SAMPLES:
                bam = [bam for bam in input if sample in bam][0]
                f.write('{0}\t{1}\n'.format(sample, str(bam)))

rule qualimap_multiqc:
    input:
        rules.write_qualimap_multiqc_file.output
    output:
        directory('{0}/qualimap/qualimap_multiqc'.format(QC_DIR))
    log: 'logs/qualimap/multiqc.log'
    conda: '../envs/qualimap.yaml'
    resources:
        time = '04:00:00'
    shell:
        """
        unset DISPLAY;
        qualimap multi-bamqc --data {{input}}\
            --paint-chromosome-limits \
            -outdir {0}/qualimap/qualimap_multiqc \
            -outformat html >> {{log}} 2>&1
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
