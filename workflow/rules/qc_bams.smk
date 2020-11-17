rule qualimap_bam_qc:
    input:
        get_bam
    output:
        temp(directory('../results/qualimap/{sample}_qualimap_bamqc'))
    log: 'logs/qualimap/{sample}_bamqc.log'
    conda: '../envs/qualimap.yaml'
    resources:
        cpus = 4,
        mem_mb = lambda wildcards, input: 4 * int(input.size_mb)
    shell:
        """
        unset DISPLAY;
        qualimap bamqc -bam {input} \
            --paint-chromosome-limits \
            --collect-overlap-pairs \
            -nt {resources.cpus} \
            -outdir ../results/qualimap/{wildcards.sample}_qualimap_bamqc \
            -outformat html \
            --java-mem-size={resources.mem_mb}M >> {log} 2>&1
        """

rule write_qualimap_multiqc_file:
    input:
        expand(rules.qualimap_bam_qc.output, sample=SAMPLES)
    output:
        '../results/program_resources/bamqcPaths_forQualimap_multiQC.txt'
    run:
        with open(output[0], 'w') as f:
            for sample in SAMPLES:
                bam = [bam for bam in input if sample in bam][0]
                f.write('{0}\t{1}\n'.format(sample, str(bam)))

rule qualimap_multiqc:
    input:
        rules.write_qualimap_multiqc_file.output
    output:
        directory('../results/qualimap/qualimap_multiqc')
    log: 'logs/qualimap/multiqc.log'
    conda: '../envs/qualimap.yaml'
    shell:
        """
        unset DISPLAY;
        qualimap multi-bamqc --data {input} \
            --paint-chromosome-limits \
            -outdir ../results/qualimap/qualimap_multiqc \
            -outformat html >> {log} 2>&1
        """

rule bamutil_validate:
    input:
        get_bam
    output:
        "../results/bamutil_validate/{sample}_validation.txt"
    log: "logs/bamutil_validate/{sample}_validation.log"
    conda: "../envs/bamutil.yaml"
    shell:
        """
        bam validate --in {input} \
            --so_coord \
            --verbose 2> {output}
        """
