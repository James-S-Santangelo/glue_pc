rule qualimap_bam_qc:
    input:
        get_bam
    output:
        temp(directory('../results/qualimap_bamqc/{sample}_qualimap_bamqc'))
    log: 'logs/qualimap_bamqc/{sample}_bamqc.log'
    conda: '../envs/qualimap.yaml'
    threads: 2
    shell:
        """
        unset DISPLAY;
        qualimap bamqc -bam {input} \
            --paint-chromosome-limits \
            --collect-overlap-pairs \
            -nt {threads} \
            -outdir ../results/qualimap_bamqc/{wildcards.sample}_qualimap_bamqc \
            -outformat html >> {log} 2>&1
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
