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
