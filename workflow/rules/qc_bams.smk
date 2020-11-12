rule qualimap_bam_qc:
    input:
        get_all_bams
    output:
        '{sample}_test.txt'
    shell:
        "echo {input} > {output}"   
