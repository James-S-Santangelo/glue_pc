rule fastqc_raw_reads:
    input:
        get_raw_reads()
    output:
        fastqc_target_files()
    conda: "envs/fastqc.yaml"
    log: "logs/fastqc_raw_reads/fastqc_raw_reads.log"
    threads: 12
    shell:
        #"touch {output}"
        """
        fastqc --threads {{threads}} --outdir {0} --noextract --quiet --dir {1} {{input}} 2> {{log}}
        """.format(RAW_READ_FASTQC_DIR, TMPDIR)

