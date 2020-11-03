rule fastqc_raw_reads:
    input:
        tmp = rules.create_tmp_dir.output,
        reads = get_raw_reads()
    output:
        protected(fastqc_target_files()),
        touch('fastqc_raw_reads.done')
    conda: "../envs/fastqc.yaml"
    log: "logs/fastqc_raw_reads/fastqc_raw_reads.log"
    threads: 12
    shell:
        #"touch {output}"
        """
        fastqc --threads {{threads}} --outdir {0} --noextract --quiet --dir {{input.tmp}} {{input.reads}} 2> {{log}}
        """.format(RAW_READ_FASTQC_DIR)

