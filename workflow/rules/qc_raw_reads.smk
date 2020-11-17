rule fastqc_raw_reads:
    input:
        tmp = rules.create_tmp_dir.output,
        reads = get_raw_reads()
    output:
        fastqc_target_files(),
        touch('../results/flag_files/fastqc_raw_reads.done')
    conda: "../envs/fastqc.yaml"
    log: "logs/fastqc_raw_reads/fastqc_raw_reads.log"
    resources: 
        cpus = lambda wildcards, input: len(input.reads),
        mem_mb = lambda wildcards, input: len(input.reads) * 300        
    shell:
        """
        fastqc --threads {{resources.cpus}} --outdir {0} --noextract --quiet --dir {{input.tmp}} {{input.reads}} 2> {{log}}
        """.format(RAW_READ_FASTQC_DIR)

