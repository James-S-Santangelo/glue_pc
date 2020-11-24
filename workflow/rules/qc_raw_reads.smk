rule fastqc_raw_reads:
    input:
        tmp = rules.create_tmp_dir.output,
        reads = get_raw_reads()
    output:
        fastqc_target_files(),
        touch('{0}/fastqc_raw_reads.done'.format(FLAG_FILES_DIR))
    conda: '../envs/fastqc.yaml'
    log: 'logs/fastqc_raw_reads/fastqc_raw_reads.log'
    resources: 
        cpus = lambda wildcards, input: len(input.reads),
        mem_mb = lambda wildcards, input: len(input.reads) * 300,
        time = '02:00:00'
    shell:
        """
        fastqc --threads {{resources.cpus}} --outdir {0}/fastqc_raw_reads --noextract --quiet --dir {{input.tmp}} {{input.reads}} 2> {{log}}
        """.format(QC_DIR)

