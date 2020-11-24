rule fastqc_trimmed_reads:
    input:
        reads = expand(['{0}/{{sample}}/{{sample}}_trimmed_1.fq.gz'.format(TRIMMED_READ_DIR), 
            '{0}/{{sample}}/{{sample}}_trimmed_2.fq.gz'.format(TRIMMED_READ_DIR)], sample=SAMPLES)
    output:
        expand(['{0}/fastqc_trimmed_reads/{{sample}}_trimmed_1{{ext}}'.format(QC_DIR), 
            '{0}/fastqc_trimmed_reads/{{sample}}_trimmed_2{{ext}}'.format(QC_DIR)],
            sample=SAMPLES, ext=['_fastqc.html', '_fastqc.zip']),
        touch('{0}/fastqc_trimmed_reads.done'.format(FLAG_FILES_DIR))
    conda: '../envs/fastqc.yaml'
    log: 'logs/fastqc_trimmed_reads/fastqc_trimmed_reads.log'
    resources:
        cpus = lambda wildcards, input: len(input.reads),
        mem_mb = lambda wildcards, input: len(input.reads) * 300,
        time = '04:00:00'
    shell:
        """
        fastqc --threads {{resources.cpus}} --outdir {0}/fastqc_trimmed_reads --noextract --quiet --dir {1} {{input.reads}} 2> {{log}}
        """.format(QC_DIR, TMPDIR)


