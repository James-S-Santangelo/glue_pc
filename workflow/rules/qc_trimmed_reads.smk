rule fastqc_trimmed_reads:
    input:
        reads = expand(['{0}/{{sample}}/{{sample}}_trimmed_1.fq.gz'.format(TRIMMED_READ_DIR), 
            '{0}/{{sample}}/{{sample}}_trimmed_2.fq.gz'.format(TRIMMED_READ_DIR)], sample=SAMPLES)
    output:
        protected(expand(['{0}/{{sample}}_trimmed_1{{ext}}'.format(TRIMMED_READ_FASTQC_DIR), 
            '{0}/{{sample}}_trimmed_2{{ext}}'.format(TRIMMED_READ_FASTQC_DIR)],
            sample=SAMPLES, ext=['_fastqc.html', '_fastqc.zip'])),
        temp(touch('fastqc_trimmed_reads.done'))
    conda: "../envs/fastqc.yaml"
    log: "logs/fastqc_trimmed_reads/fastqc_trimmed_reads.log"
    threads: 12
    shell:
        #"touch {output}"
        """
        fastqc --threads {{threads}} --outdir {0} --noextract --quiet --dir {1} {{input.reads}} 2> {{log}}
        """.format(TRIMMED_READ_FASTQC_DIR, TMPDIR)


