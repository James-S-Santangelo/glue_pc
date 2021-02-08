rule create_bam_list_forSpeciesID:
    input:
        expand(rules.samtools_markdup.output.bam, sample=SAMPLES)
    output:
        '{0}/allSamples_bams_forSpeciesID.list'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/create_bam_list/create_bam_list.log'
    run:
        import os
        with open(output[0], 'w') as f:
            for bam in input:
                sample = os.path.basename(bam).split('_merged')[0]
                f.write('{0}\n'.format(bam))

rule bcftools_chloroplast_gene_variants:
    input:
        rules.create_bam_list_forSpeciesID.output
    output:
        '{0}/{{gene}}/allSamples_{{gene}}.vcf.gz'.format(SPECIES_ID_DIR)
    log: 'logs/chloroplast_gene_variants/chloroplast_{gene}_variants.log'
    conda: '../envs/species_id.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time = '00:30:00'
    params:
        region = lambda w: 'VCDJ01010680.1:6656-8083' if w.gene == 'rbcl' else 'VCDJ01010680.1:3872-5392'
    shell:
        """
        ( bcftools mpileup \
            --output-type u \
            --fasta-ref {0} \
            --count-orphans \
            --regions {{params.region}} \
            --bam-list {{input}} | \
            bcftools call \
            --output-type z \
            --variants-only \
            --multiallelic-caller \
            --ploidy 1 \
            --output {{output}} ) 2> {{log}}
        """.format(REFERENCE_GENOME)

rule chloroplast_gene_fasta:
    input:
       REFERENCE_GENOME
    output:
       '{0}/{{gene}}/{{gene}}.fna'.format(SPECIES_ID_DIR)
    log: 'logs/chloroplast_gene_fasta/{gene}_fasta.log'
    conda: '../envs/species_id.yaml'
    params:
        region = lambda w: 'VCDJ01010680.1:6656-8083' if w.gene == 'rbcl' else 'VCDJ01010680.1:3872-5392'
    shell:
        """
        samtools faidx {input} {params.region} > {output} 2> {log}
        """

rule index_chloroplast_gene_vcf:
    input:
        rules.bcftools_chloroplast_gene_variants.output
    output:
        '{0}/{{gene}}/allSamples_{{gene}}.vcf.gz.tbi'.format(SPECIES_ID_DIR)
    log: 'logs/index_chloroplast_gene_vcf/index_{gene}_vcf.log'
    conda: '../envs/species_id.yaml'
    shell:
        """
        tabix {input} 2> {log}
        """

rule chloroplast_gene_consensus:
    input:
        ref = rules.chloroplast_gene_fasta.output,
        vcf = rules.bcftools_chloroplast_gene_variants.output,
        idx = rules.index_chloroplast_gene_vcf.output
    output:
        temp('{0}/{{gene}}/consensus_fasta/{{sample}}_{{gene}}.fna'.format(SPECIES_ID_DIR))
    log: 'logs/{gene}_consensus/{sample}_{gene}_consensus.log'
    conda: '../envs/species_id.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1000,
        time = '00:30:00'
    shell:
        """
        bcftools consensus --fasta-ref {input.ref} \
            --sample {wildcards.sample} \
            {input.vcf} > {output} 2> {log}
        """

rule concat_fasta:
    input:
        get_fastas_to_concat
    output:
        '{0}/{{gene}}/consensus_fasta/allSamples_{{gene}}.fasta'.format(SPECIES_ID_DIR)
    log: 'logs/concat_fasta/{gene}_concat_fastas.log'
    run:
        import os
        with open(output[0], 'w') as fout:
            for fasta in input:
                sample = os.path.basename(fasta).split('_{0}'.format(wildcards.gene))[0]
                with open(fasta, 'r') as fin:
                    lines = fin.readlines()
                    seq = ''.join(line.strip() for line in lines[1:])
                    fout.write('>{0};{1}\n{2}\n'.format(sample, wildcards.gene, seq))

rule download_nr_database:
    output:
        '{0}/ncbi_nr_database/nr.gz'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/download_nr_database/download_nr_database.log'
    params:
        out = '{0}/ncbi_nr_database/'.format(PROGRAM_RESOURCE_DIR)
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz -P {params.out} 2> {log}
        """

rule species_id_done:
    input:
        expand(rules.concat_fasta.output, gene = ['rbcl','matk']),
        '{0}/ncbi_nr_database/nr.gz'.format(PROGRAM_RESOURCE_DIR)
    output:
        '{0}/species_id.done'.format(SPECIES_ID_DIR)
    shell:
        """
        touch {output}
        """

