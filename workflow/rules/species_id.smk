rule bcftools_chloroplast_gene_variants:
    input:
        rules.create_bam_list.output
    output:
        '{0}/{{gene}}/allSamples_{{gene}}.vcf.gz'.format(SPECIES_ID_DIR)
    log: 'logs/chloroplast_gene_variants/chloroplast_{gene}_variants.log'
    conda: '../envs/species_id.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
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
       '{0}/{{gene}}/gene.fna'.format(SPECIES_ID_DIR)
    log: 'logs/chloroplast_gene_fasta/{gene}_fasta.log'
    conda: '../envs/species_id.yaml'
    params:
        region = lambda w: 'VCDJ01010680.1:6656-8083' if w.gene == 'rbcl' else 'VCDJ01010680.1:3872-5392'
    shell:
        """
        samtools faidx {input} {params.region} > {output} 2> {log}
        """

rule reheader_chloroplast_gene_vcf:
    input:
        rules.bcftools_chloroplast_gene_variants.output
    output:
        '{0}/{{gene}}/allSamples_{{gene}}_reheader.vcf.gz'.format(SPECIES_ID_DIR)
    log: 'logs/reheader_chloroplast_gene_vcf/reheader_{gene}_vcf.log'
    conda: '../envs/species_id.yaml'
    shell:
        """
        cat {0} | tail -n +2 | cut -f1 > tmp.samples;
        bcftools reheader --samples tmp.samples {{input}} > {{output}} 2> {{log}};
        rm tmp.samples
        """.format(config['samples'])

rule index_chloroplast_gene_vcf:
    input:
        rules.reheader_chloroplast_gene_vcf.output
    output:
        '{0}/{{gene}}/allSamples_{{gene}}_reheader.vcf.gz.tbi'.format(SPECIES_ID_DIR)
    log: 'logs/index_chloroplast_gene_vcf/index_{gene}_vcf.log'
    conda: '../envs/species_id.yaml'
    shell:
        """
        tabix {input} 2> {log}
        """

rule chloroplast_gene_consensus:
    input:
        ref = rules.chloroplast_gene_fasta.output,
        vcf = rules.reheader_chloroplast_gene_vcf.output,
        idx = rules.index_chloroplast_gene_vcf.output
    output:
        '{0}/{{gene}}/consensus_fasta/{{sample}}_{{gene}}.fna'.format(SPECIES_ID_DIR)
    log: 'logs/{gene}_consensus/{sample}_{gene}_consensus.log'
    conda: '../envs/species_id.yaml'
    shell:
        """
        bcftools consensus --fasta-ref {input.ref} \
            --sample {wildcards.sample} \
            {input.vcf} > {output} 2> {log}
        """
