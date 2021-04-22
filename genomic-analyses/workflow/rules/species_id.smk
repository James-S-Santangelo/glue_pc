# Rules to generate per-sample consensus sequences of matK and rbcL chloroplast genes. 
# Sequences are blasted against NCBI nt database to identify species.
# Used during QC to confirm species identity. Likely not actually needed since Qualimap 
# alignment error rate seems to catch mistaken species. 

rule create_bam_list_forSpeciesID:
    """
    Create text files with paths to all BAM files
    """
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
    """
    Call variants in matK and rbcL genes using bcftools. Generates single VCF for each gene. 
    """
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
    """
    Extract FASTA sequence for matK and rbcL genes from reference genome using samtools faidx.
    One FASTA for each gene.
    """
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
    """
    Index VCF files with variants in each gene. Required for generating consensus
    """
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
    """
    Use genotype calls in VCFs to generate per-sample matK and rbcL consensus sequences. 
    Outputs two FASTAs per sample (one for each gene)
    """
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
    """
    Concatenate per-sample consensus sequences into single FASTA file. One FASTA for each gene.
    """
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

rule download_nt_database:
    """
    Download NCBI "nt" database. 

    TODO: Checkpoint this since database will likely grow in the future. Currently expects a fixed
    number of database files (specified in configfile)
    """
    output:
        file_db = temp(expand('{0}/ncbi_nt_database/nt.{{num}}.{{ext}}'.format(PROGRAM_RESOURCE_DIR), num=NT_DB_FILE, ext=NT_DB_FILE_EXT)), 
        db = temp(expand('{0}/ncbi_nt_database/nt.{{ext}}'.format(PROGRAM_RESOURCE_DIR), ext=NT_DB_FINAL_EXT)),
        taxon = temp(expand('{0}/ncbi_nt_database/taxdb.{{ext}}'.format(PROGRAM_RESOURCE_DIR), ext=['btd', 'bti'])),
        done = '{0}/ncbi_nt_database/nt_db_download.done'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/download_nt_database/download_nt_database.log'
    conda: '../envs/species_id.yaml'
    params:
        out = '{0}/ncbi_nt_database/'.format(PROGRAM_RESOURCE_DIR)
    shell:
        """
        update_blastdb.pl --decompress nt 2> {log}
        mv nt* {params.out};
        touch {output.done}
        """

rule blast_chloroplast_genes:
    """
    BLAST matK and rbcL FASTA sequences against local NCBI "nt" database using blastn. 
    Outputs single tab-separated text files for each gene.
    """
    input:
        fasta = rules.concat_fasta.output,
        done = rules.download_nt_database.output.done
    output:
        '{0}/{{gene}}/{{gene}}_blast_results.txt'.format(SPECIES_ID_DIR)
    log: 'logs/blast_chloroplast_genes/{gene}_blast.log'
    conda: '../envs/species_id.yaml'
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '03:00:00'
    params:
        db = '{0}/ncbi_nt_database/nt'.format(PROGRAM_RESOURCE_DIR),
        outfmt = "'6 qseqid qlen sseqid slen evalue bitscore lenth pident qcovhsp ssciname scomname'"
    shell:
        """
        blastn -query {input.fasta} \
            -db {params.db} \
            -out {output} \
            -num_threads {threads} \
            -max_hsps 1 \
            -max_target_seqs 1 \
            -outfmt {params.outfmt} 2> {log}
        """

rule species_id_done:
    """
    Collect final output files and write flag file signalling successful completion of species ID
    """
    input:
        expand(rules.concat_fasta.output, gene = ['rbcl','matk']),
        expand(rules.blast_chloroplast_genes.output, gene = ['rbcl','matk'])
    output:
        '{0}/species_id.done'.format(SPECIES_ID_DIR)
    shell:
        """
        touch {output}
        """

