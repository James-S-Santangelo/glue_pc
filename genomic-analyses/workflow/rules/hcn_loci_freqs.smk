# Rules to get allele frequencies at Ac and Li loci from read count data

rule read_count_data:
    """
    Calculate number of reads overlapping target region in Ac or Li locus for each sample
    """
    input:
        expand(rules.samtools_markdup.output.bam, sample = SAMPLES)
    output:
        '{0}/{{gene}}_read_counts.txt'.format(HCN_LOCI_DIR)
    log: 'logs/read_count_data/{gene}_counts.log'
    conda: '../envs/hcn_loci_freqs.yaml'
    params: 
        region = lambda w: 'CM019108.1:30218214-30229250 CM019108.1:30230911-30247247' if w.gene == 'li' else 'CM019103.1:19559221-19573344' 
    shell:
        """
        ( for bam in {input}; do
            count=( $( samtools view $bam {params.region} | wc -l ) )
            printf "%s\t%s\n" $bam $count >> {output}
        done ) 2> {log}
        """

rule calculate_hcn_loci_frequencies:
    """
    Estimate per-sample genotype likelihoods and per-city deletion frequencies from read counts
    at Ac or Li locus. Writes two dataframes to disk per locus.
    """
    input:
        multiqc = rules.multiqc.output,
        counts = rules.read_count_data.output
    output:
        freqs = '{0}/{{gene}}_freqs.txt'.format(HCN_LOCI_DIR),
        likes = '{0}/{{gene}}_GLs.txt'.format(HCN_LOCI_DIR)
    log: 'logs/calculate_hcn_loci_frequencies/{gene}_freqs.log'
    conda: '../envs/hcn_loci_freqs.yaml'
    script:
        "../scripts/python/hcn_loci_GLs_freqs.py"

rule hcn_loci_freq_done:
    """
    Generate empty flag file to signal successful completion of GL and deletion frequency estimation
    """
    input:
        expand(rules.calculate_hcn_loci_frequencies.output, gene=['li', 'ac'])
    output:
        '{0}/hcn_loci_freq_done'.format(HCN_LOCI_DIR)
    shell:
        """
        touch {output}
        """

rule hcn_loci_notebook:
    input:
        rules.hcn_loci_freq_done.output
    output:
        '{0}/hcn_loci_notebook.done'.format(HCN_LOCI_DIR)
    conda: '../envs/notebooks.yaml'
    notebook:
        "../notebooks/hcn_loci.r.ipynb"
