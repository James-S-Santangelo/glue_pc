# Rules to get allele frequencies at Ac and Li loci from read count data

rule read_count_data:
    input:
        expand(rules.samtools_markdup.output.bam, sample = SAMPLES)
    output:
        '{0}/{{gene}}_read_counts.txt'.format(HCN_LOCI_DIR)
    log: 'logs/read_count_data/{gene}_counts.log'
    conda: '../envs/mapping.yaml'
    params: 
        region = lambda w: 'CM019108.1:30218214-30229250 CM019108.1:30230911-30247247' if w.gene == 'li' else 'CM019103.1:19559221-19573344' 
    shell:
        """
        ( for bam in {input}; do
            count=( $( samtools view $bam {params.region} | wc -l ) )
            printf "%s\t%s\n" $bam $count >> {output}
        done ) 2> {log}
        """

rule hcn_loci_freq_done:
    input:
        expand(rules.read_count_data.output, gene=['ac', 'li'])
    output:
        '{0}/hcn_loci_freq_done'.format(HCN_LOCI_DIR)
    shell:
        """
        touch {output}
        """
