rule create_regions_equal_coverage:
    input: 
        qc_done = rules.multiqc.output,
        bam = get_representative_bam,
        ref_index = '{0}.fai'.format(REFERENCE_GENOME),
    output:
        temp('{0}/{{chrom}}_forFreebayes.regions'.format(PROGRAM_RESOURCE_DIR))
    log: 'logs/create_regions_equal_cov/{chrom}_create_regions_equal_cov.log'
    conda: '../envs/variant_calling.yaml'
    threads: 8
    resources:
        time = '06:00:00'
    shell:
        """
        ( samtools view --threads {{threads}} -b -s 0.20 {{input.bam}} {{wildcards.chrom}} |\
            bamtools coverage |\
            coverage_to_regions.py {{input.ref_index}} {0} > {{output}} ) 2> {{log}} 
        """.format(NUM_REGIONS_PER_CHROM)

rule region_files_forFreebayes:
    input:
        rules.create_regions_equal_coverage.output
    output:
        '{0}/{{chrom}}_{{node}}_forFreebayes.regions'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/regions_files_forFreebayes/{chrom}_{node}_forFreebayes.log'
    run:
        with open(input[0], 'r') as fin:
            inpath = os.path.dirname(input[0])
            if NODES_PER_CHROM == 1:
                out = '{0}/{1}_node1_forFreebayes.regions'.format(inpath, wildcards.chrom)
                with open(out, 'w') as f:
                    f.write('worked')

rule create_bam_list:
    input:
        expand(rules.samtools_markdup.output.bam, sample=SAMPLES)
    output:
        '{0}/freebayes_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/create_bam_list/create_bam_list.log'
    run:
        with open(output[0], 'w') as f:
            for bam in input:
                f.write('{0}\n'.format(bam))

# rule freebayes_call_variants:
#     input:
#         bams = rules.create_bam_list.output,
#         regions = rules.create_regions_equal_coverage.output
#         #regions = '../resources/{chrom}_forFreebayes.regions'
#     output:
#         temp('{0}/vcf/{{chrom}}/{{chrom}}_allSamples.vcf'.format(VARIANT_DIR))
#     log: 'logs/freebayes/{chrom}_freebayes.log'
#     conda: '../envs/variant_calling.yaml'
#     threads: 5
#     resources:
#         time = '24:00:00'
#     shell:
#         """
#         ( freebayes-parallel {{input.regions}} {{threads}} \
#             --fasta-reference {0} \
#             --bam-list {{input.bams}} \
#             --use-best-n-alleles 4 \
#             --report-monomorphic \
#             --max-complex-gap 1 \
#             --haplotype-length 1 \
#             --genotype-qualities > {{output}} ) 2> {{log}}
#         """.format(REFERENCE_GENOME)
# 
# rule bgzip_vcf:
#     input:
#         rules.freebayes_call_variants.output
#     output:
#         temp('{0}/vcf/{{chrom}}/{{chrom}}_allSamples.vcf.gz'.format(VARIANT_DIR))
#     log: 'logs/bgzip/{chrom}_bgzip.log'
#     conda: '../envs/variant_calling.yaml',
#     threads: 8
#     resources:
#         time = '01:00:00'
#     shell:
#         """
#         bgzip -@ {threads} {input}
#         """
# 
# rule bcftools_sort:
#     input:
#         rules.bgzip_vcf.output
#     output:
#         '{0}/vcf/{{chrom}}/{{chrom}}_allSamples_sorted.vcf.gz'.format(VARIANT_DIR)
#     log: 'logs/bcftools_sort/{chrom}_bcftools_sort.log'
#     conda: '../envs/variant_calling.yaml'
#     resources:
#         time = '12:00:00'
#     shell:
#         """
#         mkdir {0}/{{wildcards.chrom}};
#         bcftools sort -O z -o {{output}} -T {0}/{{wildcards.chrom}} {{input}} 2> {{log}}
#         """.format(TMPDIR)
# 
# rule bcftools_split_variants:
#     input:
#         vcf = rules.bcftools_sort.output
#     output:
#         '{0}/vcf/{{chrom}}/{{chrom}}_allSamples_{{site_type}}_sorted.vcf.gz'.format(VARIANT_DIR)
#     log: 'logs/bcftools_split_variants/{chrom}_bcftools_split_variants_{site_type}.log'
#     conda: '../envs/variant_calling.yaml'
#     wildcard_constraints:
#         site_type='snps|indels|invariant|mnps|other'
#     threads: 8
#     resources:
#         time = '12:00:00'
#     shell:
#         """
#         if [ {{wildcards.site_type}} = 'invariant' ]; then
#             bcftools view --threads {{threads}} -O z --include 'N_ALT = 0' {{input}} > {{output}} 2> {{log}}
#         elif [ {{wildcards.site_type}} = 'snps' ]; then
#             ( bcftools view --threads {{threads}} -O v --types {{wildcards.site_type}} {{input}} |\
#             vcfallelicprimitives --keep-info --keep-geno |\
#             bcftools view --threads {{threads}} --types {{wildcards.site_type}} --min-alleles 2 --max-alleles 2 |\
#             bcftools sort -O z -T {0}/{{wildcards.chrom}} -o {{output}} ) 2> {{log}}
#         else
#             bcftools view --threads {{threads}} -O z --types {{wildcards.site_type}} {{input}} > {{output}} 2> {{log}}
#         fi
#         """.format(TMPDIR)
# 
# rule tabix_vcf:
#     input: get_tabix_files
#     output:
#         '{0}/vcf/{{chrom}}/{{chrom}}_allSamples_{{site_type}}_sorted.vcf.gz.tbi'.format(VARIANT_DIR)
#     log: 'logs/tabix/{chrom}_tabix_{site_type}.log'
#     conda: '../envs/variant_calling.yaml'
#     shell:
#         """
#         tabix {input}
#         """
# 
# rule vcf_to_zarr:
#     input:
#         rules.bcftools_split_variants.output
#     output:
#         directory('{0}/zarr_db/{{chrom}}/{{chrom}}_allSamples_{{site_type}}_sorted.zarr'.format(VARIANT_DIR))
#     log: 'logs/vcf_to_zarr/{chrom}_vcf_to_zarr_{site_type}.log'
#     conda: '../envs/variant_calling.yaml'
#     wildcard_constraints:
#         site_type='snps|invariant'
#     threads: 1
#     resources:
#         time = '01:00:00'
#     script:
#         "../scripts/python/vcf_to_zarr.py"
