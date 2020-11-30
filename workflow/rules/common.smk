# Python functions used throughout snakemake workflow

def create_raw_read_dict(RAW_READ_DIR, SAMPLES):
    raw_read_dict = {}
    for sample in SAMPLES:
        R1 = glob.glob('{0}/{1}/{1}_*_1.fq.gz'.format(RAW_READ_DIR, sample))[0]
        R2 = glob.glob('{0}/{1}/{1}_*_2.fq.gz'.format(RAW_READ_DIR, sample))[0]
        raw_read_dict[sample] = {'R1': R1, 'R2': R2}
    return raw_read_dict

def get_representative_bam(wildcards):
    bam_index_files = expand(rules.index_bam.output, sample=SAMPLES)
    for i in bam_index_files:
        if REPRESENTATIVE_SAMPLE in i:
            bam = os.path.splitext(i)[0]
    return bam

def get_tabix_files(wildcards):
    if wildcards.site_type == 'allSites':
        vcf_in = rules.bcftools_sort.output
    else:
        vcf_in = rules.bcftools_split_variants.output
    return vcf_in
