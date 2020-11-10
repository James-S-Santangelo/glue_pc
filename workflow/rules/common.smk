# Python functions used throughout snakemake workflow

def get_raw_reads():
    all_raw_reads = []
    reads = [1, 2]
    for sample in SAMPLES:
        for read in reads:
            read_path = str(glob.glob('{0}/{1}/{1}_*_{2}.fq.gz'.format(RAW_READ_DIR, sample, read))[0])
            all_raw_reads.append(read_path)
    return all_raw_reads

def get_basename(filename):
    basename = os.path.basename(filename).split(os.extsep)[0]
    return basename 

def fastqc_target_files():
    html_files = ['{0}/{1}{2}'.format(RAW_READ_FASTQC_DIR, get_basename(f), '_fastqc.html') for f in get_raw_reads()]
    zip_files = ['{0}/{1}{2}'.format(RAW_READ_FASTQC_DIR, get_basename(f), '_fastqc.zip') for f in get_raw_reads()]
    return html_files + zip_files

def get_representative_bam(wildcards):
    bam_index = checkpoints.index_bam.get(sample = REPRESENTATIVE_SAMPLE).output[0]
    bam = os.path.splitext(bam_index)[0]
    return bam

def get_tabix_files(wildcards):
    if wildcards.site_type == 'allSites':
        vcf_in = rules.bcftools_sort.output
    else:
        vcf_in = rules.bcftools_split_variants.output
    return vcf_in
