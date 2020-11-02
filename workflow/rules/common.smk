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

