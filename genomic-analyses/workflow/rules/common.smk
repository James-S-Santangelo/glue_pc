# Python functions used throughout snakemake workflow

def get_all_bams(BAM_DIR):
    """
    Returns list with paths to GLUE bams 
    """
    bams = expand(BAM_DIR + '/{sample}_{site}.bam', sample=SAMPLES, site='4fold') 
    return bams

def get_bed(wildcards):
    """
    Get correct BED file for conversion to ANGSD sites format
    """
    bed = expand(rules.get_fourfold_zerofold.output, site=['4fold', '0fold'])
    bed_out = [x for x in bed if wildcards.site in os.path.basename(x)]
    return bed_out

def get_angsd_gl_toConcat(wildcards):
    """
    Returns list with correct genotype likelihood files for concatenation, depending on combination
    of "sample_set", "site", and "maf" wildcard values
    """
    out = expand(rules.angsd_gl_allSamples.output.gls, chrom=CHROMOSOMES, site=wildcards.site)
    return out

def get_angsd_maf_toConcat(wildcards):
    """
    Returns list with correct minor allele frequency files for concatenation, depending on combination
    of "sample_set", "site", and "maf" wildcard values
    """
    out = expand(rules.angsd_gl_allSamples.output.mafs, chrom=CHROMOSOMES, site=wildcards.site)
    return out

def get_files_for_saf_estimation_byCity_byHabitat(wildcards):
    """
    Get files to estimate SAF likelihhods for urban and rural habitats by city.
    """
    sites_idx = expand(rules.angsd_index_random_degen_sites.output.idx, site=wildcards.site)
    sites = expand(rules.select_random_degenerate_sites.output, site=wildcards.site)
    ref = rules.glue_dnaSeqQC_unzip_reference.output
    bams = expand(rules.create_bam_list_byCity_byHabitat.output, city=wildcards.city, habitat=wildcards.habitat, site = wildcards.site)
    return { 'bams' : bams, 'sites_idx' : sites_idx , 'sites' : sites, 'ref' : ref }

def get_files_for_alleleFreq_estimation_byCity_byHabitat(wildcards):
    """
    Get files to estimate Allele Frequencies for urban and rural habitats by city.
    """
    sites_idx = expand(rules.angsd_index_city_snps.output.idx, site=wildcards.site, city=wildcards.city)
    sites = expand(rules.snps_forAlleleFreqs_byCity_byHabitat.output, site=wildcards.site, city=wildcards.city)
    ref = rules.glue_dnaSeqQC_unzip_reference.output
    bams = expand(rules.create_bam_list_byCity_byHabitat.output, city=wildcards.city, habitat=wildcards.habitat, site = wildcards.site)
    chroms = config['chromosomes']
    return { 'bams' : bams, 'sites_idx' : sites_idx , 'sites' : sites, 'ref' : ref }

def get_habitat_saf_files_byCity(wildcards):
    """
    Returns list with 4fold urban and rural SAF files by city
    """
    city_saf_files = expand(rules.angsd_saf_likelihood_byCity_byHabitat.output.saf_idx, city=wildcards.city, habitat=HABITATS, site=wildcards.site)
    return city_saf_files

def get_urban_rural_bam_lists(wildcards):
    """
    Collect files with paths to urban and rural bams by City. Return as dictionary. 
    """
    urban = expand(rules.create_bam_list_byCity_byHabitat.output, city=wildcards.city, habitat='u', site=wildcards.site)[0]
    rural = expand(rules.create_bam_list_byCity_byHabitat.output, city=wildcards.city, habitat='r', site=wildcards.site)[0]
    return { 'urban_bams' : urban, 'rural_bams' : rural }

def get_files_for_permuted_saf_estimation(wildcards):
    """
    Get files to estimate SAF likelihoods for permuted versions of "urban" and "rural" populations
    """
    sites_idx = expand(rules.angsd_index_random_degen_sites.output.idx, site=wildcards.site)
    sites = expand(rules.select_random_degenerate_sites.output, site=wildcards.site)
    ref = rules.glue_dnaSeqQC_unzip_reference.output
    if wildcards.habitat == 'u':
        bams = expand(rules.create_random_bam_list_byCity_byHabitat.output.urban, city=wildcards.city, seed=wildcards.seed, site=wildcards.site)
    elif wildcards.habitat == 'r':
        bams = expand(rules.create_random_bam_list_byCity_byHabitat.output.rural, city=wildcards.city, seed=wildcards.seed, site=wildcards.site)
    return { 'bams' : bams, 'sites_idx' : sites_idx , 'sites' : sites, 'ref' : ref }

def get_habitat_saf_files_byCity_permuted(wildcards):
    """
    Returns list with 4fold urban and rural SAF files by city
    """
    city_saf_files = expand(rules.angsd_permuted_saf_likelihood_byCity_byHabitat.output.saf_idx, city=wildcards.city, habitat=HABITATS, site=wildcards.site, seed=wildcards.seed)
    return city_saf_files

def get_bamLists_toConcat(wildcards):
    """
    Collect text files with paths to urban and rural bams by city
    """
    all_bam_lists = expand(rules.create_bam_list_byCity_byHabitat.output, city = wildcards.city, habitat = HABITATS, site = wildcards.site)
    return all_bam_lists

def get_bams_for_read_counts(wildcards):
    """
    Returns the correct GLUE or Toronto BAM file
    """
    tor_bams = expand(rules.glue_dnaSeqQC_downsample_toronto_bam.output, sample=TOR_SAMPLES)
    glue_bams = expand(rules.glue_dnaSeqQC_samtools_markdup.output.bam, sample=SAMPLES)
    glue_bams = [bam for bam in glue_bams if not os.path.basename(bam).startswith('s_')]
    return tor_bams + glue_bams

def get_files_for_saf_estimation_snps_hcn_chroms(wildcards):
    """
    Get files to estimate SAF likelihhods for urban and rural habitats by city.
    """
    if wildcards.gene == 'li':
        sites_idx = expand(rules.index_chromosomal_angsd_sites.output.idx, chrom='CM019108.1', site='4fold')
        sites = expand(rules.split_angsd_sites_byChrom.output, chrom='CM019108.1', site='4fold')
    elif wildcards.gene == 'ac':
        sites_idx = expand(rules.index_chromosomal_angsd_sites.output.idx, chrom='CM019103.1', site='4fold')
        sites = expand(rules.split_angsd_sites_byChrom.output, chrom='CM019103.1', site='4fold')
    ref = rules.glue_dnaSeqQC_unzip_reference.output
    bams = expand(rules.create_bam_list_byCity_byHabitat.output, city=wildcards.city, habitat=wildcards.habitat, site=wildcards.site)
    return { 'bams' : bams, 'sites_idx' : sites_idx , 'sites' : sites, 'ref' : ref }

def get_habitat_saf_files_byCity_hcn_chroms(wildcards):
    """
    Returns list with 4fold urban and rural SAF files by city
    """
    saf_files = expand(rules.angsd_saf_likelihood_snps_hcn_chroms.output.saf_idx, city=wildcards.city, habitat=HABITATS, site=wildcards.site, gene=wildcards.gene)
    return saf_files

def get_sites_for_hcn_chrom_analysis(wildcards):
    """
    Return sites file for degenerate sites along the correct chromosome, based on 'gene' wildcard
    """
    if wildcards.gene == 'li':
        sites = expand(rules.split_angsd_sites_byChrom.output, chrom = 'CM019108.1', site=wildcards.site)
    elif wildcards.gene == 'ac':
        sites = expand(rules.split_angsd_sites_byChrom.output, chrom = 'CM019103.1', site=wildcards.site)
    return sites
