# Python functions used throughout snakemake workflow

def get_subset_bams_degeneracy_input(wildcards):
    """
    Returns the correct GLUE or Toronto BAM file
    """
    all_degen_bed_files = expand(rules.get_fourfold_zerofold.output, site=['0fold', '4fold'])
    regions = [x for x in all_degen_bed_files if wildcards.site in os.path.basename(x)]
    if wildcards.sample.startswith('s_'):
        bam = expand(rules.glue_dnaSeqQC_downsample_toronto_bam.output, sample=wildcards.sample)
        idx = expand(rules.glue_dnaSeqQC_index_toronto_bam.output, sample = wildcards.sample)
    else:
        bam = expand(rules.glue_dnaSeqQC_samtools_markdup.output.bam, sample=wildcards.sample)
        idx = expand(rules.glue_dnaSeqQC_index_bam.output, sample = wildcards.sample)
    return { 'bam' : bam, 'idx' : idx, 'regions' : regions }

def get_all_bams(wildcards):
    """
    Returns list with paths to GLUE bams 
    """
    bams = expand(rules.subset_bams_degeneracy.output, sample=SAMPLES, site=wildcards.site)
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
    sites_idx = expand(rules.angsd_index_sites_allChroms.output.idx, site=wildcards.site)
    sites = expand(rules.convert_sites_for_angsd.output, site=wildcards.site)
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
    all_bam_lists = expand(rules.create_bam_list_byCity_byHabitat.output, city = wildcards.city, habitat = ['u', 'r'])
    return all_bam_lists

def get_files_saf_estimation_byPopulation(wildcards):
    """
    Collect files required to estimate SAF likelihoods by population in each city. Population ID extracted
    from files created by Checkpoint.
    """
    bams = '{0}/bam_lists/by_city/{{city}}/by_pop/{{city}}_{{popu}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    sites_idx = expand(rules.angsd_index_degenerate.output.idx, chrom='CM019101.1', site='4fold')
    sites = expand(rules.split_angsd_sites_byChrom.output, chrom='CM019101.1', site='4fold') 
    ref = REFERENCE_GENOME
    return { 'bams' : bams, 'sites_idx' : sites_idx, 'sites' : sites, 'ref' : ref }
    
def aggregate_input_theta(wildcards):
    """
    Collect population ID ('popu') wildcard values from checkpoint.
    Collects output for estimation of thetas
    """
    checkpoint_output = checkpoints.populations_byCity_byHabitat.get(**wildcards).output[0]
    pops = glob_wildcards(os.path.join(checkpoint_output, '{{city}}_{{popu}}_bams.list'.format(PROGRAM_RESOURCE_DIR))).popu
    return expand('{0}/summary_stats/thetas/by_city/{{city}}/by_pop/{{city}}_{{popu}}_{{site}}.thetas.idx.pestPG'.format(ANGSD_DIR), city=wildcards.city, popu=pops, site='4fold')

def get_population_saf_files_byCity(wildcards):
    """
    Get SAF files for two populations for which to estimate joint SFS.
    """
    checkpoint_output = checkpoints.populations_byCity_byHabitat.get(**wildcards).output[0]
    pops = glob_wildcards(os.path.join(checkpoint_output, '{{city}}_{{popu}}_bams.list'.format(PROGRAM_RESOURCE_DIR))).popu
    all_saf_files = expand(rules.angsd_saf_likelihood_byCity_byPopulation.output.saf_idx, city = wildcards.city, popu=pops, site='4fold')
    pop1 = wildcards.pop_comb.split('_')[0]
    pop2 = wildcards.pop_comb.split('_')[1]
    saf1 = [x for x in all_saf_files if '_{0}_'.format(pop1) in os.path.basename(x)]
    saf2 = [x for x in all_saf_files if '_{0}_'.format(pop2) in os.path.basename(x)]
    return saf1 + saf2

def get_population_saf_and_sfs_files_byCity(wildcards):
    """
    Get SAF and SFS files for two populations for which Fst should be estimated. 
    """
    checkpoint_output = checkpoints.populations_byCity_byHabitat.get(**wildcards).output[0]
    pops = glob_wildcards(os.path.join(checkpoint_output, '{{city}}_{{popu}}_bams.list'.format(PROGRAM_RESOURCE_DIR))).popu
    all_saf_files = expand(rules.angsd_saf_likelihood_byCity_byPopulation.output.saf_idx, city = wildcards.city, popu=pops, site='4fold')
    pop1 = wildcards.pop_comb.split('_')[0]
    pop2 = wildcards.pop_comb.split('_')[1]
    saf1 = [x for x in all_saf_files if '_{0}_'.format(pop1) in os.path.basename(x)]
    saf2 = [x for x in all_saf_files if '_{0}_'.format(pop2) in os.path.basename(x)]
    saf_files = saf1 + saf2
    sfs = expand(rules.angsd_estimate_joint_sfs_populations.output, city = wildcards.city, site='4fold', pop_comb=wildcards.pop_comb)
    return { 'saf_files' : saf_files, 'sfs' : sfs }


def aggregate_input_fst(wildcards):
    """
    Collect population ID ('popu') wildcard values from checkpoint.
    Collects output for estimation of fst
    """
    checkpoint_output = checkpoints.populations_byCity_byHabitat.get(**wildcards).output[0]
    pops = glob_wildcards(os.path.join(checkpoint_output, '{{city}}_{{popu}}_bams.list'.format(PROGRAM_RESOURCE_DIR))).popu
    pop_combinations = [c[0] + '_' + c[1] for c in list(itertools.combinations(pops, 2))]
    return expand('{0}/summary_stats/fst/fst1/{{city}}/pairwise/{{city}}_{{site}}_{{pop_comb}}_readable.fst'.format(ANGSD_DIR), city=wildcards.city, pop_comb=pop_combinations, site='4fold')

def get_ngsadmix_logfiles_byCity(wildcards):
    """
    Get NGSadmix logfiles for each city. Used to generate input file for CLUMPAK
    """
    return expand(rules.ngsadmix.output.lf, city=wildcards.city, k=NGSADMIX_K, seed=NGSADMIX_SEEDS, site='4fold', maf='0.05')

def get_files_for_saf_estimation_snps_hcn_chroms(wildcards):
    """
    Get files to estimate SAF likelihhods for urban and rural habitats by city.
    """
    if wildcards.gene == 'li':
        sites_idx = expand(rules.angsd_index_degenerate.output.idx, chrom='CM019108.1', site='4fold')
        sites = expand(rules.split_angsd_sites_byChrom.output, chrom='CM019108.1', site='4fold')
    elif wildcards.gene == 'ac':
        sites_idx = expand(rules.angsd_index_degenerate.output.idx, chrom='CM019103.1', site='4fold')
        sites = expand(rules.split_angsd_sites_byChrom.output, chrom='CM019103.1', site='4fold')
    ref = REFERENCE_GENOME
    bams = expand(rules.create_bam_list_byCity_byHabitat.output, city=wildcards.city, habitat=wildcards.habitat)
    return { 'bams' : bams, 'sites_idx' : sites_idx , 'sites' : sites, 'ref' : ref }

def get_habitat_saf_files_byCity_hcn_chroms(wildcards):
    """
    Returns list with 4fold urban and rural SAF files by city
    """
    saf_files = expand(rules.angsd_saf_likelihood_snps_hcn_chroms.output.saf_idx, city=wildcards.city, habitat=HABITATS, site=['4fold'], gene=wildcards.gene)
    return saf_files
