# Python functions used throughout snakemake workflow

def create_raw_read_dict(RAW_READ_DIR, SAMPLES):
    raw_read_dict = {}
    for sample in SAMPLES:
        R1 = glob.glob('{0}/{1}/{1}_*_1.fq.gz'.format(RAW_READ_DIR, sample))[0]
        R2 = glob.glob('{0}/{1}/{1}_*_2.fq.gz'.format(RAW_READ_DIR, sample))[0]
        raw_read_dict[sample] = {'R1': R1, 'R2': R2}
    return raw_read_dict

def get_fastas_to_concat(wildcards):
    if wildcards.gene == 'rbcl':
        return expand(rules.chloroplast_gene_consensus.output, sample=SAMPLES, gene='rbcl')
    elif wildcards.gene == 'matk':
        return expand(rules.chloroplast_gene_consensus.output, sample=SAMPLES, gene='matk')

def get_toronto_bam(wildcards):
    bam = glob.glob('{0}/{1}_*.bam'.format(TOR_BAMS, wildcards.tor_sample))
    return bam

def get_bed_to_subset(wildcards):
    all_bed_files = rules.get_fourfold_zerofold.output
    bed = [bed for bed in all_bed_files if wildcards.site in os.path.basename(bed)]
    return bed

def get_bams_for_angsd_gls(wildcards):
    if wildcards.sample_set == 'highQualSamples':
        return rules.create_bam_list_highQualSamples.output
    elif wildcards.sample_set == 'finalSamples_relatedRemoved':
        return rules.create_bam_list_allFinalSamples.output

def angsd_sfs_input(wildcards):
    saf_idx = rules.angsd_saf_likelihood_allSites.output.saf_idx
    sites_idx = rules.angsd_index_sites.output.idx
    if wildcards.site == 'allSites':
        sites = rules.extract_angsd_allSites.output
    else:
        sites = rules.split_angsd_sites_byChrom.output
    return { 'saf_idx' : saf_idx, 'sites_idx' : sites_idx, 'sites' : sites }

def angsd_estimate_thetas_input(wildcards):
    saf_idx = rules.angsd_saf_likelihood_allSites.output.saf_idx
    sfs = rules.angsd_estimate_sfs.output
    sites_idx = rules.angsd_index_sites.output.idx
    if wildcards.site == 'allSites':
        sites = rules.extract_angsd_allSites.output
    else:
        sites = rules.split_angsd_sites_byChrom.output
    return { 'saf_idx' : saf_idx, 'sfs' : sfs, 'sites_idx' : sites_idx, 'sites' : sites }

def get_angsd_stats_toConcat(wildcards):
    if wildcards.site == '0fold':
        return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='0fold', sample_set='finalSamples_relatedRemoved')
    elif wildcards.site == '4fold':
        return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='4fold', sample_set='finalSamples_relatedRemoved')
    else:
        return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='allSites', sample_set='finalSamples_relatedRemoved')

def get_angsd_sfs_toConcat(wildcards):
    if wildcards.site == '0fold':
        return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='0fold', sample_set='finalSamples_relatedRemoved')
    elif wildcards.site == '4fold':
        return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='4fold', sample_set='finalSamples_relatedRemoved')
    else:
        return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='allSites', sample_set='finalSamples_relatedRemoved')

def get_sites_for_angsd_index(wildcards):
    if wildcards.site == 'allSites':
        return rules.extract_angsd_allSites.output
    elif wildcards.site == '0fold' or wildcards.site == '4fold':
        return rules.split_angsd_sites_byChrom.output

def get_angsd_gl_toConcat(wildcards):
    if wildcards.site == '0fold' and wildcards.maf == '0.05' and wildcards.sample_set == 'highQualSamples':
        return expand(rules.subset_angsd_gl.output, site='0fold', maf='0.05', chrom=CHROMOSOMES, sample_set='highQualSamples')
    elif wildcards.site == '0fold' and wildcards.maf == '0.05' and wildcards.sample_set == 'finalSamples_relatedRemoved':
        return expand(rules.subset_angsd_gl.output, site='0fold', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_relatedRemoved')
    elif wildcards.site == '4fold' and wildcards.maf == '0.05' and wildcards.sample_set == 'highQualSamples':
        return expand(rules.subset_angsd_gl.output, site='4fold', maf='0.05', chrom=CHROMOSOMES, sample_set='highQualSamples')
    elif wildcards.site == '4fold' and wildcards.maf == '0.05' and wildcards.sample_set == 'finalSamples_relatedRemoved':
        return expand(rules.subset_angsd_gl.output, site='4fold', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_relatedRemoved')
    elif wildcards.site == 'allSites' and wildcards.maf == '0.05' and wildcards.sample_set == 'highQualSamples':
        return expand(rules.angsd_gl_allSites.output.gls, maf='0.05', chrom=CHROMOSOMES, sample_set='highQualSamples')
    elif wildcards.site == 'allSites' and wildcards.maf == '0.05' and wildcards.sample_set == 'finalSamples_relatedRemoved':
        return expand(rules.angsd_gl_allSites.output.gls, maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_relatedRemoved')

def get_angsd_maf_toConcat(wildcards):
    if wildcards.site == '0fold' and wildcards.maf == '0.05' and wildcards.sample_set == 'highQualSamples':
        return expand(rules.subset_angsd_maf.output, site='0fold', maf='0.05', chrom=CHROMOSOMES, sample_set='highQualSamples')
    elif wildcards.site == '0fold' and wildcards.maf == '0.05' and wildcards.sample_set == 'finalSamples_relatedRemoved':
        return expand(rules.subset_angsd_maf.output, site='0fold', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_relatedRemoved')
    elif wildcards.site == '4fold' and wildcards.maf == '0.05' and wildcards.sample_set == 'highQualSamples':
        return expand(rules.subset_angsd_maf.output, site='4fold', maf='0.05', chrom=CHROMOSOMES, sample_set='highQualSamples')
    elif wildcards.site == '4fold' and wildcards.maf == '0.05' and wildcards.sample_set == 'finalSamples_relatedRemoved':
        return expand(rules.subset_angsd_maf.output, site='4fold', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_relatedRemoved')
    elif wildcards.site == 'allSites' and wildcards.maf == '0.05' and wildcards.sample_set == 'highQualSamples':
        return expand(rules.angsd_gl_allSites.output.mafs, maf='0.05', chrom=CHROMOSOMES, sample_set='highQualSamples')
    elif wildcards.site == 'allSites' and wildcards.maf == '0.05' and wildcards.sample_set == 'finalSamples_relatedRemoved':
        return expand(rules.angsd_gl_allSites.output.mafs, maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_relatedRemoved')
