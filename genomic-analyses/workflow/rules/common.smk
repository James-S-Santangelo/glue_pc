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

def get_bams_for_angsd(wildcards):
    if wildcards.sample_set == 'highErrorRemoved':
        return rules.create_bam_list_highErrorRemoved.output
    elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
        return rules.create_bam_list_finalSamples_lowCovRemoved.output

def angsd_sfs_input(wildcards):
    saf_idx = rules.angsd_saf_likelihood_allSites.output.saf_idx
    if wildcards.site == 'allSites':
        sites = rules.extract_angsd_allSites.output.sites
        idx = rules.angsd_index_allSites.output.idx
    else:
        sites = rules.split_angsd_sites_byChrom.output.sites
        idx = rules.angsd_index_degenerate.output.idx
    return { 'saf_idx' : saf_idx, 'sites' : sites, 'sites_idx' : idx }

def angsd_estimate_thetas_input(wildcards):
    saf_idx = rules.angsd_saf_likelihood_allSites.output.saf_idx
    sfs = rules.angsd_estimate_sfs.output
    if wildcards.site == 'allSites':
        sites = rules.extract_angsd_allSites.output.sites
        idx = rules.angsd_index_allSites.output.idx
    else:
        sites = rules.split_angsd_sites_byChrom.output.sites
        idx = rules.angsd_index_degenerate.output.idx
    return { 'saf_idx' : saf_idx, 'sfs' : sfs, 'sites' : sites, 'sites_idx' : idx }

def get_angsd_stats_toConcat(wildcards):
    if wildcards.sample_set == 'highErrorRemoved':
        if wildcards.site == '0fold':
            return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='0fold', sample_set='highErrorRemoved')
        elif wildcards.site == '4fold':
            return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='4fold', sample_set='highErrorRemoved')
        else:
            return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='allSites', sample_set='highErrorRemoved')
    elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
        if wildcards.site == '0fold':
            return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='0fold', sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == '4fold':
            return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='4fold', sample_set='finalSamples_lowCovRemoved')
        else:
            return expand(rules.angsd_diversity_neutrality_stats.output, chrom=CHROMOSOMES, site='allSites', sample_set='finalSamples_lowCovRemoved')

def get_angsd_sfs_toConcat(wildcards):
    if wildcards.sample_set == 'highErrorRemoved':
        if wildcards.site == '0fold':
            return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='0fold', sample_set='highErrorRemoved')
        elif wildcards.site == '4fold':
            return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='4fold', sample_set='highErrorRemoved')
        else:
            return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='allSites', sample_set='highErrorRemoved')
    elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
        if wildcards.site == '0fold':
            return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='0fold', sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == '4fold':
            return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='4fold', sample_set='finalSamples_lowCovRemoved')
        else:
            return expand(rules.angsd_estimate_sfs.output, chrom=CHROMOSOMES, site='allSites', sample_set='finalSamples_lowCovRemoved')

def get_angsd_gl_toConcat(wildcards):
    if wildcards.maf == '0.005':
        if wildcards.site == '0fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_gl.output, site='0fold', maf='0.005', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_gl.output, site='0fold', maf='0.005', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == '4fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_gl.output, site='4fold', maf='0.005', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_gl.output, site='4fold', maf='0.005', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == 'allSites':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.angsd_gl_allSites.output.gls, site='allSites', maf='0.005', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.angsd_gl_allSites.output.gls, site='allSites', maf='0.005', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
    elif wildcards.maf == '0.01':
        if wildcards.site == '0fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_gl.output, site='0fold', maf='0.01', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_gl.output, site='0fold', maf='0.01', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == '4fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_gl.output, site='4fold', maf='0.01', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_gl.output, site='4fold', maf='0.01', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == 'allSites':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.angsd_gl_allSites.output.gls, site='allSites', maf='0.01', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.angsd_gl_allSites.output.gls, site='allSites', maf='0.01', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
    elif wildcards.maf == '0.05':
        if wildcards.site == '0fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_gl.output, site='0fold', maf='0.05', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_gl.output, site='0fold', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == '4fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_gl.output, site='4fold', maf='0.05', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_gl.output, site='4fold', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == 'allSites':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.angsd_gl_allSites.output.gls, site='allSites', maf='0.05', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.angsd_gl_allSites.output.gls, site='allSites', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
    return out

def get_angsd_maf_toConcat(wildcards):
    if wildcards.maf == '0.005':
        if wildcards.site == '0fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_maf.output, site='0fold', maf='0.005', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_maf.output, site='0fold', maf='0.005', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == '4fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_maf.output, site='4fold', maf='0.005', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_maf.output, site='4fold', maf='0.005', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == 'allSites':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.angsd_gl_allSites.output.mafs, site='allSites', maf='0.005', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.angsd_gl_allSites.output.mafs, site='allSites', maf='0.005', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
    elif wildcards.maf == '0.01':
        if wildcards.site == '0fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_maf.output, site='0fold', maf='0.01', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_maf.output, site='0fold', maf='0.01', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == '4fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_maf.output, site='4fold', maf='0.01', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_maf.output, site='4fold', maf='0.01', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == 'allSites':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.angsd_gl_allSites.output.mafs, site='allSites', maf='0.01', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.angsd_gl_allSites.output.mafs, site='allSites', maf='0.01', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
    elif wildcards.maf == '0.05':
        if wildcards.site == '0fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_maf.output, site='0fold', maf='0.05', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_maf.output, site='0fold', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == '4fold':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.subset_angsd_maf.output, site='4fold', maf='0.05', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.subset_angsd_maf.output, site='4fold', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
        elif wildcards.site == 'allSites':
            if wildcards.sample_set == 'highErrorRemoved':
                out = expand(rules.angsd_gl_allSites.output.mafs, site='allSites', maf='0.05', chrom=CHROMOSOMES, sample_set='highErrorRemoved')
            elif wildcards.sample_set == 'finalSamples_lowCovRemoved':
                out = expand(rules.angsd_gl_allSites.output.mafs, site='allSites', maf='0.05', chrom=CHROMOSOMES, sample_set='finalSamples_lowCovRemoved')
    return out

def get_habitat_saf_files_byCity(wildcards):
    all_saf_files = expand(rules.angsd_saf_likelihood_byCity_byHabitat.output.saf_idx, city=CITIES, habitat=HABITATS, site=['4fold'])
    city_saf_files = [x for x in all_saf_files if wildcards.city in x and wildcards.site in x]
    return city_saf_files

def get_toronto_bam_pi_fst_test(wildcards):
    bam = glob.glob('{0}/{1}_*.bam'.format(TOR_BAMS, wildcards.tor_test_sample))
    return bam

def toronto_pi_fst_test_saf_input(wildcards):
    if wildcards.group == 'highCov':
        bams = rules.urban_rural_toronto_bam_lists.output
    elif wildcards.group == 'lowCov':
        bams = rules.urban_rural_toronto_downsampled_bam_lists.output
    sites = expand(rules.split_angsd_sites_byChrom.output, chrom=CHROMOSOMES, site=['0fold', '4fold'])
    sites = [x for x in sites if 'CM019101.1' in x and '4fold' in x]
    sites_idx = expand(rules.angsd_index_degenerate.output.idx, chrom=CHROMOSOMES, site=['0fold', '4fold'])
    sites_idx = [x for x in sites if 'CM019101.1' in x and '4fold' in x]
    return {'bams' : bams, 'sites' : sites, 'sites_idx' : sites_idx}

def get_saf_joint_sfs_pi_fst_test(wildcards):
    if wildcards.group == 'highCov':
        out = expand(rules.angsd_saf_likelihood_toronto_pi_fst_test.output.saf_idx, group = 'highCov', habitat = ['urban', 'rural'], ss = '7')
    elif wildcards.group == 'lowCov':
        all_safs = expand(rules.angsd_saf_likelihood_toronto_pi_fst_test.output.saf_idx, group = 'lowCov', habitat = ['urban', 'rural'], ss = SS_PI_FST_TEST)
        sample_sizes = wildcards.joint_sfs_ss
        urban_ss = int(sample_sizes.split('_')[0][1])
        rural_ss = int(sample_sizes.split('_')[1][1])
        urban_saf = expand(rules.angsd_saf_likelihood_toronto_pi_fst_test.output.saf_idx, group = 'lowCov', habitat = 'urban', ss = urban_ss)
        rural_saf = expand(rules.angsd_saf_likelihood_toronto_pi_fst_test.output.saf_idx, group = 'lowCov', habitat = 'rural', ss = rural_ss)
        out = urban_saf + rural_saf
    return out

