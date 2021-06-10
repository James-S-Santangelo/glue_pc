# Script to estimate genotype likelihoods and allele frequencies at Ac and Li from read counts
# Writes two dataframes to disk per gene: (1) Dataframe with per-sample normalized genotype likelihoods
# (2) Dataframe with per-city frequency of deletion
#
# Author: Rob W. Ness
# Adapted by: James S. Santangelo

###############
#### SETUP ####
###############

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq, least_squares
from scipy.stats import nbinom

###################
#### FUNCTIONS ####
###################

def convert_params(mean, var):
    """
    Convert mean/dispersion parameterization of a negative binomial to the ones scipy supports

    See https://mathworld.wolfram.com/NegativeBinomialDistribution.html
    """
    p = mean/var
    n = mean*p/(1.0 - p)
    return n, p

def PMF_genotype(read_count, coverage_fraction, total_reads, region_size, genome_size, error_rate, dispersion_factor, scale_factor):
    '''
    Calculate expected number of samples with a given genotype, for a given read_count
    coverage_fraction=0, 0.5 or 1 for the different number of present alleles
    '''
    prob_hit = scale_factor * ((coverage_fraction*(region_size/genome_size)) + error_rate)
    mean = total_reads*prob_hit
    var = mean * dispersion_factor
    n,p=convert_params(mean, var)
    rv = nbinom(n, p)
    prob =rv.pmf(read_count)
    return prob

def fraction_from_read_counts(read_count, coeffs):
    #read_count = x
    # y is fraction of population with that read_count
    dispersion_factor, p_a, scale_factor, error_rate = coeffs
    p_aa = p_a**2
    p_Aa = 2*(p_a)*(1-p_a)
    p_AA = (1-p_a)**2
    expected_aas_w_this_read_count = (p_aa*sample_size) * PMF_genotype(read_count,0.0, total_reads, region_size, genome_size, error_rate, dispersion_factor, scale_factor)
    expected_Aas_w_this_read_count = (p_Aa*sample_size) * PMF_genotype(read_count,0.5, total_reads, region_size, genome_size, error_rate, dispersion_factor, scale_factor)
    expected_AAs_w_this_read_count = (p_AA*sample_size) * PMF_genotype(read_count,1.0, total_reads, region_size, genome_size, error_rate, dispersion_factor, scale_factor)
    return sum([expected_aas_w_this_read_count,expected_Aas_w_this_read_count,expected_AAs_w_this_read_count])
        
def residuals(coeffs, y_data, x_values):
    """
    coeffs - parameters to be estimated
    x_values = array of possible read_counts
    y_data. = observed counts of each read count
    """
    return y_data - fraction_from_read_counts(x_values, coeffs)

def p_aa(focal_read_count,total_reads, region_size, genome_size, error_rate, dispersion_factor, scale_factor):
    prob_hit= scale_factor* (0*(region_size/genome_size) + error_rate)
    mean=total_reads*prob_hit
    var = mean * dispersion_factor
    n,p=convert_params(mean, var)
    rv = nbinom(n, p)
    prob =rv.pmf(focal_read_count)
    return prob

def p_Aa(focal_read_count,total_reads, region_size, genome_size, error_rate, dispersion_factor, scale_factor):
    prob_hit= scale_factor* (0.5*(region_size/genome_size) + error_rate)
    mean=total_reads*prob_hit
    var = mean * dispersion_factor
    n,p=convert_params(mean, var)
    rv = nbinom(n, p)
    prob =rv.pmf(focal_read_count)
    return prob

def p_AA(focal_read_count,total_reads, region_size, genome_size, error_rate, dispersion_factor, scale_factor):
    prob_hit= scale_factor* (1.0*(region_size/genome_size) + error_rate)
    mean=total_reads*prob_hit
    var = mean * dispersion_factor
    n,p=convert_params(mean, var)
    rv = nbinom(n, p)
    prob =rv.pmf(focal_read_count)
    return prob

def gt_likelihoods(focal_read_count,total_reads, region_size, genome_size, error_rate, dispersion_factor, scale_factor):
    aa=p_aa(focal_read_count,total_reads, region_size, genome_size, error_rate, dispersion_factor, scale_factor)
    Aa=p_Aa(focal_read_count,total_reads, region_size, genome_size, error_rate, dispersion_factor, scale_factor)
    AA=p_AA(focal_read_count,total_reads, region_size, genome_size, error_rate, dispersion_factor, scale_factor)
    return [aa,Aa,AA]

##################
#### ANALYSIS ####
##################

with open(snakemake.log[0], 'w') as f:
    sys.stderr = sys.stdout = f

    try:
        # Create dictionary with info for all samples 
        sample_info = {}
        for l in open(snakemake.input.counts[0], 'r'):
            b, c = l.strip().split()
            s = b.split("/")[-1].split("_merged")[0]
            city = "_".join(s.split("_")[:-2])
            population = city+ "_" + s.split("_")[-2]
            sample_info[s] = {'bam' : b, 'read_count' : int(c), 'city' : city, 'population' : population}
            
        # Path to multiQC files with mapping statistics
        qual_datafile = '../results/qc/multiqc/multiqc_data/multiqc_bamtools_stats_bamtools_stats.txt'

        # Add number of mapped and duplicate reads to sample dictionary
        for l in open(qual_datafile).readlines()[1:]:
            s = l.split("\t")[0].split("|")[-1].strip()
            mapped_reads = float(l.split("\t")[2])
            duplicates = float(l.split("\t")[10])
            try:
                sample_info[s]['mapped_reads'] = mapped_reads
                sample_info[s]['duplicates'] = duplicates
            except KeyError:
                print(s)
        
        # Create dataframe from sample dictionary
        df = pd.DataFrame.from_dict(sample_info, orient = 'index')
        df['read_count_normalized'] = df.read_count/(df.mapped_reads-df.duplicates)
        
        # Define constants
        # TODO Automate region size calculation from Snakemake parameters 
        total_reads = df.mapped_reads.mean()
        if snakemake.wildcards.gene == 'li':
            region_size = (30229250 - 30218214 + 1) + (30247247 - 30230911 + 1)
        else:
            region_size = 19573344 - 19559221 + 1
        genome_size = 1e9 * 0.75  # Only 75% of genome is mappable
        sample_size = df.shape[0]

        # Define coefficients and optimize fit by least squares
        coeffs0 = np.array([15, 0.65, 2, 2.5e-6], dtype='float')
        coeff_bounds = ([1,0,0,0], [np.inf,1, np.inf,1])
        largest_read_count = int(round(max(df.mapped_reads.mean() * df.read_count_normalized)))
        x_values = np.array(range(largest_read_count), dtype="float")
        y_data, bins = plt.hist(df.mapped_reads.mean() * df.read_count_normalized, bins = range(0, largest_read_count + 1))[:2]
        fitted_coeffs = least_squares(residuals, coeffs0, args=(y_data, x_values), bounds=coeff_bounds)
        
        # Fitted parameters
        dispersion_factor, p_a, scale_factor, error_rate = fitted_coeffs.x
        
        # Estimate genotype likelihoods from standardized read counts
        df['read_count_standardized'] = df.read_count_normalized * total_reads
        df['read_count_standardized'] = df['read_count_standardized'].astype(int)
        df['l_aa'] = p_aa(df.read_count_standardized, total_reads, region_size, genome_size, error_rate, dispersion_factor, scale_factor)
        df['l_Aa'] = p_Aa(df.read_count_standardized, total_reads, region_size, genome_size, error_rate, dispersion_factor, scale_factor)
        df['l_AA'] = p_AA(df.read_count_standardized, total_reads, region_size, genome_size, error_rate, dispersion_factor, scale_factor)
        
        # Normalize genotype likelihoods so they sum to 1 
        df["l_aa_norm"] = df.l_aa / sum([df.l_aa,  df.l_Aa, df.l_AA])
        df["l_Aa_norm"] = df.l_Aa / sum([df.l_aa,  df.l_Aa, df.l_AA])
        df["l_AA_norm"] = df.l_AA / sum([df.l_aa,  df.l_Aa, df.l_AA]) 
        
        # Write genotype likelihoods to disk
        df.to_csv(snakemake.output.likes, sep = '\t', index_label = False) 
       
        # Calculate allele frequencies from genotpe likelihoods
        with open(snakemake.output.freqs, 'w') as freqs_out:
            freqs_out.write('city\tp\n')  # Header
            for city in df.city.unique():
                num_aa = sum(df[df.city==city].l_aa_norm)
                num_Aa = sum(df[df.city==city].l_Aa_norm)
                num_AA = sum(df[df.city==city].l_AA_norm)
                p = (num_aa + (0.5 * num_Aa)) / (num_aa + num_Aa + num_AA)
                if city == 's': city = 'Toronto'
                freqs_out.write('{0}\t{1}\n'.format(city, p))

    except Exception as e:
        print(e) 
