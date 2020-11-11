# Script to convert SNP VCF to Zarr database

import argparse
import sys
import os
from itertools import repeat

# Temporaty fix for know issue while importing allel
# https://github.com/cggh/scikit-allel/issues/285
os.environ["NUMEXPR_MAX_THREADS"]="272"

import pandas as pd
import multiprocessing as mp
 
import zarr
import numcodecs
import allel

def chroms_to_list(chromosome_file):
    """Converts text file with chromosome names to list

    Args:
        chromosomes (str): Path to text file with chromsomes

    Returns: 
        chrom_list (:obj:`list` of :obj:`str`): List with
        chromosomes names as elements.
    """

    chrom_list = pd.read_table(chromosome_file, header=None).iloc[:,0].tolist()
    
    return chrom_list


def vcf_to_zarr(vcf_in, outpath, chrom):
    """Convert on-disk VCF to on-disk Zarr database using
    scikit-allele and zarr modules

    Zarr database written to same directory as input VCF
    
    Args:
        vcf_in (str): Path to input VCF on disk
        outpath (str): Path to derectory to output Zarr db
        chrom (str): Chromosome for which Zarr database should be created

    Returns:
        None
    """
    # Rename 'numalt' field. Required by Zarr to distinguish `NUMALT` from `numalt`
    # `numalt` is automatically computed by scikit-allel
    rename_dict = {'variants/numalt':'variants/numalt_sci'}

    # Use vcf_to_zarr function from scikit-allel to create zarr database
    # Currently optimized for biallelic SNP VCF but easy to extend functionality
    allel.vcf_to_zarr(
            input=vcf_in,
            output=outpath,
            overwrite=True,
            group=chrom,
            rename_fields=rename_dict,
            fields='*',
            alt_number=1,
            region=chrom,
            compressor=numcodecs.Blosc(cname='zstd', clevel=1, shuffle=False)
            )


def main():

    vcf_in = snakemake.input[0]
    outpath = snakemake.output[0]
    chrom_path = str(snakemake.config['chromosomes'])
    
    chrom_list = chroms_to_list(chrom_path)

    iterator = zip(repeat(vcf_in), repeat(outpath), chrom_list)

    with mp.Pool(processes = len(chrom_list)) as pool:
        pool.starmap(vcf_to_zarr, iterator)
   

with open(snakemake.log[0], 'w') as f:
    sys.stderr = sys.stdout = f
    main()
