# # Import Eucalyptus data to scikit-allel
#
# This notebook does the basic import of the acanthophis-derived varants for
# Helen's Decra sample set into zarr for use with scikit allel.


import numpy as np
import allel
import zarr
import sys



# First step: set up the fields we wish to import. We want most bits of the
# FORMAT (per-genotype) fileds, and a few relevant columns from the INFO/table
# (per-variant) part of the VCF.
#
# CHROM/POS/REF/ALT/QUAL are the vcf columns, variants/AD & DP are the INFO
# fields, per variant depth and allelic depth. calldata/GT is genotype, DP is
# per-genotype depth, AD is per allle read depth, PL is phread-scaled
# likelihood.

fields=['samples', 'variants/CHROM', 'variants/POS', 'variants/REF',
        'variants/ALT', 'variants/AD', 'variants/QUAL', 'calldata/GT',
        'calldata/AD']


# Now convert the vcf to a zarr array. We save each reference's calls in their
# own group (e.g. 'egrandis'), so we can have all three in the one zarr.
# 
# This takes ages, and I do it in parallel with concurrent.futures


from concurrent.futures import ProcessPoolExecutor, as_completed


with ProcessPoolExecutor(3) as exc:
    jobs = []

    jobs.append(exc.submit(
        allel.vcf_to_zarr,
        "data/variants-acanthophis/enhanced_filter/mpileup~bwa~Egrandis_phytozome13_v2~enhanced_filter.vcf.gz",
        "data/variants-acanthophis/mpileup~bwa~HBDecra.zarr",
        group="egrandis",
        compressor=zarr.Blosc(cname="zstd", clevel=5, shuffle=1),
        fields=fields,
        chunk_width=1000,
        chunk_length=100000,
        #log=sys.stderr,
    ))

    jobs.append(exc.submit(
        allel.vcf_to_zarr,
        "data/variants-acanthophis/enhanced_filter/mpileup~bwa~Emelliodora_GCA_004368105.2~enhanced_filter.vcf.gz",
        "data/variants-acanthophis/mpileup~bwa~HBDecra.zarr",
        group="emelliodora",
        compressor=zarr.Blosc(cname="zstd", clevel=5, shuffle=1),
        fields=fields,
        chunk_width=1000,
        chunk_length=100000,
        #log=sys.stderr,
    ))

    jobs.append(exc.submit(
        allel.vcf_to_zarr,
        "data/variants-acanthophis/enhanced_filter/mpileup~bwa~Esideroxylon_GCA_014182405.1~enhanced_filter.vcf.gz",
        "data/variants-acanthophis/mpileup~bwa~HBDecra.zarr",
        group="esideroxylon",
        compressor=zarr.Blosc(cname="zstd", clevel=5, shuffle=1),
        fields=fields,
        chunk_width=1000,
        chunk_length=100000,
        #log=sys.stderr,
    ))

    for j in as_completed(jobs):
        print("Done one")
