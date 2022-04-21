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

fields=['samples', 'variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT',
        'variants/AD', 'variants/QUAL', 'calldata/GT', 'calldata/AD']


# Now convert the vcf to a zarr array. We save each reference's calls in their
# own group (e.g. 'egrandis'), so we can have all three in the one zarr.
# 
# This takes ages.


allel.vcf_to_zarr(
    "data/variants-acanthophis/enhanced_filter/mpileup~bwa~Egrandis_phytozome13_v2~enhanced_filter.vcf.gz",
    "data/variants-acanthophis/smalltest.zarr",
    group="egrandis",
    compressor=zarr.Blosc(cname="zstd", clevel=5, shuffle=1),
    fields=fields,
    chunk_width=1000,
    chunk_length=100000,
    log=sys.stderr,
    region='Chr01:1-20000000',
    overwrite=True,
)


za = zarr.open("data/variants-acanthophis/smalltest.zarr", mode='r')
za.tree()
egr = za["egrandis"]
egr.tree()

egr_idx = allel.SortedMultiIndex(egr["variants/CHROM"], egr["variants/POS"])
egr_snp = allel.GenotypeChunkedArray(egr["calldata/GT"])

egr_snp

egr_gn = egr_snp.to_n_alt()
dist = allel.pairwise_distance(egr_gn, metric='euclidean')
