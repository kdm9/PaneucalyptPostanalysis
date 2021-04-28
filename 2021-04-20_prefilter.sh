# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     text_representation:
#       extension: .sh
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Bash
#     language: bash
#     name: bash
# ---

# # Pre-filter variants
#
# Here we pre-filter the new acanthophis variants with some slightly more
# enhanced filtering than is done in Acanthophis. This is motivated by the
# fairly basic fact that these raw variant calls are *enormous*: nearly a
# terabyte each. This process should throw out the long tail of either very
# rare or erroneous variants, while keeping nearly all variants useful for any
# purpose. For nearly all forms of subsequent inference one would want to
# filter again in a more application specific fashion, but that should be
# possible from these variant files without referencing the massive bcfs again.
#
# To minimise repetion, here are the filtering steps as a bash function:

function filter_variants() {
    input="$1"
    output="$2"

    pv -n -i 20 --name $(basename "$input") "$input" \
    | bcftools filter                       \
        --threads ${PBS_NCPUS:-4}           \
        -Ou                                 \
        --IndelGap 10                       \
        --SnpGap 10                         \
        -e 'F_MISSING > 0.8 | QUAL < 100'   \
    | bcftools norm -m -snps                \
        --do-not-normalize                  \
        --threads ${PBS_NCPUS:-4}           \
        -Ou                                 \
    | bcftools filter                       \
        --threads ${PBS_NCPUS:-4}           \
        -Ou                                 \
        -e "INFO/AD < 20 | INFO/AC < 15"    \
    | bcftools norm                         \
        --threads ${PBS_NCPUS:-4}           \
        --do-not-normalize                  \
        -m +snps                            \
        -Oz                                 \
    | tee "${output}.tmp"                   \
    | bcftools index                        \
        --force                             \
        --csi                               \
        -o "${output}.csi.tmp"

    mv "${output}.tm"p "${output}"
    mv "${output}.csi.tmp" "${output}.csi"
}


# In short, we:
#
# - Filter variants within 10bp of an indel
# - remove variants with either >80% missing data, or quality < 100 
#   ($p > 10^{-10}$)
# - Split multialleic sites into multiple bialleic sites
# - Filter each allele independently, removing alleles with fewer than 20
#   supporting reads, and fewer than 15 calls (2 for hom, 1 for het) across all
#   samples
# - Re-join any remaining split mulitalleic sites into single vcf entries
# - On-the-fly, create an index and write out both the index and vcf.gz


# # Actual execution
#
# The above was just a bash function, now we need to actually run this. Here we
# run it for Helen's decra samples aligned against the three genomes
# (mel/sid/grandis).


mkdir -p data/variants-acanthophis

export -f filter_variants

# A parallel loop that operates over the three genomes
parallel --lb filter_variants \
        data/variants-acanthophis/orig/mpileup~bwa~{}~HBDecra~filtered-default.vcf.gz \
        data/variants-acanthophis/enhanced_filter/mpileup~bwa~{}~enhanced_filter.vcf.gz \
    ::: Esideroxylon_GCA_014182405.1  Egrandis_phytozome13_v2 Emelliodora_GCA_004368105.2
