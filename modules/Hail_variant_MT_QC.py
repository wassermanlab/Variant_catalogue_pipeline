#!/usr/bin/env python
# coding: utf-8

# In[47]:


pon_predictions_table='/mnt/common/SILENT/Act3/MT_references/pon_mt_trna_predictions_08_27_2020.txt'
artifact_prone_sites_bed = '/mnt/common/SILENT/Act3/MT_references/artifact_prone_sites.bed'
GRCh38_MT_local_fasta='/mnt/common/SILENT/Act3/MT_references/Homo_sapiens_assembly38.chrM.fasta'
GRCh38_MT_local_fai='/mnt/common/SILENT/Act3/MT_references/Homo_sapiens_assembly38.chrM.fasta.fai'
mitotip_predictions_table='/mnt/common/SILENT/Act3/MT_references/mitotip_scores_08_27_2020.txt'

#Created through the nextflow pipeline
import sys

MT_Step1_input_tsv=sys.argv[1]
MT_Step2_participant_data=sys.argv[2]
MT_participants_to_subset =sys.argv[3]
MT_Step3_participant_data=sys.argv[4]

#Created through the nextflow pipeline
#MT_Step1_input_tsv='MT_Step1_input_tsv.tsv'
#MT_Step2_participant_data= 'MT_Step2_participant_data.tsv'
#MT_participants_to_subset = 'MT_participants_to_subset.txt'
#MT_Step3_participant_data= 'work/2c/aa6a7391f6b2ee9ce2bb9b85be0e3b/MT_Step3_participant_data.tsv'

#In the step3 participant data, the median wgs coverage is replaced by the mean wgs coverage


# Hail and plot initialisation 

# In[2]:


import hail as hl
from hail.plot import output_notebook, show
from hail.utils.java import info
from typing import Dict
from hail.genetics import ReferenceGenome

hl.init()
output_notebook()


# In[3]:


from gnomad.utils.annotations import age_hists_expr
from gnomad.utils.reference_genome import add_reference_sequence
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import vep_struct_to_csq
#from gnomad_qc.v3.resources.meta import meta  # pylint: disable=import-error
from gnomad.resources.grch38.gnomad import POPS
from gnomad.resources.grch38.reference_data import dbsnp, _import_dbsnp


# In[4]:


import argparse
import logging
import math
import os

from pprint import pprint
from bokeh.models import Span
hl.plot.output_notebook()
from bokeh.models import Range1d
from bokeh.plotting import figure, output_file, show, save

import pandas as pd
from typing import Tuple
import string

import re

from collections import Counter
from textwrap import dedent

from os.path import dirname

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("Annotate coverage")
logger.setLevel(logging.INFO)


# gnomAD pipeline:
# 
# https://github.com/broadinstitute/gnomad-mitochondria
# 
# some are available here : https://github.com/broadinstitute/gnomad-mitochondria/tree/main/gnomad_mitochondria/resources
# 
# In gnomAD, split in several steps and several python script, here one python script / step
# 
# 1. Overview of the pipeline
# 2. Functions
# 3. Main scripts
# 
# Step 1: annotate_coverage.py
# 
# combine the per base coverage files across many samples into a mt(MatrixTable), ht(HailTable), and tsv file, and will also calculate the following aggregate statistics per base:
# 
# - Mean coverage
# - Median coverage
# - Fraction of samples with > 100x coverage
# - Fraction of samples with > 10000x coverage
# 
# Step 2:combine_vcfs.py
# 
# takes individual sample vcfs and combines them into one vcf/mt
# 
# Step 3:add_annotations.py
# 
# Annotations to the mt/vcf are added

# **Functions**

# In[5]:


def multi_way_union_mts(mts: list, temp_dir: str, chunk_size: int) -> hl.MatrixTable:
    """
    Hierarchically join together MatrixTables in the provided list.
    :param mts: List of MatrixTables to join together
    :param temp_dir: Path to temporary directory for intermediate results
    :param chunk_size: Number of MatrixTables to join per chunk (the number of individual VCFs that should be combined at a time)
    :return: Joined MatrixTable
    """
    # Convert the MatrixTables to tables where entries are an array of structs
    staging = [mt.localize_entries("__entries", "__cols") for mt in mts]
    stage = 0
    while len(staging) > 1:
        # Calculate the number of jobs to run based on the chunk size
        n_jobs = int(math.ceil(len(staging) / chunk_size))
        info(f"multi_way_union_mts: stage {stage}: {n_jobs} total jobs")
        next_stage = []

        for i in range(n_jobs):
            # Grab just the tables for the given job
            to_merge = staging[chunk_size * i : chunk_size * (i + 1)]
            info(
                f"multi_way_union_mts: stage {stage} / job {i}: merging {len(to_merge)} inputs"
            )

            # Multiway zip join will produce an __entries annotation, which is an array where each element is a struct containing the __entries annotation (array of structs) for that sample
            merged = hl.Table.multi_way_zip_join(to_merge, "__entries", "__cols")
            # Flatten __entries while taking into account different entry lengths at different samples/variants (samples lacking a variant will be NA)
            merged = merged.annotate(
                __entries=hl.flatten(
                    hl.range(hl.len(merged.__entries)).map(
                        # Coalesce will return the first non-missing argument, so if the entry info is not missing, use that info, but if it is missing, create an entries struct with the correct element type for each null entry annotation (such as int32 for DP)
                        lambda i: hl.coalesce(
                            merged.__entries[i].__entries,
                            hl.range(hl.len(merged.__cols[i].__cols)).map(
                                lambda j: hl.null(
                                    merged.__entries.__entries.dtype.element_type.element_type
                                )
                            ),
                        )
                    )
                )
            )

            # Flatten col annotation from array<struct{__cols: array<struct{s: str}>} to array<struct{s: str}>
            merged = merged.annotate_globals(
                __cols=hl.flatten(merged.__cols.map(lambda x: x.__cols))
            )

            next_stage.append(
                merged.checkpoint(
                    os.path.join(temp_dir, f"stage_{stage}_job_{i}.ht"), overwrite=True
                )
            )
        info(f"Completed stage {stage}")
        stage += 1
        staging.clear()
        staging.extend(next_stage)

    # Unlocalize the entries, and unfilter the filtered entries and populate fields with missing values
    return (
        staging[0]
        ._unlocalize_entries("__entries", "__cols", list(mts[0].col_key))
        .unfilter_entries()
    )


# In[6]:


def collect_vcf_paths(
    participant_data: str, vcf_col_name: str, participants_to_subset: str = None,
) -> Dict[str, str]:
    """
    Create dictionary of VCF paths for only the samples specified in participants_to_subset.
    .. note::
        Participant data should be a tab-delimited file with at minimum columns for:
        - 'entity:participant_id': sample name with prohibited characters replaced with underscores
        - 's': sample name
        - path to the Mutect2 VCF output, where name of this column is supplied to the `vcf_col_name` parameter
    :param participant_data: Participant data (the downloaded data tab from Terra)
    :param vcf_col_name: Name of column that contains VCF output
    :param participants_to_subset: Path to file of participant_ids to which the data should be subset
    :return: Dictionary with sample name as key and path to VCF as value
    """
    vcf_paths = {}
    # Load in data from Terra
    participant_ht = hl.import_table(participant_data)

    # Remove participants that don't have VCF output
    participant_ht.filter(participant_ht[vcf_col_name] != "")

    # Subset participants if specified
    if participants_to_subset:
        participants_of_interest = hl.import_table(
            participants_to_subset
        ).participant.collect()
        participant_ht = participant_ht.filter(
            hl.literal(participants_of_interest).contains(
                participant_ht["entity:participant_id"]
            )
        )

    # Add the vcf path to a dictionary with sample name as key
    df = participant_ht.to_pandas()

    for _, row in df.iterrows():
        vcf_paths[row["s"]] = row[vcf_col_name]

    return vcf_paths


# In[7]:


def multi_way_union_mts(mts: list, temp_dir: str, chunk_size: int) -> hl.MatrixTable:
    """
    Hierarchically join together MatrixTables in the provided list.
    :param mts: List of MatrixTables to join together
    :param temp_dir: Path to temporary directory for intermediate results
    :param chunk_size: Number of MatrixTables to join per chunk (the number of individual VCFs that should be combined at a time)
    :return: Joined MatrixTable
    """
    # Convert the MatrixTables to tables where entries are an array of structs
    staging = [mt.localize_entries("__entries", "__cols") for mt in mts]
    stage = 0
    while len(staging) > 1:
        # Calculate the number of jobs to run based on the chunk size
        n_jobs = int(math.ceil(len(staging) / chunk_size))
        info(f"multi_way_union_mts: stage {stage}: {n_jobs} total jobs")
        next_stage = []

        for i in range(n_jobs):
            # Grab just the tables for the given job
            to_merge = staging[chunk_size * i : chunk_size * (i + 1)]
            info(
                f"multi_way_union_mts: stage {stage} / job {i}: merging {len(to_merge)} inputs"
            )

            # Multiway zip join will produce an __entries annotation, which is an array where each element is a struct containing the __entries annotation (array of structs) for that sample
            merged = hl.Table.multi_way_zip_join(to_merge, "__entries", "__cols")
            # Flatten __entries while taking into account different entry lengths at different samples/variants (samples lacking a variant will be NA)
            merged = merged.annotate(
                __entries=hl.flatten(
                    hl.range(hl.len(merged.__entries)).map(
                        # Coalesce will return the first non-missing argument, so if the entry info is not missing, use that info, but if it is missing, create an entries struct with the correct element type for each null entry annotation (such as int32 for DP)
                        lambda i: hl.coalesce(
                            merged.__entries[i].__entries,
                            hl.range(hl.len(merged.__cols[i].__cols)).map(
                                lambda j: hl.null(
                                    merged.__entries.__entries.dtype.element_type.element_type
                                )
                            ),
                        )
                    )
                )
            )

            # Flatten col annotation from array<struct{__cols: array<struct{s: str}>} to array<struct{s: str}>
            merged = merged.annotate_globals(
                __cols=hl.flatten(merged.__cols.map(lambda x: x.__cols))
            )

            next_stage.append(
                merged.checkpoint(
                    os.path.join(temp_dir, f"stage_{stage}_job_{i}.ht"), overwrite=True
                )
            )
        info(f"Completed stage {stage}")
        stage += 1
        staging.clear()
        staging.extend(next_stage)

    # Unlocalize the entries, and unfilter the filtered entries and populate fields with missing values
    return (
        staging[0]
        ._unlocalize_entries("__entries", "__cols", list(mts[0].col_key))
        .unfilter_entries()
    )


# In[8]:


def join_mitochondria_vcfs_into_mt(
    vcf_paths: Dict[str, str], temp_dir: str, chunk_size: int = 100
) -> hl.MatrixTable:
    """
    Reformat and join individual mitochondrial VCFs into one MatrixTable.
    :param vcf_paths: Dictionary of samples to combine (sample as key, path to VCF as value)
    :param temp_dir: Path to temporary directory for intermediate results
    :param chunk_size: Number of MatrixTables to join per chunk (the number of individual VCFs that should be combined at a time)
    :return: Joined MatrixTable of samples given in vcf_paths dictionary
    """
    mt_list = []
    for sample, vcf_path in vcf_paths.items():
        try:
            mt = hl.import_vcf(vcf_path, reference_genome="GRCh38", force_bgz = True)
        except Exception as e:
            raise ValueError(
                f"vcf path {vcf_path} does not exist for sample {sample}"
            ) from e

        # Because the vcfs are split, there is only one AF value, although misinterpreted as an array because Number=A in VCF header
        # Second value of MMQ is the value of the mapping quality for the alternate allele
        # Add FT annotation for sample genotype filters (pull these from filters annotations of the single-sample VCFs)
        mt = mt.select_entries("DP", HL=mt.AF[0])
        mt = mt.annotate_entries(
            MQ=hl.float(mt.info["MMQ"][1]),
            TLOD=mt.info["TLOD"][0],
            FT=hl.if_else(hl.len(mt.filters) == 0, {"PASS"}, mt.filters),
        )
        # Use GRCh37 reference as most external resources added in downstream scripts use GRCh37 contig names
        # (although note that the actual sequences of the mitochondria in both GRCh37 and GRCh38 are the same)
        mt = mt.key_rows_by(
            locus=hl.locus("MT", mt.locus.position, reference_genome="GRCh37"),
            alleles=mt.alleles,
        )
        mt = mt.key_cols_by(s=sample)
        mt = mt.select_rows()
        mt_list.append(mt)

    combined_mt = multi_way_union_mts(mt_list, temp_dir, chunk_size)
    
    return combined_mt


# In[9]:


def remove_genotype_filters(
    mt: hl.MatrixTable,
    filters_to_remove: set = {
        "possible_numt",
        "mt_many_low_hets",
        "FAIL",
        "blacklisted_site",
    },
) -> hl.MatrixTable:
    """
    Remove unneeded sample-level genotype filters (in FT field of the VCF) specified by the filters_to_remove parameter.
    By default, remove the 'possible_numt', 'mt_many_low_hets', and 'FAIL' filters because these filters were found to have low performance.
    Also remove the 'blacklisted_site' filter because this filter did not always behave as expected in early GATK versions. This filter can be reimplemented with the apply_mito_artifact_filter function.
    :param mt: MatrixTable containing genotype filters in the FT field of the VCF that should be removed
    :param filters_to_remove: List of genptype filters (in FT field of VCF) that should be removed from the entries
    :return: MatrixTable with specific genotype filters (in FT field of VCF) removed
    """
    mt = mt.annotate_entries(FT=mt.FT.difference(filters_to_remove))

    # If no filters exist after removing those specified above, set the FT field to PASS
    mt = mt.annotate_entries(FT=hl.if_else(hl.len(mt.FT) == 0, {"PASS"}, mt.FT))
    
    return mt


# In[10]:


def determine_hom_refs(
    mt: hl.MatrixTable, coverage_mt_path: str, minimum_homref_coverage: int = 100
) -> hl.MatrixTable:
    """
    Use coverage to distinguish between homref and missing sites.
    :param mt: MatrixTable from initial multi-sample merging, without homref sites determined
    :param coverage_mt_path: MatrixTable of sample level coverage at each position (per-sample and per-base; can be generated by running annotate_coverage.py)
    :param minimum_homref_coverage: Minimum depth of coverage required to call a genotype homoplasmic reference rather than missing
    :return: MatrixTable with missing genotypes converted to homref depending on coverage
    """
    # Convert coverage to build GRCh37 to match contig names
    # Note: the mitochondrial reference genome is the same for GRCh38 and GRCh37
    coverages = hl.read_matrix_table(coverage_mt_path)
    coverages = coverages.key_rows_by(
        locus=hl.locus("MT", coverages.locus.position, reference_genome="GRCh37")
    )

    mt = mt.annotate_entries(
        DP=hl.if_else(hl.is_missing(mt.HL), coverages[mt.locus, mt.s].coverage, mt.DP)
    )

    hom_ref_expr = hl.is_missing(mt.HL) & (mt.DP > minimum_homref_coverage)

    mt = mt.annotate_entries(
        HL=hl.if_else(hom_ref_expr, 0.0, mt.HL),
        FT=hl.if_else(hom_ref_expr, {"PASS"}, mt.FT),
        DP=hl.if_else(
            hl.is_missing(mt.HL) & (mt.DP <= minimum_homref_coverage),
            hl.null(hl.tint32),
            mt.DP,
        ),
    )
  
    return mt


# In[11]:


def apply_mito_artifact_filter(
    mt: hl.MatrixTable, artifact_prone_sites_path: str,
) -> hl.MatrixTable:
    """
    Add in artifact_prone_site filter.
    :param mt: MatrixTable to be annotated with artifact_prone_sites filter
    :param artifact_prone_sites_path: Path to BED file of artifact_prone_sites to flag in the filters column
    :return: MatrixTable with artifact_prone_sites filter
    """
    # Apply "artifact_prone_site" filter to any SNP or deletion that spans a known problematic site
    bed = hl.import_bed(artifact_prone_sites_path)
    bed = bed.annotate(target="artifact")


    # Create a region annotation containing the interval that the variant overlaps (for SNP will be one position, but will be longer for deletions based on the length of the deletion)
    mt = mt.annotate_rows(
        region=hl.interval(
            hl.locus("MT", mt.locus.position, reference_genome="GRCh37"),
            hl.locus(
                "MT",
                mt.locus.position + hl.len(mt.alleles[0]) - 1,
                reference_genome="GRCh37",
            ),
            includes_end=True,
        )
    )

    
    # Annotate if the start of the variant overlaps an interval in the bed file
    mt = mt.annotate_rows(start_overlaps=bed.index(mt.region.start, all_matches=True))

    # Annotate if the end of the variant overlaps an interval in the bed file
    mt = mt.annotate_rows(end_overlaps=bed.index(mt.region.end, all_matches=True))

    # Create struct containing locus and allele (need to the check if any position of the allele overlaps an artifact-prone site, not just the locus)
    mt_temp = mt.annotate_rows(variant=hl.struct(locus=mt.locus, alleles=mt.alleles))
    mt_temp = mt_temp.key_rows_by(mt_temp.region)

    # Need to account for cases where the start and end of the variant interval don't fall within a bed interval, but start before and after the interval (the bed interval falls completely within the variant interval)
    bed_temp = bed.annotate(
        contained_mt_alleles=mt_temp.index_rows(
            bed.interval.start, all_matches=True
        ).variant
    )

    # Explode so that each allele is on its own row and create locus and allele annotations
    bed_temp = bed_temp.explode(bed_temp.contained_mt_alleles).rename(
        {"contained_mt_alleles": "contained_mt_allele"}
    )
    bed_temp = bed_temp.annotate(
        locus=bed_temp.contained_mt_allele.locus,
        alleles=bed_temp.contained_mt_allele.alleles,
    )
    bed_temp = bed_temp.key_by(bed_temp.locus, bed_temp.alleles)

    # Annotate back onto the original mt cases where the bed interval falls completely within the variant interval
    mt = mt.annotate_rows(start_and_end_span=bed_temp[mt.locus, mt.alleles].target)

    # Add artifact-prone site filter to any SNP/deletion that starts within, ends within, or completely overlaps an artifact-prone site
    
    mt = mt.annotate_rows(
        filters=hl.if_else(
            (hl.len(mt.start_overlaps) > 0)
            | (hl.len(mt.end_overlaps) > 0)
            | (hl.is_defined(mt.start_and_end_span)),
            {"artifact_prone_site"},
            {"PASS"},
        )
    )

    mt = mt.drop("start_overlaps", "end_overlaps", "start_and_end_span") 
    
    return mt


# In[12]:


def add_descriptions(
    input_mt: hl.MatrixTable,
    min_hom_threshold: float = 0.95,
    vaf_filter_threshold: float = 0.01,
    min_het_threshold: float = 0.10,
) -> hl.MatrixTable:
    """
    Add descriptions of annotations to globals.
    :param input_mt: MatrixTable
    :param min_hom_threshold: minimum cutoff to define a variant as homoplasmic
    :param min_het_threshold: minimum cutoff to define a variant as a PASS heteroplasmic variant, genotypes below this threshold will count towards the below_min_het_threshold
    :param vaf_filter_threshold: should match vaf_filter_threshold supplied to Mutect2, variants below this value will be set to homoplasmic reference after calculating the common_low_heteroplasmy filter    :return: MatrixTable with descriptions of annotations stored in  globals
    :rtype: MatrixTable
    """
    # pull out variables that are needed in description dictionaries
    #hap_order = hl.eval(input_mt.hap_order)
    #population_order = hl.eval(input_mt.pop_order)
    hl_hist_bin_edges = input_mt.hl_hist.bin_edges.take(1)[0]

    global_annotation_dict = hl.struct(
        vep_version=hl.struct(Description="VEP version"),
        dbsnp_version=hl.struct(Description="Version of dbSNP used"),
        hap_order=hl.struct(
            Description="The order in which haplogroups are reported for haplogroup-related annotations"
        ),
        dp_hist_all_variants_bin_freq=hl.struct(
            Description="Histogram values for depth (DP) across all variants"
        ),
        dp_hist_all_variants_n_larger=hl.struct(
            Description="Number of values greater than largest histogram bin for depth (DP) across all variants"
        ),
        dp_hist_all_variants_bin_edges=hl.struct(
            Description="Values of bin edges for histogram of depth (DP) across all variants"
        ),
        mq_hist_all_variants_bin_freq=hl.struct(
            Description="Histogram values for mapping quality (MQ) across all variants"
        ),
        mq_hist_all_variants_n_larger=hl.struct(
            Description="Number of values greater than largest histogram bin for mapping quality (MQ) across all variants"
        ),
        mq_hist_all_variants_bin_edges=hl.struct(
            Description="Values of bin edges for histogram of mapping quality (MQ) across all variants"
        ),
        tlod_hist_all_variants_bin_freq=hl.struct(
            Description="Histogram values for tumor log odds (TLOD) across all variants"
        ),
        tlod_hist_all_variants_n_larger=hl.struct(
            Description="Number of values greater than largest histogram bin for tumor log odds (TLOD) across all variants"
        ),
        tlod_hist_all_variants_bin_edges=hl.struct(
            Description="Values of bin edges for histogram of tumor log odds (TLOD) across all variants"
        ),
    )

    col_annotation_dict = hl.struct(
        s=hl.struct(Description="Sample ID"),
        participant_id=hl.struct(
            Description="Participant ID which is used on the terra platform, not always the same as sample ID"
        ),
        contamination=hl.struct(
            Description="Estimate of mitochondrial contamination that is output by Haplocheck, reported as a proportion"
        ),
        freemix_percentage=hl.struct(
            Description="Estimate of nuclear contamination that is output by verifyBamID, reported as a percentage"
        ),
        major_haplogroup=hl.struct(
            Description="The major haplogroup that is output by Haplogrep"
        ),
        wgs_mean_coverage=hl.struct(
            Description="The mean depth of coverage on the autosomes for whole genome sequencing data which is output by Picard’s CollectWgsMetrics tool"
        ),
        mt_mean_coverage=hl.struct(
            Description="Mean mitochondrial depth of coverage that is output by MuTect2"
        ),
        release=hl.struct(Description="True if the sample is in the gnomAD release"),
        hard_filters=hl.struct(
            Description="A set containing any of the hard filters that were applied to the sample for QC on the nuclear genome"
        ),
        research_project=hl.struct(
            Description="Description of the research project to which the sample belongs"
        ),
        project_id=hl.struct(
            Description="The Project ID for the sample, typically a RP or G project for internal Broad samples or a short description for external samples"
        ),
        product=hl.struct(Description="Sequencing platform"),
        sample_pi=hl.struct(
            Description="Principal investigator of the research project of the sample"
        ),
        sex_karyotype=hl.struct(
            Description="Sample’s sex karyotype (combined X and Y karyotype)"
        ),
        broad_external=hl.struct(
            Description='Whether the sample was sequenced internally ("broad") or externally ("external")'
        ),
        mito_cn=hl.struct(
            Description="The estimated mitochondrial copy number of the sample, calculated as 2*mt_mean_coverage/wgs_mean_coverage"
        ),
        over_85_mean=hl.struct(
            Description="Mean heteroplasmy level (restricted to heteroplasmy levels >= 0.85) for the sample for PASS variants that are haplogroup-defining variants and are not at an artifact-prone site"
        ),
        over_85_count=hl.struct(
            Description="Number of variants with heteroplasmy levels >= 0.85 for the sample for PASS variants that are haplogroup-defining variants and are not at an artifact-prone site"
        ),
        bt_85_and_99_mean=hl.struct(
            Description="Mean heteroplasmy level (restricted to heteroplasmy levels >= 0.85 and <= 0.998) for the sample for PASS variants that are haplogroup-defining variants and are not at an artifact-prone site"
        ),
        bt_85_and_99_count=hl.struct(
            Description="Number of variants with heteroplasmy levels >= 0.85 and <= 0.998 for the sample for PASS variants that are haplogroup-defining variants and are not at an artifact-prone site"
        ),
        contam_high_het=hl.struct(
            Description="Internal estimate of contamination. It is defined for each sample as one minus the mean heteroplasmy of any haplogroup-defining variants observed with heteroplasmy 85-99.8 percent alternate alleles if at least 3 such variants are present; otherwise estimated contamination is defined as one minus the mean heteroplasmy of haplogroup-defining variants with heteroplasmy 85-100 percent"
        ),
        indel_pos_counter=hl.struct(
            Description="Dictionary containing the count of the number of indels at each position that contains at least one indel for the sample"
        ),
    )

    row_annotation_dict = hl.struct(
        locus=hl.struct(
            Description="Hail LocusExpression containing contig and position"
        ),
        alleles=hl.struct(
            Description="Alternate allele (multiallelic sites are split)"
        ),
        rsid=hl.struct(Description="The RSID obtained from dbSNP"),
        filters=hl.struct(
            Description="Site or allele-specific filters applied to the variant"
        ),
        variant_collapsed=hl.struct(Description="Variant in format of RefPosAlt"),
        pon_mt_trna_prediction=hl.struct(
            Description="tRNA pathogenicity classification from PON-mt-tRNA"
        ),
        pon_ml_probability_of_pathogenicity=hl.struct(
            Description="tRNA ML_probability_of_pathogenicity from PON-mt-tRNA"
        ),
        mitotip_score=hl.struct(Description="MitoTip raw score"),
        mitotip_trna_prediction=hl.struct(Description="MitoTip score interpretation"),
        region=hl.struct(
            Description="Region (control, non-control, or NA) as obtained from the file supplied in add_variant_context"
        ),
        variant_context=hl.struct(
            Description="Variant in format of REF>ALT with strand (heavy or light) appended"
        ),
        vep=hl.struct(Description="Annotations from VEP"),
        common_low_heteroplasmy=hl.struct(
            Description=f"Flag, true if the overall allele frequency is > 0.001 for samples with a heteroplasmy > 0 and < 0.50 (includes variants < {vaf_filter_threshold} heteroplasmy which are subsequently filtered), is evaluated before other sample-level filters are applied"
        ),
        base_qual_hist=hl.struct(
            Description=f"Histogram of number of individuals failing the base_qual filter (alternate allele median base quality) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}"
        ),
        position_hist=hl.struct(
            Description=f"Histogram of number of individuals failing the position filter (median distance of alternate variants from end of reads) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}"
        ),
        strand_bias_hist=hl.struct(
            Description=f"Histogram of number of individuals failing the strand_bias filter (evidence for alternate allele comes from one read direction only) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}"
        ),
        weak_evidence_hist=hl.struct(
            Description=f"Histogram of number of individuals failing the weak_evidence filter (mutation does not meet likelihood threshold) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}"
        ),
        contamination_hist=hl.struct(
            Description=f"Histogram of number of individuals failing the contamination filter across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}"
        ),
        heteroplasmy_below_min_het_threshold_hist=hl.struct(
            Description=f"Histogram of number of individuals with a heteroplasmy level below {min_het_threshold}, bin edges are: {hl_hist_bin_edges}"
        ),
        excluded_AC=hl.struct(
            Description="Excluded allele count (number of individuals in which the variant was filtered out)"
        ),
        AN=hl.struct(
            Description="Overall allele number (number of samples with non-missing genotype)"
        ),
        AC_hom=hl.struct(
            Description=f"Allele count restricted to variants with a heteroplasmy level >= {min_hom_threshold}"
        ),
        AC_het=hl.struct(
            Description=f"Allele count restricted to variants with a heteroplasmy level >= 0.10 and < {min_hom_threshold}"
        ),
        hl_hist=hl.struct(
            Description="Histogram of heteroplasmy levels",
            sub_annotations=hl.struct(
                bin_edges=hl.struct(Description="Bin edges of the histogram"),
                bin_freq=hl.struct(
                    Description="Count for number of values in each bin"
                ),
                n_smaller=hl.struct(
                    Description="Number of values smaller than the start of the first bin"
                ),
                n_larger=hl.struct(
                    Description="Number of values larger than the end of the last bin"
                ),
            ),
        ),
        dp_hist_all=hl.struct(
            Description="Histogram of dp values for all individuals",
            sub_annotations=hl.struct(
                bin_edges=hl.struct(Description="Bin edges of the histogram"),
                bin_freq=hl.struct(
                    Description="Count for number of values in each bin"
                ),
                n_smaller=hl.struct(
                    Description="Number of values smaller than the start of the first bin"
                ),
                n_larger=hl.struct(
                    Description="Number of values larger than the end of the last bin"
                ),
            ),
        ),
        dp_hist_alt=hl.struct(
            Description="Histogram of dp values for individuals with the alternative allele",
            sub_annotations=hl.struct(
                bin_edges=hl.struct(Description="Bin edges of the histogram"),
                bin_freq=hl.struct(
                    Description="Count for number of values in each bin"
                ),
                n_smaller=hl.struct(
                    Description="Number of values smaller than the start of the first bin"
                ),
                n_larger=hl.struct(
                    Description="Number of values larger than the end of the last bin"
                ),
            ),
        ),
        dp_mean=hl.struct(Description="Mean depth across all individuals for the site"),
        mq_mean=hl.struct(
            Description="Mean MMQ (median mapping quality) across individuals with a variant for the site"
        ),
        tlod_mean=hl.struct(
            Description="Mean TLOD (Log 10 likelihood ratio score of variant existing versus not existing) across individuals with a variant for the site"
        ),
        AF_hom=hl.struct(
            Description=f"Allele frequency restricted to variants with a heteroplasmy level >= {min_hom_threshold}"
        ),
        AF_het=hl.struct(
            Description=f"Allele frequency restricted to variants with a heteroplasmy level >= 0.10 and < {min_hom_threshold}"
        ),
        max_hl=hl.struct(
            Description="Maximum heteroplasmy level observed among all samples with the variant"
        ),
        faf_hapmax_hom=hl.struct(
            Description="Maximum filtering allele frequency across haplogroups restricted to homoplasmic variants"
        ),
    )

    # add descriptions to the matrix table
    input_mt = input_mt.annotate_globals(
        global_annotation_descriptions=hl.literal(global_annotation_dict),
        col_annotation_descriptions=hl.literal(col_annotation_dict),
        row_annotation_descriptions=hl.literal(row_annotation_dict),
    )

    return input_mt


# In[13]:


def adjust_descriptions(input_ht: hl.Table) -> hl.Table:
    """
    Remove descriptions for globals that are no longer present in the ht.
    :param input_ht: Hail Table of variants
    :return: Hail Table with globals descriptions edited to remove descriptions for annotations no longer present in the table
    :rtype: hl.Table
    """
    row_descriptions = [x for x in input_ht.row_annotation_descriptions.keys()]
    row_fields = list(input_ht.row)
    rows_descriptions_to_drop = set(row_descriptions) - set(row_fields)

    global_descriptions = [x for x in input_ht.global_annotation_descriptions.keys()]
    gobal_fields = list(input_ht.globals)
    global_descriptions_to_drop = set(global_descriptions) - set(gobal_fields)

    input_ht = input_ht.annotate_globals(
        row_annotation_descriptions=input_ht.row_annotation_descriptions.drop(
            *rows_descriptions_to_drop
        ),
        global_annotation_descriptions=input_ht.global_annotation_descriptions.drop(
            *global_descriptions_to_drop
        ),
    )
    input_ht = input_ht.drop("col_annotation_descriptions")

    return input_ht


# In[14]:


def add_genotype(mt_path: str, min_hom_threshold: float = 0.95) -> hl.MatrixTable:
    """
    Add in genotype annotation based on heteroplasmy level.
    If the heteroplasmy level is above the min_hom_threshold, set the genotype to 1/1.
    If the heteroplasmy level is less than the min_hom_threshold, but greater than 0, set the genotype to 0/1.
    Otherwise set the genotype to 0/0.
    :param mt_path: Path to the MatrixTable (this MatrixTable can be generated by running combine_vcfs.py)
    :param min_hom_threshold: Minimum heteroplasmy level to define a variant as homoplasmic
    :return: MatrixTable with GT field added
    """
    logger.info("Reading in MT...")
    mt = hl.read_matrix_table(mt_path)

    # Add in genotype (GT) based on min_hom_threshold
    mt = mt.annotate_entries(
        GT=(
            hl.case()
            .when((mt.HL < min_hom_threshold) & (mt.HL > 0.0), hl.parse_call("0/1"))
            .when(mt.HL >= min_hom_threshold, hl.parse_call("1/1"))
            .when(mt.HL == 0, hl.parse_call("0/0"))
            .default(hl.missing(hl.tcall))
        ),
    )

    return mt


# In[15]:


def add_gnomad_metadata(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Add select gnomAD metadata to the MatrixTable.
    :param input_mt: MatrixTable
    :return: MatrixTable with select gnomAD metadata added
    """
    # TODO: Add option here to accomodate non-gnomAD metadata
    genome_meta_ht = meta.versions["3.1"].ht()

    genome_meta_struct = genome_meta_ht[input_mt.s]

    input_mt = input_mt.annotate_cols(
        release=genome_meta_struct.release,
        hard_filters=genome_meta_struct.sample_filters.hard_filters,
        research_project=genome_meta_struct.project_meta.research_project,
        project_id=genome_meta_struct.project_meta.project_id,
        product=genome_meta_struct.project_meta.product,
        sample_pi=genome_meta_struct.project_meta.sample_pi,
        sex_karyotype=genome_meta_struct.sex_imputation.sex_karyotype,
        age=hl.if_else(
            hl.is_defined(genome_meta_struct.project_meta.age),
            genome_meta_struct.project_meta.age,
            genome_meta_struct.project_meta.age_alt,
        ),
        broad_external=genome_meta_struct.project_meta.broad_external,
        pop=genome_meta_struct.population_inference.pop,
    )

    return input_mt


# In[16]:


def filter_by_copy_number(
    input_mt: hl.MatrixTable, keep_all_samples: bool = False
) -> hl.MatrixTable:
    """
    Calculate the mitochondrial copy number based on mean mitochondrial coverage and median nuclear coverage. Filter out samples with more extreme copy numbers.
    Note that median and mean coverage for mitochondria are very similar. Mean mitochondria coverage was used based on metrics available at the time, but releases will switch to using median mitochondria coverage.
    :param hl.MatrixTable input_mt: MatrixTable
    :param keep_all_samples: If True, keep all samples (calculate mitochondrial copy number, but do not filter any samples based on this metric)
    :return: MatrixTable filtered to samples with a copy number of at least 50 and less than 500, number samples below 50 removed, number samples above 500 removed
    """
    # Calculate mitochondrial copy number, if median autosomal coverage is not present default to a wgs_median_coverage of 30x
    input_mt = input_mt.annotate_cols(
        mito_cn=(2 * input_mt.mt_mean_coverage)/ hl.if_else(
            hl.is_missing(input_mt.wgs_mean_coverage),
            30,
            input_mt.wgs_mean_coverage,
        )
    )
    n_removed_below_cn = input_mt.aggregate_cols(
        hl.agg.count_where(input_mt.mito_cn < 50)
    )
    n_removed_above_cn = input_mt.aggregate_cols(
        hl.agg.count_where(input_mt.mito_cn > 500)
    )

    if not keep_all_samples:
        # Remove samples with a mitochondrial copy number below 50 or greater than 500
        input_mt = input_mt.filter_cols(
            (input_mt.mito_cn >= 50) & (input_mt.mito_cn <= 500)
        )
    input_mt = input_mt.filter_rows(hl.agg.any(input_mt.HL > 0))

    return input_mt, n_removed_below_cn, n_removed_above_cn


# In[17]:


def filter_by_contamination(
    input_mt: hl.MatrixTable, output_dir: str, keep_all_samples: bool = False
) -> hl.MatrixTable:
    """
    Solenne : Removed the freemix variable as VerifyBamID is not performed as well as haplotype variables
    
    Calculate contamination based on internal algorithm and filter out samples with contamination above 2%.
    Contamination takes into account:
    a) mitochondria contamination output by HaploCheck
    b) nuclear contamination (freemix) output by VerifyBamID
    c) an internal algorithm with utilizes the PASS haplogroup-defining variants which should be homoplasmic (100% alternate alleles), but in contaminated samples show multiple alleles with heteroplasmy 85-99.8%
    :param input_mt: MatrixTable
    :param output_dir: Output directory to which results should be written
    :param keep_all_samples: If True, keep all samples (calculate contamination, but do not filter any samples based on this metric)
    :return: MatrixTable filtered to samples without contamination, number of contaminated samples removed
    """
    # Generate expression for genotypes with >= 85% heteroplasmy and no FT filters at haplogroup-defining sites that are not filtered as artifact-prone sites
    over_85_expr = (
        (input_mt.HL >= 0.85)
        & (input_mt.FT == {"PASS"})
        #& input_mt.hap_defining_variant
        #& ~hl.str(input_mt.filters).contains("artifact_prone_site")
    )

    input_mt = input_mt.annotate_cols(
        over_85_mean=hl.agg.filter(over_85_expr, hl.agg.mean(input_mt.HL)),
        over_85_count=hl.agg.filter(
            over_85_expr, hl.agg.count_where(hl.is_defined(input_mt.HL))
        ),
        bt_85_and_99_mean=hl.agg.filter(
            over_85_expr & (input_mt.HL <= 0.998), hl.agg.mean(input_mt.HL)
        ),
        bt_85_and_99_count=hl.agg.filter(
            over_85_expr & (input_mt.HL <= 0.998),
            hl.agg.count_where(hl.is_defined(input_mt.HL)),
        ),
    )

    input_mt = input_mt.annotate_cols(
        contam_high_het=hl.if_else(
            input_mt.bt_85_and_99_count >= 3,
            1 - input_mt.bt_85_and_99_mean,
            1 - input_mt.over_85_mean,
        )
    )

    # If contam_high_het is nan, set to 0 (to avoid filtering out missing values which would be more common with haplogroups closer to the reference haplogroup)
    input_mt = input_mt.annotate_cols(
        contam_high_het=hl.if_else(
            hl.is_nan(input_mt.contam_high_het), 0, input_mt.contam_high_het
        )
    )

    # Find samples on border of .02 that may flip between < 0.02 and > 0.02 from issues with floating point precision and mark these samples for removal
    epsilon = 0.000001
    border_samples = input_mt.aggregate_cols(
        hl.agg.filter(
            (input_mt.contam_high_het > (0.02 - epsilon))
            & (input_mt.contam_high_het < (0.02 + epsilon)),
            hl.agg.collect((input_mt.s)),
        )
    )

    border_samples = (
        hl.literal(border_samples) if border_samples else hl.empty_array(hl.tstr)
    )

    # Add annotation to keep only samples with a contamination less than 2%
    input_mt = input_mt.annotate_cols(
        keep=(input_mt.contamination < 0.02)
        #& (input_mt.freemix_percentage < 2)
        & (input_mt.contam_high_het < 0.02)
        & ~border_samples.contains(input_mt.s)
    )
    # Save sample contamination information to separate file
    n_contaminated = input_mt.aggregate_cols(hl.agg.count_where(~input_mt.keep))

    sample_data = input_mt.select_cols(
        "contamination",
        #"freemix_percentage",
        "contam_high_het",
        "over_85_mean",
        "over_85_count",
        "bt_85_and_99_mean",
        "bt_85_and_99_count",
        "keep",
    )
    data_export = sample_data.cols()
    data_export.export(f"{output_dir}/sample_contamination.tsv")

    if not keep_all_samples:
        logger.info(
            "Removing %d samples with contamination above 2 percent", n_contaminated
        )
        input_mt = input_mt.filter_cols(input_mt.keep)
    input_mt = input_mt.drop("keep")

    input_mt = input_mt.filter_rows(hl.agg.any(input_mt.HL > 0))

    return input_mt, n_contaminated


# In[18]:


def add_terra_metadata(
    input_mt: hl.MatrixTable, participant_data: str
) -> hl.MatrixTable:
    """
    Add Terra metadata to the MatrixTable.
    Solenne changes : The VerifyBamID is not used, so freemix is not calculated.
    To do : Create the metadata table in Hail and input it (currently done manually to test)
    
    The participant_data file can be obtained by downloading the participant data after running Mutect2 in Terra. This file should contain the following columns:
        - entity:participant_id: Participant ID uploaded to Terra by user
        - s: Sample ID uploaded to Terra by user
        - contamination: Output by Mutect2, gives the estimate of mitochondrial contamination
        - wgs_median_coverage: Uploaded to Terra by user, can be calculated with Picard's CollectWgsMetrics
        - mt_mean_coverage: Output by Mutect2, gives the mean mitochondrial coverage
        
    Removed by Solenne
        - freemix_percentage: Uploaded to Terra by user, can be calculated with VerifyBamID
        - major_haplogroup: Output by Mutect2 which utilizes Haplogrep    
    
    :param input_mt: MatrixTable
    :param participant_data: Path to metadata file downloaded from Terra
    :return: MatrixTable with Terra metadata annotations added
    """
    # Add haplogroup and Mutect2/Terra output annotations
    ht = hl.import_table(
        participant_data,
        types={
            "contamination": hl.tfloat64,
            #"freemix_percentage": hl.tfloat64,
            "mt_mean_coverage": hl.tfloat64,
            "wgs_mean_coverage": hl.tfloat64,
        },
        missing="",
    ).key_by("s")
    ht = ht.rename({"entity:participant_id": "participant_id"})

    ht = ht.select(
        "participant_id",
        "contamination",
        #"freemix_percentage",
        #"major_haplogroup",
        "wgs_mean_coverage",
        "mt_mean_coverage",
    )

    input_mt = input_mt.annotate_cols(**ht[input_mt.s])

       
    # Annotate the high level haplogroup by taking the first letter, with the exception of H and L haplogroups which are more commonly referred to using the first two letters
    #input_mt = input_mt.annotate_cols(
    #    hap=hl.if_else(
    #        input_mt.major_haplogroup.startswith("HV")
    #        | input_mt.major_haplogroup.startswith("L"),
    #        input_mt.major_haplogroup[0:2],
    #        input_mt.major_haplogroup[0],
    #    )
    #)

    return input_mt


# In[19]:


def get_indel_expr(input_mt: hl.MatrixTable) -> hl.expr.BooleanExpression:
    """
    Generate expression for filtering to indels that should be used to evaluate indel stacks.
    To be considered a variant to be used to evaluate indel stacks, the variant should:
    a) be an indel
    b) have a heteroplasmy level >= 0.01 and <= 0.95
    c) have a PASS genotype
    :param input_mt: MatrixTable
    :return: Expression to be used for determining if a variant is an indel that should to be used to evaluate indel stacks
    """
    indel_expr = (
        hl.is_indel(input_mt.alleles[0], input_mt.alleles[1])
        & (input_mt.HL <= 0.95)
        & (input_mt.HL >= 0.01)
        & (input_mt.FT == {"PASS"})
    )

    return indel_expr


# In[20]:


def generate_expressions(
    input_mt: hl.MatrixTable, min_hom_threshold: float = 0.95
) -> hl.MatrixTable:
    """
    Create expressions to use for annotating the MatrixTable.
    The expressions include AC, AN, AF, filtering allele frequency (FAF) split by homplasmic/heteroplasmic, haplgroup, and population.
    Also includes calcuations of mean DP, MQ, and TLOD.
    :param input_mt: MatrixTable
    :param min_hom_threshold: Minimum heteroplasmy level to define a variant as homoplasmic
    :return: Tuple of hail expressions
    """
    # Calculate AC and AN
    AC = hl.agg.count_where((input_mt.HL > 0.0))
    AN = hl.agg.count_where(hl.is_defined(input_mt.HL))
    # Note: if AN is zero, AFs will evaluate to NaN, which may need to be converted to zero for downstream tools
    AF = AC / AN

    # Calculate AC for het and hom variants, and histogram for HL
    AC_hom = hl.agg.count_where(input_mt.HL >= min_hom_threshold)
    AC_het = hl.agg.count_where((input_mt.HL < min_hom_threshold) & (input_mt.HL > 0.0))
    HL_hist = hl.agg.filter(input_mt.HL > 0, hl.agg.hist(input_mt.HL, 0, 1, 10))
    DP_hist_alt = hl.agg.filter(
        input_mt.GT.is_non_ref(), hl.agg.hist(input_mt.DP, 0, 2000, 10)
    )
    DP_hist_all = hl.agg.hist(input_mt.DP, 0, 2000, 10)
    DP_mean = hl.agg.mean(input_mt.DP)
    MQ_mean = hl.agg.mean(input_mt.MQ)
    TLOD_mean = hl.agg.mean(input_mt.TLOD)

    # Calculate AF
    # Note: if AN is zero, AFs will evaluate to NaN, which may need to be converted to zero for downstream tools
    AF_hom = AC_hom / AN
    AF_het = AC_het / AN

    # Calculate max individual heteroplasmy
    max_HL = hl.agg.max(input_mt.HL)

    # Haplogroup annotations
    #pre_hap_AC = hl.agg.group_by(input_mt.hap, AC)
    #pre_hap_AN = hl.agg.group_by(input_mt.hap, AN)
    #pre_hap_AF = hl.agg.group_by(input_mt.hap, AF)
    #pre_hap_AC_het = hl.agg.group_by(input_mt.hap, AC_het)
    #pre_hap_AC_hom = hl.agg.group_by(input_mt.hap, AC_hom)
    #pre_hap_AF_hom = hl.agg.group_by(input_mt.hap, AF_hom)
    #pre_hap_AF_het = hl.agg.group_by(input_mt.hap, AF_het)
    #pre_hap_HL_hist = hl.agg.group_by(input_mt.hap, HL_hist.bin_freq)
    #pre_hap_FAF = hl.agg.group_by(
    #    input_mt.hap,
    #    hl.experimental.filtering_allele_frequency(hl.int32(AC), hl.int32(AN), 0.95),
    #)
    #pre_hap_FAF_hom = hl.agg.group_by(
    #    input_mt.hap,
    #    hl.experimental.filtering_allele_frequency(
    #        hl.int32(AC_hom), hl.int32(AN), 0.95
    #    ),
    #)

    # population annotations
    #pre_pop_AC = hl.agg.group_by(input_mt.pop, AC)
    #pre_pop_AN = hl.agg.group_by(input_mt.pop, AN)
    #pre_pop_AF = hl.agg.group_by(input_mt.pop, AF)
    #pre_pop_AC_het = hl.agg.group_by(input_mt.pop, AC_het)
    #pre_pop_AC_hom = hl.agg.group_by(input_mt.pop, AC_hom)
    #pre_pop_AF_hom = hl.agg.group_by(input_mt.pop, AF_hom)
    #pre_pop_AF_het = hl.agg.group_by(input_mt.pop, AF_het)
    #pre_pop_HL_hist = hl.agg.group_by(input_mt.pop, HL_hist.bin_freq)

    return hl.struct(
        AC=AC,
        AN=AN,
        AF=AF,
        AC_hom=AC_hom,
        AC_het=AC_het,
        hl_hist=HL_hist,
        dp_hist_all=DP_hist_all,
        dp_hist_alt=DP_hist_alt,
        dp_mean=DP_mean,
        mq_mean=MQ_mean,
        tlod_mean=TLOD_mean,
        AF_hom=AF_hom,
        AF_het=AF_het,
        max_hl=max_HL,
        #pre_hap_AC=pre_hap_AC,
        #pre_hap_AN=pre_hap_AN,
        #pre_hap_AF=pre_hap_AF,
        #pre_hap_AC_het=pre_hap_AC_het,
        #pre_hap_AF_het=pre_hap_AF_het,
        #pre_hap_AC_hom=pre_hap_AC_hom,
        #pre_hap_AF_hom=pre_hap_AF_hom,
        #pre_hap_hl_hist=pre_hap_HL_hist,
        #pre_hap_faf=pre_hap_FAF,
        #pre_hap_faf_hom=pre_hap_FAF_hom,
        #pre_pop_AN=pre_pop_AN,
        #pre_pop_AC_het=pre_pop_AC_het,
        #pre_pop_AF_het=pre_pop_AF_het,
        #pre_pop_AC_hom=pre_pop_AC_hom,
        #pre_pop_AF_hom=pre_pop_AF_hom,
        #pre_pop_hl_hist=pre_pop_HL_hist,
    )


# In[21]:


def add_quality_histograms(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Add histogram annotations for quality metrics to the MatrixTable.
    :param input_mt: MatrixTable
    :return: MatrixTable annotated with quality metric histograms
    """
    # Generate histogram for site quality metrics across all variants
    # TODO: decide on bin edges
    dp_hist_all_variants = input_mt.aggregate_rows(
        hl.agg.hist(input_mt.dp_mean, 0, 4000, 40)
    )
    input_mt = input_mt.annotate_globals(
        dp_hist_all_variants_bin_freq=dp_hist_all_variants.bin_freq,
        dp_hist_all_variants_n_larger=dp_hist_all_variants.n_larger,
        dp_hist_all_variants_bin_edges=dp_hist_all_variants.bin_edges,
    )

    mq_hist_all_variants = input_mt.aggregate_rows(
        hl.agg.hist(input_mt.mq_mean, 0, 80, 40)
    )  # is 80 the actual max value here?
    input_mt = input_mt.annotate_globals(
        mq_hist_all_variants_bin_freq=mq_hist_all_variants.bin_freq,
        mq_hist_all_variants_n_larger=mq_hist_all_variants.n_larger,
        mq_hist_all_variants_bin_edges=mq_hist_all_variants.bin_edges,
    )

    tlod_hist_all_variants = input_mt.aggregate_rows(
        hl.agg.hist(input_mt.tlod_mean, 0, 40000, 40)
    )
    input_mt = input_mt.annotate_globals(
        tlod_hist_all_variants_bin_freq=tlod_hist_all_variants.bin_freq,
        tlod_hist_all_variants_n_larger=tlod_hist_all_variants.n_larger,
        tlod_hist_all_variants_bin_edges=tlod_hist_all_variants.bin_edges,
    )

    # Generate histogram for overall age distribution
    #age_hist_all_samples = input_mt.aggregate_cols(
    #    hl.agg.hist(input_mt.age, 30, 80, 10)
    #)
    #input_mt = input_mt.annotate_globals(
    ##    age_hist_all_samples_bin_freq=age_hist_all_samples.bin_freq,
    #    age_hist_all_samples_n_larger=age_hist_all_samples.n_larger,
    #    age_hist_all_samples_n_smaller=age_hist_all_samples.n_smaller,
    #    age_hist_all_samples_bin_edges=age_hist_all_samples.bin_edges,
    #)

    # Add age histograms per variant type (heteroplasmic or homoplasmic)
    #age_data = age_hists_expr(True, input_mt.GT, input_mt.age)
    #input_mt = input_mt.annotate_rows(
    #    age_hist_hom=age_data.age_hist_hom, age_hist_het=age_data.age_hist_het
    #)

    return input_mt


# In[22]:


def apply_common_low_het_flag(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Apply the common_low_heteroplasmy flag to the MatrixTable.
    The common_low_heteroplasmy flag marks variants where the overall frequency is > 0.001 for samples with a heteroplasmy level > 0 and < 0.50 and either "low_allele_frac" or "PASS" for the genotype filter
    NOTE: The "low_allele_frac" is applied by Mutect2 to variants with a heteroplasmy level below the supplied vaf_filter_threshold
    :param input_mt: MatrixTable
    :return: MatrixTable with the common_low_heteroplasmy flag added
    """
    input_mt = input_mt.annotate_rows(
        AC_mid_het=hl.agg.count_where(
            (input_mt.HL < 0.50)
            & (input_mt.HL > 0.0)
            & ((input_mt.FT == {"PASS"}) | (input_mt.FT == {"low_allele_frac"}))
        )
    )
    input_mt = input_mt.annotate_rows(
        AF_mid_het=input_mt.AC_mid_het
        / hl.agg.count_where(
            hl.is_defined(input_mt.HL)
            & ((input_mt.FT == {"PASS"}) | (input_mt.FT == {"low_allele_frac"}))
        )
    )
    input_mt = input_mt.annotate_rows(
        common_low_heteroplasmy=input_mt.AF_mid_het > 0.001
    )

    return input_mt


# In[23]:


def remove_low_allele_frac_genotypes(
    input_mt: hl.MatrixTable, vaf_filter_threshold: float = 0.01
) -> hl.MatrixTable:
    """
    Remove low_allele_frac genotypes and sets the call to homoplasmic reference.
    NOTE: vaf_filter_threshold should match what was supplied to the vaf_filter_threshold when running Mutect2, variants below this value will be set to homoplasmic reference after calculating the common_low_heteroplasmy filter, Mutect2 will have flagged these variants as "low_allele_frac"
    :param input_mt: MatrixTable
    :param vaf_filter_threshold: Should match vaf_filter_threshold supplied to Mutect2, variants below this value will be set to homoplasmic reference after calculating the common_low_heteroplasmy filter
    :return: MatrixTable with genotypes below the vaf_filter_threshold set to homoplasmic reference
    """
    # Set HL to 0 if < vaf_filter_threshold and remove variants that no longer have at least one alt call
    input_mt = input_mt.annotate_entries(
        HL=hl.if_else(
            (input_mt.HL > 0) & (input_mt.HL < vaf_filter_threshold), 0, input_mt.HL
        )
    )
    # Filter field for all variants with a heteroplasmy of 0 should be set to PASS
    # This step is needed to prevent homref calls that are filtered
    input_mt = input_mt.annotate_entries(
        FT=hl.if_else(input_mt.HL < vaf_filter_threshold, {"PASS"}, input_mt.FT)
    )
    input_mt = input_mt.annotate_entries(
        GT=hl.if_else(
            input_mt.HL < vaf_filter_threshold, hl.parse_call("0/0"), input_mt.GT
        )
    )

    # Check that variants no longer contain the "low_allele_frac" filter (vaf_filter_threshold should be set to appropriate level to remove these variants)
    laf_rows = input_mt.filter_rows(
        hl.agg.any(hl.str(input_mt.FT).contains("low_allele_frac"))
    )
    n_laf_rows = laf_rows.count_rows()
    if n_laf_rows > 0:
        sys.exit(
            "low_allele_frac filter should no longer be present after applying vaf_filter_threshold (vaf_filter_threshold should equal the vaf_filter_threshold supplied to Mutect2)"
        )
    input_mt = input_mt.filter_rows(hl.agg.any(input_mt.HL > 0))

    return input_mt


# In[24]:


def apply_indel_stack_filter(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Apply the indel_stack filter to the MatrixTable.
    The indel_stack filter marks alleles where all samples with the variant call had at least 2 different indels called at the position
    :param input_mt: MatrixTable
    :return: MatrixTable with the indel_stack filter added
    """
    # Add variant-level indel_stack at any indel allele where all samples with a variant call had at least 2 different indels called at that position
    # If any sample had a solo indel at that position, do not filter
    indel_expr = get_indel_expr(input_mt)
    input_mt = input_mt.annotate_cols(
        indel_pos_counter=hl.agg.filter(
            indel_expr, hl.agg.counter(input_mt.locus.position)
        )
    )
    indel_expr = get_indel_expr(input_mt)
    input_mt = input_mt.annotate_entries(
        indel_occurences=(
            hl.case()
            .when(
                (
                    indel_expr
                    & (input_mt.indel_pos_counter.get(input_mt.locus.position) >= 2)
                ),
                "stack",
            )
            .when(
                (
                    indel_expr
                    & (input_mt.indel_pos_counter.get(input_mt.locus.position) == 1)
                ),
                "solo",
            )
            .or_missing()
        )
    )

    
    # If stack is true and solo is false, the indel is stack only and should be filtered out
    input_mt = input_mt.annotate_rows(
        filters=hl.if_else(
            hl.agg.any(input_mt.indel_occurences == "stack")
            & ~hl.agg.any(input_mt.indel_occurences == "solo"),
            input_mt.filters.add("indel_stack"),
            input_mt.filters,
        )
    )

    return input_mt


# In[25]:


def filter_genotypes_below_min_het_threshold(
    input_mt: hl.MatrixTable, min_het_threshold: float = 0.10
) -> hl.MatrixTable:
    """
    Filter out genotypes with a heteroplasmy below the min_het_threshold.
    This filter is a genotype level filter to remove variants with a heteroplasmy level below the specified min_het_threshold
    NOTE: Should later parameterize this function to allow other heteroplasmy cutoffs?
    :param input_mt: MatrixTable
    :param min_het_threshold: Minimum heteroplasmy level to define a variant as a PASS heteroplasmic variant, genotypes below this threshold will count towards the heteroplasmy_below_min_het_threshold filter and be set to missing
    :return: MatrixTable with the heteroplasmy_below_min_het_threshold in the FT field added where applicable
    """
    input_mt = input_mt.annotate_entries(
        FT=hl.if_else(
            (input_mt.HL < min_het_threshold) & (input_mt.GT.is_het()),
            input_mt.FT.add("heteroplasmy_below_min_het_threshold"),
            input_mt.FT,
        )
    )

    # Remove "PASS" from FT if it's not the only filter
    input_mt = input_mt.annotate_entries(
        FT=hl.if_else(input_mt.FT != {"PASS"}, input_mt.FT.remove("PASS"), input_mt.FT)
    )

    return input_mt


# In[26]:


def apply_npg_filter(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Apply the npg filter to the MatrixTable.
    The npg (no pass genotypes) filter marks sites that don't have at least one pass alt call
    :param input_mt: MatrixTable
    :return: MatrixTable with the npg filter added
    """
    input_mt = input_mt.annotate_rows(
        filters=hl.if_else(
            ~(hl.agg.any((input_mt.HL > 0.0) & (input_mt.FT == {"PASS"}))),
            input_mt.filters.add("npg"),
            input_mt.filters,
        )
    )

    return input_mt


# In[27]:


def generate_filter_histogram(
    input_mt: hl.MatrixTable, filter_name: str
) -> hl.ArrayExpression:
    """
    Generate histogram for number of indiviudals with the specified sample-level filter at different heteroplasmy levels.
    :param input_mt: MatrixTable
    :param filter_name: Name of sample-filter for which to generate a histogram
    :return: Histogram containing the counts of individuals with a variant filtered by the specified filter name across binned heteroplasmy levels
    """
    filter_histogram = hl.agg.filter(
        hl.str(input_mt.FT).contains(filter_name), hl.agg.hist(input_mt.HL, 0, 1, 10)
    ).bin_freq

    return filter_histogram


# In[28]:


def add_filter_annotations(
    input_mt: hl.MatrixTable,
    vaf_filter_threshold: float = 0.01,
    min_het_threshold: float = 0.10,
) -> hl.MatrixTable:
    """
    Generate histogram for number of individuals with the specified sample-level filter at different heteroplasmy levels.
    :param input_mt: MatrixTable
    :param vaf_filter_threshold: Should match vaf_filter_threshold supplied to Mutect2, variants below this value will be set to homoplasmic reference after calculating the common_low_heteroplasmy filter
    :param min_het_threshold: Minimum heteroplasmy level to define a variant as a PASS heteroplasmic variant, genotypes below this threshold will count towards the heteroplasmy_below_min_het_threshold filter and be set to missing
    :return: MatrixTable with added annotations for sample and variant level filters and number of genotypes with heteroplasmy_below_min_het_threshold
    """
    # TODO: pull these from header instead?
    filters = [
        "base_qual",
        "position",
        "strand_bias",
        "weak_evidence",
        "contamination",
        "heteroplasmy_below_min_het_threshold",
    ]

    logger.info("Applying common low heteroplasmy flag...")
    input_mt = apply_common_low_het_flag(input_mt)

    logger.info("Removing low_allele_frac genotypes...")
    input_mt = remove_low_allele_frac_genotypes(input_mt, vaf_filter_threshold)

    logger.info("Applying indel_stack filter...")
    input_mt = apply_indel_stack_filter(input_mt)

    logger.info(
        "Filtering genotypes below with heteroplasmy below the min_het_threshold..."
    )
    input_mt = filter_genotypes_below_min_het_threshold(input_mt, min_het_threshold)
    n_het_below_min_het_threshold = input_mt.aggregate_entries(
        hl.agg.count_where(
            hl.str(input_mt.FT).contains("heteroplasmy_below_min_het_threshold")
        )
    )

    logger.info("Applying npg filter...")
    input_mt = apply_npg_filter(input_mt)

    logger.info("Generating filter histograms and calculating excluded_AC...")
    for i in filters:
        annotation_name = i + "_hist"
        input_mt = input_mt.annotate_rows(
            **{annotation_name: generate_filter_histogram(input_mt, i)}
        )
    input_mt = input_mt.annotate_rows(
        excluded_AC=hl.agg.count_where(input_mt.FT != {"PASS"})
    )

    # Remove "PASS" from filters column if it's not the only filter
    input_mt = input_mt.annotate_rows(
        filters=hl.if_else(
            input_mt.filters != {"PASS"},
            input_mt.filters.remove("PASS"),
            input_mt.filters,
        )
    )

    return input_mt, n_het_below_min_het_threshold


# In[29]:


def filter_genotypes(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Set all genotype field values to missing if the variant is not "PASS" for that sample.
    :param input_mt: MatrixTable
    :return: MatrixTable with filtered genotype fields set to missing
    """
    pass_expr = input_mt.FT == {"PASS"}

    input_mt = input_mt.annotate_entries(
        GT=hl.or_missing(pass_expr, input_mt.GT),
        DP=hl.or_missing(pass_expr, input_mt.DP),
        HL=hl.or_missing(pass_expr, input_mt.HL),
        FT=hl.or_missing(pass_expr, input_mt.FT),
        MQ=hl.or_missing(pass_expr, input_mt.MQ),
        TLOD=hl.or_missing(pass_expr, input_mt.TLOD),
    )

    return input_mt


# In[30]:


def add_sample_annotations(
    input_mt: hl.MatrixTable, min_hom_threshold: float = 0.95
) -> hl.MatrixTable:
    """
    Add sample annotations to the MatrixTable.
    These sample annotations include the callrate, number of heteroplasmic/homoplasmic SNPs/indels, and the number of singletons.
    Filtered variants (such as artifact_prone_site, npg, and indel_stack) are excluded from these calculations.
    :param input_mt: MatrixTable
    :param min_hom_threshold: Minimum heteroplasmy level to define a variant as homoplasmic
    :return: MatrixTable with sample annotations added
    """
    # Count number of variants
    num_rows = input_mt.count_rows()

    # Add sample qc annotations
    filter_expr = hl.len(input_mt.filters) == 0
    input_mt = input_mt.annotate_cols(
        callrate=hl.agg.filter(
            filter_expr, (hl.agg.count_where(hl.is_defined(input_mt.HL))) / num_rows
        ),
        n_singletons_het=hl.agg.filter(
            filter_expr,
            hl.agg.count_where(
                (input_mt.AC_het == 1)
                & ((input_mt.HL < min_hom_threshold) & (input_mt.HL > 0.0))
            ),
        ),
        n_singletons_hom=hl.agg.filter(
            filter_expr,
            hl.agg.count_where(
                (input_mt.AC_hom == 1) & (input_mt.HL >= min_hom_threshold)
            ),
        ),
        n_snp_het=hl.agg.filter(
            filter_expr,
            hl.agg.count_where(
                (input_mt.HL < min_hom_threshold)
                & (input_mt.HL > 0.0)
                & hl.is_snp(input_mt.alleles[0], input_mt.alleles[1])
            ),
        ),
        n_snp_hom=hl.agg.filter(
            filter_expr,
            hl.agg.count_where(
                (input_mt.HL >= min_hom_threshold)
                & hl.is_snp(input_mt.alleles[0], input_mt.alleles[1])
            ),
        ),
        n_indel_het=hl.agg.filter(
            filter_expr,
            hl.agg.count_where(
                (input_mt.HL < min_hom_threshold)
                & (input_mt.HL > 0.0)
                & (~hl.is_snp(input_mt.alleles[0], input_mt.alleles[1]))
            ),
        ),
        n_indel_hom=hl.agg.filter(
            filter_expr,
            hl.agg.count_where(
                (input_mt.HL >= min_hom_threshold)
                & (~hl.is_snp(input_mt.alleles[0], input_mt.alleles[1]))
            ),
        ),
    )

    return input_mt


# In[31]:


def add_vep(input_mt: hl.MatrixTable, run_vep: bool, vep_output: str) -> hl.MatrixTable:
    """
    Add vep annotations to the MatrixTable.
    :param input_mt: MatrixTable
    :param run_vep: Whether or not to run vep
    :param vep_output: Path to the MatrixTable output vep results (either the existing results or where to ouput new vep results)
    :return: MatrixTable with vep annotations
    """
    if run_vep:
        vep_mt = hl.vep(input_mt)
        vep_mt = vep_mt.checkpoint(vep_output, overwrite=True)
    else:
        vep_mt = hl.read_matrix_table(vep_output)

    input_mt = input_mt.annotate_rows(
        vep=vep_mt.index_rows(input_mt.locus, input_mt.alleles).vep
    )
    # TODO: get vep version directly from config file
    input_mt = input_mt.annotate_globals(vep_version="v101")

    # If only filter is END_TRUNC, change lof for LC to HC and remove the END_TRUNC filter
    # Remove SINGLE_EXON flags because all exons are single exon in the mitochondria
    input_mt = input_mt.annotate_rows(
        vep=input_mt.vep.annotate(
            transcript_consequences=input_mt.vep.transcript_consequences.map(
                lambda x: x.annotate(
                    lof=hl.if_else(x.lof_filter == "END_TRUNC", "HC", x.lof),
                    lof_filter=hl.if_else(
                        x.lof_filter == "END_TRUNC", hl.missing(hl.tstr), x.lof_filter
                    ),
                    lof_flags=hl.if_else(
                        x.lof_flags == "SINGLE_EXON", hl.missing(hl.tstr), x.lof_flags
                    ),
                )
            )
        )
    )

    end_trunc_count = input_mt.filter_rows(
        hl.str(input_mt.vep.transcript_consequences[0].lof_filter).contains("END_TRUNC")
    ).count_rows()
    if end_trunc_count > 0:
        sys.exit(
            f"END_TRUNC filter should no longer be present but was found for {end_trunc_count} variants"
        )

    single_exon_count = input_mt.filter_rows(
        hl.str(input_mt.vep.transcript_consequences[0].lof_flags).contains(
            "SINGLE_EXON"
        )
    ).count_rows()
    if single_exon_count > 0:
        sys.exit(
            f"SINGLE_EXON flag should no longer be present but was found for {single_exon_count} variants"
        )

    return input_mt


# In[32]:


def add_rsids(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Add rsid annotations to the MatrixTable.
    :param input_mt: MatrixTable
    :return: MatrixTable with rsid annotations added
    """
    dbsnp_import_args = dbsnp.versions["b154"].import_args
    # Replace the contig recoding with just the chrM mapping
    dbsnp_import_args.update({"contig_recoding": {"NC_012920.1": "chrM"}})
    dbsnp_ht = _import_dbsnp(**dbsnp_import_args)

    input_mt = input_mt.annotate_rows(
        rsid=dbsnp_ht[input_mt.locus, input_mt.alleles].rsid
    )
    input_mt = input_mt.annotate_globals(dbsnp_version="b154")

    return input_mt


# In[71]:


def export_simplified_variants(input_ht: hl.Table, output_dir: str) -> None:
    """
    Export a text file containing only several high-level variant annotations.
    :param input_ht: Hail Table of variants
    :param output_dir: Output directory to which results should be output
    :return: None
    """
    reduced_ht = (
        input_ht.key_by(
            chromosome=input_ht.locus.contig,
            position=input_ht.locus.position,
            ref=input_ht.alleles[0],
            alt=input_ht.alleles[1],
        )
        .select("filters", "AC_hom", "AC_het", "AF_hom", "AF_het", "AN", "max_hl", "hl_hist")
        .rename({"max_hl": "max_observed_heteroplasmy", "hl_hist":"heteroplasmy_histogram"})
    )
    reduced_ht = reduced_ht.annotate(
        filters=hl.if_else(
            hl.len(reduced_ht.filters) == 0,
            "PASS",
            hl.str(",").join(hl.array(reduced_ht.filters)),
        )
    )

    reduced_ht.export(f"{output_dir}/reduced_annotations.txt")


# In[34]:


def generate_output_paths(
    output_dir: str, file_name: str, subset_name: str, extension: str
) -> list:
    """
    Generate output paths for results files based on the given output directory, file name, subset name, and extension.
    :param output_dir: Output directory to which results should be output
    :param file_name: Name of the file, preceeds subset_name
    :param subset_name: Name that should be appended to output file names
    :param extension: Extension for the output file
    :return: Path for the file
    """
    # set up output paths for callset
    file_path = f"{output_dir}/{file_name}{subset_name}.{extension}"

    return file_path


# In[35]:


def report_stats(
    input_mt: hl.MatrixTable,
    output_dir: str,
    pass_only: bool,
    n_samples_below_cn: int,
    n_samples_above_cn: int,
    n_samples_contam: int,
    n_het_below_min_het_threshold: int,
    min_het_threshold: float = 0.10,
    min_hom_threshold: float = 0.95,
) -> None:
    """
    Generate output report with basic stats.
    :param input_mt: MatrixTable
    :param output_dir: Output directory to which results should be output
    :param pass_only: Whether or not directory should be filtered to pass_only variants
    :param n_samples_below_cn: Number of samples removed because mitochondrial number is less than 50
    :param n_samples_above_cn: Number of samples removed because mitochondrial number is above 500
    :param n_samples_contam: Number of samples removed because of contamination
    :param n_het_below_min_het_threshold: Number of genotypes filtered because the heteroplasmy levels was below the min_het_threshold
    :param min_het_threshold: Minimum heteroplasmy level to define a variant as a PASS heteroplasmic variant, genotypes below this threshold will count towards the heteroplasmy_below_min_het_threshold filter and be set to missing
    :param min_hom_threshold: Minimum heteroplasmy level to define a variant as homoplasmic
    :return: None
    """
    if pass_only:
        suffix = "_pass"
        input_mt = input_mt.filter_rows(hl.len(input_mt.filters) == 0)
    else:
        suffix = ""
    out_stats = hl.hadoop_open(f"{output_dir}/stats{suffix}.txt", "w")

    if pass_only:
        out_stats.write("Below metrics are for PASS-only variants\n\n")

    # Report numbers of filtered samples/genotypes
    out_stats.write(
        f"Number of samples removed because contamination above 2%: {n_samples_contam}\n"
    )
    out_stats.write(
        f"Number of samples removed because mitochondrial copy number below 50: {n_samples_below_cn}\n"
    )
    out_stats.write(
        f"Number of samples removed because mitochondrial copy number above 500: {n_samples_above_cn}\n"
    )
    out_stats.write(
        f'Number of genotypes filtered because "heteroplasmy_below_min_het_threshold": {n_het_below_min_het_threshold}\n\n'
    )

    # Count variant, samples, bases
    unique_variants, samples = input_mt.count()
    out_stats.write(f"Number of samples: {samples}\n")
    out_stats.write(f"Number of unique variants: {unique_variants}\n")
    bases_w_variant = len(set(input_mt.locus.position.collect()))
    out_stats.write(f"Number of bases with variation: {bases_w_variant}\n\n")

    # Count number of filters
    for filter_name, filter_count in Counter(
        [i for sublist in input_mt.filters.collect() for i in sublist]
    ).items():
        out_stats.write(
            f'Number of variants with "{filter_name}" filter: {filter_count} variants\n'
        )

    # Calculate row stats
    row_stats = input_mt.aggregate_rows(
        hl.struct(
            common_low_het_count=hl.agg.count_where(input_mt.common_low_heteroplasmy),
            het_only_sites=hl.agg.count_where(
                (input_mt.AC_het > 0) & (input_mt.AC_hom == 0)
            ),
            hom_only_sites=hl.agg.count_where(
                (input_mt.AC_hom > 0) & (input_mt.AC_het == 0)
            ),
            het_and_hom_sites=hl.agg.count_where(
                (input_mt.AC_hom > 0) & (input_mt.AC_het > 0)
            ),
            #hap_defining_sites=hl.agg.count_where(input_mt.hap_defining_variant),
            snps=hl.agg.count_where(
                hl.is_snp(input_mt.alleles[0], input_mt.alleles[1])
            ),
            indels=hl.agg.count_where(
                hl.is_indel(input_mt.alleles[0], input_mt.alleles[1])
            ),
            transitions=hl.agg.count_where(
                hl.is_transition(input_mt.alleles[0], input_mt.alleles[1])
            ),
            transversions=hl.agg.count_where(
                hl.is_transversion(input_mt.alleles[0], input_mt.alleles[1])
            ),
        )
    )

    # Calculate col stats
    ##col_stats = input_mt.aggregate_cols(
    #    hl.struct(
    #        unique_haplogroups=hl.len(hl.agg.collect_as_set(input_mt.major_haplogroup)),
    #        unique_top_level_haplogroups=hl.len(hl.agg.collect_as_set(input_mt.hap)),
    #    )
    #)

    # Calculate entry stats
    entry_stats = input_mt.aggregate_entries(
        hl.struct(
            total_variants=hl.agg.count_where(input_mt.HL > 0),
            total_hom_variants=hl.agg.count_where(input_mt.HL >= min_hom_threshold),
            total_het_variants=hl.agg.count_where(
                (input_mt.HL < min_hom_threshold) & (input_mt.HL >= min_het_threshold)
            ),
            min_hl=hl.agg.filter(input_mt.HL > 0, hl.agg.min(input_mt.HL)),
            max_hl=hl.agg.filter(input_mt.HL > 0, hl.agg.max(input_mt.HL)),
        )
    )

    
    # Count number of flags
    out_stats.write(
        f'Number of variants with "common_low_heteroplasmy" flag: {row_stats["common_low_het_count"]} variants\n\n'
    )

    # Count variants
    out_stats.write(f'Total number of variants: {entry_stats["total_variants"]}\n')

    # Count of homoplasmic/heteroplasmic variants
    out_stats.write(
        f'Number of homoplasmic-only sites: {row_stats["hom_only_sites"]}\n'
    )
    out_stats.write(
        f'Number of heteroplasmic-only sites: {row_stats["het_only_sites"]}\n'
    )
    out_stats.write(f'Number of het and hom sites: {row_stats["het_and_hom_sites"]}\n')

    # To add when more samples are available
    #percent_hom = round(
    #    entry_stats["total_hom_variants"] / entry_stats["total_variants"], 2
    #)
    #percent_het = round(
    #    entry_stats["total_het_variants"] / entry_stats["total_variants"], 2
    #)
    out_stats.write(
        f'Total number of homoplasmic variants: {entry_stats["total_hom_variants"]}\n'
    )
    #out_stats.write(f"Percent homoplasmic variants: {percent_hom}\n")
    out_stats.write(
        f'Total number of heteroplasmic variants: {entry_stats["total_het_variants"]}\n'
    )
    #out_stats.write(f"Percent heteroplasmic variants: {percent_het}\n\n")

    out_stats.write(f'Minimum heteroplasmy detected: {entry_stats["min_hl"]}\n')
    out_stats.write(f'Maximum heteroplasmy detected: {entry_stats["max_hl"]}\n\n')

    # Count number of snps and indels
    out_stats.write(f'Number of SNPs: {row_stats["snps"]}\n')
    out_stats.write(f'Number of indels: {row_stats["indels"]}\n')

    # Count number of transitions and transversions
    out_stats.write(f'Number of transitions: {row_stats["transitions"]}\n')
    out_stats.write(f'Number of transversions: {row_stats["transversions"]}\n')
    #out_stats.write(
    #    f'Number of haplogroup defining variants: {row_stats["hap_defining_sites"]}\n\n'
    #)

    # Count number of haplogroups
    #out_stats.write(
    #    f'Number of unique haplogroups: {col_stats["unique_haplogroups"]}\n'
    #)
    #out_stats.write(
    #    f'Number of top-level haplogroups: {col_stats["unique_top_level_haplogroups"]}\n'
    #)

    out_stats.close()


# In[36]:


def change_to_grch38_chrm(input_mt: hl.MatrixTable) -> None:
    """
    Change build to GRCh38 and filters reference genome to chrM.
    :param input_mt: MatrixTable
    :return: MatrixTable with GRCh38 reference genome subsetted to just chrM (so that extraneous contigs will be excluded from VCF output)
    """
    ref = hl.get_reference("GRCh38")
    my_ref = hl.ReferenceGenome(
        "GRCh38_chrM", contigs=["chrM"], lengths={"chrM": ref.lengths["chrM"]}
    )
    assert "chrM" in ref.contigs
    input_mt = input_mt.key_rows_by(
        locus=hl.locus("chrM", input_mt.locus.position, reference_genome="GRCh38_chrM"),
        alleles=input_mt.alleles,
    )

    return input_mt


# In[37]:


def format_vcf(
    input_mt: hl.MatrixTable,
    output_dir: str,
    min_hom_threshold: float = 0.95,
    vaf_filter_threshold: float = 0.01,
    min_het_threshold: float = 0.10,
) -> dict:
    """
    Generate dictionary for VCF header annotations.
    :param input_mt: MatrixTable
    :param min_hom_threshold: Minimum heteroplasmy level to define a variant as homoplasmic
    :param vaf_filter_threshold: Should match vaf_filter_threshold supplied to Mutect2, variants below this value will be set to homoplasmic reference after calculating the common_low_heteroplasmy filter
    :param min_het_threshold: Minimum heteroplasmy level to define a variant as a PASS heteroplasmic variant, genotypes below this threshold will count towards the heteroplasmy_below_min_het_threshold filter and be set to missing
    :param output_dir: Output directory to which appended header info should be written
    :return: MatrixTable with VCF annotations in the info field and dictionary of filter, info, and format fields to be output in the VCF header; path of VCF headers to append
    """
    input_mt = change_to_grch38_chrm(input_mt)

    #haplogroup_order = hl.eval(input_mt.hap_order)
    #population_order = hl.eval(input_mt.pop_order)

    #age_hist_hom_bin_edges = input_mt.age_hist_hom.bin_edges.take(1)[0]
    #age_hist_het_bin_edges = input_mt.age_hist_het.bin_edges.take(1)[0]
    hl_hist_bin_edges = input_mt.hl_hist.bin_edges.take(1)[0]
    dp_hist_all_bin_edges = input_mt.dp_hist_all.bin_edges.take(1)[0]
    dp_hist_alt_bin_edges = input_mt.dp_hist_alt.bin_edges.take(1)[0]

    input_mt = input_mt.annotate_rows(
        hl_hist=input_mt.hl_hist.bin_freq,
        #age_hist_hom_bin_freq=input_mt.age_hist_hom.bin_freq,
        ##age_hist_hom_n_smaller=input_mt.age_hist_hom.n_smaller,
        #age_hist_hom_n_larger=input_mt.age_hist_hom.n_larger,
        #age_hist_het_bin_freq=input_mt.age_hist_het.bin_freq,
        #age_hist_het_n_smaller=input_mt.age_hist_het.n_smaller,
        #age_hist_het_n_larger=input_mt.age_hist_het.n_larger,
        dp_hist_all_n_larger=input_mt.dp_hist_all.n_larger,
        dp_hist_alt_n_larger=input_mt.dp_hist_alt.n_larger,
        dp_hist_all_bin_freq=input_mt.dp_hist_all.bin_freq,
        dp_hist_alt_bin_freq=input_mt.dp_hist_alt.bin_freq,
    )

    #To add when VEP is installed through Hail
    #input_mt = input_mt.annotate_rows(vep=vep_struct_to_csq(input_mt.vep))

    # Get length of annotations to use in Number fields in the VCF where necessary
    #len_hap_hl_hist = len(input_mt.hap_hl_hist.take(1)[0])
    #len_pop_hl_hist = len(input_mt.pop_hl_hist.take(1)[0])

    # Output appended header info to file
    #vcf_header_file = output_dir + "/extra_fields_for_header.tsv"
    #appended_vcf_header = dedent(
    #    f"""
    ##VEP version={hl.eval(input_mt.vep_version)}
    ##dbSNP version={hl.eval(input_mt.dbsnp_version)}
    ##age distributions=bin_edges:{hl.eval(input_mt.age_hist_all_samples_bin_edges)}, bin_freq:{hl.eval(input_mt.age_hist_all_samples_bin_freq)}, n_smaller:{hl.eval(input_mt.age_hist_all_samples_n_smaller)}, n_larger:{hl.eval(input_mt.age_hist_all_samples_n_larger)}
    #"""
    #)
    #with hl.hadoop_open(vcf_header_file, "w") as out:
    #    out.write(appended_vcf_header)

    # Drop intermediate annotations
    input_mt = input_mt.drop(
        "region",
        #"variant_context",
        #"age_hist_het",
        #"age_hist_hom",
        "dp_hist_all",
        "dp_hist_alt",
    )

    ht = input_mt.rows()
    # Move row annotations into info struct
    input_mt = input_mt.annotate_rows(info=hl.struct())
    #input_mt = input_mt.annotate_rows(
    #    info=input_mt.info.annotate(**ht[input_mt.row_key]).drop("rsid")
    #).select_rows(
    #    "rsid", "filters", "info"
    #)  # create info annotation

    # Convert "rsid" array to str for VCF output
    #input_mt = input_mt.annotate_rows(rsid=hl.str(";").join(input_mt.rsid))

    # Convert "," to "|" for array annotations
    for key, value in input_mt.row_value.info.items():
        if str(value).startswith("<Array"):
            if str(value.dtype).startswith("array<array"):
                # If value is an array of arrays, only replace the commas within each individual array
                input_mt = input_mt.annotate_rows(
                    info=input_mt.info.annotate(
                        **{
                            key: hl.map(
                                lambda x: hl.delimit(x, delimiter="|"),
                                input_mt.info[key],
                            )
                        }
                    )
                )
            else:
                input_mt = input_mt.annotate_rows(
                    info=input_mt.info.annotate(
                        **{key: hl.delimit(input_mt.info[key], delimiter="|")}
                    )
                )

    meta_dict = {
        "filter": {
            "artifact_prone_site": {
                "Description": "Variant overlaps site that is commonly reported in literature to be artifact prone"
            },
            "npg": {
                "Description": "No-pass-genotypes site (no individuals were PASS for the variant)"
            },
            "indel_stack": {
                "Description": "Allele where all samples with the variant call had at least 2 different heteroplasmic indels called at the position"
            },
        },
        "info": {
            "variant_collapsed": {
                "Description": "Variant in format of RefPosAlt",
                "Number": "1",
                "Type": "String",
            },
            "hap_defining_variant": {
                "Description": "Present if variant is present as a haplogroup defining variant in PhyloTree build 17",
                "Number": "0",
                "Type": "Flag",
            },
            "common_low_heteroplasmy": {
                "Description": f"Present if variant is found at an overall frequency of .001 across all samples with a heteroplasmy level > 0 and < 0.50 (includes variants <{vaf_filter_threshold} heteroplasmy which are subsequently filtered)",
                "Number": "0",
                "Type": "Flag",
            },
            "AN": {
                "Description": "Overall allele number (number of samples with non-missing genotype)",
                "Number": "1",
                "Type": "Integer",
            },
            "AC_hom": {
                "Description": f"Allele count restricted to variants with a heteroplasmy level >= {min_hom_threshold}",
                "Number": "1",
                "Type": "Integer",
            },  # should put in threshold variable
            "AC_het": {
                "Description": f"Allele count restricted to variants with a heteroplasmy level >= 0.10 and < {min_hom_threshold}",
                "Number": "1",
                "Type": "Integer",
            },
            
            "AF_hom": {
                "Description": f"Allele frequency restricted to variants with a heteroplasmy level >= {min_hom_threshold}",
                "Number": "1",
                "Type": "Float",
            },
            "AF_het": {
                "Description": f"Allele frequency restricted to variants with a heteroplasmy level >= 0.10 and < {min_hom_threshold}",
                "Number": "1",
                "Type": "Float",
            },
            "max_hl": {
                "Description": "Maximum heteroplasmy level observed among all samples for that variant",
                "Number": "1",
                "Type": "Float",
            },
            "bin_edges_hl_hist": {
                "Description": "Bin edges for histogram of heteroplasmy levels",
                "Number": "1",
                "Type": "String",
            },
            "dp_hist_all_bin_freq": {
                "Description": f"Histogram of dp values for all individuals; bin edges are: {dp_hist_all_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "dp_hist_alt_bin_freq": {
                "Description": f"Histogram of dp values for individuals with the alternative allele; bin edges are: {dp_hist_alt_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "dp_mean": {
                "Description": "Mean depth across all individuals for the site",
                "Number": "1",
                "Type": "Float",
            },
            "mq_mean": {
                "Description": "Mean MMQ (median mapping quality) across individuals with a variant for the site",
                "Number": "1",
                "Type": "Float",
            },
            "tlod_mean": {
                "Description": "Mean TLOD (Log 10 likelihood ratio score of variant existing versus not existing) across individuals with a variant for the site",
                "Number": "1",
                "Type": "Float",
            },
            "pon_mt_trna_prediction": {
                "Description": "tRNA pathogenicity classification from PON-mt-tRNA",
                "Number": "1",
                "Type": "String",
            },
            "pon_ml_probability_of_pathogenicity": {
                "Description": "tRNA ML_probability_of_pathogenicity from PON-mt-tRNA",
                "Number": "1",
                "Type": "Float",
            },
            "mitotip_score": {
                "Description": "MitoTip raw score",
                "Number": "1",
                "Type": "Float",
            },
            "mitotip_trna_prediction": {
                "Description": "MitoTip score interpretation",
                "Number": "1",
                "Type": "String",
            },
            "vep": {
                "Description": "Consequence annotations from Ensembl VEP; note that the SINGLE_EXON flag and END_TRUNC filters have been removed from the LOFTEE annotations to avoid misinterpretation in context of the mitochondrial genome. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info",
                "Number": ".",
                "Type": "String",
            },
            "filters": {
                "Description": "Site-level filters",
                "Number": ".",
                "Type": "String",
            },
            "base_qual_hist": {
                "Description": f"Histogram of number of individuals failing the base_qual filter (alternate allele median base quality) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "heteroplasmy_below_min_het_threshold_hist": {
                "Description": f"Histogram of number of individuals with a heteroplasmy level below {min_het_threshold}, bin edges are: {hl_hist_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "position_hist": {
                "Description": f"Histogram of number of individuals failing the position filter (median distance of alternate variants from end of reads) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "strand_bias_hist": {
                "Description": f"Histogram of number of individuals failing the strand_bias filter (evidence for alternate allele comes from one read direction only) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "weak_evidence_hist": {
                "Description": f"Histogram of number of individuals failing the weak_evidence filter (mutation does not meet likelihood threshold) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "contamination_hist": {
                "Description": f"Histogram of number of individuals failing the contamination filter across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "excluded_AC": {
                "Description": "Excluded allele count (number of individuals in which the variant was filtered out)",
                "Number": "1",
                "Type": "String",
            },
        },
        "format": {
            "GT": {
                "Description": f"Genotype, 1/1 if heteroplasmy level >= {min_hom_threshold}, and 0/1 if heteroplasmy level < {min_hom_threshold}",
                "Number": "1",
                "Type": "String",
            },
            "DP": {
                "Description": "Depth of coverage",
                "Number": "1",
                "Type": "Integer",
            },
            "FT": {
                "Description": "Sample-level filters",
                "Number": ".",
                "Type": "String",
            },
            "HL": {"Description": "Heteroplasmy level", "Number": "1", "Type": "Float"},
            "MQ": {"Description": "Mapping quality", "Number": "1", "Type": "Float"},
            "TLOD": {
                "Description": "Log 10 likelihood ratio score of variant existing versus not existing",
                "Number": "1",
                "Type": "Float",
            },
        },
    }

    return input_mt, meta_dict,
    #return input_mt, meta_dict, vcf_header_file


# In[38]:


def add_trna_predictions(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Add tRNA predictions on pathogenicity from PON-mt-tRNA and MitoTIP to the MatrixTable.
    :param input_mt: MatrixTable
    :return: MatrixTable with tRNA predictions of pathogenicity added
    """
    
    input_mt = input_mt.annotate_rows(
        variant_collapsed=input_mt.alleles[0]
        + hl.str(input_mt.locus.position)
        + input_mt.alleles[1]
    )
    
    # Add PON-mt-tRNA predictions
    pon_predictions = hl.import_table(pon_predictions_table)

    # If reference allele from fasta doesn't match Reference_nucleotide, PON-mt-tRNA is reporting the allele of opposite strand and need to get reverse complement for ref and alt
    hl.ReferenceGenome.from_fasta_file('GRCh38_MT_local',GRCh38_MT_local_fasta,GRCh38_MT_local_fai)
    
    #add_reference_sequence(hl.get_reference('GRCh38'))
    pon_predictions = pon_predictions.annotate(
        ref=hl.get_sequence(
            "MT", hl.int(pon_predictions.mtDNA_position), reference_genome="GRCh38_MT_local"
        )
    )
    pon_predictions = pon_predictions.annotate(
        alt=hl.if_else(
            pon_predictions.Reference_nucleotide == pon_predictions.ref,
            pon_predictions.New_nucleotide,
            hl.reverse_complement(pon_predictions.New_nucleotide),
        )
    )
    pon_predictions = pon_predictions.key_by(
        variant_id=pon_predictions.ref
        + hl.str(pon_predictions.mtDNA_position)
        + pon_predictions.alt
    )
    input_mt = input_mt.annotate_rows(
        pon_mt_trna_prediction=pon_predictions[input_mt.variant_collapsed]
        .Classification.lower()
        .replace(" ", "_"),
        pon_ml_probability_of_pathogenicity=hl.float(
            pon_predictions[input_mt.variant_collapsed].ML_probability_of_pathogenicity
        ),
    )

    # Add MitoTIP predictions
    mitotip_predictions = hl.import_table(mitotip_predictions_table)
    mitotip_predictions = mitotip_predictions.key_by(
        variant_id=mitotip_predictions.rCRS
        + hl.str(mitotip_predictions.Position)
        + mitotip_predictions.Alt
    )
    input_mt = input_mt.annotate_rows(
        mitotip_score=hl.float(
            mitotip_predictions[input_mt.variant_collapsed].MitoTIP_Score
        )
    )
    # Set pathogenicity based on MitoTIP scores, classifications obtained from MitoTIP's website
    input_mt = input_mt.annotate_rows(
        mitotip_trna_prediction=(
            hl.case()
            .when(input_mt.mitotip_score > 16.25, "likely_pathogenic")
            .when(
                (input_mt.mitotip_score <= 16.25) & (input_mt.mitotip_score > 12.66),
                "possibly_pathogenic",
            )
            .when(
                (input_mt.mitotip_score <= 12.66) & (input_mt.mitotip_score >= 8.44),
                "possibly_benign",
            )
            .when((input_mt.mitotip_score < 8.44), "likely_benign")
            .or_missing()
        )
    )

    return input_mt


# **Main scripts**

# In[39]:


def main_step1(args):  # noqa: D103
    input_tsv = args.input_tsv
    output_ht = args.output_ht
    temp_dir = args.temp_dir
    chunk_size = args.chunk_size
    overwrite = args.overwrite

    if args.overwrite == False and hl.hadoop_exists(output_ht):
        logger.warning(
            "Overwrite is set to False but file already exists at %s, script will run but output will not be written",
            output_ht,
        )
    # Ensure that user supplied ht extension for output_ht
    if not output_ht.endswith(".ht"):
        sys.exit("Path supplied as output_ht must end with .ht extension")

    mt_list = []
    logger.info(
        "Reading in individual coverage files as matrix tables and adding to a list of matrix tables..."
    )
    with hl.hadoop_open(input_tsv, "r") as f:
        next(f)
        for line in f:
            line = line.rstrip()
            items = line.split("\t")
            participant_id, base_level_coverage_metrics, sample = items[0:3]
            mt = hl.import_matrix_table(
                base_level_coverage_metrics,
                delimiter="\t",
                row_fields={"chrom": hl.tstr, "pos": hl.tint, "target": hl.tstr},
                row_key=["chrom", "pos"],
            ).drop("target")
            mt = mt.rename({"x": "coverage"})
            mt = mt.key_cols_by(s=sample)
            mt_list.append(mt)

    logger.info("Joining individual coverage mts...")
    out_dir = dirname(output_ht)

    cov_mt = multi_way_union_mts(mt_list, temp_dir, chunk_size)
    n_samples = cov_mt.count_cols()

    logger.info("Adding coverage annotations...")
    # Calculate the mean and median coverage as well the fraction of samples above 100x or 1000x coverage at each base
    cov_mt = cov_mt.annotate_rows(
        locus=hl.locus(cov_mt.chrom, cov_mt.pos, reference_genome="GRCh38"),
        mean=hl.float(hl.agg.mean(cov_mt.coverage)),
        median=hl.median(hl.agg.collect(cov_mt.coverage)),
        over_100=hl.float((hl.agg.count_where(cov_mt.coverage > 100) / n_samples)),
        over_1000=hl.float((hl.agg.count_where(cov_mt.coverage > 1000) / n_samples)),
    )

    cov_mt = cov_mt.key_rows_by("locus").drop("chrom", "pos")

    output_mt = re.sub(r"\.ht$", ".mt", output_ht)
    output_tsv = re.sub(r"\.ht$", ".tsv", output_ht)
    output_samples = re.sub(r"\.ht$", "_sample_level.txt", output_ht)

    logger.info("Writing sample level coverage...")
    sample_mt = cov_mt.key_rows_by(pos=cov_mt.locus.position)
    sample_mt.coverage.export(output_samples)

    logger.info("Writing coverage mt and ht...")
    cov_mt.write(output_mt, overwrite=overwrite)
    cov_ht = cov_mt.rows()
    cov_ht = cov_ht.checkpoint(output_ht, overwrite=overwrite)
    cov_ht.export(output_tsv)


# In[40]:


args_step1 = argparse.Namespace(input_tsv=  MT_Step1_input_tsv,
                          output_ht=  'MT_Step1_output.ht',
                          temp_dir= 'MT_Step1_temp_dir',
                          chunk_size= 100,
                          overwrite = True)

main_step1(args_step1)


# In[41]:


META_DICT = {
    "filter": {
        "artifact_prone_site": {
            "Description": "Variant overlaps an artifact-prone site"
        }
    },
    "format": {
        "DP": {"Description": "Depth of coverage", "Number": "1", "Type": "Integer"},
        "FT": {
            "Description": "Sample-level genotype filters",
            "Number": ".",
            "Type": "String",
        },
        "HL": {"Description": "Heteroplasmy level", "Number": "1", "Type": "Float"},
        "MQ": {"Description": "Mapping quality", "Number": "1", "Type": "Float"},
        "TLOD": {
            "Description": "Log 10 likelihood ratio score of variant existing versus not existing",
            "Number": "1",
            "Type": "Float",
        },
    },
}


# In[42]:


def main_step2(args):  # noqa: D103
    participant_data = args.participant_data
    coverage_mt_path = args.coverage_mt_path
    vcf_col_name = args.vcf_col_name
    artifact_prone_sites_path = args.artifact_prone_sites_path
    output_bucket = args.output_bucket
    file_name = args.file_name
    minimum_homref_coverage = args.minimum_homref_coverage
    chunk_size = args.chunk_size

    output_path_mt = f"{output_bucket}/raw_combined.mt"

    if args.overwrite == False and hl.hadoop_exists(output_path_mt):
        logger.warning(
            "Overwrite is set to False but file already exists at %s, script will run but output will not be written",
            output_path_mt,
        )

    logger.info("Collecting VCF paths for samples to subset...")
    vcf_paths = collect_vcf_paths(
        participant_data, vcf_col_name, args.participants_to_subset
    )

    logger.info("Combining VCFs...")
    combined_mt = join_mitochondria_vcfs_into_mt(vcf_paths, args.temp_dir, chunk_size)
    combined_mt = combined_mt.checkpoint(output_path_mt, overwrite=args.overwrite)

    logger.info("Removing select sample-level filters...")
    combined_mt = remove_genotype_filters(combined_mt)

    logger.info("Determining homoplasmic reference sites...")
    combined_mt = determine_hom_refs(
        combined_mt, coverage_mt_path, minimum_homref_coverage
    )

    logger.info("Applying artifact_prone_site fiter...")
    combined_mt = apply_mito_artifact_filter(combined_mt, artifact_prone_sites_path)

    logger.info("Writing combined MT and VCF...")
    # Set the file names for output files
    out_vcf = f"{output_bucket}/{file_name}.vcf.bgz"
    out_mt = f"{output_bucket}/{file_name}.mt"

    combined_mt = combined_mt.checkpoint(out_mt, overwrite=args.overwrite)
    # For the VCF output, join FT values by semicolon
    combined_mt = combined_mt.annotate_entries(
        FT=hl.str(";").join(hl.array(combined_mt.FT))
    )
    hl.export_vcf(combined_mt, out_vcf, metadata=META_DICT)


# In[44]:


args_step2 = argparse.Namespace(minimum_homref_coverage= 5,
                                participant_data= MT_Step2_participant_data,
                                artifact_prone_sites_path= artifact_prone_sites_bed,
                                participants_to_subset = MT_participants_to_subset,
                                coverage_mt_path= 'MT_Step1_output.mt/',
                                vcf_col_name= 'vcf_path',
                                output_bucket = 'MT_Step2_output',
                                file_name= 'MT_Step2_output',
                                temp_dir = 'MT_Step2_temp_dir',
                                chunk_size = 100,
                                overwrite = True)

main_step2(args_step2)


# In[45]:


def main_step3(args):  # noqa: D103
    mt_path = args.mt_path
    output_dir = args.output_dir
    participant_data = args.participant_data
    vep_results = args.vep_results
    min_hom_threshold = args.min_hom_threshold
    vaf_filter_threshold = args.vaf_filter_threshold
    min_het_threshold = args.min_het_threshold
    gnomad_subset = args.subset_to_gnomad_release
    keep_all_samples = args.keep_all_samples
    run_vep = args.run_vep

    logger.info("Cutoff for homoplasmic variants is set to %.2f...", min_hom_threshold)

    # Define mt path, output directory, subset name
    subset_name = ""

    logger.info("Adding genotype annotation...")
    mt = add_genotype(mt_path, min_hom_threshold)

    logger.info("Adding annotations from Terra...")
    mt = add_terra_metadata(mt, participant_data)

    #logger.info("Annotating haplogroup-defining variants...")
    #No haplogroup in the IBVL
    #mt = add_hap_defining(mt)

    logger.info("Annotating tRNA predictions...")
    mt = add_trna_predictions(mt)

    #No gnomAD data for now
    #If 'subset-to-gnomad-release' is set, 'age' and 'pop' are added by the add_gnomad_metadata function.
    # If 'subset-to-gnomad-release' is not set, the user should include an 'age' and 'pop' column in the file supplied to `participant-data`.
    #if gnomad_subset:
    #    logger.info("Adding gnomAD metadata sample annotations...")
    #    mt = add_gnomad_metadata(mt)
    #else:
    #    logger.info("Adding age and pop annotations...")
    #    mt = add_age_and_pop(mt, participant_data)

    #logger.info("Adding variant context annotations...")
    #Require download of a resource file : Not applicable
    #mt = add_variant_context(mt)

    # If specified, subet to only the gnomAD samples in the current release
    #if gnomad_subset:
    #    logger.warning("Subsetting results to gnomAD release samples...")
    #    subset_name = "_gnomad"

        # Subset to release samples and filter out rows that no longer have at least one alt call
        #mt = mt.filter_cols(mt.release)  # Filter to cols where release is true
        #mt = mt.filter_rows(hl.agg.any(mt.HL > 0))

    
    logger.info("Checking for samples with low/high mitochondrial copy number...")
    mt, n_removed_below_cn, n_removed_above_cn = filter_by_copy_number(
        mt, keep_all_samples
    )
    
    logger.info("Checking for contaminated samples...")
    mt, n_contaminated = filter_by_contamination(mt, output_dir, keep_all_samples)

    logger.info("Switch build and checkpoint...")
    # Switch build 37 to build 38
    mt = mt.key_rows_by(
        locus=hl.locus("chrM", mt.locus.position, reference_genome="GRCh38"),
        alleles=mt.alleles,
    )
    mt = mt.checkpoint(f"{output_dir}/prior_to_vep.mt", overwrite=args.overwrite)

    #logger.info("Adding vep annotations...")
    #TO DO ONCE VEP is set up in Hail
    #mt = add_vep(mt, run_vep, vep_results)

    #logger.info("Adding dbsnp annotations...")
    #Will be done through VEP as gnomad module need to download a specific file (inaccessible to GPCC)
   #mt = add_rsids(mt)

    logger.info("Setting up output paths...")
    annotated_mt_path = generate_output_paths(
        output_dir, "annotated_combined", subset_name, "mt"
    )
    sites_ht_path = generate_output_paths(
        output_dir, "combined_sites_only", subset_name, "ht"
    )
    sites_txt_path = generate_output_paths(
        output_dir, "combined_sites_only", subset_name, "txt"
    )
    sites_vcf_path = generate_output_paths(
        output_dir, "combined_sites_only", subset_name, "vcf.bgz"
    )
    samples_txt_path = generate_output_paths(
        output_dir, "sample_annotations", subset_name, "txt"
    )
    samples_vcf_path = generate_output_paths(
        output_dir, "sample_vcf", subset_name, "vcf.bgz"
    )

    logger.info("Results will be output to the following files:")
    print(
        "\n".join(
            [
                annotated_mt_path,
                sites_ht_path,
                sites_txt_path,
                sites_vcf_path,
                samples_txt_path,
                samples_vcf_path,
            ]
        )
    )

    logger.info("Annotating MT...")
    mt, n_het_below_min_het_threshold = add_filter_annotations(
    mt, vaf_filter_threshold, min_het_threshold
    )

    mt = mt.checkpoint(
        f"{output_dir}/prior_to_filter_genotypes.mt", overwrite=args.overwrite
    )

    mt = filter_genotypes(mt)
    # Add variant annotations such as AC, AF, and AN
    mt = mt.annotate_rows(**dict(generate_expressions(mt, min_hom_threshold)))
    # Checkpoint to help avoid Hail errors from large queries
    mt = mt.checkpoint(f"{output_dir}/temp.mt", overwrite=args.overwrite)
    mt = add_quality_histograms(mt)
    #mt = add_annotations_by_hap_and_pop(mt) --> Removed because haplotype and pop are not used in the IBVL
    
    mt = add_descriptions(
        mt, min_hom_threshold, vaf_filter_threshold, min_het_threshold
    )
    mt = mt.checkpoint(
        annotated_mt_path, overwrite=args.overwrite
    )  # Full matrix table for internal use

    logger.info("Generating summary statistics reports...")
    report_stats(
        mt,
        output_dir,
        False,
        n_removed_below_cn,
        n_removed_above_cn,
        n_contaminated,
        n_het_below_min_het_threshold,
    )
    report_stats(
        mt,
        output_dir,
        True,
        n_removed_below_cn,
        n_removed_above_cn,
        n_contaminated,
        n_het_below_min_het_threshold,
    )
    
    logger.info("Writing ht...")
    variant_ht = mt.rows()
    variant_ht = variant_ht.drop("region")
    #variant_ht = variant_ht.drop("region", "variant_context")
    #To add when gnomAD works through Hail
    #variant_ht = adjust_descriptions(variant_ht)
    variant_ht.export(sites_txt_path)  # Sites-only txt file for external use
    variant_ht.write(
        sites_ht_path, overwrite=args.overwrite
    )  # Sites-only ht for external use

    logger.info("Writing sample annotations...")
    mt = add_sample_annotations(mt, min_hom_threshold)
    sample_ht = mt.cols()
    #sample_ht.group_by(sample_ht.hap).aggregate(n=hl.agg.count()).export(
    #    f"{output_dir}/haplogroup_counts.txt"
    #)  # Counts of top level haplogroups
    sample_ht.export(samples_txt_path)  # Sample annotations txt file for internal use#

    logger.info("Formatting and writing VCF...")
    rows_ht = mt.rows()
    export_simplified_variants(rows_ht, output_dir)
    vcf_mt, vcf_meta = format_vcf(mt, output_dir, min_hom_threshold)
    #vcf_mt, vcf_meta, vcf_header_file = format_vcf(mt, output_dir, min_hom_threshold)
    hl.export_vcf(
        vcf_mt,
        samples_vcf_path,
        metadata=vcf_meta,
        #append_to_header=vcf_header_file,
        tabix=True,
    )  # Full VCF for internal use
    vcf_variant_ht = vcf_mt.rows()
    rows_mt = hl.MatrixTable.from_rows_table(vcf_variant_ht).key_cols_by(s="foo")
    hl.export_vcf(
        rows_mt,
        sites_vcf_path,
        metadata=vcf_meta,
    #    append_to_header=vcf_header_file,
        tabix=True,
    )  # Sites-only VCF for external use

    logger.info("All annotation steps are completed")


# In[ ]:


args_step3 = argparse.Namespace(mt_path = 'MT_Step2_output/MT_Step2_output.mt',
                                participant_data = MT_Step3_participant_data,
                                output_dir = 'MT_Step3_output_dir',
                                vep_results = 'MT_Step3_vep_results',
                                min_hom_threshold = 0.95,
                                vaf_filter_threshold = 0.01,
                                min_het_threshold = 0.10,
                                subset_to_gnomad_release = True,
                                keep_all_samples = True,
                                run_vep = True,
                                overwrite = True)
    
main_step3(args_step3)


# In[ ]:




