#!/usr/bin/env python
# coding: utf-8

# temp dir setting
import sys
temp_directory=sys.argv[3]

# Hail and plot initialisation
import hail as hl
from hail.plot import output_notebook, show
hl.init(tmp_dir=temp_directory)
output_notebook()

from hail.plot import show
from pprint import pprint
from bokeh.models import Span
hl.plot.output_notebook()
from bokeh.models import Range1d
from bokeh.plotting import figure, output_file, show, save

import pandas as pd
import os
from typing import Tuple
import string

from typing import Optional, Dict, List, Union

#Created through the nextflow pipeline

# Phil Richmond added 2023-01-30
# Create ref genome
# Based on - https://hail.is/docs/0.2/genetics/hail.genetics.ReferenceGenome.html#hail.genetics.ReferenceGenome.from_fasta_file
# classmethod from_fasta_file(name, fasta_file, index_file, x_contigs=[], y_contigs=[], mt_contigs=[], par=[])
# I can't reuse GRCh38, so instead I'll make a silly 'GRCh38ForHail'
# I pass the "assembly" variables in sys.argv[4], so I'll pass the fasta in sys.argv[5], the fasta index in sys.argv[6]

# here catching the 'X', 'Y', 'MT'
try:
    RefGenome = hl.ReferenceGenome.from_fasta_file('%sForHail'%sys.argv[4],sys.argv[5],sys.argv[6], x_contigs=['X'], y_contigs=['Y'], mt_contigs=['MT'])
    hl.import_vcf(sys.argv[1],array_elements_required=False, force_bgz=True, reference_genome=RefGenome).write('filtered_samples_vcf.mt', overwrite=True)
except:
# here catching the chrX, chrY, chrM
    RefGenome = hl.ReferenceGenome.from_fasta_file('%sForHail'%sys.argv[4],sys.argv[5],sys.argv[6], x_contigs=['X'], y_contigs=['Y'], mt_contigs=['MT'])
    hl.import_vcf(sys.argv[1],array_elements_required=False, force_bgz=True, reference_genome=RefGenome).write('filtered_samples_vcf.mt', overwrite=True)

#hl.import_vcf(sys.argv[1],array_elements_required=False, force_bgz=True).write('filtered_samples_vcf.mt', overwrite=True)
sex_table = (hl.import_table(sys.argv[2], impute=True).key_by('s'))

#**Import file**
#Inport a vcf file and read it as a matrix table (mt, hail specific file type)

mt = hl.read_matrix_table('filtered_samples_vcf.mt')

#**Graph functions**
#In order to create the graph, 3 functions were needed
#- stat : To calcualte the mean, standard deviation and other metrics for each parameter
#- plot_histo : To create the histogram as expected
#- plot_sp : To create the scatter plots as expected

def stat(table):
    Mean = table[table.columns[2]]. mean()  
    StdDev = table[table.columns[2]]. std()
    Low_threashold = Mean - 3*StdDev
    High_threashold = Mean + 3*StdDev
    min_graph = table[table.columns[2]]. min() - 3*StdDev
    max_graph = table[table.columns[2]]. max() + 3*StdDev
    return Mean, StdDev, Low_threashold, High_threashold, min_graph, max_graph

def plot_histo (table_plot, mt_plot, variable) :
    output_file(filename=os.path.join(("variant_QC_"+variable+".html")), title="Variant QC HTML file")
    p = hl.plot.histogram(mt_plot,
                      range = (stat(table_plot) [4], stat(table_plot) [5]),
                      bins = 60,
                      legend=variable,
                      title="Red lines are Mean +/- 3xStdDev")
    annot = Span(dimension="height",location=stat(table_plot) [2],line_dash='dashed', line_width=3,line_color="red")
    p.add_layout(annot)
    annot2 = Span(dimension="height",location=stat(table_plot) [3],line_dash='dashed', line_width=3,line_color="red")
    p.add_layout(annot2)
    p.yaxis.axis_label = 'Count'
    return save(p)

def plot_sp (table_x_axis, mt_x_axis, table_y_axis, mt_y_axis, x_variable, y_variable) :
    output_file(filename=os.path.join(("variant_QC_"+x_variable+"_"+y_variable+".html")), title="Variant QC HTML file")
    p = hl.plot.scatter(x=mt_x_axis,
                   y=mt_y_axis,
                  xlabel=x_variable,
                  ylabel=y_variable,
                title="Red lines are Mean +/- 3xStdDev",
                  size=5)
    annot = Span(dimension="height",location=stat(table_x_axis) [2],line_dash='dashed', line_width=3,line_color="red")
    annot2 = Span(dimension="height",location=stat(table_x_axis) [3],line_dash='dashed', line_width=3,line_color="red")
    annot3 = Span(dimension="width",location=stat(table_y_axis) [2],line_dash='dashed', line_width=3,line_color="red")
    annot4 = Span(dimension="width",location=stat(table_y_axis) [3],line_dash='dashed', line_width=3,line_color="red")
    p.add_layout(annot)
    p.add_layout(annot2)
    p.add_layout(annot3)
    p.add_layout(annot4)
    p.x_range=Range1d(stat(table_x_axis) [4], stat(table_x_axis) [5])
    p.y_range=Range1d(stat(table_y_axis) [4], stat(table_y_axis) [5])
    return save(p)


# **Variant QC**

mt = hl.variant_qc(mt)

# List of variables for which we will create a table, calculate the standard deviation (StdDev) and the mean (Mean) for sample QC:
# - DP (mt_sample_qc.variant_qc.dp_stats.mean)
# - QG (mt_sample_qc.vaiant_qc.gq_stats.mean)
# - call_rate (mt_sample_qc.variant_qc.call_rate)
# - AN (mt_sample_qc.variant_qc.AN)
# - n_not_called (mt_sample_qc.variant_qc.n_not_called)
# - p_value_hwe (mt_sample_qc.variant_qc.p_value_hwe)
# - het_freq_hwe (mt_sample_qc.variant_qc.het_freq_hwe)
# - n_het (mt_sample_qc.variant_qc.n_het)

mt.variant_qc.dp_stats.mean.export('DP_SNV.tsv')
mt.variant_qc.gq_stats.mean.export('GQ_SNV.tsv')
mt.variant_qc.call_rate.export('call_rate_SNV.tsv')
mt.variant_qc.AN.export('AN_SNV.tsv')
mt.variant_qc.n_not_called.export('n_not_called_SNV.tsv')
mt.variant_qc.p_value_hwe.export('p_value_hwe_SNV.tsv')
mt.variant_qc.het_freq_hwe.export('het_freq_hwe_SNV.tsv')
mt.variant_qc.n_het.export('n_het_SNV.tsv')

DP_SNV_table=pd.read_table('DP_SNV.tsv')
GQ_SNV_table=pd.read_table('GQ_SNV.tsv')
call_rate_SNV_table=pd.read_table('call_rate_SNV.tsv')
AN_SNV_table=pd.read_table('AN_SNV.tsv')
n_not_called_SNV_table=pd.read_table('n_not_called_SNV.tsv')
p_value_hwe_SNV_table=pd.read_table('p_value_hwe_SNV.tsv')
het_freq_hwe_SNV_table=pd.read_table('het_freq_hwe_SNV.tsv')
n_het_SNV_table=pd.read_table('n_het_SNV.tsv')

DP_SNV_table.rename(columns = {DP_SNV_table.columns[2]:'DP'}, inplace = True)
GQ_SNV_table.rename(columns = {GQ_SNV_table.columns[2]:'GQ'}, inplace = True)
call_rate_SNV_table.rename(columns = {call_rate_SNV_table.columns[2]:'call_rate'}, inplace = True)
AN_SNV_table.rename(columns = {AN_SNV_table.columns[2]:"AN"}, inplace = True)
n_not_called_SNV_table.rename(columns = {n_not_called_SNV_table.columns[2]:"n_not_called"}, inplace = True)
p_value_hwe_SNV_table.rename(columns = {p_value_hwe_SNV_table.columns[2]:"p_value_hwe"}, inplace = True)
het_freq_hwe_SNV_table.rename(columns = {het_freq_hwe_SNV_table.columns[2]:"het_freq_hwe"}, inplace = True)
n_het_SNV_table.rename(columns = {n_het_SNV_table.columns[2]:"n_het"}, inplace = True)


plot_histo(DP_SNV_table,
           mt.variant_qc.dp_stats.mean,
           'Mean Depth per variant - Unfiltered SNV')

plot_histo(GQ_SNV_table,
           mt.variant_qc.gq_stats.mean,
           'Mean Genotype Quality per variant - Unfiltered SNV')

plot_histo(call_rate_SNV_table,
           mt.variant_qc.call_rate,
           'Call rate per variant - Unfiltered SNV')

plot_histo(AN_SNV_table,
           mt.variant_qc.AN,
           'Allele number per variant - Unfiltered SNV')

plot_histo(n_not_called_SNV_table,
           mt.variant_qc.n_not_called,
           'Number of samples with a missing genotype per variant - Unfiltered SNV')

plot_histo(p_value_hwe_SNV_table,
           mt.variant_qc.p_value_hwe,
           'p-value from two-sided test of Hardy-Weinberg equilibrium per variant - Unfiltered SNV')

plot_sp (het_freq_hwe_SNV_table,
         mt.variant_qc.het_freq_hwe,
         n_het_SNV_table,
         mt.variant_qc.n_het,
         'Expected frequency of heterozygous samples under Hardy-Weinberg equilibrium',
         'Number of heterozygous samples per variant - Unfiltered SNV')


#**Variant filtering**

#1. Autosome and sexual chromosomes only
#Keep only variants on chr 1-22 + X, Y (remove MT and GL)

#2. Length
#Remove small insertions / deletions (indel) of length > 50bp 
#Goal : Avoid overlap with the SV pipeline

#3. Variant quality
#Filter the variants based on the threasholds repersented on the figures 
#Low_threashold = Mean - 3*StdDev = stat(table) [2]
#High_threashold = Mean + 3*StdDev = stat(table) [3]

#Filters :
#- Mean DP per variant lower than the low threshoold
#- Mean Genotype Quality per variant lower than the low threshoold
#- Call rate per variant lower than the low threshoold
#- Allele Number (AN) per variant lower than the low threshoold
#- Number of samples with missing genotype per variant lower than the low threshoold
#Not implemented :  Hardy–Weinberg values

intervals = [hl.parse_locus_interval(x) for x in ['X', 'Y', '1-22']]
SNV_mt_var_filtered = hl.filter_intervals(mt, intervals, keep=True)

SNV_mt_var_filtered = SNV_mt_var_filtered.filter_rows(
    (SNV_mt_var_filtered.variant_qc.dp_stats.mean > stat(DP_SNV_table) [2]) &
    (SNV_mt_var_filtered.variant_qc.gq_stats.mean > stat(GQ_SNV_table) [2]) &
    (SNV_mt_var_filtered.variant_qc.call_rate > stat(call_rate_SNV_table) [2]) &
    (SNV_mt_var_filtered.variant_qc.AN > stat(AN_SNV_table) [2]) &
    (SNV_mt_var_filtered.variant_qc.n_not_called > stat(n_not_called_SNV_table) [2]) &
    (hl.len(SNV_mt_var_filtered.alleles[0]) < 50) &
    (hl.len(SNV_mt_var_filtered.alleles[1]) < 50)
)

#**QC Report**
#Write the report of the number of filtered out variants and the reason they were filtered out

def calc_removed_variant(mt, mt_var, stat_table) :
    input_mt = mt.annotate_rows(
        keep=(mt_var > stat_table [2]))
    n_removed = input_mt.aggregate_rows(hl.agg.count_where(~input_mt.keep))
    
    return n_removed

def report_stats():
    """
    Generate output report with basic stats.
    """
    out_stats = hl.hadoop_open(f"SNV_indel_QC_report.txt", "w")
    # Report numbers of filtered SNV/indels
    out_stats.write(
        f"Number of SNV/indels not located on autosomes or sexual chromosomes : {n_non_chr}\n"
        f"Number of SNV/indels removed because of deletion superior to 50bp: {n_large_del}\n"
        f"Number of SNV/indels removed because of insertion superior to 50bp: {n_large_ins}\n"
        f"Number of SNV/indels removed because of depth metrics: {DP_var_removed}\n"
        f"Number of SNV/indels removed because of genotype quality metrics: {GQ_var_removed}\n"
        f"Number of SNV/indels removed because of call rate metrics: {CR_var_removed}\n"
        f"Number of SNV/indels removed because of allele number (AN): {n_AN_removed}\n"
        f"Number of SNV/indels removed because of number of not called: {n_not_called_removed}\n"
        f"Total number of SNV/indels removed : {n_var_removed}\n"
        f"Percentage of the SNV/indels filtered out: {perc_removed_variants}\n"
    )
    out_stats.close()
    
n_non_chr = mt.count()[0] - hl.filter_intervals(mt, intervals, keep=True).count()[0]

n_large_del = mt.filter_rows(hl.len(mt.alleles[0]) > 50).count()[0]
n_large_ins = mt.filter_rows(hl.len(mt.alleles[1]) > 50).count()[0]

DP_var_removed = calc_removed_variant(mt, mt.variant_qc.dp_stats.mean, stat(DP_SNV_table))
GQ_var_removed = calc_removed_variant(mt, mt.variant_qc.gq_stats.mean, stat(GQ_SNV_table))
CR_var_removed = calc_removed_variant(mt, mt.variant_qc.call_rate, stat(call_rate_SNV_table))
n_AN_removed = calc_removed_variant(mt, mt.variant_qc.AN, stat(AN_SNV_table))
n_not_called_removed = calc_removed_variant(mt, mt.variant_qc.n_not_called, stat(n_not_called_SNV_table))
n_var_removed = (mt.count()[0]-SNV_mt_var_filtered.count()[0])
perc_removed_variants = (mt.count()[0]-SNV_mt_var_filtered.count()[0])/mt.count()[0] * 100


report_stats()


#**Calculate sex specific frequencies**
#Sex is defined using F-stat in Hail_sample_QC or file with sample sex can be loaded by user
#Calculate AF, AC, AN and number of homozygotes
#Code adapted from gnomAD : https://github.com/broadinstitute/gnomad_qc/blob/main/gnomad_qc/v3/annotations/generate_freq_data.py

SNV_mt_var_filtered = SNV_mt_var_filtered.annotate_cols(**sex_table[SNV_mt_var_filtered.s])

def annotate_freq(
    mt: hl.MatrixTable,
    sex_expr: Optional[hl.expr.StringExpression] = None,
    pop_expr: Optional[hl.expr.StringExpression] = None,
    subpop_expr: Optional[hl.expr.StringExpression] = None,
    additional_strata_expr: Optional[Dict[str, hl.expr.StringExpression]] = None,
    downsamplings: Optional[List[int]] = None,
) -> hl.MatrixTable:
    """
    Adapted from gnomAD to remove adj and population specific
    Annotate `mt` with stratified allele frequencies.

    The output Matrix table will include:
        - row annotation `freq` containing the stratified allele frequencies
        - global annotation `freq_meta` with metadata
        - global annotation `freq_sample_count` with sample count information

    .. note::

        Currently this only supports bi-allelic sites.
        The input `mt` needs to have the following entry fields:
        - GT: a CallExpression containing the genotype
        - adj: a BooleanExpression containing whether the genotype is of high quality or not.
        All expressions arguments need to be expression on the input `mt`.

    .. rubric:: `freq` row annotation

    The `freq` row annotation is an Array of Struct, with each Struct containing the following fields:

        - AC: int32
        - AF: float64
        - AN: int32
        - homozygote_count: int32

    Each element of the array corresponds to a stratification of the data, and the metadata about these annotations is
    stored in the globals.

    .. rubric:: Global `freq_meta` metadata annotation

    The global annotation `freq_meta` is added to the input `mt`. It is a list of dict.
    Each element of the list contains metadata on a frequency stratification and the index in the list corresponds
    to the index of that frequency stratification in the `freq` row annotation.

    .. rubric:: Global `freq_sample_count` annotation

    The global annotation `freq_sample_count` is added to the input `mt`. This is a sample count per sample grouping
    defined in the `freq_meta` global annotation.

    .. rubric:: The `downsamplings` parameter

    If the `downsamplings` parameter is used, frequencies will be computed for all samples and by population
    (if `pop_expr` is specified) by downsampling the number of samples without replacement to each of the numbers specified in the
    `downsamplings` array, provided that there are enough samples in the dataset.
    In addition, if `pop_expr` is specified, a downsampling to each of the exact number of samples present in each population is added.
    Note that samples are randomly sampled only once, meaning that the lower downsamplings are subsets of the higher ones.

    :param mt: Input MatrixTable
    :param sex_expr: When specified, frequencies are stratified by sex. If `pop_expr` is also specified, then a pop/sex stratifiction is added.
    :param pop_expr: When specified, frequencies are stratified by population. If `sex_expr` is also specified, then a pop/sex stratifiction is added.
    :param subpop_expr: When specified, frequencies are stratified by sub-continental population. Note that `pop_expr` is required as well when using this option.
    :param additional_strata_expr: When specified, frequencies are stratified by the given additional strata found in the dict. This can e.g. be used to stratify by platform.
    :param downsamplings: When specified, frequencies are computed by downsampling the data to the number of samples given in the list. Note that if `pop_expr` is specified, downsamplings by population is also computed.
    :return: MatrixTable with `freq` annotation
    """

    if additional_strata_expr is None:
        additional_strata_expr = {}

    _freq_meta_expr = hl.struct(**additional_strata_expr)
    if sex_expr is not None:
        _freq_meta_expr = _freq_meta_expr.annotate(sex=sex_expr)
    
    # Annotate cols with provided cuts
    mt = mt.annotate_cols(_freq_meta=_freq_meta_expr)

    # Get counters for sex if set
    cut_dict = {
        cut: hl.agg.filter(
            hl.is_defined(mt._freq_meta[cut]), hl.agg.counter(mt._freq_meta[cut])
        )
        for cut in mt._freq_meta
        if cut != "subpop"
    }

    cut_data = mt.aggregate_cols(hl.struct(**cut_dict))
    sample_group_filters = []


    # Add all desired strata, starting with the full set and ending with downsamplings (if any)
    sample_group_filters = (
        [({}, True)]
        + [({"pop": pop}, mt._freq_meta.pop == pop) for pop in cut_data.get("pop", {})]
        + [({"sex": sex}, mt._freq_meta.sex == sex) for sex in cut_data.get("sex", {})]
        + [
            (
                {"pop": pop, "sex": sex},
                (mt._freq_meta.sex == sex) & (mt._freq_meta.pop == pop),
            )
            for sex in cut_data.get("sex", {})
            for pop in cut_data.get("pop", {})
        ]
        + [
            (
                {"subpop": subpop.subpop, "pop": subpop.pop},
                (mt._freq_meta.pop == subpop.pop)
                & (mt._freq_meta.subpop == subpop.subpop),
            )
            for subpop in cut_data.get("subpop", {})
        ]
        + [
            ({strata: str(s_value)}, mt._freq_meta[strata] == s_value)
            for strata in additional_strata_expr
            for s_value in cut_data.get(strata, {})
        ]
        + sample_group_filters
    )

    freq_sample_count = mt.aggregate_cols(
        [hl.agg.count_where(x[1]) for x in sample_group_filters]
    )

    # Annotate columns with group_membership
    mt = mt.annotate_cols(group_membership=[x[1] for x in sample_group_filters])

    # Create frequency expression array from the sample groups
    # Adding sample_group_filters_range_array to reduce memory usage in this array_agg
    mt = mt.annotate_rows(
        sample_group_filters_range_array=hl.range(len(sample_group_filters))
    )
    freq_expr = hl.agg.array_agg(
        lambda i: hl.agg.filter(
            mt.group_membership[i], hl.agg.call_stats(mt.GT, mt.alleles)
        ),
        mt.sample_group_filters_range_array,
    )

    # Select non-ref allele (assumes bi-allelic)
    freq_expr = freq_expr.map(
        lambda cs: cs.annotate(
            AC=cs.AC[1],
            AF=cs.AF[
                1
            ],  # TODO This is NA in case AC and AN are 0 -- should we set it to 0?
            homozygote_count=cs.homozygote_count[1],
        )
    )

    # Return MT with freq row annotation
    return mt.annotate_rows(freq=freq_expr).drop("_freq_meta")


SNV_mt_var_filtered = annotate_freq(
                SNV_mt_var_filtered,
                sex_expr=SNV_mt_var_filtered.sex,
            )

SNV_mt_var_filtered = SNV_mt_var_filtered.annotate_rows(
    info = SNV_mt_var_filtered.info.annotate(AC_tot_XX_XY=SNV_mt_var_filtered.freq.AC,
                                             AF_tot_XX_XY=SNV_mt_var_filtered.freq.AF,
                                             AN_tot_XX_XY=SNV_mt_var_filtered.freq.AN,
                                             hom_tot_XX_XY=SNV_mt_var_filtered.freq.homozygote_count)
                     )

SNV_mt_var_filtered = SNV_mt_var_filtered.annotate_rows(info=SNV_mt_var_filtered.info.drop('AF', "AC", "AN", "AQ"))

SNV_mt_var_filtered_no_geno = SNV_mt_var_filtered.rows()


#**Export files of interest**
#- Variants passing QC, with variant frequencies per sex (total, XX, XY) and individual genotypes
#        File name : SNV_filtered_with_geno.vcf.bgz
#- Variants passing QC, with variant frequencies per sex (total, XX, XY), without individual genotypes
#        File name : SNV_filtered_frequ_only.vcf.bgz


hl.export_vcf(SNV_mt_var_filtered, 'SNV_filtered_with_geno.vcf.bgz', tabix=True)

hl.export_vcf(SNV_mt_var_filtered_no_geno, 'SNV_filtered_frequ_only.vcf.bgz', tabix=True)
