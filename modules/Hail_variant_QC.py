#!/usr/bin/env python
# coding: utf-8

# In[1]:


import hail as hl
from hail.plot import output_notebook, show
from bokeh.models import Span
from bokeh.models import Range1d
from bokeh.plotting import  output_file, save
import pandas as pd
import os
import sys
from typing import Optional, Dict, List
hl.plot.output_notebook()


# In[2]:


temp_directory=sys.argv[3]
genome=sys.argv[4] # assembly - GRCh37 or GRCh38
ref_fasta=sys.argv[5]
ref_fasta_index=sys.argv[6]
chr = sys.argv[7]
vcf_file=sys.argv[1]
sex_table=sys.argv[2]


# In[3]:

# In[4]:


# Hail and plot initialisation
# Configure Spark properties
spark_conf = {
    'spark.driver.memory': '8g'  # Set the driver memory, e.g., to 8 GB
}


# In[5]:


# Initialize Hail with custom Spark configuration
hl.init(master='local[*]', spark_conf=spark_conf)
output_notebook()


# In[6]:


# import vcf - sys_argv[1] is output from Hail_sample_QC.py
# created through the nextflow pipeline
# import with recoding for GRCh38 data as data follows GRCh37 labelling
# (no 'chr' prefix) and Hail requires the 'chr' prefix for GRCh38 data
if genome == 'GRCh37':
    mt = hl.import_vcf(vcf_file,
                array_elements_required=False,
                force_bgz=True,
                reference_genome=genome)

elif genome == 'GRCh38':

    recode = {"MT":"chrM", **{f"{i}":f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}}

    mt = hl.import_vcf(vcf_file,
                array_elements_required=False,
                force_bgz=True,
                reference_genome=genome,
                contig_recoding=recode)  
else:
    raise Exception("Must use either GRCh37 or GRCh38!")


# In[7]:


# Sex table (sys_argv[2] is the sample sex table created in Hail_sample_QC.py)
sex_table = (hl.import_table(sex_table, impute=True).key_by('s'))


# In[8]:


#Handle interval
if (chr == "autosomal"):
    #recode for hail functions to use hail GRCh38, select chromosomes 1-22
    if genome == "GRCh37":
        contigs = [f"{i}" for i in (list(range(1, 23)))]
    elif genome =="GRCh38":
        contigs = [f"chr{i}" for i in (list(range(1, 23)))]
    else:
        raise ValueError("must use a valid assembly - GRCh37 or GRCh38")

    intervals = [hl.parse_locus_interval(x, reference_genome=genome) for x in contigs]
    mt = hl.filter_intervals(mt, intervals, keep=True)
    
elif ((chr == "X") | (chr == "x") |(chr == "chrX") | (chr == "chrx")):
    if genome == "GRCh37":
        contigs = 'X'
    elif genome =="GRCh38":
        contigs = 'chrX'
        
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval(contigs, reference_genome=genome)], keep=True)

elif ((chr == "Y") | (chr == "y") | (chr == "chrY") | (chr == "chry")):
    if genome == "GRCh37":
        contigs = 'Y'
    elif genome =="GRCh38":
        contigs = 'chrY'
        
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval(contigs, reference_genome=genome)], keep=True)

else: 
    raise Exception("Invalid chromosome interval/contig in hail_variant_qc.py:" + str(chr))
    


# hl.import_vcf(sys.argv[1],array_elements_required=False, force_bgz=True, reference_genome=genome).write('filtered_samples_vcf.mt', overwrite=True)

# In[9]:


#-----------------XX Samples -------------------------------------------------------
# set genotypes on the Y chromosome to missing for all XX samples
# do prior to filtering samples for interval, so both X and Y are present
# done for Y chromosome QC only
#-----------------------------------------------------------------------------------
#make a Y contig
if((chr == "Y") | (chr == "y") | (chr == "chrY") | (chr == "chry")):
    if(genome == "GRCh37"):
        y_locus = hl.locus("Y", 1, reference_genome=genome)
    if(genome == "GRCh38"):
        y_locus = hl.locus("chrY", 1, reference_genome=genome)
        
    GT_count_before_XX = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
    XX_filter_count = mt.aggregate_entries(hl.agg.count_where((sex_table[mt.s].sex == "XX") & (mt.locus.contig == y_locus.contig) & hl.is_defined(mt.GT)))
   
    mt = mt.transmute_entries(GT = hl.if_else(((sex_table[mt.s].sex == "XX") & (mt.locus.contig == y_locus.contig)), hl.missing(hl.tcall), mt.GT))
   
    GT_count_after_XX = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
else: 
    XX_filter_count = 0
    GT_count_before_XX = 0
    GT_count_after_XX = 0


# In[10]:


#------------------variant QC (hail) & caching on disk--------------------------------------
mt = hl.variant_qc(mt)
mt.write(f"SNV_variant_QC_{chr}.mt", overwrite=True) 
mt = hl.read_matrix_table(f"SNV_variant_QC_{chr}.mt")


# **Import file**
# Inport a vcf file and read it as a matrix table (mt, hail specific file type)

# **Graph functions**
# In order to create the graph, 3 functions were needed
# - stat : To calcualte the mean, standard deviation and other metrics for each parameter
# - plot_histo : To create the histogram as expected
# - plot_sp : To create the scatter plots as expected

# In[11]:


def stat(table):
    Mean = table[table.columns[2]]. mean()  
    StdDev = table[table.columns[2]]. std()
    Low_threashold = Mean - 3*StdDev
    High_threashold = Mean + 3*StdDev
    min_graph = table[table.columns[2]]. min() - 3*StdDev
    max_graph = table[table.columns[2]]. max() + 3*StdDev
    return Mean, StdDev, Low_threashold, High_threashold, min_graph, max_graph


# In[12]:


def plot_histo (table_plot, mt_plot, variable) :
    output_file(filename=os.path.join(("variant_QC_"+variable+".html")), title=f"Variant QC HTML file {chr}")
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


# In[13]:


def plot_sp (table_x_axis, mt_x_axis, table_y_axis, mt_y_axis, x_variable, y_variable) :
    output_file(filename=os.path.join(("variant_QC_"+x_variable+"_"+y_variable+".html")), title=f"Variant QC HTML file {chr}")
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


# In[14]:


#Calculate AF, AC, AN and number of homozygotes
#Code adapted from gnomAD : https://github.com/broadinstitute/gnomad_qc/blob/main/gnomad_qc/v3/annotations/generate_freq_data.py


# In[15]:


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


# **Variant QC**

# List of variables for which we will create a table, calculate the standard deviation (StdDev) and the mean (Mean) for sample QC:
# - DP (mt_sample_qc.variant_qc.dp_stats.mean)
# - QG (mt_sample_qc.vaiant_qc.gq_stats.mean)
# - call_rate (mt_sample_qc.variant_qc.call_rate)
# - AN (mt_sample_qc.variant_qc.AN)
# - n_not_called (mt_sample_qc.variant_qc.n_not_called)
# - p_value_hwe (mt_sample_qc.variant_qc.p_value_hwe)
# - het_freq_hwe (mt_sample_qc.variant_qc.het_freq_hwe)
# - n_het (mt_sample_qc.variant_qc.n_het)

# In[16]:


mt.variant_qc.dp_stats.mean.export(f"DP_SNV_{chr}.tsv")
mt.variant_qc.gq_stats.mean.export(f"GQ_SNV_{chr}.tsv")
mt.variant_qc.call_rate.export(f"Call_Rate_SNV_{chr}.tsv")
mt.variant_qc.AN.export(f"AN_SNV_{chr}.tsv")
mt.variant_qc.n_not_called.export(f"n_notcalled_SNV_{chr}.tsv")
mt.variant_qc.p_value_hwe.export(f"p_value_HWE_SNV_{chr}.tsv")
mt.variant_qc.het_freq_hwe.export(f"het_freq_HWE_SNV_{chr}.tsv")
mt.variant_qc.n_het.export(f"n_het_SNV_{chr}.tsv")


# In[17]:


DP_SNV_table=pd.read_table(f"DP_SNV_{chr}.tsv")
GQ_SNV_table=pd.read_table(f"GQ_SNV_{chr}.tsv")
call_rate_SNV_table=pd.read_table(f"Call_Rate_SNV_{chr}.tsv")
AN_SNV_table=pd.read_table(f"AN_SNV_{chr}.tsv")
n_not_called_SNV_table=pd.read_table(f"n_notcalled_SNV_{chr}.tsv")
p_value_hwe_SNV_table=pd.read_table(f"p_value_HWE_SNV_{chr}.tsv")
het_freq_hwe_SNV_table=pd.read_table(f"het_freq_HWE_SNV_{chr}.tsv")
n_het_SNV_table=pd.read_table(f"n_het_SNV_{chr}.tsv")


# In[18]:


DP_SNV_table.rename(columns = {DP_SNV_table.columns[2]:'DP'}, inplace = True)
GQ_SNV_table.rename(columns = {GQ_SNV_table.columns[2]:'GQ'}, inplace = True)
call_rate_SNV_table.rename(columns = {call_rate_SNV_table.columns[2]:'call_rate'}, inplace = True)
AN_SNV_table.rename(columns = {AN_SNV_table.columns[2]:"AN"}, inplace = True)
n_not_called_SNV_table.rename(columns = {n_not_called_SNV_table.columns[2]:"n_not_called"}, inplace = True)
p_value_hwe_SNV_table.rename(columns = {p_value_hwe_SNV_table.columns[2]:"p_value_hwe"}, inplace = True)
het_freq_hwe_SNV_table.rename(columns = {het_freq_hwe_SNV_table.columns[2]:"het_freq_hwe"}, inplace = True)
n_het_SNV_table.rename(columns = {n_het_SNV_table.columns[2]:"n_het"}, inplace = True)


# In[19]:


plot_histo(DP_SNV_table,
           mt.variant_qc.dp_stats.mean,
           f"Mean Depth per variant {chr} - Unfiltered SNV")


# In[20]:


plot_histo(GQ_SNV_table,
           mt.variant_qc.gq_stats.mean,
           f"Mean Genotype Quality per variant {chr} - Unfiltered SNV")


# In[21]:


plot_histo(call_rate_SNV_table,
           mt.variant_qc.call_rate,
           f"Call rate per variant {chr} - Unfiltered SNV")


# In[22]:


plot_histo(AN_SNV_table,
           mt.variant_qc.AN,
           f"Allele number per variant {chr}- Unfiltered SNV")


# In[23]:


plot_histo(n_not_called_SNV_table,
           mt.variant_qc.n_not_called,
           f"Number of samples with a missing genotype per variant {chr} - Unfiltered SNV")


# In[24]:


plot_histo(p_value_hwe_SNV_table,
           mt.variant_qc.p_value_hwe,
           f"p-value from two-sided test of Hardy-Weinberg equilibrium per variant {chr}- Unfiltered SNV")


# In[25]:


plot_sp (het_freq_hwe_SNV_table,
         mt.variant_qc.het_freq_hwe,
         n_het_SNV_table,
         mt.variant_qc.n_het,
         f"Expected frequency of heterozygous samples under Hardy-Weinberg equilibrium {chr}",
         f"Number of heterozygous samples per variant {chr} - Unfiltered SNV")


# In[26]:


#---------------------- Variant Metrics and Filtering -----------------------------------------
# The previous method was deprecated and removed. It used 3-standard deviations for 
# hard filtering variants based on all of the hail variant_qc() metrics which was 
# arbitrary and not aligned with the actual quality of the sites. The plotting and 
# distributions has been kept, and the only filtering is for indels > 50bp, which
# are now also saved in their own vcf file for easy reference.
#
# The following will be reported: GQ within certain bands (10, 15, 20)
#
# The following will be filtered:
#  - AC0 : any site where no genotype is of high quality (GQ>=20, DP>=10 and allele balance > 0.2 for heterozygous)
#  - no variant: any site with an allele count (AC) of 0 - a result from sample filtering
#  - XX genotypes : any XX sample's genotypes on the Y chromosome
#


# In[27]:


# Large indels (> 50BP)

#make mt to export to vcf
large_indels_mt = mt.filter_rows((hl.len(mt.alleles[0]) > 50 ) | (hl.len(mt.alleles[1]) > 50), keep=True)

# note this 50 bp is the interval of the indel, rather than the size of the 
# number of bases inserted or deleted
n_large_del = mt.aggregate_rows(hl.agg.count_where((hl.len(mt.alleles[0]) > 50)))
n_large_ins = mt.aggregate_rows(hl.agg.count_where((hl.len(mt.alleles[1]) > 50)))

#count before
count_GT_before_indel_filter = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
count_rows_before_indel_filter =  mt.aggregate_rows(hl.agg.count())

#FILTER LARGE INDELS
mt = mt.filter_rows(((hl.len(mt.alleles[1]) > 50) | (hl.len(mt.alleles[0]) > 50)), keep=False)

#count after
total_GTs_after_indel = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
total_rows_after_indel = mt.aggregate_rows(hl.agg.count())

#count removed
count_indels_removed = count_rows_before_indel_filter - total_rows_after_indel
percent_indels_removed = (count_indels_removed / count_rows_before_indel_filter) * 100
count_indels_GT_removed = count_GT_before_indel_filter - total_GTs_after_indel


# In[28]:


#-------------- AC0 FILTER ----------------------------------------------------------
# No sample had a high quality genotype at this variant site  
# high quality definition: (GQ>=20, DP>=10 and allele balance > 0.2 for heterozygous)
# mark as "AC0" in filter, which will be filtered out in file for downstream processing
#---------------------------------------------------------------------------------------

#count before
total_GTs_before_AC0 = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
total_rows_before_AC0 = mt.aggregate_rows(hl.agg.count())

#FILTER
#mt = mt.filter_rows(mt.filters.contains("AC0"), keep=False)
mt = mt.filter_rows(~(hl.agg.any((mt.GQ>=20) & (mt.DP>=10) 
                              & (hl.if_else(mt.GT.is_het(),
                                (hl.if_else(((((mt.AD[1])/mt.DP) > 0.2) & (((mt.AD[1])/mt.DP) < 0.8)) ,True, False)),
                                True)))), keep=False)

#count after
total_GTs_after_AC0 = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
total_rows_after_AC0 = mt.aggregate_rows(hl.agg.count())

AC0_count = total_rows_before_AC0 - total_rows_after_AC0
AC0_percent = (AC0_count/count_rows_before_indel_filter)  * 100
AC0_GT_count = total_GTs_before_AC0 - total_GTs_after_AC0


# **---------------Calculate sex strata and stratified AC, AN, AF---------------------------------**
# #Sex is defined using F-stat in Hail_sample_QC or file with sample sex can be loaded by user

# In[29]:


SNV_mt_var_filtered = mt.annotate_cols(**sex_table[mt.s]).checkpoint(f"SNV_sex_annotation_{chr}.mt")


# In[30]:


SNV_mt_var_filtered = annotate_freq(
                SNV_mt_var_filtered,
                sex_expr=SNV_mt_var_filtered.sex,
            )


# In[31]:


SNV_mt_var_filtered = SNV_mt_var_filtered.annotate_rows(
    info = SNV_mt_var_filtered.info.annotate(AC_tot_XX_XY=SNV_mt_var_filtered.freq.AC,
                                             AF_tot_XX_XY=SNV_mt_var_filtered.freq.AF,
                                             AN_tot_XX_XY=SNV_mt_var_filtered.freq.AN,
                                             hom_tot_XX_XY=SNV_mt_var_filtered.freq.homozygote_count)
                     )


# In[32]:


SNV_mt_var_filtered = SNV_mt_var_filtered.annotate_rows(info=SNV_mt_var_filtered.info.drop('AF', "AC", "AN", "AQ"))


# In[33]:


#----------------NO VARIANT----------------------------------------------------------------
# This is caused by the filtering of samples (previous script) and genotypes.
# Remove any variants where the AC count became 0
#---------------------------------------------------------------------------------------

total_GTs_before_novariant = SNV_mt_var_filtered.aggregate_entries(hl.agg.count_where(hl.is_defined(SNV_mt_var_filtered.GT)))
total_rows_before_novariant = SNV_mt_var_filtered.aggregate_rows(hl.agg.count())

#filter the sites with no alt call (after stratification step, index 0 is alt allele for whole population)
SNV_mt_var_filtered = SNV_mt_var_filtered.filter_rows(SNV_mt_var_filtered.info.AC_tot_XX_XY[0] == 0, keep=False)

total_GTs_after_novariant = SNV_mt_var_filtered.aggregate_entries(hl.agg.count_where(hl.is_defined(SNV_mt_var_filtered.GT)))
total_rows_after_novariant = SNV_mt_var_filtered.aggregate_rows(hl.agg.count())

novar_count = total_rows_before_novariant - total_rows_after_novariant
novar_percent = (novar_count/count_rows_before_indel_filter * 100)
novar_GT_count = total_GTs_before_novariant - total_GTs_after_novariant


# In[34]:


#---------------GENOTYPE QUALITY METRICS -------------------------------------------
# --> Output metrics on genotype quality.
# Currently there is no filtering being done at the GT level.
# Only the AC0 filter accounts for it at the SITE level
#-------------------------------------------------------------------------------------
total_GTs = SNV_mt_var_filtered.aggregate_entries(hl.agg.count_where(hl.is_defined(SNV_mt_var_filtered.GT)))

GQ_10= mt.aggregate_entries(hl.agg.count_where((hl.is_defined(mt.GT)) & (mt.GQ < 10 )))
GQ_15= mt.aggregate_entries(hl.agg.count_where((hl.is_defined(mt.GT)) & (mt.GQ < 15 )))
GQ_20= mt.aggregate_entries(hl.agg.count_where((hl.is_defined(mt.GT)) & (mt.GQ < 20 )))


p20 = (GQ_20/total_GTs) * 100
p15 = (GQ_15/total_GTs) * 100
p10 = (GQ_10/total_GTs) * 100


# In[35]:


SNV_mt_var_filtered_no_geno = SNV_mt_var_filtered.rows()


# In[36]:


total_filtered = novar_count + AC0_count + count_indels_removed 
percent_filtered = (total_filtered/count_rows_before_indel_filter) * 100
total_after_strata = SNV_mt_var_filtered.aggregate_rows(hl.agg.count())

#----------------------REPORT-------------------------------------------------
out_stats = hl.hadoop_open(f"SNV_variant_QC_report_{chr}.txt", "w")
out_stats.write(
       
        f"\n\n-----------------------VARIANT QC REPORT--------------------\n"
         f"------------------------------------------------------------\n\n\n"


        f"\n------------------------------------------------------------\n"
        f"\nSUMMARY\n"
        f"\n (whole dataset)\n"
        f"------------------------------------------------------------\n"
        f"\nTotal sites before filtering: {count_rows_before_indel_filter}\n"
        f"\nTotal called genotypes (before filters): {count_GT_before_indel_filter}\n\n\n"

    
        f"\n------------------------------------------------------------\n"
        f"Indels > 50 bp\n"
        f"Indels > 50bp are filtered (but saved in a separate vcf) to"
        f"\navoid overlap with the results of the SV pipeline\n"
        f"------------------------------------------------------------\n"
        f"Number of indels filtered because of deletion larger than 50bp: {n_large_del}\n"
        f"Number of indels filtered because of insertion larger than 50bp: {n_large_ins}\n\n"

        f"Combined large indels sites filtered: {count_indels_removed}\n"
        f"Combined percent large indels sites filtered: {percent_indels_removed}\n"
        f"Number of large indels genotypes filtered: {count_indels_GT_removed}\n\n"

        f"Number of rows before indel filter: {count_rows_before_indel_filter }\n"
        f"Number of GTs before indel filter: {count_GT_before_indel_filter}\n"
        f"Number of rows after indel filter: {total_rows_after_indel}\n"
        f"Number of GTs after indel filter: {total_GTs_after_indel}\n\n\n"



        f"\n------------------------------------------------------------\n"
        f"\n\nGENOTYPE QUALITIES\n" 
        f"------------------------------------------------------------\n"
        f"Number of genotypes with GQ < 20: {GQ_20},  percent: {p20} \n"
        f"Number of genotypes with GQ < 15: {GQ_15},  percent: {p15} \n"
        f"Number of genotypes with GQ < 10: {GQ_10},  percent: {p10}\n\n\n"


        f"\n------------------------------------------------------------\n"
        f"Sites with AC0\n"
        f"-No sample had a high quality genotype at this variant site\n"
        f"-High quality definition: (GQ>=20, DP>=10 and allele balance > 0.2 for heterozygotes)\n"
        f"------------------------------------------------------------\n"
        f"\nAC0 sites (filtered): {AC0_count}\n"
        f"\nAC0 percentage of sites (filtered): {AC0_percent}\n"
        f"\nAC0 genotypes filtered: {AC0_GT_count}\n\n\n"

        f"Number of genotypes before AC0 filter: {total_GTs_before_AC0}\n"
        f"Number of rows before AC0 filter: {total_rows_before_AC0}\n\n"

        f"Number of genotypes after AC0 filter: {total_GTs_after_AC0}\n"
        f"Number of rows after AC0 filter: {total_rows_after_AC0}\n\n\n"

        f"------------------------------------------------------------\n"
        f"Sites with No Variant\n"
        f"Mostly from rare variants that belonged only to filtered samples,"
        f"so AC becomes 0\n"
        f"------------------------------------------------------------\n"
        f"\nSites filtered with no variant: {novar_count}\n"
        f"\nPercentage of sites with no variant: {novar_percent}\n"
        f"\nNumber of genotypes filtered with no variant: {novar_GT_count}\n\n"

        
        f"Number of genotypes before no-variant filter: {total_GTs_before_novariant}\n"
        f"Number of rows before no-variant filter: {total_rows_before_novariant}\n\n\n"


        f"Number of genotypes after no-variant filter: {total_GTs_after_novariant}\n"
        f"Number of rows after no-variant filter: {total_rows_after_novariant}\n\n\n"

        f"\n------------------------------------------------------------\n"
        f"\nXX Filtering (Y-chromosome only)\n"
        f"------------------------------------------------------------\n"

        f"\nGenotypes on Y chromosome before XX filtering: {GT_count_before_XX}\n"
        f"\nGenotypes filtered for XX samples on Y chromosome: {XX_filter_count}\n"
        f"\nGenotypes on Y chromosome after XX filtering: {GT_count_after_XX}\n"

        f"\n------------------------------------------------------------\n"
        f"FILTERING SUMMARY (combined)\n"
        f"Combined number of variants filered for: no variant, indel > 50bp, or for AC0\n"
        f"------------------------------------------------------------\n"
        f"\nTotal variants (sites) filtered: {total_filtered}\n"
        f"\nTotal sites before filtering: {count_rows_before_indel_filter}\n"
        f"\nTotal stratified sites after filtering: {total_after_strata}\n"
        f"\nPercentage filtered: {percent_filtered}\n\n\n"

    )
out_stats.close()


# **Export files of interest**
# - Variants passing QC, with variant frequencies per sex (total, XX, XY) and individual genotypes
#        File name : SNV_filtered_with_geno.vcf.bgz
# - Variants passing QC, with variant frequencies per sex (total, XX, XY), without individual genotypes
#        File name : SNV_filtered_frequ_only.vcf.bgz
# - A VCF with all the indels > 50bp that pass QC
# 
# - ** There will be one autosomal file (with genotypes)

# In[37]:


#create an autosome file for autosomal variants (mostly for analysis, does not get used downstream in pipeline)
if(chr == "autosomal"):
    hl.export_vcf(SNV_mt_var_filtered, f'SNV_autosomal.vcf.bgz', tabix=True)

#split by chromosome for downstream use in pipeline, if autosomal, otherwise it will be 
#already split for chrX and chrY
    for c in contigs :
        interval = [hl.parse_locus_interval(c, reference_genome=genome)]
        contig_mt = hl.filter_intervals(SNV_mt_var_filtered, interval, keep=True)
        contig_mt_no_geno = hl.filter_intervals(SNV_mt_var_filtered_no_geno, interval, keep=True)
        hl.export_vcf(contig_mt, f'SNV_filtered_with_geno_{c}.vcf.bgz', tabix=True)
        hl.export_vcf(contig_mt_no_geno, f'SNV_filtered_frequ_only_{c}.vcf.bgz',tabix=True)

else:
    hl.export_vcf(SNV_mt_var_filtered, f'SNV_filtered_with_geno_{chr}.vcf.bgz', tabix=True)
    hl.export_vcf(SNV_mt_var_filtered_no_geno, f'SNV_filtered_frequ_only_{chr}.vcf.bgz',
     tabix=True)


# In[38]:


#save file with indels >50 bp that passed quality control
hl.export_vcf(large_indels_mt, f'SNV_large_indels_{chr}.vcf.bgz', tabix=True)


# In[ ]:

