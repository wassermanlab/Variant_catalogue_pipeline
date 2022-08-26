#!/usr/bin/env python
# coding: utf-8

# In[1]:


import hail as hl
from hail.plot import output_notebook, show
hl.init()
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
import sys
import string

from typing import Optional, Dict, List, Union


# #Created through the nextflow pipeline
 hl.import_vcf(sys.argv[1],array_elements_required=False, force_bgz=True).write('MEI_vcf.mt', overwrite=True)
 sex_table = (hl.import_table(sys.argv[2], impute=True).key_by('s'))


# In[2]:


#vcf_path = '/mnt/scratch/SILENT/Act3/Processed/Individual/GRCh37/Batch_DryRun_alpha/Version_0.0.1/MEI/'
#hl.import_vcf(os.path.join(vcf_path,'MEI_Version_0.0.1.vcf.gz'),
#              array_elements_required=False, force_bgz=True).write('MEI_vcf.mt', overwrite=True)
#sex_table = (hl.import_table('/mnt/scratch/SILENT/Act3/Processed/Workflow/Version_0.0.1/IBVL_pipeline/work/1e/6f3f64c0158271c52b26ca865a75d9/filtered_samples_sex.tsv', impute=True)
#         .key_by('s'))


# **Import file**
# 
# Inport a vcf file and read it as a matrix table (mt, hail specific file type)

# In[3]:


mt = hl.read_matrix_table('MEI_vcf.mt')


# In[4]:


mt = mt.annotate_cols(sample_name=mt.s.split('_')[0])


# In[5]:


mt = mt.annotate_cols(**sex_table[mt.sample_name])


# **Graph functions**
# 
# In order to create the graph, 3 functions were needed
# - stat : To calcualte the mean, standard deviation and other metrics for each parameter
# - plot_histo : To create the histogram as expected
# - plot_sp : To create the scatter plots as expected

# In[6]:


def stat(table):
    Mean = table[table.columns[2]]. mean()  
    StdDev = table[table.columns[2]]. std()
    Low_threashold = Mean - 3*StdDev
    High_threashold = Mean + 3*StdDev
    min_graph = table[table.columns[2]]. min() - 3*StdDev
    max_graph = table[table.columns[2]]. max() + 3*StdDev
    return Mean, StdDev, Low_threashold, High_threashold, min_graph, max_graph


# In[7]:


def plot_histo (table_plot, mt_plot, variable) :
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
    return show(p)


# In[8]:


#    output_file(filename=os.path.join("plot_hail", "Variant_QC",("variant_QC_"+variable+".html")), title="Variant QC HTML file")

    #output_file(filename=os.path.join("plot_hail", "Variant_QC",("variant_QC_"+x_variable+"_"+y_variable+".html")), title="Variant QC HTML file")



# In[9]:


def plot_sp (table_x_axis, mt_x_axis, table_y_axis, mt_y_axis, x_variable, y_variable) :
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
    return show(p)


# **SV QC**
# 
# SV : Exploration of the info column
# 
# - ASSESS (table)
# - internal (table)
# - SVTYPE (table)
# - SVLEN (histo)
# - LP (histo)
# - RP (histo)
# - SR (histo)
# - LP/RP (plot)
# 
# All info available in the vcf : 
# 
# ASSESS Number=1 Type=Integer Description= Provides information on evidence availible to decide insertion site.0 = No overlapping reads at site;1 = Imprecise breakpoint due to greater than expected distance between evidence;2 = discordant pair evidence only -- No split read information;3 = left side TSD evidence only;4 = right side TSD evidence only;5 = TSD decided with split reads  highest possible quality 
# 
# TSD Number=1 Type=String Description= Precise Target Site Duplication for bases  if unknown  value will be NULL 
# 
# INTERNAL Number=2 Type=String Description= If insertion internal or close to a gene  listed here followed by a discriptor of the location in the gene (either INTRON  EXON_   5_UTR  3_UTR  PROMOTER  or TERMINATOR). If multiple genes intersected  will be seperated by '\|' 
# 
# SVTYPE Number=1 Type=String Description= Type of structural variant 
# 
# SVLEN Number=1 Type=Integer Description= Difference in length between REF and ALT alleles; If unknown  will be -1 
# 
# MEINFO Number=4 Type=String Description= Mobile element info of the form NAME START END POLARITY; If START or END is unknown  will be -1; If
# 
# POLARITY is unknown  will be 'null' 
# 
# DIFF Number=. Type=String Description= Coverage and Differences in relation to the ALU reference. Form is %2XCoverage:Differences  with differences delimited by ' ' 
# 
# LP Number=1 Type=Integer Description= Total number of discordant pairs supporting the left side of the breakpont 
# 
# RP Number=1 Type=Integer Description= Total number of discordant pairs supporting the right side of the breakpont 
# 
# RA Number=1 Type=Float Description= Ratio between LP and RP  reported as log2(LP / RP) 
# 
# PRIOR Number=1 Type=String Description= True if this site was not discovered in this dataset  but was included on a provided priors list 
# 
# SR Number=1 Type=Integer Description= Total number of SRs at the estimated breakpoint for this site. Recomended to filter sites with <= 2 SRs 

# In[10]:


assess_table = mt.group_rows_by(mt.info.ASSESS).aggregate(n=hl.agg.count())


# In[11]:


len_hist = mt.aggregate_entries(hl.expr.aggregators.hist(mt.info.SVLEN,
            mt.aggregate_rows(hl.agg.min(mt.info.SVLEN)),
            mt.aggregate_rows(hl.agg.max(mt.info.SVLEN)), 100))
p = hl.plot.histogram(len_hist, legend='length (bp)', title='Length Histogram')
show(p)


# In[12]:


LP_hist = mt.aggregate_entries(hl.expr.aggregators.hist(mt.info.LP,
            mt.aggregate_rows(hl.agg.min(mt.info.LP)),
            mt.aggregate_rows(hl.agg.max(mt.info.LP)), 100))
p = hl.plot.histogram(LP_hist, legend='LP',
                      title='Total number of discordant pairs supporting the left side of the breakpoint')
show(p)


# In[13]:


RP_hist = mt.aggregate_entries(hl.expr.aggregators.hist(mt.info.RP,
            mt.aggregate_rows(hl.agg.min(mt.info.RP)),
            mt.aggregate_rows(hl.agg.max(mt.info.RP)), 100))
p = hl.plot.histogram(RP_hist, legend='RP',
                      title='Total number of discordant pairs supporting the right side of the breakpoint')
show(p)


# In[14]:


RA_hist = mt.aggregate_entries(hl.expr.aggregators.hist(mt.info.RA,
            mt.aggregate_rows(hl.agg.min(mt.info.RA)),
            mt.aggregate_rows(hl.agg.max(mt.info.RA)), 100))
p = hl.plot.histogram(RA_hist, legend='RA',
                      title='Ratio between LP and RP reported as log2(LP / RP)')
show(p)


# In[15]:


p = hl.plot.scatter(x=mt.info.RP,
                   y=mt.info.LP,
                  xlabel="RP : Total number of discordant pairs supporting the right side of the breakpoint",
                  ylabel="LP : Total number of discordant pairs supporting the left side of the breakpoint",
                title="RP/LP")
show(p)


# In[16]:


SR_hist = mt.aggregate_entries(hl.expr.aggregators.hist(mt.info.SR,
            mt.aggregate_rows(hl.agg.min(mt.info.SR)),
            mt.aggregate_rows(hl.agg.max(mt.info.SR)), 100))

p = hl.plot.histogram(SR_hist, legend='SR',
                      title='Total number of SRs at the estimated breakpoint for this site.')
annot = Span(dimension="height",location=2,line_dash='dashed', line_width=3,line_color="red")
p.add_layout(annot)
show(p)


# **MEI Hail QC**

# In[17]:


mt = hl.variant_qc(mt)


# In[18]:


#Replace NA and {} by "none" in the filter column

mt = mt.annotate_rows(filters = hl.if_else(hl.len(mt.filters) > 0,
                            mt.filters,
                            hl.set(["none"])))

mt = mt.annotate_rows(filters =hl.if_else(hl.is_defined(mt.filters),
                                                mt.filters,
                                                hl.set(["none"])))


# List of variant_qc variable of interest:
# - DP (mt.variant_qc.dp_stats.mean)
# - AN (mt.variant_qc.AN)
# - call_rate (mt.variant_qc.call_rate)
# - n_called (mt.variant_qc.n_called)
# - n_not_called (mt.variant_qc.n_not_called)
# - n_het (mt.variant_qc.n_het)
# - p_value_hwe (mt.variant_qc.p_value_hwe)
# - het_freq_hwe (mt.variant_qc.het_freq_hwe)
# 
# Other available
# 
# - AC
# - AF
# - Homozygote_count
# - n_called
# - n_non_ref
# - p_value_excess_het
# 
# Empty fields :
# - n_filtered (mt.variant_qc.n_filtered)

# In[19]:


mt.variant_qc.dp_stats.mean.export('DP_MEI.tsv')


# In[20]:


mt.variant_qc.AN.export('AN_MEI.tsv')
mt.variant_qc.call_rate.export('call_rate_MEI.tsv')
mt.variant_qc.n_called.export('n_called_MEI.tsv')
mt.variant_qc.n_not_called.export('n_not_called_MEI.tsv')
mt.variant_qc.n_het.export('n_het_MEI.tsv')
mt.variant_qc.p_value_hwe.export('p_value_hwe_MEI.tsv')
mt.variant_qc.het_freq_hwe.export('het_freq_hwe_MEI.tsv')


# In[21]:


DP_table=pd.read_table('DP_MEI.tsv')


# In[22]:


AN_table=pd.read_table('AN_SV.tsv')
call_rate_table=pd.read_table('call_rate_MEI.tsv')
n_called_table=pd.read_table('n_called_MEI.tsv')
n_not_called_table=pd.read_table('n_not_called_MEI.tsv')
n_het_table=pd.read_table('n_het_MEI.tsv')
p_value_hwe_table=pd.read_table('p_value_hwe_MEI.tsv')
het_freq_hwe_table=pd.read_table('het_freq_hwe_MEI.tsv')


# In[23]:


DP_table.rename(columns = {DP_table.columns[2]:'DP'}, inplace = True)
AN_table.rename(columns = {AN_table.columns[2]:"AN"}, inplace = True)
call_rate_table.rename(columns = {call_rate_table.columns[2]:'call_rate'}, inplace = True)
n_called_table.rename(columns = {n_called_table.columns[2]:"n_called"}, inplace = True)
n_not_called_table.rename(columns = {n_not_called_table.columns[2]:"n_not_called"}, inplace = True)
n_het_table.rename(columns = {n_het_table.columns[2]:"n_het"}, inplace = True)
p_value_hwe_table.rename(columns = {p_value_hwe_table.columns[2]:"p_value_hwe"}, inplace = True)
het_freq_hwe_table.rename(columns = {het_freq_hwe_table.columns[2]:"het_freq_hwe"}, inplace = True)


# In[24]:


plot_histo(DP_table,
           mt.variant_qc.dp_stats.mean,
           'Mean Depth per variant')


# In[25]:


plot_histo(AN_table,
           mt.variant_qc.AN,
           'Allele number per variant')


# In[26]:


plot_histo(call_rate_table,
           mt.variant_qc.call_rate,
           'Call rate per variant')


# In[27]:


plot_histo(n_called_table,
           mt.variant_qc.n_called,
           'Number of samples with a called genotype per variant')


# In[28]:


plot_histo(n_not_called_table,
           mt.variant_qc.n_not_called,
           'Number of samples with a missing genotype per variant')


# In[29]:


plot_histo(p_value_hwe_table,
           mt.variant_qc.p_value_hwe,
           'p-value from two-sided test of Hardy-Weinberg equilibrium per variant')


# In[30]:


plot_sp (het_freq_hwe_table,
         mt.variant_qc.het_freq_hwe,
         n_het_table,
         mt.variant_qc.n_het,
         'Expected frequency of heterozygous samples under Hardy-Weinberg equilibrium',
         'Number of heterozygous samples per variant')


# **Filter the variants**
# 
# Filter the variants
# 
# **Step 1 : chr 1-22, X, Y**
# 
# Remove the variants not located on chr 1-22, X/Y
# 
# **Step 2 : mt.filters column**
# 
# ==> Remove all the MEI with a non empty filters column
# 
# **Step 3 : ASSESS**
# 
# MELT ASSESS score >3
# 
# 
# **Step 4 : SR**
# 
# split reads >2
# 
# 
# **Sources :**
# 
# 
# In MELT filters (mt.filters column), remove:
# - lc : Rewmove low complexity (lc)
# - s25 : '> 25% missing (s25)
# - rSD : LP/RP > 2.0 standard deviation (rSD)
# - hDP : more discordant pairs are also split than expected (hDP)
# ==> Remove all the MEI with a non empty filters column
# from https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-16481-5/MediaObjects/41467_2020_16481_MOESM1_ESM.pdf table 1
# 
# 
# Variants that did not pass the following criteria were filtered to get a high-quality MEI call set:
# (i) not in low complexity regions; 
# (ii) be genotyped in >25.0% of individuals; 
# (iii) split reads >2;
# (iv) MELT ASSESS score >3 (at least one-side TSD evidence)
# (v) VCF FILTER column be PASS. 
#  from https://academic.oup.com/nar/article/50/5/2493/6536893
# 
# 
# 
# 
# Possible future improvments :
# Filtering of Candidate MEIs: MEIs were required to have at least four DPs of supporting evidence during initial discovery at each site for the final call set. This provided a good balance between the false negative rate (FNR) and the FDR (Table 1, main text). Putative MEIs were filtered if they mapped within reference mobile elements of the same type as annotated by RepeatMasker v. 4.0.3 at the University of California Santa Cruz (UCSC) Genome Browser website21. To control for sequence coverage variation at candidate MEI sites, 100 bp windows flanking each MEI site were sampled for depth of coverage fluctuations. Sites that fell outside of the range of 70 to 130% sequence coverage were filtered.
# from : Supplementary material for: An integrated map of structural variation in 2,504 human genomes

# **Step 1 : chr 1-22, X, Y**

# In[31]:


intervals = [hl.parse_locus_interval(x) for x in ['X', 'Y', '1-22']]
MEI_mt_locus_filtered = hl.filter_intervals(mt, intervals, keep=True)


# **Step 2 : mt.filters column**
# 
# **Step 3 : ASSESS**
# 
# **Step 4 : SR**

# **Step 5 : Hail filters**
#     
# - AN (Allele Number) : The SV with a very low AN were not genotyped correctly in the cohort and will not give a strong frequency estimation, should be removed.
# AN threshold : 
# Max AN = number of samples*2
# Min AN --> 75% of max AN
# 
# - Not implemented : DP, quality

# **AN (Allele Number)**
# 
# The SV with a very low AN were not genotyped correctly in the cohort and will not give a strong frequency estimation, should be removed.
# 
# AN threshold : 
# 
# Max AN = number of samples*2
# 
# Min AN --> 75% of max AN

# In[32]:


MEI_mt_locus_filtered = hl.variant_qc(MEI_mt_locus_filtered)


# In[33]:


min_AN=MEI_mt_locus_filtered.count()[1]*2*0.75


# In[34]:


MEI_mt_filtered = MEI_mt_locus_filtered.filter_rows(
    (MEI_mt_locus_filtered.filters == hl.set(["none"])) &
    (MEI_mt_locus_filtered.info.ASSESS >= 3) &
    (MEI_mt_locus_filtered.info.SR >= 2) &
    ((MEI_mt_locus_filtered.variant_qc.AN) >  min_AN)
)


# **Calculation of the stats post filtering for report**

# In[35]:


n_MEI_tot=mt.count()[0]

n_non_chr = n_MEI_tot - MEI_mt_locus_filtered.count()[0]
perc_non_chr = n_non_chr/n_MEI_tot * 100

n_filters_filtered = n_MEI_tot - mt.filter_rows(mt.filters == hl.set(["none"])).count()[0]
perc_filters_filtered = n_filters_filtered / n_MEI_tot *100

n_assess_filtered = n_MEI_tot - mt.filter_rows((mt.info.ASSESS >= 3)).count()[0]
perc_assess_filtered = n_assess_filtered / n_MEI_tot *100

n_SR_filtered = n_MEI_tot - mt.filter_rows((mt.info.SR >= 2)).count()[0]
perc_SR_filtered = n_SR_filtered / n_MEI_tot *100

n_AN_filtered = n_MEI_tot - mt.filter_rows((mt.variant_qc.AN) >  min_AN).count()[0]
perc_AN_filtered = n_AN_filtered / n_MEI_tot *100

n_MEI_removed = n_MEI_tot - MEI_mt_filtered.count()[0]
perc_MEI_removed = n_MEI_removed/n_MEI_tot * 100

n_MEI_tot_filtered = MEI_mt_filtered.count()[0]


# MEI filters : info from the filters column before filtration

# In[36]:


filters_ht = mt.rows()
filters_ht = filters_ht.select ('filters')
filters_df = filters_ht.to_pandas()
filters_table = filters_df['filters'].value_counts()
filters_table.to_csv(r'MEI_filters_report.txt', sep=':')


# MEI Type (before and after filtration)

# In[37]:


type_ht = mt.select_rows(mt.info.SVTYPE).rows()
type_df = type_ht.to_pandas()
n_ALU = len(type_df[type_df['SVTYPE'] == 'ALU'])
n_LINE1 = len(type_df[type_df['SVTYPE'] == 'LINE1'])
n_SVA = len(type_df[type_df['SVTYPE'] == 'SVA'])


# In[38]:


type_filtered_ht = MEI_mt_filtered.select_rows(MEI_mt_filtered.info.SVTYPE).rows()
type_filtered_df = type_filtered_ht.to_pandas()
n_ALU_filtered = len(type_filtered_df[type_filtered_df['SVTYPE'] == 'ALU'])
n_LINE1_filtered = len(type_filtered_df[type_filtered_df['SVTYPE'] == 'LINE1'])
n_SVA_filtered = len(type_filtered_df[type_filtered_df['SVTYPE'] == 'SVA'])


# MEI location (post-filters uniquely)
# 
# For internal table, need to ignore the gene name to know number by location type
# 
# (!) Not working for EXON

# In[39]:


n_internal_intron = MEI_mt_filtered.aggregate_rows(hl.agg.count_where(MEI_mt_filtered.info.INTERNAL.contains('INTRONIC')))
n_internal_exon = MEI_mt_filtered.aggregate_rows(hl.agg.count_where(MEI_mt_filtered.info.INTERNAL.contains('EXON')))
n_internal_5_UTR = MEI_mt_filtered.aggregate_rows(hl.agg.count_where(MEI_mt_filtered.info.INTERNAL.contains('5_UTR')))
n_internal_3_UTR = MEI_mt_filtered.aggregate_rows(hl.agg.count_where(MEI_mt_filtered.info.INTERNAL.contains('3_UTR')))
n_internal_promotor = MEI_mt_filtered.aggregate_rows(hl.agg.count_where(MEI_mt_filtered.info.INTERNAL.contains('PROMOTER')))
n_internal_terminator = MEI_mt_filtered.aggregate_rows(hl.agg.count_where(MEI_mt_filtered.info.INTERNAL.contains('TERMINATOR')))
n_internal_null = MEI_mt_filtered.aggregate_rows(hl.agg.count_where(MEI_mt_filtered.info.INTERNAL.contains('null')))


# **MEI QC Report**
# 
# Write the report of the number of filtered out variants and the reason they were filtered out

# In[40]:


def report_stats():
    """
    Generate output report with basic stats.
    """
    out_stats = hl.hadoop_open(f"MEI_QC_report.txt", "w")
    # Report numbers of filtered MEI
    out_stats.write(
        f"Before fiiltration \n"
        f"Number of MEI: {n_MEI_tot}\n"
        f"     Number of ALU: {n_ALU}\n"
        f"     Number of LINE1: {n_LINE1}\n"
        f"     Number of SVA: {n_SVA}\n\n\n"
        
        f"Variant filtering\n\n"
        f"Filer 1 : Location\n"        
        f"     Number of MEI removed because of location (out of chr1-22, X, Y) : {n_non_chr}\n"
        f"     Percentage of MEI removed because of location (out of chr1-22, X, Y) : {perc_non_chr} %\n"
        f"Filer 2 : Filters column\n"        
        f"     Number of MEI removed because of filters (non PASS) : {n_filters_filtered}\n"
        f"     Percentage of MEI removed because of filters (non PASS) : {perc_filters_filtered} %\n"
        f"Filer 3 : ASSESS value (ASSESS >= 3)\n"        
        f"     Number of MEI removed because of low ASSESS : {n_assess_filtered}\n"
        f"     Percentage of MEI removed because of low ASSESS : {perc_assess_filtered} %\n"
        f"Filer 4 : SR value (SR >= 2)\n"        
        f"     Number of MEI removed because of low SR : {n_SR_filtered}\n"
        f"     Percentage of MEI removed because of low SR : {perc_SR_filtered} %\n"
        f"Filer 5 : Hail QC, AN (Allele number)\n"
        f"     AN threshold for filtering : {min_AN}\n"        
        f"     Number of MEI removed because of AN threshold : {n_AN_filtered}\n"
        f"     Percentage of MEI removed because of AN threshold : {perc_AN_filtered} %\n"
        f"Conclusion\n"
        f"     Total number of MEI removed : {n_MEI_removed}\n"
        f"     Percentage of the MEI filtered out: {perc_MEI_removed} %\n\n\n"

        f"After variant filtering\n\n"
        f"Number of MEI remaining (all types, all callers): {n_MEI_tot_filtered}\n"
        f"     Number of ALU: {n_ALU_filtered}\n"
        f"     Number of LINE1: {n_LINE1_filtered}\n"
        f"     Number of SVA: {n_SVA_filtered}\n\n\n"
        f"Localisation of the MEI\n"
        f'     Intronic = {n_internal_intron}\n'
        f'     Exonic = NOT WORKING\n'
        f'     5 UTR = {n_internal_5_UTR}\n'
        f'     3 UTR = {n_internal_3_UTR}\n'
        f'     Promotor = {n_internal_promotor}\n'
        f'     Terminator = {n_internal_terminator}\n'
        f'     Null = {n_internal_null}'        
    )
    out_stats.close()


# In[41]:


report_stats()


# **Calculate sex specific frequencies**
# 
# Sex is defined using F-stat in Hail_sample_QC or file with sample sex can be loaded by user
# 
# Calculate AF, AC, AN and number of homozygotes
# 
# Code adapted from gnomAD : https://github.com/broadinstitute/gnomad_qc/blob/main/gnomad_qc/v3/annotations/generate_freq_data.py

# In[44]:


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




# In[60]:


MEI_mt_filtered_export = annotate_freq(
                MEI_mt_filtered,
                sex_expr=MEI_mt_filtered.sex,
            )


# In[61]:


MEI_mt_filtered_export = MEI_mt_filtered_export.annotate_rows(
    info = MEI_mt_filtered_export.info.annotate(AC_tot_XX_XY=MEI_mt_filtered_export.freq.AC,
                                             AF_tot_XX_XY=MEI_mt_filtered_export.freq.AF,
                                             AN_tot_XX_XY=MEI_mt_filtered_export.freq.AN,
                                             hom_tot_XX_XY=MEI_mt_filtered_export.freq.homozygote_count)
                     )


# In[62]:


MEI_mt_filtered_export = MEI_mt_filtered_export.rows()


# **Export files of interest**
# 
# - Variants passing QC, with variant frequencies per sex (total, XX, XY) and individual genotypes
# 
#         File name : MEI_filtered_with_geno.vcf.bgz
# 
# 
# 
# - Variants passing QC, with variant frequencies per sex (total, XX, XY), without individual genotypes
# 
# 
#         File name : MEI_filtered_frequ_only.vcf.bgz

# In[49]:


hl.export_vcf(MEI_mt_filtered_export, 'MEI_filtered_with_geno.vcf.bgz', tabix=True)


# In[50]:


hl.export_vcf(MEI_mt_filtered_export_no_geno, 'MEI_filtered_frequ_only.vcf.bgz', tabix=True)


# In[ ]:





# In[ ]:





# In[ ]:




