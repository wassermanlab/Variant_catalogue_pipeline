#!/usr/bin/env python
# coding: utf-8

# Hail and plot initialisation

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

#Created through the nextflow pipeline
hl.import_vcf(sys.argv[1],array_elements_required=False, force_bgz=True).write('filtered_samples_vcf.mt', overwrite=True)
sex_table = (hl.import_table(sys.argv[2], impute=True).key_by('s'))


hl.import_vcf(sys.argv[1],
              array_elements_required=False, force_bgz=True).write('filtered_samples_vcf.mt', overwrite=True)



SNV_mt = hl.read_matrix_table('filtered_samples_vcf.mt')


# In order to create the graph, 3 functions were needed
# - stat : To calcualte the mean, standard deviation and other metrics for each parameter
# - plot_histo : To create the histogram as expected
# - plot_sp : To create the scatter plots as expected

# In[5]:


def stat(table):
    Mean = table[table.columns[2]]. mean()  
    StdDev = table[table.columns[2]]. std()
    Low_threashold = Mean - 3*StdDev
    High_threashold = Mean + 3*StdDev
    min_graph = table[table.columns[2]]. min() - 3*StdDev
    max_graph = table[table.columns[2]]. max() + 3*StdDev
    return Mean, StdDev, Low_threashold, High_threashold, min_graph, max_graph


# In[6]:


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


# In[7]:


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

# In[8]:


SNV_mt = hl.variant_qc(SNV_mt)

# List of variables for which we will create a table, calculate the standard deviation (StdDev) and the mean (Mean) for sample QC:
# - DP (mt_sample_qc.variant_qc.dp_stats.mean)
# - QG (mt_sample_qc.vaiant_qc.gq_stats.mean)
# - call_rate (mt_sample_qc.variant_qc.call_rate)
# - AN (mt_sample_qc.variant_qc.AN)
# - n_not_called (mt_sample_qc.variant_qc.n_not_called)
# - p_value_hwe (mt_sample_qc.variant_qc.p_value_hwe)
# - het_freq_hwe (mt_sample_qc.variant_qc.het_freq_hwe)
# - n_het (mt_sample_qc.variant_qc.n_het)

# In[11]:


SNV_mt.variant_qc.dp_stats.mean.export('DP_variant.tsv')
SNV_mt.variant_qc.gq_stats.mean.export('GQ_variant.tsv')
SNV_mt.variant_qc.call_rate.export('call_rate_variant.tsv')
SNV_mt.variant_qc.AN.export('AN_variant.tsv')
SNV_mt.variant_qc.n_not_called.export('n_not_called_variant.tsv')
SNV_mt.variant_qc.p_value_hwe.export('p_value_hwe_variant.tsv')
SNV_mt.variant_qc.het_freq_hwe.export('het_freq_hwe_variant.tsv')
SNV_mt.variant_qc.n_het.export('n_het_variant.tsv')


# In[13]:


DP_variant_table=pd.read_table('DP_variant.tsv')
GQ_variant_table=pd.read_table('GQ_variant.tsv')
call_rate_variant_table=pd.read_table('call_rate_variant.tsv')
AN_variant_table=pd.read_table('AN_variant.tsv')
n_not_called_variant_table=pd.read_table('n_not_called_variant.tsv')
p_value_hwe_variant_table=pd.read_table('p_value_hwe_variant.tsv')
het_freq_hwe_variant_table=pd.read_table('het_freq_hwe_variant.tsv')
n_het_variant_table=pd.read_table('n_het_variant.tsv')


# In[15]:


DP_variant_table.rename(columns = {DP_variant_table.columns[2]:'DP'}, inplace = True)
GQ_variant_table.rename(columns = {GQ_variant_table.columns[2]:'GQ'}, inplace = True)
call_rate_variant_table.rename(columns = {call_rate_variant_table.columns[2]:'call_rate'}, inplace = True)
AN_variant_table.rename(columns = {AN_variant_table.columns[2]:"AN"}, inplace = True)
n_not_called_variant_table.rename(columns = {n_not_called_variant_table.columns[2]:"n_not_called"}, inplace = True)
p_value_hwe_variant_table.rename(columns = {p_value_hwe_variant_table.columns[2]:"p_value_hwe"}, inplace = True)
het_freq_hwe_variant_table.rename(columns = {het_freq_hwe_variant_table.columns[2]:"het_freq_hwe"}, inplace = True)
n_het_variant_table.rename(columns = {n_het_variant_table.columns[2]:"n_het"}, inplace = True)


# In[17]:


plot_histo(DP_variant_table,
           SNV_mt.variant_qc.dp_stats.mean,
           'Mean Depth per variant')


# In[18]:


plot_histo(GQ_variant_table,
           SNV_mt.variant_qc.gq_stats.mean,
           'Mean Genotype Quality per variant')


# In[19]:


plot_histo(call_rate_variant_table,
           SNV_mt.variant_qc.call_rate,
           'Call rate per variant')


# In[20]:


plot_histo(AN_variant_table,
           SNV_mt.variant_qc.AN,
           'Allele number per variant')


# In[21]:


plot_histo(n_not_called_variant_table,
           SNV_mt.variant_qc.n_not_called,
           'Number of samples with a missing genotype per variant')


# In[23]:


plot_histo(p_value_hwe_variant_table,
           SNV_mt.variant_qc.p_value_hwe,
           'p-value from two-sided test of Hardy-Weinberg equilibrium per variant')


# In[24]:


plot_sp (het_freq_hwe_variant_table,
         SNV_mt.variant_qc.het_freq_hwe,
         n_het_variant_table,
         SNV_mt.variant_qc.n_het,
         'Expected frequency of heterozygous samples under Hardy-Weinberg equilibrium',
         'Number of heterozygous samples per variant')


# Filter the samples based on the threasholds repersented on the figures 
# 
# Low_threashold = Mean - 3*StdDev = stat(table) [2]
# 
# High_threashold = Mean + 3*StdDev = stat(table) [3]
# 
# Filters :
# - Mean DP per variant lower than the low threshoold
# - Mean Genotype Quality per variant lower than the low threshoold
# - Call rate per variant lower than the low threshoold
# - Allele Number (AN) per variant lower than the low threshoold
# - Number of samples with missing genotype per variant lower than the low threshoold
# 
# ?? Hardy–Weinberg values

#Filter the small insertions / deletions (indel) of length > 50bp (Should be called by the SV pipeline)

# In[25]:

SNV_mt_var_filtered = SNV_mt.filter_rows(
    (SNV_mt.variant_qc.dp_stats.mean > stat(DP_variant_table) [2]) &
    (SNV_mt.variant_qc.gq_stats.mean > stat(GQ_variant_table) [2]) &
    (SNV_mt.variant_qc.call_rate > stat(call_rate_variant_table) [2]) &
    (SNV_mt.variant_qc.AN > stat(AN_variant_table) [2]) &
    (SNV_mt.variant_qc.n_not_called > stat(n_not_called_variant_table) [2]) &
    (hl.len(SNV_mt.alleles[0]) < 50) &
    (hl.len(SNV_mt.alleles[1]) < 50)
)

# In[55]:


hl.export_vcf(SNV_mt_var_filtered, 'SNV_filtered_samples_variants.vcf.bgz', tabix=True)

#Write the report of the number of filtered out variants and the reason they were filtered out

n_large_del = SNV_mt.filter_rows(hl.len(SNV_mt.alleles[0]) > 50).count()[0]
n_large_ins = SNV_mt.filter_rows(hl.len(SNV_mt.alleles[1]) > 50).count()[0]

def calc_removed_variant(mt, mt_var, stat_table) :
    input_mt = mt.annotate_rows(
        keep=(mt_var > stat_table [2]))

    n_removed = input_mt.aggregate_rows(hl.agg.count_where(~input_mt.keep))

    return n_removed

def report_stats():
    """
    Generate output report with basic stats.
    """
    out_stats = hl.hadoop_open(f"variant_QC.txt", "w")
    # Report numbers of filtered samples
    out_stats.write(
        f"Number of variants removed because of deletion superior to 50bp: {n_large_del}\n"
        f"Number of variants removed because of insertion superior to 50bp: {n_large_ins}\n"
        f"Number of variants removed because of depth metrics: {DP_var_removed}\n"
        f"Number of variants removed because of genotype quality metrics: {GQ_var_removed}\n"
        f"Number of variants removed because of call rate metrics: {CR_var_removed}\n"
        f"Number of variants removed because of allele number (AN): {n_AN_removed}\n"
        f"Number of variants removed because of number of not called: {n_not_called_removed}\n"
        f"Total number of variants removed : {n_var_removed}\n"
        f"Percentage of the variants filtered out: {perc_removed_variants}\n"
    )
    out_stats.close()

DP_var_removed = calc_removed_variant(SNV_mt, SNV_mt.variant_qc.dp_stats.mean, stat(DP_variant_table))
GQ_var_removed = calc_removed_variant(SNV_mt, SNV_mt.variant_qc.gq_stats.mean, stat(GQ_variant_table))
CR_var_removed = calc_removed_variant(SNV_mt, SNV_mt.variant_qc.call_rate, stat(call_rate_variant_table))
n_AN_removed = calc_removed_variant(SNV_mt, SNV_mt.variant_qc.AN, stat(AN_variant_table))
n_not_called_removed = calc_removed_variant(SNV_mt, SNV_mt.variant_qc.n_not_called, stat(n_not_called_variant_table))
n_var_removed = (SNV_mt.count()[0]-SNV_mt_var_filtered.count()[0])
perc_removed_variants = (SNV_mt.count()[0]-SNV_mt_var_filtered.count()[0])/SNV_mt.count()[0] * 100

report_stats()




# **Calculate statistic**
# 
# Calculate AF, AC, AN and numb of homozygotes (Total) 

# Sex is defined using F-stat in Hail_sample_QC or file with sample sex can be loaded by user

# Initially calculation is done for total, then XX, then XY specific frequencies are added to the vcf info tab


SNV_mt_var_filtered = hl.variant_qc(SNV_mt_var_filtered)
SNV_mt_var_filtered_tot = SNV_mt_var_filtered.annotate_rows(
    info = SNV_mt_var_filtered.info.annotate(
        chrom=SNV_mt_var_filtered.rsid.split('_')[0],
        pos=SNV_mt_var_filtered.rsid.split('_')[1],
        ref=SNV_mt_var_filtered.rsid.split('_')[2],
        alt=SNV_mt_var_filtered.rsid.split('_')[3],        
        qual=SNV_mt_var_filtered.qual,
        af_tot=SNV_mt_var_filtered.variant_qc.AF[1],
        ac_tot=SNV_mt_var_filtered.variant_qc.AC[1],
        an_tot=SNV_mt_var_filtered.variant_qc.AN,
        hom_alt_tot=SNV_mt_var_filtered.variant_qc.homozygote_count[1],
    )
) 

# Calculate the variants frequency per sex
# 

SNV_mt_var_filtered_sex = SNV_mt_var_filtered_tot.annotate_cols(**sex_table[SNV_mt_var_filtered_tot.s])
SNV_mt_var_filtered_XX = SNV_mt_var_filtered_sex.filter_cols(SNV_mt_var_filtered_sex.sex == 'XX')
SNV_mt_var_filtered_XX = hl.variant_qc(SNV_mt_var_filtered_XX)
SNV_mt_var_filtered_XX = SNV_mt_var_filtered_XX.annotate_rows(
    info = SNV_mt_var_filtered_XX.info.annotate(
        af_xx=SNV_mt_var_filtered_XX.variant_qc.AF[1],
        ac_xx=SNV_mt_var_filtered_XX.variant_qc.AC[1],
        an_xx=SNV_mt_var_filtered_XX.variant_qc.AN,
        hom_alt_xx=SNV_mt_var_filtered_XX.variant_qc.homozygote_count[1],
    )
)

SNV_mt_var_filtered_XY = SNV_mt_var_filtered_sex.filter_cols(SNV_mt_var_filtered_sex.sex == 'XY')
SNV_mt_var_filtered_XY = hl.variant_qc(SNV_mt_var_filtered_XY)
SNV_mt_var_filtered_XY = SNV_mt_var_filtered_XY.annotate_rows(
    info = SNV_mt_var_filtered_XY.info.annotate(
        af_xy=SNV_mt_var_filtered_XY.variant_qc.AF[1],
        ac_xy=SNV_mt_var_filtered_XY.variant_qc.AC[1],
        an_xy=SNV_mt_var_filtered_XY.variant_qc.AN,
        hom_alt_xy=SNV_mt_var_filtered_XY.variant_qc.homozygote_count[1],
    )
)

# **Save version without individual genotype**
# 
# As the frequencies are calculated using hail, it is not necessary to export the individual genotype for the last part of the pipeline (Annotation)

#Vcf : saved separately for total, XX and XY For tsv : Saved as one file with all the info

# Exporting the vcf with the varaint frequencies (No individual genotype) and xx samples frequencies (No XY for now)

SNV_mt_var_filtered_XX_export = SNV_mt_var_filtered_XX.drop('variant_qc')
SNV_mt_var_filtered_XX_export = SNV_mt_var_filtered_XX_export.annotate_rows(info=SNV_mt_var_filtered_XX_export.info.drop('AF', "AQ", "AC", "AN", "chrom", "pos", "ref", "alt"))
SNV_mt_var_filtered_XX_export = SNV_mt_var_filtered_XX_export.rows()
hl.export_vcf(SNV_mt_var_filtered_XX_export, 'SNV_filtered_frequ_total_xx.vcf.bgz', tabix=True)


# Exporting the tsv file with total, XX and XY frequencies
SNV_mt_var_filtered_tot_info = SNV_mt_var_filtered_tot.select_rows(
    SNV_mt_var_filtered_tot.info.chrom,
    SNV_mt_var_filtered_tot.info.pos,
    SNV_mt_var_filtered_tot.info.ref,
    SNV_mt_var_filtered_tot.info.alt,
    SNV_mt_var_filtered_tot.qual,
    SNV_mt_var_filtered_tot.rsid,
    SNV_mt_var_filtered_tot.info.af_tot,
    SNV_mt_var_filtered_tot.info.ac_tot,
    SNV_mt_var_filtered_tot.info.an_tot,
    SNV_mt_var_filtered_tot.info.hom_alt_tot,
).rows()

SNV_mt_var_filtered_XX_info = SNV_mt_var_filtered_XX.select_rows(
    SNV_mt_var_filtered_XX.info.af_xx,
    SNV_mt_var_filtered_XX.info.ac_xx,
    SNV_mt_var_filtered_XX.info.an_xx,
    SNV_mt_var_filtered_XX.info.hom_alt_xx,
).rows()

SNV_mt_var_filtered_XY_info = SNV_mt_var_filtered_XY.select_rows(
    SNV_mt_var_filtered_XY.info.af_xy,
    SNV_mt_var_filtered_XY.info.ac_xy,
    SNV_mt_var_filtered_XY.info.an_xy,
    SNV_mt_var_filtered_XY.info.hom_alt_xy,
).rows()

SNV_mt_var_filtered_tot_XX_info = SNV_mt_var_filtered_tot_info.join(SNV_mt_var_filtered_XX_info, how='left')
SNV_mt_var_filtered_tot_XX_XY_info = SNV_mt_var_filtered_tot_XX_info.join(SNV_mt_var_filtered_XY_info, how='left')
SNV_mt_var_filtered_tot_XX_XY_info.export('SNV_filtered_tot_XX_XY.tsv.bgz')








