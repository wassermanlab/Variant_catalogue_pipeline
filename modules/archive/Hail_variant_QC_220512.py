#!/usr/bin/env python
# coding: utf-8

# Hail and plot initialisation

# In[1]:


import hail as hl
from hail.plot import output_notebook, show
hl.init()
output_notebook()

import sys

# In[2]:


from hail.plot import show
from pprint import pprint
from bokeh.models import Span
hl.plot.output_notebook()
from bokeh.models import Range1d
from bokeh.plotting import figure, output_file, show, save

import pandas as pd
import os


# Import a vcf file and read it as a matrix table (mt, hail specific file type)
# For specific on how to look at the mt file, refer to the bottom of this Jupyter notebook)
# 
# For now : Import only the SNV vcf file, following the sample filtering step from previously
# Later : Load in the 5 vcf and apply the sample filtering on all of them

# In[ ]:


hl.import_vcf(sys.argv[1],
              array_elements_required=False, force_bgz=True).write('filtered_samples_vcf.mt', overwrite=True)

#hl.import_vcf('vcf_to_try_hail/filtered_samples.vcf.bgz'),
#              array_elements_required=False, force_bgz=True).write('filtered_samples_vcf.mt', overwrite=True)

# vcf_path = '/mnt/scratch/SILENT/Act3/Processed/Individual/GRCh37/Batch_DryRun/Run_20220426/vcf_pre_hail/'
# 
# hl.import_vcf(os.path.join(vcf_path,'DeepVariant_GLnexus_Run_20220426.vcf.gz'),
#               array_elements_required=False, force_bgz=True).write('hail/SNV_vcf.mt', overwrite=True)
# 
# hl.import_vcf(os.path.join(vcf_path,'MEI_Run_20220426.vcf.gz'),
#               array_elements_required=False, force_bgz=True).write('hail/MEI_vcf.mt', overwrite=True)
#               
# hl.import_vcf(os.path.join(vcf_path,'STR_Run_20220426.vcf.gz'),
#               array_elements_required=False, force_bgz=True).write('hail/STR_vcf.mt', overwrite=True)
#               
# hl.import_vcf(os.path.join(vcf_path,'SV_Run_20220426.vcf.gz'),
#               array_elements_required=False, force_bgz=True).write('hail/SV_vcf.mt', overwrite=True)

# For MT, need to find a different technique as it is not a diploiod genome
# 
# 
# Error summary: VCFParseError: ploidy > 2 not supported

# hl.import_vcf(os.path.join(vcf_path,'MT_Run_20220426.vcf.gz'), 
#               array_elements_required=False, force_bgz=True,
#               reference_genome='GRCh38').write('hail/MT_vcf.mt', overwrite=True)

# In[4]:


SNV_mt = hl.read_matrix_table('filtered_samples_vcf.mt')


# MEI_mt = hl.read_matrix_table('hail/MEI_vcf.mt')
# STR_mt = hl.read_matrix_table('hail/STR_vcf.mt')
# SV_mt = hl.read_matrix_table('hail/SV_vcf.mt')

# The vcf should be merged into one vcf to avoid redundancy (Possible calling overlap for indel witht he deepvaraint and SV pipeline)

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
    output_file(filename=("variant_QC_"+variable+".html"), title="Variant QC HTML file")
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
    output_file(filename=("variant_QC_"+x_variable+"_"+y_variable+".html"), title="Variant QC HTML file")
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


SNV_mt_variant_qc = hl.variant_qc(SNV_mt)


# In[9]:


SNV_mt_variant_qc.s


# In[10]:


SNV_mt_variant_qc.describe()


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


SNV_mt_variant_qc.variant_qc.dp_stats.mean.export('vcf_to_try_hail/11samples/DP_variant.tsv')


# In[12]:


SNV_mt_variant_qc.variant_qc.gq_stats.mean.export('vcf_to_try_hail/11samples/GQ_variant.tsv')
SNV_mt_variant_qc.variant_qc.call_rate.export('vcf_to_try_hail/11samples/call_rate_variant.tsv')
SNV_mt_variant_qc.variant_qc.AN.export('vcf_to_try_hail/11samples/AN_variant.tsv')
SNV_mt_variant_qc.variant_qc.n_not_called.export('vcf_to_try_hail/11samples/n_not_called_variant.tsv')
SNV_mt_variant_qc.variant_qc.p_value_hwe.export('vcf_to_try_hail/11samples/p_value_hwe_variant.tsv')
SNV_mt_variant_qc.variant_qc.het_freq_hwe.export('vcf_to_try_hail/11samples/het_freq_hwe_variant.tsv')
SNV_mt_variant_qc.variant_qc.n_het.export('vcf_to_try_hail/11samples/n_het_variant.tsv')


# In[13]:


DP_variant_table=pd.read_table('vcf_to_try_hail/11samples/DP_variant.tsv')


# In[14]:


GQ_variant_table=pd.read_table('vcf_to_try_hail/11samples/GQ_variant.tsv')
call_rate_variant_table=pd.read_table('vcf_to_try_hail/11samples/call_rate_variant.tsv')
AN_variant_table=pd.read_table('vcf_to_try_hail/11samples/AN_variant.tsv')
n_not_called_variant_table=pd.read_table('vcf_to_try_hail/11samples/n_not_called_variant.tsv')
p_value_hwe_variant_table=pd.read_table('vcf_to_try_hail/11samples/p_value_hwe_variant.tsv')
het_freq_hwe_variant_table=pd.read_table('vcf_to_try_hail/11samples/het_freq_hwe_variant.tsv')
n_het_variant_table=pd.read_table('vcf_to_try_hail/11samples/n_het_variant.tsv')


# In[15]:


DP_variant_table.rename(columns = {DP_variant_table.columns[2]:'DP'}, inplace = True)
GQ_variant_table.rename(columns = {GQ_variant_table.columns[2]:'GQ'}, inplace = True)
call_rate_variant_table.rename(columns = {call_rate_variant_table.columns[2]:'call_rate'}, inplace = True)
AN_variant_table.rename(columns = {AN_variant_table.columns[2]:"AN"}, inplace = True)
n_not_called_variant_table.rename(columns = {n_not_called_variant_table.columns[2]:"n_not_called"}, inplace = True)
p_value_hwe_variant_table.rename(columns = {p_value_hwe_variant_table.columns[2]:"p_value_hwe"}, inplace = True)
het_freq_hwe_variant_table.rename(columns = {het_freq_hwe_variant_table.columns[2]:"het_freq_hwe"}, inplace = True)
n_het_variant_table.rename(columns = {n_het_variant_table.columns[2]:"n_het"}, inplace = True)


# In[16]:


plot_histo(DP_variant_table,
           SNV_mt_variant_qc.variant_qc.dp_stats.mean,
           'Mean Depth per variant')


# In[17]:


plot_histo(GQ_variant_table,
           SNV_mt_variant_qc.variant_qc.gq_stats.mean,
           'Mean Genotype Quality per variant')


# In[18]:


plot_histo(call_rate_variant_table,
           SNV_mt_variant_qc.variant_qc.call_rate,
           'Call rate per variant')


# In[19]:


plot_histo(AN_variant_table,
           SNV_mt_variant_qc.variant_qc.AN,
           'Allele number per variant')


# In[20]:


plot_histo(n_not_called_variant_table,
           SNV_mt_variant_qc.variant_qc.n_not_called,
           'Number of samples with a missing genotype per variant')


# In[21]:


plot_histo(p_value_hwe_variant_table,
           SNV_mt_variant_qc.variant_qc.p_value_hwe,
           'p-value from two-sided test of Hardy-Weinberg equilibrium per variant')


# In[22]:


plot_sp (het_freq_hwe_variant_table,
         SNV_mt_variant_qc.variant_qc.het_freq_hwe,
         n_het_variant_table,
         SNV_mt_variant_qc.variant_qc.n_het,
         'Expected frequency of heterozygous samples under Hardy-Weinberg equilibrium',
         'Number of heterozygous samples per variant')


# In[ ]:





# In[ ]:





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

# In[28]:


filtered_SNV_mt_sample_qc_variants = SNV_mt_variant_qc.filter_rows((SNV_mt_variant_qc.variant_qc.dp_stats.mean > stat(DP_variant_table) [2]) &
                                                                 (SNV_mt_variant_qc.variant_qc.gq_stats.mean > stat(GQ_variant_table) [2]) &
                                                                 (SNV_mt_variant_qc.variant_qc.call_rate > stat(call_rate_variant_table) [2]) &
                                                                 (SNV_mt_variant_qc.variant_qc.AN > stat(AN_variant_table) [2]) &
                                                                 (SNV_mt_variant_qc.variant_qc.n_not_called > stat(n_not_called_variant_table) [2]) 
                                                                )


# In[ ]:


hl.export_vcf(filtered_SNV_mt_sample_qc_variants, 'filtered_samples_variants.vcf.bgz')


# In[32]:


SNV_mt_variant_qc.count()


# In[30]:


filtered_SNV_mt_sample_qc_variants.count()


# **Percentage of varaints removed by the filtering**

# In[37]:


perc_removed_varaints = (SNV_mt_variant_qc.count()[0]-filtered_SNV_mt_sample_qc_variants.count()[0])/SNV_mt_variant_qc.count()[0] * 100


# In[44]:


print("%.2f %% of the varaints were filtered out." % perc_removed_varaints)


# In[42]:





# In[ ]:




