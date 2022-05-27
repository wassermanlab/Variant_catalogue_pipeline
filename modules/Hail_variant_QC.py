#!/usr/bin/env python
# coding: utf-8

# Hail and plot initialisation

# In[1]:


import hail as hl
from hail.plot import output_notebook, show
hl.init()
output_notebook()


# In[2]:


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


# Import a vcf file and read it as a matrix table (mt, hail specific file type)
# For specific on how to look at the mt file, refer to the bottom of this Jupyter notebook)
# 
# Currently : Import only the SNV vcf file, following the sample filtering step from previously

# **Mohammed part**
# 
# Load several vcf and identify and remove the variants potentially overlapping
# 
# To be included

# In[ ]:


hl.import_vcf(sys.argv[1],
              array_elements_required=False, force_bgz=True).write('filtered_samples_vcf.mt', overwrite=True)


# vcf_path = '/mnt/scratch/SILENT/Act3/Processed/Individual/GRCh37/Batch_DryRun/Run_20220426/SNV/'
# 
# hl.import_vcf(os.path.join(vcf_path,'DeepVariant_GLnexus_Run_20220426.vcf.gz'),
#               array_elements_required=False, force_bgz=True).write('filtered_samples_vcf', overwrite=True)

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


# In[9]:


SNV_mt.s


# In[10]:


SNV_mt.describe()


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


SNV_mt.variant_qc.dp_stats.mean.export('vcf_to_try_hail/11samples/DP_variant.tsv')


# In[12]:


SNV_mt.variant_qc.gq_stats.mean.export('vcf_to_try_hail/11samples/GQ_variant.tsv')
SNV_mt.variant_qc.call_rate.export('vcf_to_try_hail/11samples/call_rate_variant.tsv')
SNV_mt.variant_qc.AN.export('vcf_to_try_hail/11samples/AN_variant.tsv')
SNV_mt.variant_qc.n_not_called.export('vcf_to_try_hail/11samples/n_not_called_variant.tsv')
SNV_mt.variant_qc.p_value_hwe.export('vcf_to_try_hail/11samples/p_value_hwe_variant.tsv')
SNV_mt.variant_qc.het_freq_hwe.export('vcf_to_try_hail/11samples/het_freq_hwe_variant.tsv')
SNV_mt.variant_qc.n_het.export('vcf_to_try_hail/11samples/n_het_variant.tsv')


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

# In[25]:


SNV_mt_var_filtered = SNV_mt.filter_rows((SNV_mt.variant_qc.dp_stats.mean > stat(DP_variant_table) [2]) &
                                         (SNV_mt.variant_qc.gq_stats.mean > stat(GQ_variant_table) [2]) &
                                         (SNV_mt.variant_qc.call_rate > stat(call_rate_variant_table) [2]) &
                                         (SNV_mt.variant_qc.AN > stat(AN_variant_table) [2]) &
                                         (SNV_mt.variant_qc.n_not_called > stat(n_not_called_variant_table) [2])
                                        )


# In[55]:


hl.export_vcf(SNV_mt_var_filtered, 'SNV_filtered_samples_variants.vcf.bgz', tabix=True)


# In[29]:


SNV_mt.count()


# In[28]:


SNV_mt_var_filtered.count()


# **Percentage of variants removed by the filtering**

# In[30]:


perc_removed_variants = (SNV_mt.count()[0]-SNV_mt_var_filtered.count()[0])/SNV_mt.count()[0] * 100


# In[32]:


print("%.2f %% of the variants were filtered out." % perc_removed_variants)


# **Calculate statistic**
# 
# To save time in R, calculate AF, AC, AN and numb of homozygotes (Total) 
# 
# Column order and name expected for Oracle Apex output: variant, af_total, af_xx, af_xy, ac_total, ac_xx, ac_xy, an_total, an_xx, an_xy, hom_alt_total, hom_alt_xx, hom_alt_xy, quality
# 
# Initially calculation is done for total, then XX, then XY specific frequencies are added to the vcf info tab

# In[106]:


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
# It was considered to use hail sex inputation to define the sex but it relies only on the F-Coeff
# 
# Using a separate script (relyine on mosdepth output and plink) allow to rely on both the F_coeff and normalized coverage on the sexual chromosomes for better imputation.
# 
# Steps :
# - Import the file with sex
# -  Merge the sex table with the matrix table
# - Calculate the info (AC, AF, AN, numb of hom) per sex

# In[ ]:


sex_table = (hl.import_table(sys.argv[2], impute=True)
         .key_by('sample'))


# sex_table = (hl.import_table('/mnt/scratch/SILENT/Act3/Processed/Individual/GRCh37/Batch_DryRun/Run_20220426/QC/Aggregated/R/QC_sample.tsv', impute=True)
#          .key_by('sample'))

# In[121]:


SNV_mt_var_filtered_sex = SNV_mt_var_filtered.annotate_cols(**sex_table[SNV_mt_var_filtered.s])


# In[122]:


SNV_mt_var_filtered_XX = SNV_mt_var_filtered_sex.filter_cols(SNV_mt_var_filtered_sex.Sex == 'XX')


# In[123]:


SNV_mt_var_filtered_XX = hl.variant_qc(SNV_mt_var_filtered_XX)


# In[124]:


SNV_mt_var_filtered_XX = SNV_mt_var_filtered_XX.annotate_rows(
    info = SNV_mt_var_filtered_XX.info.annotate(
        af_xx=SNV_mt_var_filtered_XX.variant_qc.AF[1],
        ac_xx=SNV_mt_var_filtered_XX.variant_qc.AC[1],
        an_xx=SNV_mt_var_filtered_XX.variant_qc.AN,
        hom_alt_xx=SNV_mt_var_filtered_XX.variant_qc.homozygote_count[1],
    )
) 


# In[125]:


SNV_mt_var_filtered_XY = SNV_mt_var_filtered_sex.filter_cols(SNV_mt_var_filtered_sex.Sex == 'XY')


# In[126]:


SNV_mt_var_filtered_XY = hl.variant_qc(SNV_mt_var_filtered_XY)


# In[127]:


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
# As the frequencies are calculated using hail, it is not necessary to export the individual genotype for the last part of the pipeline (Annotation and data organization)
# 
# Aggregate the info columns with frequenies for total, XX and XY

# In[186]:


SNV_mt_var_filtered_tot_info = SNV_mt_var_filtered_tot.select_rows(
    SNV_mt_var_filtered_tot.info.chrom,
    SNV_mt_var_filtered_tot.info.pos,
    SNV_mt_var_filtered_tot.info.ref,
    SNV_mt_var_filtered_tot.info.alt,
    SNV_mt_var_filtered_tot.rsid,
    SNV_mt_var_filtered_tot.qual,
    SNV_mt_var_filtered_tot.info.af_tot,
    SNV_mt_var_filtered_tot.info.ac_tot,
    SNV_mt_var_filtered_tot.info.an_tot,
    SNV_mt_var_filtered_tot.info.hom_alt_tot,
).rows()


# In[187]:


SNV_mt_var_filtered_XX_info = SNV_mt_var_filtered_XX.select_rows(
    SNV_mt_var_filtered_XX.info.af_xx,
    SNV_mt_var_filtered_XX.info.ac_xx,
    SNV_mt_var_filtered_XX.info.an_xx,
    SNV_mt_var_filtered_XX.info.hom_alt_xx,
).rows()


# In[188]:


SNV_mt_var_filtered_XY_info = SNV_mt_var_filtered_XY.select_rows(
    SNV_mt_var_filtered_XY.info.af_xy,
    SNV_mt_var_filtered_XY.info.ac_xy,
    SNV_mt_var_filtered_XY.info.an_xy,
    SNV_mt_var_filtered_XY.info.hom_alt_xy,
).rows()


# In[189]:


SNV_mt_var_filtered_tot_XX_info = SNV_mt_var_filtered_tot_info.join(SNV_mt_var_filtered_XX_info, how='left')


# In[190]:


SNV_mt_var_filtered_tot_XX_XY_info = SNV_mt_var_filtered_tot_XX_info.join(SNV_mt_var_filtered_XY_info, how='left')


# In[192]:


SNV_mt_var_filtered_tot_XX_XY_info.export('SNV_mt_var_filtered_tot_XX_XY_info.tsv')


# In[ ]:





# In[ ]:





# In[ ]:



