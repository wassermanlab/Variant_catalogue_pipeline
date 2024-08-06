#!/usr/bin/env python
# coding: utf-8

# Hail and plot initialisation

# Last Updated: June 1, 2024 - Stephanie Petrone
#  - Update: fixing issue with contig naming discrepency for GRCh38 files 
#   and Hail's expected labeling conventions ('chr' prefix for GRCh38)
#  - Update: sex imputation cut-off (f-stat) changed from 0.2 to 0.5 as XX samples
#    were being missed. Histogram added to visualize during run for QC check. 
# In[91]:

import sys
temp_directory=sys.argv[2]
genome=sys.argv[3] #either GRCh37 or GRCh38 (params.assembly)
ref_fasta=sys.argv[4]
ref_fasta_index=sys.argv[5]

import hail as hl
from hail.plot import output_notebook, show
hl.init(tmp_dir=temp_directory)
output_notebook()


# In[ ]:


from hail.plot import show
from pprint import pprint
from bokeh.models import Span
hl.plot.output_notebook()
from bokeh.models import Range1d
from bokeh.plotting import figure, output_file, show, save
from bokeh.io import export_png

import pandas as pd
import os


# Import a vcf file and read it as a matrix table (mt, hail specific file type)
# For specific on how to look at the mt file, refer to the bottom of this Jupyter notebook)
# Use different function calls for GRCh37 and GRCh38 to accounting
# for different contig labeling and Hail's requirements
# (Hail requires GRCh37 to have no 'chr' prefix and GRCh38 to use it
if genome == 'GRCh37':
    hl.import_vcf(sys.argv[1],
        array_elements_required=False, 
        force_bgz=True, 
        reference_genome=genome).write('SNV_vcf.mt', overwrite=True)

elif genome == 'GRCh38':
    # add extra re-coding step, as reference genome and vcf files
    # for GRCh38 data are using GRCh37 labelling (no 'chr' prefix)
    # and Hail requires that they use the 'chr' prefix
    recode = {"MT":"chrM", **{f"{i}":f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}} 

    hl.import_vcf(sys.argv[1], 
        array_elements_required=False, 
        force_bgz=True, 
        reference_genome=genome,
        contig_recoding=recode).write('SNV_vcf.mt', overwrite=True)

else:
    raise Exception("Must use GRCh37 or GRCh38 assembly!")

# In[ ]:


mt = hl.read_matrix_table('SNV_vcf.mt')


# The vcf should be merged into one vcf to avoid redundancy (Possible calling overlap for indel witht he deepvaraint and SV pipeline)

# In order to create the graph, 3 functions were needed
# - stat : To calcualte the mean, standard deviation and other metrics for each parameter
# - plot_histo : To create the histogram as expected
# - plot_sp : To create the scatter plots as expected

# In[ ]:


def stat(table):
    Mean = table[table.columns[1]]. mean()  
    StdDev = table[table.columns[1]]. std()
    Low_threashold = Mean - 3*StdDev
    High_threashold = Mean + 3*StdDev
    min_graph = table[table.columns[1]]. min() - 3*StdDev
    max_graph = table[table.columns[1]]. max() + 3*StdDev
    return Mean, StdDev, Low_threashold, High_threashold, min_graph, max_graph


# In[106]:


def plot_histo (table_plot, mt_plot, variable) :
    output_file(filename=("sample_QC_"+variable+".html"), title="Sample QC HTML file")
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


# In[107]:


def plot_sp (table_x_axis, mt_x_axis, table_y_axis, mt_y_axis, x_variable, y_variable) :
    output_file(filename=("sample_QC_"+x_variable+"_"+y_variable+".html"), title="Sample QC HTML file")
    p = hl.plot.scatter(x=mt_x_axis,
                   y=mt_y_axis,
                  xlabel=x_variable,
                  ylabel=y_variable,
                title="Red lines are Mean +/- 3xStdDev",
                  hover_fields={'ID':mt.s},
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


# **Generate the sample quality control metrics using hail**

# In[ ]:


mt = hl.sample_qc(mt)


# List the sample quality control metric that were generated by hail


# Create plots for sample QC
# 
# Following the lab meeting, it was decided to define the threasholds based on the standard deviations and no hard filters. This will allow to identify outliers without relying on gnomAD values.
# 
# The Standard deviation is available thourgh the summarize function, but I cannot figure out how to extract the value. I can't figure ouot how to create a pd DataFrame from Hail matric neither (to_Pandas() is not working). Solution : Save the file with values as a csv and reopen it, calculate the Std Dev based on this.

# List of variables for which we will create a table, calculate the standard deviation (StdDev) and the mean (Mean) for sample QC:
# - DP (mt_sample_qc.sample_qc.dp_stats.mean)
# - QG (mt_sample_qc.sample_qc.gq_stats.mean)
# - call_rate (mt_sample_qc.sample_qc.call_rate)
# - r_het_hom_var (mt_sample_qc.sample_qc.r_het_hom_var)
# - n_het (mt_sample_qc.sample_qc.n_het)
# - n_hom_var (mt_sample_qc.sample_qc.n_hom_var)
# - n_snp (mt_sample_qc.sample_qc.n_snp)
# - n_singleton (mt_sample_qc.sample_qc.n_singleton)
# - r_insertion_deletion (mt_sample_qc.sample_qc.r_insertion_deletion)
# - n_insertion (mt_sample_qc.sample_qc.n_insertion)
# - n_deletion (mt_sample_qc.sample_qc.n_deletion)
# - r_ti_tv (mt_sample_qc.sample_qc.r_ti_tv)
# - n_transition (mt_sample_qc.sample_qc.n_transition)
# - n_transversion (mt_sample_qc.sample_qc.n_transversion)

# Save the values as table

# In[ ]:


mt.sample_qc.dp_stats.mean.export('DP.tsv')
mt.sample_qc.gq_stats.mean.export('GQ.tsv')
mt.sample_qc.call_rate.export('call_rate.tsv')
mt.sample_qc.r_het_hom_var.export('r_het_hom_var.tsv')
mt.sample_qc.n_het.export('n_het.tsv')
mt.sample_qc.n_hom_var.export('n_hom_var.tsv')
mt.sample_qc.n_snp.export('n_snp.tsv')
mt.sample_qc.n_singleton.export('n_singleton.tsv')
mt.sample_qc.r_insertion_deletion.export('r_insertion_deletion.tsv')
mt.sample_qc.n_insertion.export('n_insertion.tsv')
mt.sample_qc.n_deletion.export('n_deletion.tsv')
mt.sample_qc.r_ti_tv.export('r_ti_tv.tsv')
mt.sample_qc.n_transition.export('n_transition.tsv')
mt.sample_qc.n_transversion.export('n_transversion.tsv')

# Open the tables as data frame

# In[ ]:


DP_table=pd.read_table('DP.tsv')
GQ_table=pd.read_table('GQ.tsv')
call_rate_table=pd.read_table('call_rate.tsv')
r_het_hom_var_table=pd.read_table('r_het_hom_var.tsv')
n_het_table=pd.read_table('n_het.tsv')
n_hom_var_table=pd.read_table('n_hom_var.tsv')
n_snp_table=pd.read_table('n_snp.tsv')
n_singleton_table=pd.read_table('n_singleton.tsv')
r_insertion_deletion_table=pd.read_table('r_insertion_deletion.tsv')
n_insertion_table=pd.read_table('n_insertion.tsv')
n_deletion_table=pd.read_table('n_deletion.tsv')
r_ti_tv_table=pd.read_table('r_ti_tv.tsv')
n_transition_table=pd.read_table('n_transition.tsv')
n_transversion_table=pd.read_table('n_transversion.tsv')


# Rename the column of the tables

# In[ ]:


DP_table.rename(columns = {DP_table.columns[1]:'DP'}, inplace = True)
GQ_table.rename(columns = {GQ_table.columns[1]:'GQ'}, inplace = True)
call_rate_table.rename(columns = {call_rate_table.columns[1]:'call_rate'}, inplace = True)
r_het_hom_var_table.rename(columns = {r_het_hom_var_table.columns[1]:'r_het_hom_var'}, inplace = True)
n_het_table.rename(columns = {n_het_table.columns[1]:'n_het'}, inplace = True)
n_hom_var_table.rename(columns = {n_hom_var_table.columns[1]:'n_hom_var'}, inplace = True)
n_snp_table.rename(columns = {n_snp_table.columns[1]:'n_snp'}, inplace = True)
n_singleton_table.rename(columns = {n_singleton_table.columns[1]:'n_singleton'}, inplace = True)
r_insertion_deletion_table.rename(columns = {r_insertion_deletion_table.columns[1]:'r_insertion_deletion'}, inplace = True)
n_insertion_table.rename(columns = {n_insertion_table.columns[1]:'n_insertion'}, inplace = True)
n_deletion_table.rename(columns = {n_deletion_table.columns[1]:'n_deletion'}, inplace = True)
r_ti_tv_table.rename(columns = {r_ti_tv_table.columns[1]:'r_ti_tv'}, inplace = True)
n_transition_table.rename(columns = {n_transition_table.columns[1]:'n_transition'}, inplace = True)
n_transversion_table.rename(columns = {n_transversion_table.columns[1]:'n_transversion'}, inplace = True)


# Create the graphs

# In[108]:


plot_histo(DP_table, mt.sample_qc.dp_stats.mean, 'Mean Depth per sample')


# In[109]:


plot_histo(GQ_table,
           mt.sample_qc.gq_stats.mean,
           'Mean Genotype quality per sample')


# In[110]:


plot_histo(call_rate_table,
           mt.sample_qc.call_rate,
           'Call Rate per sample')


# In[111]:


plot_histo(r_het_hom_var_table,
           mt.sample_qc.r_het_hom_var,
           'Ratio heterozygous to homozygous variants per sample')


# In[112]:


plot_sp (n_het_table,
         mt.sample_qc.n_het,
         n_hom_var_table,
         mt.sample_qc.n_hom_var,
         'Number of Heterozygous Variants',
         'Number of homozygous variants')


# In[113]:


plot_histo(n_snp_table,
           mt.sample_qc.n_snp,
           'Number of SNPs per sample')


# In[115]:


plot_histo(n_singleton_table,
           mt.sample_qc.n_singleton,
           'Number of singletons per sample')


# In[116]:


plot_sp (n_insertion_table,
         mt.sample_qc.n_insertion,
         n_deletion_table,
         mt.sample_qc.n_deletion,
         'Number of insertions',
         'Number of deletions')


# In[117]:


plot_histo(r_insertion_deletion_table,
           mt.sample_qc.r_insertion_deletion,
           'Ratio insertions to deletions per sample')


# In[118]:


plot_sp (n_transition_table,
         mt.sample_qc.n_transition,
         n_transversion_table,
         mt.sample_qc.n_transversion,
         'Number of transitions',
         'Number of transversions')


# In[119]:


plot_histo(r_ti_tv_table,
           mt.sample_qc.r_ti_tv,
           'Ratio transitions to transversions per sample')


## In[ ]:
#
#
#
#
#
## In[ ]:





# **Filter the samples based on the threasholds repersented on the figures**
# 
# On the test with 11 samples, no samples shoould be removed
# 
# Low_threashold = Mean - 3*StdDev = stat(table) [2]
# 
# High_threashold = Mean + 3*StdDev = stat(table) [3]
# 
# Filters :
# - Mean DP per sample lower than the low threshoold
# - Mean Genotype Quality per sample lower than the low threshoold
# - Mean call rate per sample lower than the low threshoold
# - Ratio heterozygous to homozygous variants per sample lower than the low threshoold or higher than high threshold
# - Number of SNPs per sample lower than the low threshoold or higher than high threshold
# - Number of singletons per sample lower than the low threshoold or higher than high threshold
# - ?? Ratio insertions to deletion per sample lower than the low threshoold or higher than high threshold
# - Ratio transition to transversions per sample lower than the low threshoold or higher than high threshold

# In[120]:


filtered_mt = mt.filter_cols(
    (stat(DP_table) [3] > mt.sample_qc.dp_stats.mean)  &
    (mt.sample_qc.dp_stats.mean > stat(DP_table) [2]) &
    (stat(GQ_table) [3] > mt.sample_qc.gq_stats.mean) &
    (mt.sample_qc.gq_stats.mean > stat(GQ_table) [2]) &
    (stat(call_rate_table) [3] > mt.sample_qc.call_rate) &
    (mt.sample_qc.call_rate > stat(call_rate_table) [2]) &
    (stat(r_het_hom_var_table) [3] > mt.sample_qc.r_het_hom_var) &
    (mt.sample_qc.r_het_hom_var > stat(r_het_hom_var_table) [2]) &
    (stat(n_snp_table) [3] > mt.sample_qc.n_snp) &
    (mt.sample_qc.n_snp > stat(n_snp_table) [2]) &
    (stat(n_singleton_table) [3] > mt.sample_qc.n_singleton) &
    (mt.sample_qc.n_singleton > stat(n_singleton_table) [2]) &
    (stat(r_insertion_deletion_table) [3] > mt.sample_qc.r_insertion_deletion) &
    (mt.sample_qc.r_insertion_deletion > stat(r_insertion_deletion_table) [2]) &
    (stat(r_ti_tv_table) [3] > mt.sample_qc.r_ti_tv) &
    (mt.sample_qc.r_ti_tv> stat(r_ti_tv_table) [2])
)


# In[124]:

hl.export_vcf(filtered_mt, 'filtered_samples.vcf.bgz', tabix = True)

# Write the report of the number of filtered out samples and the reason they were filtered out

def calc_removed_samples(mt, mt_var, stat_table) :
    # Save sample genotype quality metrics information to separate file
    input_mt = mt.annotate_cols(
        keep=(mt_var > stat_table [2]) &
            (mt_var < stat_table [3])
    )

    n_removed = input_mt.aggregate_cols(hl.agg.count_where(~input_mt.keep))
    
    return n_removed


def report_stats():
    """
    Generate output report with basic stats.
    """
    out_stats = hl.hadoop_open(f"sample_QC.txt", "w")
    # Report numbers of filtered samples
    out_stats.write(
        f"Number of samples removed because of depth metrics: {DP_removed}\n"
        f"Number of samples removed because of genotype quality metrics: {GQ_removed}\n"
        f"Number of samples removed because of call rate metrics: {CR_removed}\n"
        f"Number of samples removed because of ratio heterozygous over homozygous: {r_het_hom_removed}\n"
        f"Number of samples removed because of number of snps: {n_snps_removed}\n"
        f"Number of samples removed because of number of singletons: {n_singletons_removed}\n"
        f"Number of samples removed because of ratio insertions over deletions: {r_ins_del_removed}\n"
        f"Number of samples removed because of ratio transversions / transitions: {r_ti_tv_removed}\n"
        f"Percentage of the samples filtered out: {perc_removed_samples}\n"
    )
    out_stats.close()


DP_removed = calc_removed_samples(mt, mt.sample_qc.dp_stats.mean, stat(DP_table))
GQ_removed = calc_removed_samples(mt, mt.sample_qc.gq_stats.mean, stat(GQ_table))
CR_removed = calc_removed_samples(mt, mt.sample_qc.call_rate, stat(call_rate_table))
r_het_hom_removed = calc_removed_samples(mt, mt.sample_qc.r_het_hom_var, stat(r_het_hom_var_table))
n_snps_removed = calc_removed_samples(mt, mt.sample_qc.n_snp, stat(n_snp_table))
n_singletons_removed = calc_removed_samples(mt, mt.sample_qc.n_singleton, stat(n_singleton_table))
r_ins_del_removed = calc_removed_samples(mt, mt.sample_qc.r_insertion_deletion, stat(r_insertion_deletion_table))
r_ti_tv_removed = calc_removed_samples(mt, mt.sample_qc.r_ti_tv, stat(r_ti_tv_table))
perc_removed_samples = (mt.count()[0]-filtered_mt.count()[0])/mt.count()[0] * 100

report_stats()


# Impute sex based on F-stat (only for sample who passed QC)

# Using gnomAD hard filters :
#- Ambiguous sex: fell outside of:
#- XY: F-stat > 0.8
#- XX: F-stat < 0.5 - this value is determined by visually
#                     looking at the distribution

imputed_sex_filtered_samples = hl.impute_sex(filtered_mt.GT)
imputed_sex_filtered_samples = imputed_sex_filtered_samples.annotate(
        sex=hl.if_else(imputed_sex_filtered_samples.f_stat < 0.5,
                       "XX",
                       (hl.if_else(imputed_sex_filtered_samples.f_stat > 0.8,"XY", "ambiguous"))
                )
        )
filtered_samples_sex=imputed_sex_filtered_samples.select("sex")
filtered_samples_sex.export('filtered_samples_sex.tsv')



# Create a histogram for showing distribution for XX and XY 
# calculated with the hl.impute_sex() function (f-stat)
pl = hl.plot.histogram(imputed_sex_filtered_samples.f_stat, title="Sample Sex Distribution")
pl.yaxis.axis_label = 'Count'
pl.xaxis.axis_label = 'F-stat'
annot = Span(dimension="height", location=0.5 ,line_dash='dashed', line_width=3,line_color="red")
pl.add_layout(annot)

output_file(filename=("impute_sex_distribution.html"))
save(pl)
