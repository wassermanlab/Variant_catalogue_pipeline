#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Hail and plot initialisation


# In[1]:


temp_directory="./tmp"
import sys
import hail as hl
from hail.plot import output_notebook, show
from pprint import pprint
from bokeh.models import Span
hl.plot.output_notebook()
from bokeh.models import Range1d
from bokeh.plotting import figure, output_file, show, save
from bokeh.io import export_png
import pandas as pd
import os

# Hail and plot initialisation
# Configure Spark properties
spark_conf = {
    'spark.driver.memory': '8g'  # Set the driver memory, e.g., to 8 GB
}

# Initialize Hail with custom Spark configuration
hl.init(master='local[*]', spark_conf=spark_conf, tmp_dir=temp_directory)

output_notebook()


# In[41]:



temp_directory=sys.argv[2]
genome=sys.argv[3] #either GRCh37 or GRCh38 (params.assembly)
ref_fasta=sys.argv[4]
ref_fasta_index=sys.argv[5]
pop_vcf = sys.argv[1]

#use filter list only when it is passed as argument
if ((len(sys.argv)-1) == 6):
    filter_samples = True # flag for running filter step
    filter_table = sys.argv[6]



# Import a vcf file and read it as a matrix table (mt, hail specific file type)
# For specific on how to look at the mt file, refer to the bottom of this Jupyter notebook)
# Use different function calls for GRCh37 and GRCh38 to accounting
# for different contig labeling and Hail's requirements
# (Hail requires GRCh37 to have no 'chr' prefix and GRCh38 to use it
if genome == 'GRCh37':
    hl.import_vcf(pop_vcf,
        array_elements_required=False, 
        force_bgz=True, 
        reference_genome=genome).write('SNV_vcf.mt', overwrite=True)

elif genome == 'GRCh38':
    # add extra re-coding step, as reference genome and vcf files
    # for GRCh38 data are using ensemble labelling (no 'chr' prefix)
    # and Hail requires that they use the 'chr' prefix
    recode = {"MT":"chrM", **{f"{i}":f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}} 

    hl.import_vcf(pop_vcf, 
        array_elements_required=False, 
        force_bgz=True, 
        reference_genome=genome,
        contig_recoding=recode).write('SNV_vcf.mt', overwrite=True)

else:
    raise Exception("Must use GRCh37 or GRCh38 assembly!")


# In[141]:


mt = hl.read_matrix_table('SNV_vcf.mt')


# The vcf should be merged into one vcf to avoid redundancy (Possible calling overlap for indel witht he deepvaraint and SV pipeline)

# In order to create the graph, 3 functions were needed
# - stat : To calcualte the mean, standard deviation and other metrics for each parameter
# - plot_histo : To create the histogram as expected
# - plot_sp : To create the scatter plots as expected


# In[142]:


def stat(table):
    Mean = table[table.columns[1]]. mean()  
    StdDev = table[table.columns[1]]. std()
    Low_threashold = Mean - 3*StdDev
    High_threashold = Mean + 3*StdDev
    min_graph = table[table.columns[1]]. min() - 3*StdDev
    max_graph = table[table.columns[1]]. max() + 3*StdDev
    return Mean, StdDev, Low_threashold, High_threashold, min_graph, max_graph


# In[143]:


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


# In[6]:


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


# In[144]:


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


# In[145]:


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


# In[146]:


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


# In[147]:


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


# In[148]:


plot_histo(DP_table, mt.sample_qc.dp_stats.mean, 'Mean Depth per sample')


# In[149]:


plot_histo(GQ_table,
           mt.sample_qc.gq_stats.mean,
           'Mean Genotype quality per sample')


# In[150]:


plot_histo(call_rate_table,
           mt.sample_qc.call_rate,
           'Call Rate per sample')


# In[151]:


plot_histo(r_het_hom_var_table,
           mt.sample_qc.r_het_hom_var,
           'Ratio heterozygous to homozygous variants per sample')


# In[152]:


plot_sp (n_het_table,
         mt.sample_qc.n_het,
         n_hom_var_table,
         mt.sample_qc.n_hom_var,
         'Number of Heterozygous Variants',
         'Number of homozygous variants')


# In[153]:


plot_histo(n_snp_table,
           mt.sample_qc.n_snp,
           'Number of SNPs per sample')


# In[154]:


plot_histo(n_singleton_table,
           mt.sample_qc.n_singleton,
           'Number of singletons per sample')


# In[155]:


plot_sp (n_insertion_table,
         mt.sample_qc.n_insertion,
         n_deletion_table,
         mt.sample_qc.n_deletion,
         'Number of insertions',
         'Number of deletions')


# In[156]:


plot_histo(r_insertion_deletion_table,
           mt.sample_qc.r_insertion_deletion,
           'Ratio insertions to deletions per sample')


# In[157]:


plot_sp (n_transition_table,
         mt.sample_qc.n_transition,
         n_transversion_table,
         mt.sample_qc.n_transversion,
         'Number of transitions',
         'Number of transversions')


# In[158]:


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



# In[160]:

# **Flag the samples based on the threasholds repersented on the figures**
# (this was previously used as an automatic filter, but it was quite
# arbitrarily filtering good samples). Instead, manually look at the flagged
# samples.
# (originally Solenne's code)
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


# Return a list of flagged samples and number

def calc_flagged_samples(mt, mt_var, stat_table) :
    # Save sample genotype quality metrics information to separate file
    input_mt = mt.annotate_cols(
        keep=(mt_var > stat_table [2]) &
            (mt_var < stat_table [3])
    )

    flagged = input_mt.filter_cols(input_mt.keep == False, keep=True)
    flagged_samples = flagged.s.collect()
    
    return flagged_samples


# In[161]:


def report_stats():
    """
    Generate output report with basic stats.
    """
    out_stats = hl.hadoop_open(f"sample_QC.txt", "w")
    # Report numbers of filtered samples
    out_stats.write(
        f"Samples flagged because of depth metrics: {DP_flagged}\n"
        f"Count: {DP_count}\n\n"

        f"Samples flagged because of genotype quality metrics: {GQ_flagged}\n"
        f"Count: {GQ_count}\n\n"

        f"Samples flagged because of call rate metrics: {CR_flagged}\n"
        f"Count: {CR_count}\n\n"

        f"Samples flagged because of ratio heterozygous over homozygous: {r_het_hom_flagged}\n"
        f"Count: {r_het_hom_count}\n\n"

        f"Samples flagged because of number of snps: {n_snps_flagged}\n"
        f"Count: {n_snps_count}\n\n"

        f"Samples flagged because of number of singletons: {n_singletons_flagged}\n"
        f"Count: {n_singleton_count}\n\n"

        f"Samples flagged because of ratio insertions over deletions: {r_ins_del_flagged}\n"
        f"Count: {r_ins_del_count}\n\n"

        f"Samples flagged because of ratio transversions / transitions: {r_ti_tv_flagged}\n"
        f"Count: {r_ti_tv_count}\n\n"

        f"All samples flagged: {all_samples_flagged}\n"
        f"Count of the samples flagged: {count_flagged_samples}\n"
    )
    out_stats.close()




# In[162]:


DP_flagged = calc_flagged_samples(mt, mt.sample_qc.dp_stats.mean, stat(DP_table))
DP_count = len(DP_flagged)

GQ_flagged = calc_flagged_samples(mt, mt.sample_qc.gq_stats.mean, stat(GQ_table))
GQ_count = len(GQ_flagged)

CR_flagged = calc_flagged_samples(mt, mt.sample_qc.call_rate, stat(call_rate_table))
CR_count = len(CR_flagged)

r_het_hom_flagged = calc_flagged_samples(mt, mt.sample_qc.r_het_hom_var, stat(r_het_hom_var_table))
r_het_hom_count = len(r_het_hom_flagged)

n_snps_flagged = calc_flagged_samples(mt, mt.sample_qc.n_snp, stat(n_snp_table))
n_snps_count = len(n_snps_flagged)

n_singletons_flagged = calc_flagged_samples(mt, mt.sample_qc.n_singleton, stat(n_singleton_table))
n_singleton_count = len(n_singletons_flagged)

r_ins_del_flagged = calc_flagged_samples(mt, mt.sample_qc.r_insertion_deletion, stat(r_insertion_deletion_table))
r_ins_del_count = len(r_ins_del_flagged)

r_ti_tv_flagged = calc_flagged_samples(mt, mt.sample_qc.r_ti_tv, stat(r_ti_tv_table))
r_ti_tv_count = len(r_ti_tv_flagged)

all_samples_flagged =  list(set(DP_flagged + GQ_flagged + CR_flagged + r_het_hom_flagged + n_snps_flagged + n_singletons_flagged +
                              r_ins_del_flagged + r_ti_tv_flagged))

count_flagged_samples = len(all_samples_flagged)

report_stats()




# In[163]:


# Impute sex based on F-stat - for all samples

# Using gnomAD hard filters :
#- Ambiguous sex: fell outside of:
#- XY: F-stat > 0.8
#- XX: F-stat < 0.5 - this value is determined by visually
#                     looking at the distribution

imputed_sex_filtered_samples = hl.impute_sex(mt.GT)
imputed_sex_filtered_samples = imputed_sex_filtered_samples.annotate(
        sex=hl.if_else(imputed_sex_filtered_samples.f_stat < 0.5,
                       "XX",
                       (hl.if_else(imputed_sex_filtered_samples.f_stat > 0.8,"XY", "ambiguous"))
                )
        )
filtered_samples_sex=imputed_sex_filtered_samples.select("sex")
filtered_samples_sex.export('filtered_samples_sex.tsv')

filtered_samples_sex_fstat=imputed_sex_filtered_samples.select("sex", "f_stat")
filtered_samples_sex_fstat.export('filtered_samples_sex_f_stat.tsv')



# Create a histogram for showing distribution for XX and XY 
# calculated with the hl.impute_sex() function (f-stat)
pl = hl.plot.histogram(imputed_sex_filtered_samples.f_stat, title="Sample Sex Distribution")
pl.yaxis.axis_label = 'Count'
pl.xaxis.axis_label = 'F-stat'
annot = Span(dimension="height", location=0.5 ,line_dash='dashed', line_width=3,line_color="red")
pl.add_layout(annot)

output_file(filename=("impute_sex_distribution.html"))
save(pl)




# In[164]:


#-----------------------FILTER STEP-----------------------------------
# Filter samples if flag set (optional)
# -The filter table should be a single column table containing
# the sample IDs
#---------------------------------------------------------------------

if (filter_samples):
    # important! - hail's filter_cols does not update the info fields in the
    # row fields in the matrix table. These must be updated manually.
    
    #get a list from the table
    try:
        df = pd.read_csv(filter_table, sep='\t')

    except Exception as e:
        print(f"An error occurred opening the sample filter table: {e}")
    
    
    #Filter samples - remove any that are in the list provided
    filter_list = df['s'].tolist()
    filter_count = len(filter_list)
    mt = mt.filter_cols(hl.literal(filter_list).contains(mt.s), keep=False)
    

    out_filtered = hl.hadoop_open(f"samples_filtered.txt", "w")
    
    # Report numbers of filtered samples
    out_filtered.write(
            f"Samples filtered: {filter_list}\n"
            f"Count of samples filtered: {filter_count}\n"
        )
    out_filtered.close()


# In[168]:


#hail's filter_cols function does not adjust info fields, to adjust info fields, use the code below
# note that this is re-calculate during the next step (hail_variant_qc) during stratification
# mt = hl.variant_qc(mt)
# mt = mt.annotate_rows(info = mt.info.annotate(AC=mt.variant_qc.AC)) 
# mt = mt.annotate_rows(info = mt.info.annotate(AN=mt.variant_qc.AN)) 
# mt = mt.annotate_rows(info = mt.info.annotate(AF=mt.variant_qc.AF)) 


# In[169]:


#export file
hl.export_vcf(mt, 'filtered_samples.vcf.bgz', tabix = True)


# In[ ]:

