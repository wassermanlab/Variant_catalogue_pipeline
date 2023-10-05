#!/usr/bin/env python
# coding: utf-8

# In[1]:

import sys
temp_directory=sys.argv[3]
genome=sys.argv[4]
ref_fasta=sys.argv[5]
ref_fasta_index=sys.argv[6]

import hail as hl
from hail.plot import output_notebook, show
hl.init()
output_notebook()

from hail.plot import show
from bokeh.models import Span
hl.plot.output_notebook()
from bokeh.models import Range1d
from bokeh.plotting import output_file, show, save

import pandas as pd
import os

from typing import Optional, Dict, List


# #Created through the nextflow pipeline
# hl.import_vcf(sys.argv[1],array_elements_required=False, force_bgz=True).write('SV_vcf.mt', overwrite=True)
# Phil add 2023-09-07, define reference genome off the input fasta file, which we can pass here
# In[ ]:
try:
    hl.import_vcf(sys.argv[1], array_elements_required=False, force_bgz=True, reference_genome=genome).write('SV_vcf.mt', overwrite=True)
    referenceGenome = genome
except:
    # Phil add 2023-09-07, define reference genome off the input fasta file, which we can pass here, on the off-chance that the GRCh38 has contigs 1,2,3..X,Y,MT
    # PAR taken for GRCh38 from http://useast.ensembl.org/info/genome/genebuild/human_PARS.html
    referenceGenome = hl.genetics.ReferenceGenome.from_fasta_file("referenceGenome",ref_fasta,ref_fasta_index,x_contigs=['X'],y_contigs=['Y'],mt_contigs=['MT'],par=[('Y',10001,2781479),('X',10001,2781479),('Y',56887903,57217415),('X',155701383,156030895)])
    hl.import_vcf(sys.argv[1], array_elements_required=False, force_bgz=True, reference_genome=referenceGenome).write('SV_vcf.mt', overwrite=True)






sex_table = (hl.import_table(sys.argv[2], impute=True).key_by('s'))

# **Import file**
# 
# Inport a vcf file and read it as a matrix table (mt, hail specific file type)

# In[3]:


mt = hl.read_matrix_table('SV_vcf.mt')

# Add the sample sex info
# Only samples which passed QC are present in this file
mt = mt.annotate_cols(**sex_table[mt.s])

# Keep only samples that passed QC, and with non ambiguous sex
mt= mt.filter_cols((mt.sex != 'XX') | (mt.sex != 'XY'))


# **Graph functions**
# 
# In order to create the graph, 3 functions were needed
# - stat : To calcualte the mean, standard deviation and other metrics for each parameter
# - plot_histo : To create the histogram as expected
# - plot_sp : To create the scatter plots as expected

# In[4]:


def stat(table):
    Mean = table[table.columns[2]]. mean()  
    StdDev = table[table.columns[2]]. std()
    Low_threashold = Mean - 3*StdDev
    High_threashold = Mean + 3*StdDev
    min_graph = table[table.columns[2]]. min() - 3*StdDev
    max_graph = table[table.columns[2]]. max() + 3*StdDev
    return Mean, StdDev, Low_threashold, High_threashold, min_graph, max_graph


# In[5]:


def plot_histo (table_plot, mt_plot, variable) :
    output_file(filename=os.path.join(("SV_QC_"+variable+".html")), title="SV QC HTML file")
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
    output_file(filename=os.path.join(("SV_QC_"+x_variable+"_"+y_variable+".html")), title="SV QC HTML file")
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


# **SV QC**
# 
# SV exploration of the general info
# - qual (histo)
# - filters (table)

# In[8]:


qual_hist = mt.aggregate_entries(hl.expr.aggregators.hist(mt.qual,
            mt.aggregate_rows(hl.agg.min(mt.qual)),
            mt.aggregate_rows(hl.agg.max(mt.qual)), 100))
p = hl.plot.histogram(qual_hist, legend='Qual', title='Qual (All SV type)', log=True)
show(p)


# Replace AN and {} by "none" in the filter column

# In[9]:


mt = mt.annotate_rows(filters = hl.if_else(hl.len(mt.filters) > 0,
                            mt.filters,
                            hl.set(["none"])))


# In[10]:


mt = mt.annotate_rows(filters =hl.if_else(hl.is_defined(mt.filters),
                                                mt.filters,
                                                hl.set(["none"])))


# In[11]:


filters_ht = mt.rows()


# In[12]:


filters_table = filters_ht.group_by(filters_ht.filters).aggregate(n=hl.agg.count())


# In[13]:


filters_table.export('SV_filters.txt', delimiter=' : ')


# SV : Exploration of the info column
# 
# Info columns used to generate graphs / tables (and filter variants)
# 
# - Number of variant by SV type (table)
# - Number of variant identified by Manta / smoove / both (table)
# - Number of variant identified by Manta / smoove / both by SV type (table)
# - Length of SV (histo log)
# - Length of SV by variant type (histo log)
# 
# - SU : Number of pieces of evidence supporting the variant across all samples (histo)
# - PE : Number of paired-end reads supporting the variant across all samples (histo)
# - SR : Number of split reads supporting the variant across all samples (histo)
#  
# - SUPP : Number of samples supporting the variant  
# - VARCALLS : The number of variant calls supporting this variant  
# 
# 
# All info available in the vcf : 
# 
# IMPRECISE  Number=0  Type=Flag  Description: Imprecise structural variation  
# SVTYPE  Number=1  Type=String  Description: Type of structural variant  
# SVLEN  Number=.  Type=Integer  Description: Difference in length between REF and ALT alleles  
# END  Number=1  Type=Integer  Description: End position of the variant described in this record  
# CIPOS  Number=2  Type=Integer  Description: Confidence interval around POS  
# CIEND  Number=2  Type=Integer  Description: Confidence interval around END  
# CIGAR  Number=A  Type=String  Description: CIGAR alignment for each alternate indel allele  
# MATEID  Number=.  Type=String  Description: ID of mate breakend (NA only in vcf)
# EVENT  Number=1  Type=String  Description: ID of event associated to breakend  (NA only in vcf)
# HOMLEN  Number=.  Type=Integer  Description: Length of base pair identical homology at event breakpoints  (NA mostly in vcf)
# HOMSEQ  Number=.  Type=String  Description: Sequence of base pair identical homology at event breakpoints  
# SVINSLEN  Number=.  Type=Integer  Description: Length of insertion  (NA only in vcf)
# SVINSSEQ  Number=.  Type=String  Description: Sequence of insertion  (NA only in vcf)
# LEFT_SVINSSEQ  Number=.  Type=String  Description: Known left side of insertion for an insertion of unknown length  (NA mostly in vcf)
# RIGHT_SVINSSEQ  Number=.  Type=String  Description: Known right side of insertion for an insertion of unknown length  (NA mostly in vcf)
# INV3  Number=0  Type=Flag  Description: Inversion breakends open 3' of reported location  (FALSE mostly in vcf)
# INV5  Number=0  Type=Flag  Description: Inversion breakends open 5' of reported location  (FALSE mostly in vcf)
# BND_DEPTH  Number=1  Type=Integer  Description: Read depth at local translocation breakend  (NA mostly in vcf)
# MATE_BND_DEPTH  Number=1  Type=Integer  Description: Read depth at remote translocation mate breakend  (NA mostly in vcf)
# JUNCTION_QUAL  Number=1  Type=Integer  Description: If the SV junction is part of an EVENT (ie. a multi-adjacency variant) this field provides the QUAL value for the adjacency in question only (NA mostly in vcf) 
# 
# ##smoove_version=0.2.6
# ##smoove_count_stats=NA20769:8267920  18604928  233136  439748
# STRANDS  Number=.  Type=String  Description: Strand orientation of the adjacency in BEDPE format (DEL:+-   DUP:-+   INV:++/--)  
# CIPOS95  Number=2  Type=Integer  Description: Confidence interval (95%) around POS for imprecise variants  
# CIEND95  Number=2  Type=Integer  Description: Confidence interval (95%) around END for imprecise variants  
# SECONDARY  Number=0  Type=Flag  Description: Secondary breakend in a multi-line variants  (FALSE mostly in vcf)
# SU  Number=.  Type=Integer  Description: Number of pieces of evidence supporting the variant across all samples  
# PE  Number=.  Type=Integer  Description: Number of paired-end reads supporting the variant across all samples  
# SR  Number=.  Type=Integer  Description: Number of split reads supporting the variant across all samples  
# BD  Number=.  Type=Integer  Description: Amount of BED evidence supporting the variant across all samples  (NA mostly in vcf)
# EV  Number=.  Type=String  Description: Type of LUMPY evidence contributing to the variant call  (NA mostly in vcf)
# PRPOS  Number=.  Type=String  Description: LUMPY probability curve of the POS breakend  
# PREND  Number=.  Type=String  Description: LUMPY probability curve of the END breakend  
# SUPP_VEC  Number=1  Type=String  Description: Vector of supporting samples  
# SUPP_VEC_EXT  Number=1  Type=String  Description: Vector of supporting samples   potentially extended across multiple merges  
# SUPP  Number=1  Type=String  Description: Number of samples supporting the variant  
# SUPP_EXT  Number=1  Type=String  Description: Number of samples supporting the variant   potentially extended across multiple merges  
# IDLIST  Number=.  Type=String  Description: Variant IDs of variants merged to make this call (at most 1 per sample)  
# IDLIST_EXT  Number=.  Type=String  Description: Variant IDs of variants merged   potentially extended across multiple merges  
# SVMETHOD  Number=1  Type=String  (Jasmine only in the vcf)
# STARTVARIANCE  Number=1  Type=String  Description: Variance of start position for variants merged into this one  
# ENDVARIANCE  Number=1  Type=String  Description: Variance of end position for variants merged into this one  
# AVG_START  Number=1  Type=String  Description: Average start position for variants merged into this one  
# AVG_END  Number=1  Type=String  Description: Average end position for variants merged into this one  
# AVG_LEN  Number=1  Type=String  Description: Average length for variants merged into this one  
# PRECISE  Number=0  Type=Flag  Description: Precise structural variation  (False only in the vcf)
# ALLVARS_EXT  Number=.  Type=String  Description: A comma-separated of all variants supporting this call  
# VARCALLS  Number=1  Type=String  Description: The number of variant calls supporting this variant  
# INTRASAMPLE_IDLIST  Number=.  Type=String  Description: The IDs which were merged in the most recent round of merging  
# GRMPY_ID  Number=1  Type=String  Description: Graph ID for linking to genotypes.json.gz; matches record.graphinfo.ID in there.  
# AC  Number=A  Type=Integer  Description: Allele count in genotypes  
# AN  Number=1  Type=Integer  Description: Total number of alleles in called genotypes  

# In[14]:


DEL_mt = mt.filter_rows(mt.info.SVTYPE == "DEL")
DUP_mt = mt.filter_rows(mt.info.SVTYPE == "DUP")
INV_mt = mt.filter_rows(mt.info.SVTYPE == "INV")
INS_mt = mt.filter_rows(mt.info.SVTYPE == "INS")


# In[15]:


n_del = DEL_mt.count_rows()
n_dup = DUP_mt.count_rows()
n_ins = INS_mt.count_rows()
n_inv = INV_mt.count_rows()
n_SV_tot = mt.count_rows()


# Identification method (Manta, smoove, both)
# 
# 
# #This method indicate Jasmine, which is not the case
# SVMETHOD = extract.info(SV_vcf, element = "SVMETHOD")[i]
# 
# #This is to extract the actual method (in R)
# IDLIST = extract.info(SV_vcf, element = "IDLIST")[i]
# if (grepl("Manta", IDLIST, fixed=TRUE) ){
#     algorithm = "Manta"
# } else {
#     algorithm = "Smoove"    
# }
# 
# contains not working in Hail, as it needs to be the full length of the string
# 
# Can use info that is only provided by Smoove, such as SU
# [If SU is not empty, then it's smoove, otherwise, it's Manta)
# 
# Smoove produced a Variant Call Format (VCF) file that contains deletions, duplications, inversions, and unclassified breakends (BNDs) --> No insertions
# 

# In[16]:


n_SV_tot_smoove = mt.aggregate_rows(hl.agg.count_where(hl.len(mt.info.SU)>0))
n_SV_tot_manta = n_SV_tot-n_SV_tot_smoove

n_SV_del_smoove = DEL_mt.aggregate_rows(hl.agg.count_where(hl.len(DEL_mt.info.SU)>0))
n_SV_del_manta = n_del-n_SV_del_smoove

n_SV_dup_smoove = DUP_mt.aggregate_rows(hl.agg.count_where(hl.len(DUP_mt.info.SU)>0))
n_SV_dup_manta = n_dup-n_SV_dup_smoove

n_SV_ins_smoove = INS_mt.aggregate_rows(hl.agg.count_where(hl.len(INS_mt.info.SU)>0))
n_SV_ins_manta = n_ins-n_SV_ins_smoove

n_SV_inv_smoove = INV_mt.aggregate_rows(hl.agg.count_where(hl.len(INV_mt.info.SU)>0))
n_SV_inv_manta = n_inv-n_SV_inv_smoove


# Length of the SV

# In[17]:


min_len_tot = mt.aggregate_rows(hl.agg.min(mt.info.SVLEN[0]))
max_len_tot = mt.aggregate_rows(hl.agg.max(mt.info.SVLEN[0]))

min_len_del = DEL_mt.aggregate_rows(hl.agg.min(DEL_mt.info.SVLEN[0]))
max_len_del = DEL_mt.aggregate_rows(hl.agg.max(DEL_mt.info.SVLEN[0]))

min_len_dup = DUP_mt.aggregate_rows(hl.agg.min(DUP_mt.info.SVLEN[0]))
max_len_dup = DUP_mt.aggregate_rows(hl.agg.max(DUP_mt.info.SVLEN[0]))

min_len_ins = INS_mt.aggregate_rows(hl.agg.min(INS_mt.info.SVLEN[0]))
max_len_ins = INS_mt.aggregate_rows(hl.agg.max(INS_mt.info.SVLEN[0]))

min_len_inv = INV_mt.aggregate_rows(hl.agg.min(INV_mt.info.SVLEN[0]))
max_len_inv = INV_mt.aggregate_rows(hl.agg.max(INV_mt.info.SVLEN[0]))


# In[18]:


len_hist = mt.aggregate_entries(hl.expr.aggregators.hist(mt.info.SVLEN[0],
            min_len_tot, max_len_tot, 100))
p = hl.plot.histogram(len_hist, legend='length (bp)', title='Length Histogram (All SV type)', log= True)
show(p)


# In[19]:


len_hist_del = DEL_mt.aggregate_entries(hl.expr.aggregators.hist(DEL_mt.info.SVLEN[0],
            min_len_del, max_len_del, 100))
p = hl.plot.histogram(len_hist_del, legend='length (bp)', title='Length Histogram (Deletion)', log=True)
show(p)


# In[20]:


len_hist_dup = DUP_mt.aggregate_entries(hl.expr.aggregators.hist(DUP_mt.info.SVLEN[0], 
                                                                 min_len_dup, max_len_dup, 100))
p = hl.plot.histogram(len_hist_dup, legend='length (bp)', title='Length Histogram (Duplication)', log=True)
show(p)


# In[21]:


len_hist_ins = INS_mt.aggregate_entries(hl.expr.aggregators.hist(INS_mt.info.SVLEN[0], 
                                                                 min_len_ins, max_len_ins, 100))
p = hl.plot.histogram(len_hist_ins, legend='length (bp)', title='Length Histogram (Insertion)')
show(p)


# In[22]:


len_hist_inv = INV_mt.aggregate_entries(hl.expr.aggregators.hist(INV_mt.info.SVLEN[0], 
                                                                 min_len_ins, max_len_inv, 100))
p = hl.plot.histogram(len_hist_inv, legend='length (bp)', title='Length Histogram (Inversion)', log=True)
show(p)


# In[23]:


#Only calculated for variants identified by smoove

SU_hist_tot = mt.aggregate_entries(hl.expr.aggregators.hist(mt.info.SU[0], 
    mt.aggregate_rows(hl.agg.min(mt.info.SU[0])),
    mt.aggregate_rows(hl.agg.max(mt.info.SU[0])),
    100))

p = hl.plot.histogram(SU_hist_tot, legend='SU',
                      title='Number of pieces of evidence supporting the SV across all samples - Smoove SV, all types')
show(p)


# In[24]:


#Only calculated for variants identified by smoove

PE_hist_tot = mt.aggregate_entries(hl.expr.aggregators.hist(mt.info.PE[0], 
    mt.aggregate_rows(hl.agg.min(mt.info.PE[0])),
    mt.aggregate_rows(hl.agg.max(mt.info.PE[0])),
    100))

p = hl.plot.histogram(PE_hist_tot, legend='PE',
                      title='Number of paired-end reads supporting the SV across all samples - Smoove SV, all types')
show(p)


# In[25]:


#Only calculated for variants identified by smoove

SR_hist_tot = mt.aggregate_entries(hl.expr.aggregators.hist(mt.info.SR[0], 
    mt.aggregate_rows(hl.agg.min(mt.info.SR[0])),
    mt.aggregate_rows(hl.agg.max(mt.info.SR[0])),
    100))

p = hl.plot.histogram(SR_hist_tot, legend='SR',
                      title='Number of split reads supporting the SV across all samples - Smoove SV, all types')
show(p)


# In[26]:


SUPP = hl.int(mt.info.SUPP)
SUPP_hist_tot = mt.aggregate_entries(hl.expr.aggregators.hist(SUPP, 
    mt.aggregate_rows(hl.agg.min(SUPP)),
    mt.aggregate_rows(hl.agg.max(SUPP)),
    100))

p = hl.plot.histogram(SUPP_hist_tot, legend='SUPP',
                      title='Number of samples supporting the SV')
show(p)


# In[27]:


p = hl.plot.scatter(x=SUPP,
                   y=mt.info.AC[0],
                  xlabel="SUPP : Number of samples supporting the SV",
                  ylabel="AC : Allele count in genotypes",
                title="SUPP/AC")
show(p)


# In[28]:


VARCALLS = hl.int(mt.info.VARCALLS)
VARCALLS_hist_tot = mt.aggregate_entries(hl.expr.aggregators.hist(VARCALLS, 
    mt.aggregate_rows(hl.agg.min(VARCALLS)),
    mt.aggregate_rows(hl.agg.max(VARCALLS)),
    100))

p = hl.plot.histogram(VARCALLS_hist_tot, legend='VARCALLS',
                      title='The number of variant calls supporting this SV')
show(p)


# In[29]:


p = hl.plot.scatter(x=SUPP,
                   y=VARCALLS,
                  xlabel="SUPP : Number of samples supporting the SV",
                  ylabel="VARCALLS : The number of variant calls supporting this SV",
                title="SUPP / VARCALLS")
show(p)


# In[30]:


p = hl.plot.scatter(x=mt.info.SU[0],
                   y=mt.info.PE[0],
                  xlabel="SU : Number of pieces of evidence supporting the SV across all samples",
                  ylabel="PE : Number of paired-end reads supporting the SV across all samples",
                title="SU/PE")
show(p)


# In[31]:


p = hl.plot.scatter(x=mt.info.SU[0],
                   y=mt.info.SR[0],
                  xlabel="SU : Number of pieces of evidence supporting the SV across all samples",
                  ylabel="SR : Number of split reads supporting the SV across all samples",
                title="SU/SR")
show(p)


# In[32]:


p = hl.plot.scatter(x=mt.info.PE[0],
                   y=mt.info.SR[0],
                  xlabel="PE : Number of paired-end reads supporting the SV across all samples",
                  ylabel="SR : Number of split reads supporting the SV across all samples",
                title="PE/SR")
show(p)


# **SV Hail QC**

# In[33]:


mt = hl.variant_qc(mt)


# List of variant_qc variable of interest:
# - DP (SV_mt.variant_qc.dp_stats.mean)
# - AN (SV_mt.variant_qc.AN)
# - call_rate (SV_mt.variant_qc.call_rate)
# - n_called (SV_mt.variant_qc.n_called)
# - n_not_called (SV_mt.variant_qc.n_not_called)
# - n_het (SV_mt.variant_qc.n_het)
# - p_value_hwe (SV_mt.variant_qc.p_value_hwe)
# - het_freq_hwe (SV_mt.variant_qc.het_freq_hwe)
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
# - n_filtered (SV_mt.variant_qc.n_filtered)

# In[34]:


mt.variant_qc.dp_stats.mean.export('DP_SV.tsv')


# In[35]:


mt.variant_qc.AN.export('AN_SV.tsv')
mt.variant_qc.call_rate.export('call_rate_SV.tsv')
mt.variant_qc.n_called.export('n_called_SV.tsv')
mt.variant_qc.n_not_called.export('n_not_called_SV.tsv')
mt.variant_qc.n_het.export('n_het_SV.tsv')
mt.variant_qc.p_value_hwe.export('p_value_hwe_SV.tsv')
mt.variant_qc.het_freq_hwe.export('het_freq_hwe_SV.tsv')


# In[36]:


DP_SV_table=pd.read_table('DP_SV.tsv')


# In[37]:


AN_SV_table=pd.read_table('AN_SV.tsv')
call_rate_SV_table=pd.read_table('call_rate_SV.tsv')
n_called_SV_table=pd.read_table('n_called_SV.tsv')
n_not_called_SV_table=pd.read_table('n_not_called_SV.tsv')
n_het_SV_table=pd.read_table('n_het_SV.tsv')
p_value_hwe_SV_table=pd.read_table('p_value_hwe_SV.tsv')
het_freq_hwe_SV_table=pd.read_table('het_freq_hwe_SV.tsv')


# In[38]:


DP_SV_table.rename(columns = {DP_SV_table.columns[2]:'DP'}, inplace = True)
AN_SV_table.rename(columns = {AN_SV_table.columns[2]:"AN"}, inplace = True)
call_rate_SV_table.rename(columns = {call_rate_SV_table.columns[2]:'call_rate'}, inplace = True)
n_called_SV_table.rename(columns = {n_called_SV_table.columns[2]:"n_called"}, inplace = True)
n_not_called_SV_table.rename(columns = {n_not_called_SV_table.columns[2]:"n_not_called"}, inplace = True)
n_het_SV_table.rename(columns = {n_het_SV_table.columns[2]:"n_het"}, inplace = True)
p_value_hwe_SV_table.rename(columns = {p_value_hwe_SV_table.columns[2]:"p_value_hwe"}, inplace = True)
het_freq_hwe_SV_table.rename(columns = {het_freq_hwe_SV_table.columns[2]:"het_freq_hwe"}, inplace = True)


# In[39]:


plot_histo(DP_SV_table,
           mt.variant_qc.dp_stats.mean,
           'Mean Depth per SV')


# In[40]:


plot_histo(AN_SV_table,
           mt.variant_qc.AN,
           'Allele number per SV')


# In[41]:


plot_histo(call_rate_SV_table,
           mt.variant_qc.call_rate,
           'Call rate per SV')


# In[42]:


plot_histo(n_called_SV_table,
           mt.variant_qc.n_called,
           'Number of samples with a called genotype per SV')


# In[43]:


plot_histo(n_not_called_SV_table,
           mt.variant_qc.n_not_called,
           'Number of samples with a missing genotype per SV')


# In[44]:


plot_histo(p_value_hwe_SV_table,
           mt.variant_qc.p_value_hwe,
           'p-value from two-sided test of Hardy-Weinberg equilibrium per SV')


# In[45]:


plot_sp (het_freq_hwe_SV_table,
         mt.variant_qc.het_freq_hwe,
         n_het_SV_table,
         mt.variant_qc.n_het,
         'Expected frequency of heterozygous samples under Hardy-Weinberg equilibrium',
         'Number of heterozygous samples per SV')


# **Filter the variants**
# 
# **Step 1 : chr 1-22, X, Y**
# 
# Remove the variants not located on chr 1-22, X/Y
# 
# **Step 2 : Length**
# 
# Remove variants smaller than 50bp (Should be called as indel by the SNV pipeline)
# 
# **Step 3 : Filters**
# 
# **Step 4 : Hail variant QC**
# 
# based on the threasholds repersented on the figures 
# 
# Low_threashold = Mean - 3*StdDev = stat(table) [2]
# 
# High_threashold = Mean + 3*StdDev = stat(table) [3]
# 
# Filters :
# - Mean DP per variant lower than the low threshoold
# - Call rate per variant lower than the low threshoold
# - Allele Number (AN) per variant lower than the low threshoold
# - Number of samples with missing genotype per variant lower than the low threshoold
# 
# Not implemented :  Hardy–Weinberg values

# **Step 1 : chr 1-22, X, Y**

# In[46]:



#intervals = [hl.parse_locus_interval(x) for x in ['X', 'Y', '1-22']]

# Phil 2023-09-07: this is another place where the intervals create an issue with hard-coded expectation of contig name
try:
    contigs = referenceGenome.contigs
except:
    if genome == "GRCh37":
        contigs = [f"{i}" for i in (list(range(1, 23)) + ['X', 'Y'])]
    elif genome =="GRCh38":
        contigs = [f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])]
    else:
        raise ValueError("please enter a valid human genome assemebly value,eg GRCh37")

intervals = [hl.parse_locus_interval(x, reference_genome=referenceGenome) for x in contigs]
print(contigs)
print(intervals)

SV_mt_filtered = hl.filter_intervals(mt, intervals, keep=True)


# In[47]:


n_non_chr = mt.count()[0] - hl.filter_intervals(mt, intervals, keep=True).count()[0]
perc_non_chr = (mt.count()[0]-SV_mt_filtered.count()[0])/mt.count()[0] * 100


# **Step 2 : Length**
# 
# Remove variants smaller than 50bp (Should be called as indel by the SNV pipeline)

# In[48]:


SV_mt_len_filtered = SV_mt_filtered.filter_rows(
    (SV_mt_filtered.info.SVLEN[0] < -50) |
    (SV_mt_filtered.info.SVLEN[0] > 50) 
)


# In[49]:


n_removed_length = SV_mt_filtered.count()[0]-SV_mt_len_filtered.count()[0]
perc_SV_len_filtered = (SV_mt_filtered.count()[0]-SV_mt_len_filtered.count()[0])/SV_mt_filtered.count()[0] * 100


# **Step 3 : Filters**
# 
# Remove variants with filters tag
# 
# 1. Keep only variants with {none} tag
# 2. Generate stats

# In[50]:


SV_mt_len_filters_filtered = SV_mt_len_filtered.filter_rows(
    (SV_mt_len_filtered.filters) == hl.set(["none"]) )


# In[51]:


n_removed_filters = SV_mt_len_filtered.count()[0]-SV_mt_len_filters_filtered.count()[0]
perc_SV_len_filters_filtered = (SV_mt_len_filtered.count()[0]-SV_mt_len_filters_filtered.count()[0])/SV_mt_len_filtered.count()[0] * 100


# **Step 4 : Hail filters**
#     
# - AN (Allele Number) : The SV with a very low AN were not genotyped correctly in the cohort and will not give a strong frequency estimation, should be removed.
# AN threshold : 
# Max AN = number of samples*2
# Min AN --> 75% of max AN
# 
# - Not implemented : DP, quality

# In[52]:


SV_mt_len_filters_filtered= hl.variant_qc(SV_mt_len_filters_filtered)


# **AN (Allele Number)**
# 
# The SV with a very low AN were not genotyped correctly in the cohort and will not give a strong frequency estimation, should be removed.
# 
# AN threshold : 
# 
# Max AN = number of samples*2
# 
# Min AN --> 75% of max AN

# In[53]:


min_AN=SV_mt_len_filters_filtered.count()[1]*2*0.75


# In[54]:


SV_mt_len_filters_AN_filtered=SV_mt_len_filters_filtered.filter_rows(
    (SV_mt_len_filters_filtered.variant_qc.AN) >  min_AN)


# In[55]:


n_removed_AN = SV_mt_len_filters_filtered.count()[0]-SV_mt_len_filters_AN_filtered.count()[0]
perc_SV_len_filters_AN_filtered = (SV_mt_len_filters_filtered.count()[0]-SV_mt_len_filters_AN_filtered.count()[0])/SV_mt_len_filters_filtered.count()[0] * 100


# **Future updates**
# 
# Additional filtrations can be applied to remove false SV from the dataset

# **Calculation of the stats post filtering for report**

# In[56]:


n_SV_removed = mt.count()[0]-SV_mt_len_filters_AN_filtered.count()[0]
perc_SV_removed = (mt.count()[0]-SV_mt_len_filters_AN_filtered.count()[0])/mt.count()[0] * 100
n_var_after_filtering = SV_mt_len_filters_AN_filtered.count()[0]


# In[57]:


DEL_filtered_mt = SV_mt_len_filters_AN_filtered.filter_rows(SV_mt_len_filters_AN_filtered.info.SVTYPE == "DEL")
DUP_filtered_mt = SV_mt_len_filters_AN_filtered.filter_rows(SV_mt_len_filters_AN_filtered.info.SVTYPE == "DUP")
INV_filtered_mt = SV_mt_len_filters_AN_filtered.filter_rows(SV_mt_len_filters_AN_filtered.info.SVTYPE == "INV")
INS_filtered_mt = SV_mt_len_filters_AN_filtered.filter_rows(SV_mt_len_filters_AN_filtered.info.SVTYPE == "INS")


# In[58]:


n_del_filtered = DEL_filtered_mt.count_rows()
n_dup_filtered = DUP_filtered_mt.count_rows()
n_ins_filtered = INS_filtered_mt.count_rows()
n_inv_filtered = INV_filtered_mt.count_rows()
n_SV_tot_filtered = SV_mt_len_filters_AN_filtered.count_rows()


# In[59]:


n_del_filtered


# In[60]:


n_SV_tot_filtered


# In[61]:


n_SV_tot_filtered_smoove = SV_mt_len_filters_AN_filtered.aggregate_rows(hl.agg.count_where(hl.len(SV_mt_len_filters_AN_filtered.info.SU)>0))
n_SV_tot_filtered_manta = n_SV_tot_filtered-n_SV_tot_filtered_smoove

n_SV_del_filtered_smoove = DEL_filtered_mt.aggregate_rows(hl.agg.count_where(hl.len(DEL_filtered_mt.info.SU)>0))
n_SV_del_filtered_manta = n_del_filtered-n_SV_del_filtered_smoove

n_SV_dup_filtered_smoove = DUP_filtered_mt.aggregate_rows(hl.agg.count_where(hl.len(DUP_filtered_mt.info.SU)>0))
n_SV_dup_filtered_manta = n_dup_filtered-n_SV_dup_filtered_smoove

n_SV_ins_filtered_smoove = INS_filtered_mt.aggregate_rows(hl.agg.count_where(hl.len(INS_filtered_mt.info.SU)>0))
n_SV_ins_filtered_manta = n_ins_filtered-n_SV_ins_filtered_smoove

n_SV_inv_filtered_smoove = INV_filtered_mt.aggregate_rows(hl.agg.count_where(hl.len(INV_filtered_mt.info.SU)>0))
n_SV_inv_filtered_manta = n_inv_filtered-n_SV_inv_filtered_smoove


# In[62]:


min_len_filtered_tot = SV_mt_len_filters_AN_filtered.aggregate_rows(hl.agg.min(SV_mt_len_filters_AN_filtered.info.SVLEN[0]))
max_len_filtered_tot = SV_mt_len_filters_AN_filtered.aggregate_rows(hl.agg.max(SV_mt_len_filters_AN_filtered.info.SVLEN[0]))

min_len_filtered_del = DEL_filtered_mt.aggregate_rows(hl.agg.min(DEL_filtered_mt.info.SVLEN[0]))
max_len_filtered_del = DEL_filtered_mt.aggregate_rows(hl.agg.max(DEL_filtered_mt.info.SVLEN[0]))

min_len_filtered_dup = DUP_filtered_mt.aggregate_rows(hl.agg.min(DUP_filtered_mt.info.SVLEN[0]))
max_len_filtered_dup = DUP_filtered_mt.aggregate_rows(hl.agg.max(DUP_filtered_mt.info.SVLEN[0]))

min_len_filtered_ins = INS_filtered_mt.aggregate_rows(hl.agg.min(INS_filtered_mt.info.SVLEN[0]))
max_len_filtered_ins = INS_filtered_mt.aggregate_rows(hl.agg.max(INS_filtered_mt.info.SVLEN[0]))

min_len_filtered_inv = INV_filtered_mt.aggregate_rows(hl.agg.min(INV_filtered_mt.info.SVLEN[0]))
max_len_filtered_inv = INV_filtered_mt.aggregate_rows(hl.agg.max(INV_filtered_mt.info.SVLEN[0]))


# In[63]:


len_hist_filtered = SV_mt_len_filters_AN_filtered.aggregate_entries(hl.expr.aggregators.hist(SV_mt_len_filters_AN_filtered.info.SVLEN[0],
            min_len_filtered_tot, max_len_filtered_tot, 100))
p = hl.plot.histogram(len_hist_filtered, legend='length (bp)', title='Length Histogram (All SV type after filtration steps)', log=True)
show(p)


# In[64]:


len_hist_filtered = SV_mt_len_filters_AN_filtered.aggregate_entries(hl.expr.aggregators.hist(SV_mt_len_filters_AN_filtered.info.SVLEN[0],
            -1000, 1000, 100))
p = hl.plot.histogram(len_hist_filtered, legend='length (bp)', title='Length Histogram (All SV type after filtration steps)', log=True)
show(p)


# **SV QC Report**
# 
# Write the report of the number of filtered out variants and the reason they were filtered out

# In[65]:


def report_stats():
    """
    Generate output report with basic stats.
    """
    out_stats = hl.hadoop_open("SV_QC_report.txt", "w")
    # Report numbers of filtered SV
    out_stats.write(
        f"Prior to variant filtering\n\n"
        
        f"Number of SV identified (all types, all callers): {n_SV_tot}\n"
        f"     Variant caller\n"
        f"          Number of SV identified aith Smoove: {n_SV_tot_smoove}\n"
        f"          Number of SV identified with Manta: {n_SV_tot_manta}\n"
        f"     Length\n"
        f"          Length of the shortest SV (all types): {min_len_tot} bp\n"
        f"          Length of the longest SV (all types): {max_len_tot} bp\n"
        f"Deletions (DEL)\n"
        f"     Number of deletions (DEL) identified (all callers): {n_del}\n"
        f"          Number of deletions (DEL) identified (Manta): {n_SV_del_manta}\n"
        f"          Number of deletions (DEL) identified (Smoove): {n_SV_del_smoove}\n"
        f"          Length of the shortest deletion: {min_len_del} bp\n"
        f"          Length of the longest deletion: {max_len_del} bp\n"
        f"Duplications (DUP)\n"
        f"     Number of duplications (DUP) identified (all callers): {n_dup}\n"
        f"          Number of duplications (DUP) identified (Manta): {n_SV_dup_manta}\n"
        f"          Number of duplications (DUP) identified (Smoove): {n_SV_dup_smoove}\n"
        f"          Length of the shortest duplication: {min_len_dup} bp\n"
        f"          Length of the longest duplication: {max_len_dup} bp\n"
        f"Insertions (INS)\n"
        f"     Number of insertions (INS) identified (all callers): {n_ins}\n"
        f"          Number of insertions (INS) identified (Manta): {n_SV_ins_manta}\n"
        f"          Number of insertions (INS) identified (Smoove): {n_SV_ins_smoove}\n"
        f"          Length of the shortest insertion: {min_len_ins} bp\n"
        f"          Length of the longest insertion: {max_len_ins} bp\n"
        f"Inversions (INV)\n"
        f"     Number of inversions (INV) identified (all callers): {n_inv}\n"
        f"          Number of inversions (INV) identified (Manta): {n_SV_inv_manta}\n"
        f"          Number of inversions (INV) identified (Smoove): {n_SV_inv_smoove}\n"
        f"          Length of the shortest inversion: {min_len_inv} bp\n"
        f"          Length of the longest inversion: {max_len_inv} bp\n\n\n"
        
        f"Variant filtering\n\n"
        f"Filer 1 : Location\n"        
        f"     Number of SV removed because of location (out of chr1-22, X, Y) : {n_non_chr}\n"
        f"     Percentage of SV removed because of location (out of chr1-22, X, Y) : {perc_non_chr} %\n"
        f"Filer 2 : Length\n"        
        f"     Number of SV removed because of length (-50bp < length < 50bp) : {n_removed_length}\n"
        f"     Percentage of SV removed because of length (-50bp < length < 50bp) : {perc_SV_len_filtered} %\n"
        f"Filer 3 : Filters\n"        
        f"     Number of SV removed because of filters (non PASS) : {n_removed_filters}\n"
        f"     Percentage of SV removed because of filters (non PASS) : {perc_SV_len_filters_filtered} %\n"
        f"Filer 4 : Hail QC, AN (Allele number)\n"
        f"     AN threshold for filtering : {min_AN}\n"        
        f"     Number of SV removed because of AN threshold : {n_removed_AN}\n"
        f"     Percentage of SV removed because of AN threshold : {perc_SV_len_filters_AN_filtered} %\n"
        f"Conclusion\n"
        f"     Total number of SV removed : {n_SV_removed}\n"
        f"     Percentage of the SV filtered out: {perc_SV_removed} %\n\n\n"
        
        f"After variant filtering\n\n"
        f"Number of SV remaining (all types, all callers): {n_SV_tot_filtered}\n"
        f"     Variant caller\n"
        f"          Number of SV remaining with Smoove: {n_SV_tot_filtered_smoove}\n"
        f"          Number of SV remaining with Manta: {n_SV_tot_filtered_manta}\n"
        f"     Length\n"
        f"          Length of the shortest SV (all types): {min_len_filtered_tot} bp\n"
        f"          Length of the longest SV (all types): {max_len_filtered_tot} bp\n"
        f"Deletions (DEL)\n"
        f"     Number of deletions (DEL) remaining (all callers): {n_del_filtered}\n"
        f"          Number of deletions (DEL) remaining (Manta): {n_SV_del_filtered_manta}\n"
        f"          Number of deletions (DEL) remaining (Smoove): {n_SV_del_filtered_smoove}\n"
        f"          Length of the shortest deletion: {min_len_filtered_del} bp\n"
        f"          Length of the longest deletion: {max_len_filtered_del} bp\n"
        f"Duplications (DUP)\n"
        f"     Number of duplications (DUP) remaining (all callers): {n_dup_filtered}\n"
        f"          Number of duplications (DUP) remaining (Manta): {n_SV_dup_filtered_manta}\n"
        f"          Number of duplications (DUP) remaining (Smoove): {n_SV_dup_filtered_smoove}\n"
        f"          Length of the shortest duplication: {min_len_filtered_dup} bp\n"
        f"          Length of the longest duplication: {max_len_filtered_dup} bp\n"
        f"Insertions (INS)\n"
        f"     Number of insertions (INS) remaining (all callers): {n_ins_filtered}\n"
        f"          Number of insertions (INS) remaining (Manta): {n_SV_ins_filtered_manta}\n"
        f"          Number of insertions (INS) remaining (Smoove): {n_SV_ins_filtered_smoove}\n"
        f"          Length of the shortest insertion: {min_len_filtered_ins} bp\n"
        f"          Length of the longest insertion: {max_len_filtered_ins} bp\n"
        f"Inversions (INV)\n"
        f"     Number of inversions (INV) remaining (all callers): {n_inv_filtered}\n"
        f"          Number of inversions (INV) remaining (Manta): {n_SV_inv_filtered_manta}\n"
        f"          Number of inversions (INV) remaining (Smoove): {n_SV_inv_filtered_smoove}\n"
        f"          Length of the shortest inversion: {min_len_filtered_inv} bp\n"
        f"          Length of the longest inversion: {max_len_filtered_inv} bp\n\n\n"               
    )
    out_stats.close()


# In[66]:


report_stats()


# **Calculate sex specific frequencies**
# 
# Sex is defined using F-stat in Hail_sample_QC or file with sample sex can be loaded by user
# 
# Calculate AF, AC, AN and number of homozygotes
# 
# Code adapted from gnomAD : https://github.com/broadinstitute/gnomad_qc/blob/main/gnomad_qc/v3/annotations/generate_freq_data.py

# In[67]:


#SV_mt_filtered_export = SV_mt_len_filters_AN_filtered.annotate_cols(**sex_table[SV_mt_len_filters_AN_filtered.s])
SV_mt_filtered_export=SV_mt_len_filters_AN_filtered

# In[68]:


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

    mt.aggregate_cols(
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



# In[69]:


SV_mt_filtered_export = annotate_freq(
                SV_mt_filtered_export,
                sex_expr=SV_mt_filtered_export.sex,
            )


# In[70]:


SV_mt_filtered_export = SV_mt_filtered_export.annotate_rows(
    info = SV_mt_filtered_export.info.annotate(AC_tot_XX_XY=SV_mt_filtered_export.freq.AC,
                                             AF_tot_XX_XY=SV_mt_filtered_export.freq.AF,
                                             AN_tot_XX_XY=SV_mt_filtered_export.freq.AN,
                                             hom_tot_XX_XY=SV_mt_filtered_export.freq.homozygote_count)
                     )


# In[71]:


SV_mt_filtered_export = SV_mt_filtered_export.annotate_rows(info=SV_mt_filtered_export.info.drop("AC", "AN"))


# In[72]:


SV_mt_filtered_export_no_geno = SV_mt_filtered_export.rows()


# **Export files of interest**
# 
# - Variants passing QC, with variant frequencies per sex (total, XX, XY) and individual genotypes
# 
#         File name : SV_filtered_with_geno.vcf.bgz
# 
# 
# 
# - Variants passing QC, with variant frequencies per sex (total, XX, XY), without individual genotypes
# 
# 
#         File name : SV_filtered_frequ_only.vcf.bgz

# In[74]:


hl.export_vcf(SV_mt_filtered_export, 'SV_filtered_with_geno.vcf.bgz', tabix=True)


# In[75]:


hl.export_vcf(SV_mt_filtered_export_no_geno, 'SV_filtered_frequ_only.vcf.bgz', tabix=True)


# In[76]:


SV_mt_filtered_export_no_geno.show()


# **SV vizualisation**
# 
# To check the filters, some SV were visually checked using IGV
# 
# Below is the code to visualize SV on IGV using Jupyter Notebook

# import igv_notebook
# igv_notebook.init()

# #Get the coordinates of Del with with high AN and low DP
# low_DP_variants = SV_mt_len_filters_filtered_AN.filter_rows(
#     (SV_mt_len_filters_AN_filtered.variant_qc.dp_stats.mean < 4) &
#     (SV_mt_len_filters_AN_filtered.info.AN > 150) &
#     (SV_mt_len_filters_AN_filtered.info.SVTYPE == 'DEL')
#     )

# low_DP_variants.variant_qc.show(30)

# Sample = "NA18534"
# b = igv_notebook.Browser(
#     {
#         "genome": "hg19",
#         "locus": "chr1:913326-913416",
#         "tracks": [
#             {
#                 "name": Sample,
#                 "path": f"<Path>/BAM/{Sample}_sorted.bam",
#                 "indexPath": f"<Path>/BAM/{Sample}_sorted.bam.bai",
#                 "format": "bam",
#                 "type": "alignment"
#             }
#         ]
#     })

# In[ ]:




