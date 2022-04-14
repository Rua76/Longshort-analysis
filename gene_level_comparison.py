import pyranges as pr
import pandas as pd
import seaborn as sns; sns.set_theme(color_codes=True)
import numpy as np

from hilearn.plot import corr_plot

#reading gtf files
ref = pr.read_gtf("gencode.v39.annotation.gtf")
ref = ref.df
target = pr.read_gtf("m1.v39.genome.gtf")
target = target.df

#For gene level comparison in our study, we only consider protein_coding genes
ref_df = ref[ref['Feature'].isin(['gene'])]
ref_df = ref_df[ref_df['gene_type'].isin(['protein_coding'])]
#read target transcripts
target_df = target[target['Chromosome'].isin(list(ref['Chromosome']))]
target_df = target_df[target_df['Feature'].isin(['transcript'])]

#creating gene reference dataframe
gene_reference = pd.DataFrame(columns=['ST_chr', 'ST_gene', 'ST_transcript', 'ST_TPM', 'ENSG', 'ENST'], index=target_df.index, data='0')
gene_reference['ST_chr'] = target_df['Chromosome']
gene_reference['ST_gene'] = target_df['gene_id']
gene_reference['ST_transcript'] = target_df['transcript_id']
gene_reference['ST_TPM'] = target_df['TPM']

#define matching function
#We use a matching criteria based on Chromosome, Strand, Start and End positions
#A 500 bp shifting window is allowed for ref_gene & target_transcript matching 
#target_transcript falling into the gene region is counted as matched
def match(row):
    sub = ref_df[ref_df['Chromosome'] == row['Chromosome']]
    sub = sub[sub['Strand'] == row['Strand']]

    sub = sub[sub['Start'] - 500 < row['Start']]
    sub = sub[sub['End'] + 500 > row['End']]

    return [row.name, sub.index.values[0] if len(sub.index.values) != 0 else -1]
  
#buffering
def bufferize(target):
    buffer = []
    for index, row in target.iterrows():
        buffer.append(row)
    return buffer
  
#merging step
import multiprocessing
with multiprocessing.Pool(40) as pool:
    res = pool.map_async(match, bufferize(target_df)).get()

match_df = {}
for _res in res:
    if _res[1] != -1:
        gene_reference.at[_res[0], 'ENSG'] = ref_df.loc[_res[1]]['gene_id']
        gene_reference.at[_res[0], 'ENST'] = ref_df.loc[_res[1]]['transcript_id']

#gene filtering
gene_reference = gene_reference[gene_reference['ENSG'] != '0']

#creating csv file for storing intermediate data
gene_reference.to_csv('gene_reference.csv')

######################################################################################

#user can start here to read intermediate data, avoiding to read gtf files again
gene_reference = pd.read_csv('gene_reference.csv')

#read abundance files here
longshort = pd.read_csv('longshort.comma.csv')
fullabundance = pd.read_csv('fullAOabund.csv')
fullabundance.columns = ['transcript_id', 'gene_id', 'HAVANA_gene_id','HAVANA_trans_id','transcript_name','gene_name','length','transcript_type','_length', 'eff_length', 'onestep_est_count', 'TPM']

#In gene-level comparison, for each gene, we focus on most expressed transcripts for both standard method and our long-short hybrid method
#We select such transcript according to their TPM value
gene_reference.groupby(['ENSG'])['ST_TPM'].transform(max)
tar_idx = gene_reference.groupby(['ENSG'])['ST_TPM'].transform(max) == gene_reference['ST_TPM']
max_count_tar = gene_reference[tar_idx]

ref_idx = fullabundance.groupby(['gene_id'])['TPM'].transform(max) == fullabundance['TPM']
max_count_ref = fullabundance[ref_idx]

#assigning estimated count number to each gene
merge1 = pd.merge(left=longshort, right=max_count_tar, how='inner', left_on = ['target_id'], right_on = ['ST_transcript'])
merge2 = pd.merge(left=merge1, right=max_count_ref, how='inner', left_on = ['ENSG'], right_on = ['gene_id'])
#creating csv file for storing intermediate data
merge2.to_csv('merged_gtf_with_abundance.csv')

#plotting scatter plot with HiLearn
#Here, the estimated count for both methods are added with pseudocount of 1, and log normalized with log10
simp = merge2.loc[:,['est_counts', 'onestep_est_count']]
print (simp)
merge2_1 = pd.DataFrame(simp + 1)
merge2_1['log10_longshort_est'] = np.log10(merge2_1['est_counts'])
merge2_1['log10_onestep_est'] = np.log10(merge2_1['onestep_est_count'])
x = merge2_1.loc[:,'log10_longshort_est'].values
y = merge2_1.loc[:,'log10_onestep_est'].values
corr_plot(x, y)
#scatter plot generated

