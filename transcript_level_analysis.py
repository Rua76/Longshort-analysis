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

#selecting protein_coding transcripts from both reference gtf and target gtf
ref_trans = ref[ref['Feature'].isin(['transcript'])]
ref_trans = ref_trans[ref_trans['gene_type'].isin(['protein_coding'])]
target_df = target[target['Chromosome'].isin(list(ref['Chromosome']))]
target_df = target_df[target_df['Feature'].isin(['transcript'])]

#creating dataframe for storing transcript information
trans_mapping = pd.DataFrame(columns=['ENSG_chr', 'ENSG_gene', 'ENST_transcript', 'ENST_Start', 'ENST_End','ENST_Strand', 'ST_chr',  'ST_gene', 'ST_transcript', 'ST_Start', 'ST_End', 'ST_TPM'], index=ref.index, data='0')
trans_mapping['ENSG_chr'] = ref_trans['Chromosome']
trans_mapping['ENSG_gene'] = ref_trans['gene_id']
trans_mapping['ENST_transcript'] = ref_trans['transcript_id']
trans_mapping['ENST_Start'] = ref_trans['Start']
trans_mapping['ENST_End'] = ref_trans['End']
trans_mapping['ENST_Strand'] = ref_trans['Strand']

trans_mapping = trans_mapping.dropna(how = "any")

trans_reference = trans_mapping

#revised matching function for transcript matching
#10 bp window is allowed for matching 
def trans_match(row):
    sub = ref_trans[ref_trans['Chromosome'] == row['Chromosome']]
    sub = sub[sub['Strand'] == row['Strand']]


    sub = sub[sub['Start']  > row['Start'] - 10]
    sub = sub[sub['End'] < row['End'] + 10 ]
#row name is the index of target, sub.index.values[0] is the index of ref_trans
    return [sub.index.values[0] if len(sub.index.values) != 0 else -1, row.name]

def bufferize(target):
    buffer = []
    for index, row in target.iterrows():
        buffer.append(row)
    return buffer

import multiprocessing
with multiprocessing.Pool(40) as pool:
    res = pool.map_async(trans_match, bufferize(target_df)).get()

match_df = {}
for _res in res:
    if _res[0] != -1:

        trans_mapping.at[_res[0], 'ST_chr'] = target_df.loc[_res[1]]['Chromosome']
        trans_mapping.at[_res[0], 'ST_gene'] = target_df.loc[_res[1]]['gene_id']
        trans_mapping.at[_res[0], 'ST_transcript'] = target_df.loc[_res[1]]['transcript_id']
        trans_mapping.at[_res[0], 'ST_TPM'] = target_df.loc[_res[1]]['TPM']
        trans_mapping.at[_res[0], 'ST_Start'] = target_df.loc[_res[1]]['Start']
        trans_mapping.at[_res[0], 'ST_End'] = target_df.loc[_res[1]]['End']

#filtering transcripts
trans_mapping = trans_mapping[trans_mapping['ST_gene'] != '0']

trans_mapping.to_csv('trans_reference.csv')

trans_merging = pd.merge(trans_mapping, trans_reference, left_index=True, right_index=True)
trans_merging_simp = trans_merging[['ST_chr_x','ST_gene_x','ST_transcript_x', 'ST_Start_x', 'ST_End_x', 'ST_TPM_x','ENSG_chr_y', 'ENSG_gene_y' ,'ENST_transcript_y', 'ENST_Start_y','ENST_End_y' , 'ENST_Strand_y' ]]
trans_merging_simp.columns = ['ST_chr','ST_gene','ST_transcript', 'ST_Start', 'ST_End', 'ST_TPM','ENSG_chr', 'ENSG_gene' ,'ENST_transcript', 'ENST_Start','ENST_End' , 'ENST_Strand' ]

#grouping analysis based on the number of transcripts mapped on each gene
trans_counting = trans_merging_simp.groupby(['ENSG_gene']).count()

#save file
trans_merging_simp.to_csv('trans_merged.csv')

ref_idx = fullabundance.groupby(['gene_id'])['TPM'].transform(max) == fullabundance['TPM']
max_count_ref = fullabundance[ref_idx]
merge1 = pd.merge(left=longshort, right=trans_merging_simp, how='inner', left_on = ['target_id'], right_on = ['ST_transcript'])
trans_merged_with_abundance = pd.merge(left=merge1, right=max_count_ref, how='inner', left_on = ['ENSG_gene'], right_on = ['gene_id'])
trans_merged_with_abundance.to_csv('trans_merged_with_abundance.csv')


