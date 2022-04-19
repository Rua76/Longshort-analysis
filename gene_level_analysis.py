import pyranges as pr
import pandas as pd
import seaborn as sns; sns.set_theme(color_codes=True)
import numpy as np
import hilearn
import matplotlib.pyplot as plt

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
df_onestep = pd.read_csv(
    '/storage/yhhuang/users/yhsz/sampledata/m1_sample/onestep/abundance.tsv', sep='\t')
df_longshort = pd.read_csv(
    '/storage/yhhuang/users/yhsz/sampledata/m1_sample/longshort/abundance.tsv', sep='\t')

df_onestep.shape

df_onestep['tran_id'] = [x.split('|')[0] for x in df_onestep['target_id']]
df_onestep['gene_id'] = [x.split('|')[1] for x in df_onestep['target_id']]
df_onestep['gene_type'] = [x.split('|')[-2] for x in df_onestep['target_id']]

df_onestep = df_onestep[df_onestep['gene_type'] == 'protein_coding']

df_tmp = df_onestep[['gene_id', 'tpm']]
df_1step_ENSG = df_tmp.groupby(['gene_id']).max()
#histogram of tmp distribution
plt.hist(np.log2(df_1step_ENSG['tpm'] + 1), bins=100)
plt.show()

#abundance assignment
st_idx = hilearn.match(gene_reference['ST_transcript'].values, 
                       df_longshort['target_id'].values)
gene_reference['st_tpm_new'] = df_longshort['tpm'].values[st_idx]
df_st_ENSG = gene_reference[['ENSG', 'st_tpm_new']]
df_st_ENSG = df_st_ENSG.groupby(['ENSG']).max()

df_comb_ENSG = pd.concat([df_st_ENSG, df_1step_ENSG], axis=1, join="inner")

fig = plt.figure(dpi=80)
hilearn.plot.corr_plot(
    np.log2(df_comb_ENSG['st_tpm_new'].values + 1), 
    np.log2(df_comb_ENSG['tpm'].values + 1)
)
plt.xlabel('log2(TPM+1), StringTie gene')
plt.ylabel('log2(TPM+1), Ensembl gene')
plt.show()
