import pyranges as pr
import pandas as pd
import seaborn as sns; sns.set_theme(color_codes=True)
import numpy as np
from hilearn.plot import corr_plot

#read file generated from transcript_based_analysis.py
trans_merged_with_abundance = pd.read_csv('trans_merged_with_abundance.csv')

trans_merged_with_abundance = trans_merged_with_abundance[ ['ST_chr','ST_gene','ST_transcript', 'ST_Start', 'ST_End', 'ST_TPM','ENSG_chr', 'ENSG_gene' ,'ENST_transcript', 'ENST_Start','ENST_End' , 'ENST_Strand' , 'est_counts', 'onestep_est_count']]

#grouping
trans_counting = trans_merged_with_abundance.groupby(['ENSG_gene']).count()

#selecting genes with specific number of transcripts matched
#1-3 matching
trans_counting_1_3 = trans_counting.loc[((trans_counting['ENST_transcript'] >= 1) & (trans_counting['ENST_transcript'] <= 3))]
#3-5 matching
trans_counting_3_5 = trans_counting.loc[((trans_counting['ENST_transcript'] >= 3) & (trans_counting['ENST_transcript'] <= 5))]
#above 5 matching
trans_counting_5_ = trans_counting.loc[(trans_counting['ENST_transcript'] >5)]

#sample procedure to generate scatter plot for abundance comparison, 1-3 matching
oneto_3 = trans_merged_with_abundance.merge(trans_counting_1_3, left_on = 'ENSG_gene', right_index = True)

oneto3_idx = oneto_3.groupby(['ENSG_gene'])['ST_TPM_x'].transform(max) == oneto_3['ST_TPM_x']
oneto_3 = oneto_3[oneto3_idx]

oneto_3_s = oneto_3.loc[:,['est_counts_x', 'onestep_est_count_x']]
oneto_3_s
oneto3_1 = pd.DataFrame(oneto_3_s + 1)
oneto3_1['log10_longshort_est'] = np.log10(oneto3_1['est_counts_x'])
oneto3_1['log10_onestep_est'] = np.log10(oneto3_1['onestep_est_count_x'])
x = oneto3_1.loc[:,'log10_longshort_est'].values
y = oneto3_1.loc[:,'log10_onestep_est'].values
corr_plot(x, y)

#sample procedure to generate scatter plot for abundance comparison, 3-5 matching
three_to_5 = trans_merged_with_abundance.merge(trans_counting_3_5, left_on = 'ENSG_gene', right_index = True)
three_to_5_idx = three_to_5.groupby(['ENSG_gene'])['ST_TPM_x'].transform(max) == three_to_5['ST_TPM_x']

three_to_5s = three_to_5.loc[:,['est_counts_x', 'onestep_est_count_x']]
three_to_5_1 = pd.DataFrame(three_to_5s + 1)
three_to_5_1['log10_longshort_est'] = np.log10(three_to_5_1['est_counts_x'])
three_to_5_1['log10_onestep_est'] = np.log10(three_to_5_1['onestep_est_count_x'])
x = three_to_5_1.loc[:,'log10_longshort_est'].values
y = three_to_5_1.loc[:,'log10_onestep_est'].values
corr_plot(x, y)

#sample procedure to generate scatter plot for abundance comparison, above 5 matching
above5 = trans_merged_with_abundance.merge(trans_counting_5_, left_on = 'ENSG_gene', right_index = True)
simp = above5.loc[:,['est_counts_x', 'onestep_est_count_x']]

above5_1 = pd.DataFrame(simp + 1)
above5_1['log10_longshort_est'] = np.log10(above5_1['est_counts_x'])
above5_1['log10_onestep_est'] = np.log10(above5_1['onestep_est_count_x'])
x = above5_1.loc[:,'log10_longshort_est'].values
y = above5_1.loc[:,'log10_onestep_est'].values
corr_plot(x, y)

#sample procedure to generate scatter plot for general abundance comparison, all genes
trans_max_idx = trans_merged_with_abundance.groupby(['ENSG_gene'])['ST_TPM'].transform(max) == trans_merged_with_abundance['ST_TPM']
trans_max = trans_merged_with_abundance[trans_max_idx]

trans_simp = trans_max.loc[:,['est_counts', 'onestep_est_count']]
trans_1 = pd.DataFrame(trans_simp + 1)
trans_1['log10_longshort_est'] = np.log10(trans_1['est_counts'])
trans_1['log10_onestep_est'] = np.log10(trans_1['onestep_est_count'])
x = trans_1.loc[:,'log10_longshort_est'].values
y = trans_1.loc[:,'log10_onestep_est'].values
corr_plot(x, y)
