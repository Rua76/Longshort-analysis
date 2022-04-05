import pyranges as pr
import pandas as pd
import seaborn as sns; sns.set_theme(color_codes=True)
import numpy as np

from hilearn.plot import corr_plot

#reading gtf files
granno = pr.read_gtf("gencode.v39.annotation.gtf")
dfanno = granno.df

#change read data to dataFrame 
guided_m1 = guided_m1_gr.df

#selecting transcripts
guided_m1_trans = guided_m1[guided_m1['Feature'] == 'transcript']

#counting class codes
guided_m1_trans.groupby('class_code').count()

#reading long-short abundance file. Files containing raw count data can also be used here.
longshort = pd.read_csv('longshort.comma.csv')

#reading and adjusting onestep abundance file
onestep = pd.read_csv('onestep.sec.csv')
onestep.columns = ['transcript_id', 'gene_id', 'HAVANA_gene_id','HAVANA_trans_id','transcript_name','gene_name','length','transcript_type','onestep_est_count']

#attaching long-short estimated count to corresponding transcript by inner join ('target_id' and 'cmp_ref' are transcript id)
guided1merge = pd.merge(left=longshort, right=guided_m1_trans, how='inner', left_on = ['target_id'], right_on = ['cmp_ref'])

#attaching onestep estimated count to corresponding transcript by inner join 
guided2merge = pd.merge(left=guided1merge, right=onestep, left_on=['transcript_id'], right_on=['transcript_id'])
#Getting scatter plot
simpguided = guided2merge.loc[:,['est_counts', 'onestep_est_count']]
#adding pseudocount of 1
simpguided_1 = pd.DataFrame(simpguided + 1)
#get estimated count to log 10 scale
simpguided_1['log10_longshort_est'] = np.log10(simpguided_1['est_counts'])
simpguided_1['log10_onestep_est'] = np.log10(simpguided_1['onestep_est_count'])

x = simpguided_1.loc[:,'log10_longshort_est'].values
y = simpguided_1.loc[:,'log10_onestep_est'].values
corr_plot(x, y)


"""select highest expressed transcript for each gene and analyze the abundance""" 
#I select transcript with highest TPM for each gene, then put them into a new dataFrame
idx = transcript.groupby(['gene_id'])['TPM'].transform(max) == transcript['TPM']
primarytrans = transcript[idx]

#find the inner join of this new dataFrame and our previous count-assigned dataFrame
merged_max = pd.merge(left=primarytrans, right=guided2merge, how='inner', left_on = ['transcript_id'], right_on = ['target_id'])

#generate scatter plot agian
simpmax = merged_max.loc[:,['est_counts', 'onestep_est_count']]
simpmax_1 = pd.DataFrame(simpmax + 1)
simpmax_1['log10_longshort_est'] = np.log10(simpmax_1['est_counts'])
simpmax_1['log10_onestep_est'] = np.log10(simpmax_1['onestep_est_count'])
x = simpmax_1.loc[:,'log10_longshort_est'].values
y = simpmax_1.loc[:,'log10_onestep_est'].values
corr_plot(x, y)
