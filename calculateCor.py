# -*- coding: utf-8 -*-
import os
os.chdir('working passway')
import numpy as np, pandas as pd, copy as cp
import sys
import re
from scipy.stats import spearmanr
from pingouin import partial_corr

def sjN(x):
    x=str(x)
    if x != 'nan':
        sj=x.split(';')
        sjs=[int(x) for x in sj if int(x)>0]
        sjN=len(sjs)
    else:
        sjN=0
    return(sjN)

datDir='data passway'

sps=['hu','ma','mo']
clns=['species','tissue',
     'geneWithSjTpm1','geneWithSjTpm1_totSj','geneWithSjTpm1_totLsj','geneWithSjTpm1_totBsj','geneWithSjTpm1_totBsjRatio',
      'geneWithSjTpm1_corC_bet_br_sj','geneWithSjTpm1_corP_bet_br_sj',
      'geneWithSjTpmIntron1_corC_bet_brPerAnnIn_sj','geneWithSjTpmIntron1_corP_bet_brPerAnnIn_sj',
      'geneWithSjTpmIntron1_corC_bet_brPerObIn_sj','geneWithSjTpmIntron1_corP_bet_brPerObIn_sj',
      'geneWithSjTpm1_corC_bet_br_sj_cltAnnIn_partial','geneWithSjTpm1_corP_bet_br_sj_cltAnnIn_partial',
      'geneWithSjTpm1_corC_bet_br_sj_cltObIn_partial','geneWithSjTpm1_corP_bet_br_sj_cltObIn_partial',
      
      'geneWithSjTpm1_corC_bet_br_tpm','geneWithSjTpm1_corP_bet_br_tpm',
      'geneWithSjTpmIntron1_corC_bet_brPerAnnIn_tpm','geneWithSjTpmIntron1_corP_bet_brPerAnnIn_tpm',
      'geneWithSjTpmIntron1_corC_bet_brPerObIn_tpm','geneWithSjTpmIntron1_corP_bet_brPerObIn_tpm',
      'geneWithSjTpm1_corC_bet_br_tpm_cltAnnIn_partial','geneWithSjTpm1_corP_bet_br_tpm_cltAnnIn_partial',
      'geneWithSjTpm1_corC_bet_br_tpm_cltObIn_partial','geneWithSjTpm1_corP_bet_br_tpm_cltObIn_partial',
      
      'geneWithSjTpm1_corC_bet_br_annIn','geneWithSjTpm1_corP_bet_br_annIn',
      'geneWithSjTpmIntron1_corC_bet_brPerAnnIn_annIn','geneWithSjTpmIntron1_corP_bet_brPerAnnIn_annIn',
      
      'geneWithSjTpm1_corC_bet_br_obIn','geneWithSjTpm1_corP_bet_br_obIn',
      'geneWithSjTpmIntron1_corC_bet_brPerObIn_obIn','geneWithSjTpmIntron1_corP_bet_brPerObIn_obIn',
      
      'geneWithBsjTpm1','geneWithBsjTpm1_totSj','geneWithBsjTpm1_totLsj','geneWithBsjTpm1_totBsj','geneWithBsjTpm1_totBsjRatio',
      'geneWithBsjTpm1_corC_bet_br_sj','geneWithBsjTpm1_corP_bet_br_sj',
      'geneWithBsjTpmIntron1_corC_bet_brPerAnnIn_sj','geneWithBsjTpmIntron1_corP_bet_brPerAnnIn_sj',
      'geneWithBsjTpmIntron1_corC_bet_brPerObIn_sj','geneWithBsjTpmIntron1_corP_bet_brPerObIn_sj',
      'geneWithBsjTpm1_corC_bet_br_sj_cltAnnIn_partial','geneWithBsjTpm1_corP_bet_br_sj_cltAnnIn_partial',
      'geneWithBsjTpm1_corC_bet_br_sj_cltObIn_partial','geneWithBsjTpm1_corP_bet_br_sj_cltObIn_partial',
      
      'geneWithBsjTpm1_corC_bet_br_tpm','geneWithBsjTpm1_corP_bet_br_tpm',
      'geneWithBsjTpmIntron1_corC_bet_brPerAnnIn_tpm','geneWithBsjTpmIntron1_corP_bet_brPerAnnIn_tpm',
      'geneWithBsjTpmIntron1_corC_bet_brPerObIn_tpm','geneWithBsjTpmIntron1_corP_bet_brPerObIn_tpm',
      'geneWithBsjTpm1_corC_bet_br_tpm_cltAnnIn_partial','geneWithBsjTpm1_corP_bet_br_tpm_cltAnnIn_partial',
      'geneWithBsjTpm1_corC_bet_br_tpm_cltObIn_partial','geneWithBsjTpm1_corP_bet_br_tpm_cltObIn_partial',
      
      'geneWithBsjTpm1_corC_bet_br_annIn','geneWithBsjTpm1_corP_bet_br_annIn',
      'geneWithBsjTpmIntron1_corC_bet_brPerAnnIn_annIn','geneWithBsjTpmIntron1_corP_bet_brPerAnnIn_annIn',
      
      'geneWithBsjTpm1_corC_bet_br_obIn','geneWithBsjTpm1_corP_bet_br_obIn',
      'geneWithBsjTpmIntron1_corC_bet_brPerObIn_obIn','geneWithBsjTpmIntron1_corP_bet_brPerObIn_obIn'
      ]
cort = pd.DataFrame(columns=clns)
z = 0

for i in range(0,3,1):
	dat = pd.read_table(datDir+sps[i]+'TpmRdsSjs_addDown.xls',low_memory=False)
	proDat = dat[dat['gene_type'] == 'protein_coding'].copy()
	tisNum=int((len(proDat.columns)-11)/23)
	for j in range(11,len(proDat.columns),23):
		tisName=proDat.columns[j].split("_")[0]
		tisDat = proDat.iloc[:,[0,7,j,j+5,j+6,j+3]].copy()
		tisDat[tisName+'_totSj'] = tisDat[tisName+'_totLsj'] + tisDat[tisName+'_totBsj']
		tisDatSjTpm1 = tisDat[(tisDat[tisName+'_totSj']>=1) & (tisDat[tisName+'_allTpm']>=1)].copy() 
		tisDatSjTpm1[tisName+'_br'] = tisDatSjTpm1[tisName+'_totBsj']/tisDatSjTpm1[tisName+'_totSj']
		tisDatSjTpm1.iloc[:,1] = tisDatSjTpm1.iloc[:,1]-1
		tisDatSjTpm1.rename(columns={'exNum_lgstTpt':'intronNum_lgstTpt'},inplace=True)
		tisDatSjTpm1[tisName+'_lsjN']=tisDatSjTpm1[tisName+'_lsj'].apply(sjN)
		tisDatSjTpm1.drop([tisName+'_lsj'],axis=1,inplace=True)
		tisDatSjTpm1.loc[tisDatSjTpm1['intronNum_lgstTpt']>0,tisName+'_brPerAnnIn']=tisDatSjTpm1.loc[tisDatSjTpm1['intronNum_lgstTpt']>0,tisName+'_br']/tisDatSjTpm1.loc[tisDatSjTpm1['intronNum_lgstTpt']>0,'intronNum_lgstTpt']
		tisDatSjTpm1.loc[tisDatSjTpm1[tisName+'_lsjN']>0,tisName+'_brPerObIn']=tisDatSjTpm1.loc[tisDatSjTpm1[tisName+'_lsjN']>0,tisName+'_br']/tisDatSjTpm1.loc[tisDatSjTpm1[tisName+'_lsjN']>0,tisName+'_lsjN']
		
		result1 = [len(tisDatSjTpm1),sum(tisDatSjTpm1[tisName+'_totSj']),sum(tisDatSjTpm1[tisName+'_totLsj']),sum(tisDatSjTpm1[tisName+'_totBsj']),sum(tisDatSjTpm1[tisName+'_totBsj'])/sum(tisDatSjTpm1[tisName+'_totSj'])]+\
                   list(spearmanr(tisDatSjTpm1[tisName+'_totSj'],tisDatSjTpm1[tisName+'_br']))+\
                   list(spearmanr(tisDatSjTpm1.loc[tisDatSjTpm1['intronNum_lgstTpt'] > 0,tisName+'_totSj'],tisDatSjTpm1.loc[tisDatSjTpm1['intronNum_lgstTpt'] > 0,tisName+'_brPerAnnIn']))+\
                   list(spearmanr(tisDatSjTpm1.loc[tisDatSjTpm1[tisName+'_lsjN']>0,tisName+'_totSj'],tisDatSjTpm1.loc[tisDatSjTpm1[tisName+'_lsjN']>0,tisName+'_brPerObIn']))+\
                   list(partial_corr(data=tisDatSjTpm1,x=tisName+'_totSj',y=tisName+'_br',covar='intronNum_lgstTpt',method='spearman').iloc[0,[1,5]])+\
                   list(partial_corr(data=tisDatSjTpm1,x=tisName+'_totSj',y=tisName+'_br',covar=tisName+'_lsjN',method='spearman').iloc[0,[1,5]])+\
                   list(spearmanr(tisDatSjTpm1[tisName+'_allTpm'],tisDatSjTpm1[tisName+'_br']))+\
                   list(spearmanr(tisDatSjTpm1.loc[tisDatSjTpm1['intronNum_lgstTpt'] > 0,tisName+'_allTpm'],tisDatSjTpm1.loc[tisDatSjTpm1['intronNum_lgstTpt'] > 0,tisName+'_brPerAnnIn']))+\
                   list(spearmanr(tisDatSjTpm1.loc[tisDatSjTpm1[tisName+'_lsjN']>0,tisName+'_allTpm'],tisDatSjTpm1.loc[tisDatSjTpm1[tisName+'_lsjN']>0,tisName+'_brPerObIn']))+\
                   list(partial_corr(data=tisDatSjTpm1,x=tisName+'_allTpm',y=tisName+'_br',covar='intronNum_lgstTpt',method='spearman').iloc[0,[1,5]])+\
                   list(partial_corr(data=tisDatSjTpm1,x=tisName+'_allTpm',y=tisName+'_br',covar=tisName+'_lsjN',method='spearman').iloc[0,[1,5]])+\
                   list(spearmanr(tisDatSjTpm1['intronNum_lgstTpt'],tisDatSjTpm1[tisName+'_br']))+\
                   list(spearmanr(tisDatSjTpm1.loc[tisDatSjTpm1['intronNum_lgstTpt'] > 0,'intronNum_lgstTpt'],tisDatSjTpm1.loc[tisDatSjTpm1['intronNum_lgstTpt'] > 0,tisName+'_brPerAnnIn']))+\
                   list(spearmanr(tisDatSjTpm1[tisName+'_lsjN'],tisDatSjTpm1[tisName+'_br']))+\
                   list(spearmanr(tisDatSjTpm1.loc[tisDatSjTpm1[tisName+'_lsjN']>0,tisName+'_lsjN'],tisDatSjTpm1.loc[tisDatSjTpm1[tisName+'_lsjN']>0,tisName+'_brPerObIn']))
                   
		tisDatSjTpm1.to_csv(sps[i]+'_'+tisName+'_geneWithSjTpm1.xls',sep="\t",index=False,na_rep='NA')
		
		tisDatBsjTpm1 = tisDatSjTpm1[tisDatSjTpm1[tisName+'_totBsj']>0]
		result2 = [len(tisDatBsjTpm1),sum(tisDatBsjTpm1[tisName+'_totSj']),sum(tisDatBsjTpm1[tisName+'_totLsj']),sum(tisDatBsjTpm1[tisName+'_totBsj']),sum(tisDatBsjTpm1[tisName+'_totBsj'])/sum(tisDatBsjTpm1[tisName+'_totSj'])]+\
                   list(spearmanr(tisDatBsjTpm1[tisName+'_totSj'],tisDatBsjTpm1[tisName+'_br']))+\
                   list(spearmanr(tisDatBsjTpm1.loc[tisDatBsjTpm1['intronNum_lgstTpt'] > 0,tisName+'_totSj'],tisDatBsjTpm1.loc[tisDatBsjTpm1['intronNum_lgstTpt'] > 0,tisName+'_brPerAnnIn']))+\
                   list(spearmanr(tisDatBsjTpm1.loc[tisDatBsjTpm1[tisName+'_lsjN']>0,tisName+'_totSj'],tisDatBsjTpm1.loc[tisDatBsjTpm1[tisName+'_lsjN']>0,tisName+'_brPerObIn']))+\
                   list(partial_corr(data=tisDatBsjTpm1,x=tisName+'_totSj',y=tisName+'_br',covar='intronNum_lgstTpt',method='spearman').iloc[0,[1,5]])+\
                   list(partial_corr(data=tisDatBsjTpm1,x=tisName+'_totSj',y=tisName+'_br',covar=tisName+'_lsjN',method='spearman').iloc[0,[1,5]])+\
                   list(spearmanr(tisDatBsjTpm1[tisName+'_allTpm'],tisDatBsjTpm1[tisName+'_br']))+\
                   list(spearmanr(tisDatBsjTpm1.loc[tisDatBsjTpm1['intronNum_lgstTpt'] > 0,tisName+'_allTpm'],tisDatBsjTpm1.loc[tisDatBsjTpm1['intronNum_lgstTpt'] > 0,tisName+'_brPerAnnIn']))+\
                   list(spearmanr(tisDatBsjTpm1.loc[tisDatBsjTpm1[tisName+'_lsjN']>0,tisName+'_allTpm'],tisDatBsjTpm1.loc[tisDatBsjTpm1[tisName+'_lsjN']>0,tisName+'_brPerObIn']))+\
                   list(partial_corr(data=tisDatBsjTpm1,x=tisName+'_allTpm',y=tisName+'_br',covar='intronNum_lgstTpt',method='spearman').iloc[0,[1,5]])+\
                   list(partial_corr(data=tisDatBsjTpm1,x=tisName+'_allTpm',y=tisName+'_br',covar=tisName+'_lsjN',method='spearman').iloc[0,[1,5]])+\
                   list(spearmanr(tisDatBsjTpm1['intronNum_lgstTpt'],tisDatBsjTpm1[tisName+'_br']))+\
                   list(spearmanr(tisDatBsjTpm1.loc[tisDatBsjTpm1['intronNum_lgstTpt'] > 0,'intronNum_lgstTpt'],tisDatBsjTpm1.loc[tisDatBsjTpm1['intronNum_lgstTpt'] > 0,tisName+'_brPerAnnIn']))+\
                   list(spearmanr(tisDatBsjTpm1[tisName+'_lsjN'],tisDatBsjTpm1[tisName+'_br']))+\
                   list(spearmanr(tisDatBsjTpm1.loc[tisDatBsjTpm1[tisName+'_lsjN']>0,tisName+'_lsjN'],tisDatBsjTpm1.loc[tisDatBsjTpm1[tisName+'_lsjN']>0,tisName+'_brPerObIn']))
                   
		tisDatBsjTpm1.to_csv(sps[i]+'_'+tisName+'_geneWithBsjTpm1.xls',sep="\t",index=False,na_rep='NA')
		
		cort.loc[z] = [sps[i],tisName]+result1+result2
		z = z + 1
cort.to_csv('corrs.xls',sep="\t",index=False,na_rep='NA')