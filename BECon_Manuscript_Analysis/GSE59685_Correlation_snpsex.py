
# Read in Data
import pandas as pd
#Beta data sample names
from numpy import genfromtxt
#Beta Data beta values
WB_BLBR = genfromtxt('blood_BLBR_snpsex_Oct19.csv', delimiter=',',dtype=float, missing_values='NA')
CB_BLBR = genfromtxt('cerebellum_BLBR_snpsex_Oct19.csv', delimiter=',',dtype=float, missing_values='NA')
EC_BLBR = genfromtxt('EC_BLBR_snpsex_Oct19.csv', delimiter=',',dtype=float, missing_values='NA')
FC_BLBR = genfromtxt('FC_BLBR_snpsex_Oct19.csv', delimiter=',',dtype=float, missing_values='NA')
STG_BLBR = genfromtxt('STG_BLBR_snpsex_Oct19.csv', delimiter=',',dtype=float, missing_values='NA')


# correlation function
from scipy.stats.stats import spearmanr
import numpy as np

x=len(WB_BLBR)
# # correlation function

spearman_correlations=[]
for CpG in range(1,x):
    WB = WB_BLBR[CpG].astype(float) #shuffle
    CB = CB_BLBR[CpG].astype(float)
    EC = EC_BLBR[CpG].astype(float)
    FC = FC_BLBR[CpG].astype(float)
    STG = STG_BLBR[CpG].astype(float)


    CB=spearmanr(WB, CB)[0]
    EC=spearmanr(WB, EC)[0]
    FC=spearmanr(WB, FC)[0]
    STG=spearmanr(WB, STG)[0]
    spearman_correlations.append([CB, EC, FC, STG])

##Save
df = pd.DataFrame(spearman_correlations)
df.to_csv("correlation_snpsex_celladjusted_oct19.csv", sep='\t')
