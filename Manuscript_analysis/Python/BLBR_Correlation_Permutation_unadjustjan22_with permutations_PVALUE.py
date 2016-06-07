# Read in Data
import pandas as pd
#Beta data sample names
from numpy import genfromtxt
#Beta Data beta values
PBMC_BLBR = genfromtxt('PBMC_BLBR_unadjustjan22.csv', delimiter=',',dtype=float, missing_values='NA')
BRAIN7_BLBR = genfromtxt('BRAIN7_BLBR_unadjustjan22.csv', delimiter=',',dtype=float, missing_values='NA')
BRAIN10_BLBR = genfromtxt('BRAIN10_BLBR_unadjustjan22.csv', delimiter=',',dtype=float, missing_values='NA')
BRAIN20_BLBR = genfromtxt('BRAIN20_BLBR_unadjustjan22.csv', delimiter=',',dtype=float, missing_values='NA')

# correlation function
from scipy.stats.stats import spearmanr
import numpy as np

x=len(PBMC_BLBR)
# # correlation function
PBMC_order = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
np.random.shuffle(PBMC_order)

spearman_correlations=[]
for CpG in range(1,x):
    PBMC = PBMC_BLBR[CpG][PBMC_order].astype(float) #shuffle
    BRAIN7 = BRAIN7_BLBR[CpG][1:17].astype(float)
    BRAIN10 = BRAIN10_BLBR[CpG][1:17].astype(float)
    BRAIN20 = BRAIN20_BLBR[CpG][1:17].astype(float)
        
    brain7=spearmanr(PBMC, BRAIN7)[0]
    brain10=spearmanr(PBMC, BRAIN10)[0]
    brain20=spearmanr(PBMC, BRAIN20)[0]
    spearman_correlations.append([brain7,brain10, brain20])

##Save
df = pd.DataFrame(spearman_correlations)
df.to_csv("correlation_unadjustjan22_permutation1_PVAL.csv", sep='\t')


