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

print PBMC_BLBR[1][1:17].astype(float)
print BRAIN7_BLBR[1][1:17].astype(float)
print spearmanr(PBMC_BLBR[1][1:17].astype(float), BRAIN7_BLBR[1][1:17].astype(float))
print spearmanr(PBMC_BLBR[1].astype(float), BRAIN7_BLBR[1].astype(float))[0]

print ([spearmanr(PBMC_BLBR[1][1:17].astype(float), BRAIN7_BLBR[1][1:17].astype(float))[0],spearmanr(PBMC_BLBR[1][1:17].astype(float), BRAIN10_BLBR[1][1:17].astype(float))[0], spearmanr(PBMC_BLBR[1][1:17].astype(float), BRAIN20_BLBR[1][1:17].astype(float))[0]])

x=len(PBMC_BLBR) 	
#x=10
spearman_correlations=[]
for CpG in range(1,x):
    PBMC = PBMC_BLBR[CpG][1:17].astype(float)
    BRAIN7 = BRAIN7_BLBR[CpG][1:17].astype(float)
    BRAIN10 = BRAIN10_BLBR[CpG][1:17].astype(float)
    BRAIN20 = BRAIN20_BLBR[CpG][1:17].astype(float)
    
    brain7=spearmanr(PBMC, BRAIN7)[0]
    brain10=spearmanr(PBMC, BRAIN10)[0]
    brain20=spearmanr(PBMC, BRAIN20)[0]
    spearman_correlations.append([brain7,brain10, brain20])

print spearman_correlations[0:10]

    ##Save
df = pd.DataFrame(spearman_correlations)
print df[0:10]
#df.to_csv("correlation_sample.csv", sep='\t')
df.to_csv("correlation_unadjustjan22.csv", sep='\t')


# going to use old correlations
#spearman_correlations_random_mn=[]
#for CpG in range(1,x):
    #brain7=[]
    #brain10=[]
    #brain20=[]
    #for permute in range(1,1000):					 					##### Should be 1000
    	#print permute
        #np.random.shuffle(PBMC_order)
        #PBMC = PBMC_BLBR[CpG][PBMC_order]
        
        #BRAIN7 = BRAIN7_BLBR[CpG][BRAIN7_order]
        #BRAIN10 = BRAIN10_BLBR[CpG][BRAIN10_order]
        #BRAIN20 = BRAIN20_BLBR[CpG][BRAIN20_order]
        
        #brain7.append(spearmanr(PBMC.astype(float), BRAIN7.astype(float))[0])
        #brain10.append(spearmanr(PBMC.astype(float), BRAIN10.astype(float))[0])
        #brain20.append(spearmanr(PBMC.astype(float), BRAIN20.astype(float))[0])
        
    #spearman_correlations_random_mn.append([np.mean(brain7),np.std(brain7), 
                                            #np.mean(brain10),np.std(brain10),
                                            #np.mean(brain20), np.std(brain20)])

#dfrnd = pd.DataFrame(spearman_correlations_random_mn)
#dfrnd.to_csv("spearman_correlations_random_mn.csv", sep='\t')
#dfrnd.to_csv("spearman_correlations_random_mn_sample.csv", sep='\t')