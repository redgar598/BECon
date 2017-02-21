## BECon Script work flow

####All data objects are located on cake at `/big_data/redgar/BECon`

1. Normalization_BECon.Rmd
2. Post_Normalization_QC.Rmd
3. After_Normalization.Rmd
    * Dependent on: 
    * BLBR_Correlation_Permutation_jan22.py
    * BLBR_Correlation_Permutation_jan22_permutations.py
4. Differential_methylation_blood_brain.Rmd

7. Correlated_CpGs.Rmd

8. RefSeq_updated_june2015.R 
9. Informative_Gene_Interpretation.R
10. Informative_Gene_Expression.R
11. GSE59685_Validation_dataset.R
    * Dependent on: 
    * GSE59685_Correlation_Oct19.py
    * GSE59685_Correlation_snpsex.py
12. Informative_CpG_genomic_enrichment.Rmd
      * Dependent on: Gene_enrichment_fold_change_permutation_pvalue_function.R
13. mQTL_enrichment.R
