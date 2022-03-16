# Disentangle_AD_Age
Script for the analysis of the imbalance effect direction between Alzheimer's disease and cognitively healthy aging.

These scripts are useful to reproduce (provided the same input data) or replicate with a different dataset the analyses that were performed in the following manuscript:
The Effect of Alzheimer’s Disease-Associated Genetic Variants on Longevity. Tesi et al., 2021. Frontiers in Genetics. The manuscript is available at the following link:
https://www.frontiersin.org/articles/10.3389/fgene.2021.748781/full

# Abstract
Human longevity is influenced by the genetic risk of age-related diseases. As Alzheimer’s disease (AD) represents a common condition at old age, an interplay between genetic factors affecting AD and longevity is expected. We explored this interplay by studying the prevalence of AD-associated single-nucleotide-polymorphisms (SNPs) in cognitively healthy centenarians, and replicated findings in a parental-longevity GWAS. We found that 28/38 SNPs that increased AD-risk also associated with lower odds of longevity. For each SNP, we express the imbalance between AD- and longevity-risk as an effect-size distribution. Based on these distributions, we grouped the SNPs in three groups: 17 SNPs increased AD-risk more than they decreased longevity-risk, and were enriched for β-amyloid metabolism and immune signaling; 11 variants reported a larger longevity-effect compared to their AD-effect, were enriched for endocytosis/immune-signaling, and were previously associated with other age-related diseases. Unexpectedly, 10 variants associated with an increased risk of AD and higher odds of longevity. Altogether, we show that different AD-associated SNPs have different effects on longevity, including SNPs that may confer general neuro-protective functions against AD and other age-related diseases.

# Description of the scripts
There are two main scripts. Both scripts are written in R. The two scripts relate to:
1. The first script (BootStrap.R) is used to estimate SNP effect-sizes on longevity. To calculate effect-sizes on longevity, we used logistic regression models in a case-control setting while correcting for population stratification using PCA components (PC 1-5, arbitrarily chosen). To calculate the confidence interval, we repeated this procedure for bootstraps (B = 10,000) of the data. Additional details are available in the original publication, Methods section.
2. The second script performs all downstream analyses as well as the plots that are shown in the original publication.

# How to run the script
The scripts are written in R, therefore, as long as you have R and Rscript correctly installed in your machine, it should be fairly easy to run the scripts. You should be able to run the script with:
Rscript BootStrap.R [arg1] [arg2] ...
Rscript RotateMe.R [arg1] [arg2] ...
Please see below for an extensie description of the arguments.

# BootStrap.R
This scripts must be run with multiple arguments. The arguments are not provided here, but we offer a description of the datasets that should be used.
1. argument 1 --> 