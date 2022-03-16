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
- Rscript BootStrap.R [arg1] [arg2] ...
- Rscript RotateMe.R [arg1] [arg2] ...
Please see below for an extensive description of the arguments.

# BootStrap.R
This scripts must be run with multiple arguments. Example files of some arguments are provided here. If not provided, we describe the datasets that should be used.
Please note that the order of the arguments does matter.
1. argument 1 --> this is the path to the folder containing the .pvar files. .pvar files are generated with PLINK, and they contain SNP information (chromosome, position, alleles). For additional information about .pvar files check this link --> https://www.cog-genomics.org/plink/2.0/formats#pvar. An example of this argument is /path/to/folder/ . Make sure the name of the genotype files is chr[1-22].dose.pvar.
2. argument 2 --> tab-delinited file specifying the SNPs (in our case, AD-associated SNPs) of interest. We provide an example file (example_datasets/example_file_arg2.txt). In general, a tab-delimited file without header containing chromosome_number, position, effect_allele, other_allele, effect_size, standard_error, p_value and variant_identifier is required. Please note that the order of the columns does matter.
3. argument 3 --> tab-delimited file specifying the phenotypes of your genotype data. See example_datasets/example_file_arg3.txt as an example dataset. This should be compatible with PLINK2. A two-column file with header (IID, PHENO1) should be provided. The column IID contains the name of the samples as present in your genotype data. The column PHENO1 should contains either 1 or 2, where 1 are controls and 2 are cases. Since we looked into longevity, our cases were long-lived individuals. Additional information --> https://www.cog-genomics.org/plink/2.0/assoc
4. argument 4 --> tab-delimited file specifying the covariates of your genotype data. See example_datasets/example_file_arg4.txt as an example dataset. This should be compatible with PLINK2. A multi-column file with header should be provided. The first column should report the sample identifiers as present in the genotype data. Other columns represent the covariate name to be used in the model. As of now, only 5 covariates are used, namely PC1, PC2, PC3, PC4, PC5. If you would like to modify or use other covariates, you should also adapt the code in the script in the function BootAssoc().
5. argument 5 --> tab-delimited file reporting SNPs from which a random sample of N=1000 SNPs will be drawn. See example_datasets/example_file_arg5.txt as an example dataset. For consistency, the source (that is, the study) used for this file and the [argument 2] should be the same. Make sure that columns chromosome (CHR), position (BP), effect-allele (A1), other-allele (A2), effect-size (beta), standard error (SE), p-value (P) and snp-identifier (RS) are present. The order of the columns does not matter, but please make sure the column names are spelled correctly --> "CHR", "BP", "A1", "A2", "beta", "SE", "P", "RS".
5. argument 6 --> this is the path to the output folder that will contain results. A valid path could be /path/to/output/files/
The script produces several outputs:
1. snps_information.txt --> information about the considered SNPs (as defined in argument 2)
2. snps_information_random.txt --> information about the considered random SNPs (as defined in argument 5)
3. snps_association_random.txt --> information about the association of the random SNPs (as defined in argument 5)
4. Bootstrap_survival_10k_ADsnps.txt/RData --> main sampling output of the script in tab-delimited and R format
5. Bootstrap_survival_10k_ADsnps_random.txt/RData --> main sampling output of random SNPs in tab-delimited and R format

# RotateMe.R
This is the main downstream analysis script. The arguments for this file are:
1. argument 1 --> output directory. Please create this output directory. Please create an additional directory named 'INPUTS' inside this directory. The 'INPUTS' directory will contain all input data for this script.
2. argument 2 --> input SNPs data. Please place this file in 'INPUTS' directory. An example file is shown in INPUTS/inputSNPs.txt

Additional files that should be copied in INPUTS folder:
1. snps_association_random.txt --> this is produced by 'BootStrap.R' [not provided]
2. Bootstrap_survival_10k_ADsnps.RData --> this is produced by 'BootStrap.R' [not provided]
3. Bootstrap_survival_10k_ADsnps_random.RData --> this is produced by 'BootStrap.R' [not provided]
4. inputSNPs.txt file --> see example INPUTS/inputSNPs.txt
5. myAssoc_inputSNPs.txt --> see example INPUTS/myAssoc_inputSNPs.txt [This file contains associations (logistic regression), of the input SNPs. This does not include a bootstrapping procedure.]
6. RNAexpr_fromBrain.tab --> RNA expression of cell-types in brain. See example file INPUTS/RNAexpr_fromBrain.tab
7. semSim_out_LIN.txt --> This file can be obtained running snpxplorer functional annotation analysis at https://snpxplorer.net [provided at INPUTS/semSim_out_LIN.txt]
8. Results_enrichment_sampling_clusPro.RData --> This file can be obtained running snpxplorer functional annotation analysis at https://snpxplorer.net [provided at INPUTS/Results_enrichment_sampling_clusPro.RData]
9. snp_annotation.txt --> This file can be obtained running snpxplorer functional annotation analysis at https://snpxplorer.net [provided at INPUTS/snp_annotation.txt]
10. snp_annotation_geneList.txt --> This file can be obtained running snpxplorer functional annotation analysis at https://snpxplorer.net [provided at INPUTS/snp_annotation_geneList.txt]
11. analysis_GWAS_context.RData --> This file contains data used to draw the forest plot (figure 3) [provided at INPUTS/analysis_GWAS_context.RData]
