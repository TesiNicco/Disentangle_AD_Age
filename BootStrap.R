# LIBRARIES
    message("** Checking and Loading packages Packages")
    list.of.packages <- c("data.table", "stringr", "parallel", "lme4", "rlist")
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)
    suppressPackageStartupMessages({
    library(data.table)
    library(stringr)
    library(parallel)
    library(lme4)
    library(rlist)
    })
    args = commandArgs(trailingOnly=TRUE)

# FUNCTIONS
    ## function to grep ids of the variants of interest
    MatchIds <- function(i, inp, path_pvar_files){
        ## grep variants in the chromosome of interest
        sb <- inp[which(inp$chr == i), ]

        ## temporary output is a list
        snp_list = list()

        ## grep the position in the pvar file
        fpath <- paste(path_pvar_files, "chr", i, ".dose.pvar", sep="")
        for (snp in 1:nrow(sb)){
            cmd = paste("grep -w", sb$pos[snp], fpath, sep=" ")
            res <- system(cmd, intern=TRUE)
            res <- as.data.frame(str_split_fixed(res, "\t", 7))
            snp_list[[snp]] <- res
        }

        ## merge results
        all <- rbindlist(snp_list)
        colnames(all) <- c("chr", "pos", "id", "ref", "alt", "qc", "info")

        return(all)
    }

    ## function to check quality of each variant: maf>1% and R2>0.80
    CheckQuality <- function(all_chroms){
        ## grep info field and split
        all_info <- all_chroms$info
        all_info <- as.data.frame(str_split_fixed(all_info, ";", 5))

        ## add frequency and R2 as separate columns
        all_chroms$maf <- NA
        all_chroms$r2 <- NA
        all_chroms$type <- NA

        ## check all variants in a loop
        for (snp in 1:nrow(all_info)){
            if (all_info$V5[snp] == ""){        # if true, than it is an imputed variant
                r2_info <- as.data.frame(str_split_fixed(all_info$V4[snp], "=", 2))         # take r2 value
                maf_field <- as.data.frame(str_split_fixed(all_info$V1[snp], "=", 2))       # take maf value
                if (as.numeric(as.character(maf_field$V2)) >= 0.5) { all_chroms$maf[snp] <- 1 - as.numeric(as.character(maf_field$V2)) } else {all_chroms$maf[snp] <- as.numeric(as.character(maf_field$V2))}
                all_chroms$r2[snp] <- as.numeric(as.character(r2_info$V2))
                all_chroms$type[snp] <- "imputed"
            } else {
                maf_field <- as.data.frame(str_split_fixed(all_info$V2[snp], "=", 2))       # take maf value
                if (as.numeric(as.character(maf_field$V2)) >= 0.5) { all_chroms$maf[snp] <- 1 - as.numeric(as.character(maf_field$V2)) } else {all_chroms$maf[snp] <- as.numeric(as.character(maf_field$V2))}
                all_chroms$r2[snp] <- 1
                all_chroms$type[snp] <- "genotyped"
            }
        }

        ## filters: maf>1% and R2>0.8
        #all_chroms <- all_chroms[which(all_chroms$maf >= 0.01 & all_chroms$r2 >= 0.80), ]
        all_chroms <- all_chroms[which(all_chroms$r2 >= 0.80), ]

        return(all_chroms)
    }

    ## function to generate dosages files for these variants
    GetDosages <- function(i, clean_snps, path_pvar_files){
        ## get snps for each chromosome
        sb <- clean_snps[which(clean_snps$chr == i), ]

        ## write temporary output for the variants to keep
        write.table(sb$id, paste("chr", i, "tmp", sep="_"), quote=F, row.names=F, col.names=F)

        ## plink2 command
        fpath <- paste(path_pvar_files, "chr", i, ".dose", sep="")
        cmd = paste("plink2 --pfile ", fpath, " --extract chr_", i, "_tmp --export A --out chr_", i, "_dos", sep="")
        system(cmd, ignore.stdout=T)

        ## read file
        d <- fread(paste("chr_", i, "_dos.raw", sep=""), h=T)

        ## final clean up on file and temporary files
        d$FID <- NULL
        d$MAT <- NULL
        d$PAT <- NULL
        d$SEX <- NULL
        d$PHENOTYPE <- NULL
        cmd <- paste("rm chr_", i, "_*", sep="")
        system(cmd)

        return(d)
    }

    ## function to create bootstrap datasets
    createBootstrapDatasets <- function(i, dos){
        boot <- dos[sample(x = 1:nrow(dos), size = nrow(dos), replace=TRUE, prob=NULL), ]
        return(boot)
    }

    ## function to perform association of the bootstrapped datasets
    BootAssoc <- function(i, boot_cases, boot_controls, covar){
        ## put together cases and controls
        all <- rbind(boot_cases[[i]], boot_controls[[i]])
        all$pheno <- c(rep(1, nrow(boot_cases[[i]])), rep(0, nrow(boot_controls[[i]])))

        ## add covariates
        all <- merge(all, covar, by="IID")
        colnames(all) <- c("IID", "dos", "pheno", "pc1", "pc2", "pc3", "pc4", "pc5")
        
        ## logistic model
        m <- glm(pheno ~ dos + pc1 + pc2 + pc3 + pc4 + pc5, data=all, family="binomial")

        ## extract effect-size
        beta <- as.numeric(m$coefficients[2])
        
        return(beta)
    }

    ## function to run bootstrap
    Bootstrap <- function(i, dosages, pheno, covar, n_boot){
        print(i)
        ## get variant of interest
        sb <- dosages[, ..i]
        sb$IID <- dosages$IID

        ## define cases and controls
        cases <- pheno[which(pheno$PHENO1 == 2), ]
        controls <- pheno[which(pheno$PHENO1 == 1), ]

        ## create bootstrap sets
        boot_cases <- lapply(1:n_boot, createBootstrapDatasets, dos=sb[which(sb$IID %in% cases$IID), ])
        boot_controls <- lapply(1:n_boot, createBootstrapDatasets, dos=sb[which(sb$IID %in% controls$IID), ])

        ## perform association
        boot_assoc <- lapply(1:n_boot, BootAssoc, boot_cases=boot_cases, boot_controls=boot_controls, covar=covar)

        ## condense representation
        boot_assoc <- unlist(boot_assoc)

        return(boot_assoc)
    }

## MAIN
    ## define input SNP file
    path_pvar_files = args[1]       # this is the path to the folder containing the .pvar files -- For example, /path/to/folder/ -- Make sure the name of the file is chr[1-22].dose.pvar
    ## define AD snps file --> description in github readme (/home/nfs/ntesi/bulkALL/niccolo/lasa/nicco/new_batch2018/ADsnps/20200407_rotation_paper/INPUTS/ADsnps_GRACE.txt)
    fname <- args[2]
    ## define phenotype file --> description in github readme (/home/nfs/ntesi/bulkALL/niccolo/lasa/nicco/new_batch2018/ADsnps/20200407_rotation_paper/INPUTS/20200317_phenotypes_subsetAge.txt)
    pheno <- args[2]
    ## define covariate file --> description in github readme (/home/nfs/ntesi/bulkALL/niccolo/lasa/nicco/new_batch2018/ADsnps/20200407_rotation_paper/INPUTS/20200123_covariates_kept_Age_gwas.txt)
    covar <- args[3]
    ## define number of cores to use for parallelization
    n_cores <- args[4]
    ## define input dataset to get random snps as well --> description in github readme (/home/nfs/ntesi/bulkALL/niccolo/lasa/nicco/new_batch2018/collaboration_SPAIN/20190926_newPheno_replication/lastVersion_plots/metaGWAS_repli16dbs_20190930.1tbl.rs.gz)
    metaAnal_path <- args[5]
    ## define main directory for output files --> For example, /path/to/output/folder/ (/home/nfs/ntesi/bulkALL/niccolo/lasa/nicco/new_batch2018/ADsnps/20200407_rotation_paper/)
    MAIN = args[6]

    ## 1. read input AD snps
    message("** Read input AD SNPs")
    inp <- fread(fname, h=F)
    colnames(inp) <- c("chr", "pos", "a1", "a2", "beta", "se", "p", "rsid")
    
    ## 2. read input meta-analysis to select also random snps
    message("** Read random SNPs")
    metaAnal <- fread(metaAnal_path, h=T)
    
    ## 2.1 randomly take 1000 SNPs
    metaAnal = metaAnal[!is.na(metaAnal$RS),]
    rand_inp = metaAnal[sample(x = seq(1, nrow(metaAnal)), size = 2000, replace = F), c("CHR", "BP", "A1", "A2", "beta", "SE", "P", "RS")]
    colnames(rand_inp) <- c("chr", "pos", "a1", "a2", "beta", "se", "p", "rsid")
    rand_inp = rand_inp[with(rand_inp, order(chr, pos)), ]
    
    ## 3. extract list of chromosomes to look at -- for AD snps and random snps
    chr_list <- unique(inp$chr)
    chr_list_rand <- unique(rand_inp$chr)
    
    ## 4. run function to match ids in parallel for each chromosome -- for AD snps and random snps
    message("** Finding SNPs of interest")
    all_chroms <- rbindlist(mclapply(chr_list, MatchIds, inp=inp, path_pvar_files = path_pvar_files, mc.cores=n_cores))
    all_chroms_rand <- rbindlist(mclapply(chr_list_rand, MatchIds, inp=rand_inp, path_pvar_files = path_pvar_files, mc.cores=n_cores))
    
    ## 5. run function to check quality -- mach_r2>0.8 -- for AD snps and random snps
    message("** Cleaning set of SNPs")
    clean_snps <- CheckQuality(all_chroms)
    clean_snps_rand <- CheckQuality(all_chroms_rand)
    
    ## 5.1 make sure that there are no duplicates
    dups_rand = clean_snps_rand[duplicated(clean_snps_rand$id),]
    clean_snps_rand = clean_snps_rand[which(!(clean_snps_rand$id %in% dups_rand$id)),]
    clean_snps_rand = clean_snps_rand[which(clean_snps_rand$maf >= 0.01),]
    
    ## 6. save variant information -- for AD snps and random snps
    system(paste("mkdir ", MAIN, "RESULTS_with_random", sep=""))
    write.table(clean_snps, paste0(MAIN, "RESULTS_with_random/snps_information.txt"), quote=F, row.names=F, sep="\t")
    write.table(clean_snps_rand, paste0(MAIN, "RESULTS_with_random/snps_information_random.txt"), quote=F, row.names=F, sep="\t")

    ## 6.1 besides saving the snps, also save the actual associations of the random snps
    clean_snps_rand$SNP = paste(clean_snps_rand$chr, clean_snps_rand$pos, sep =":")
    assoc_rand = merge(clean_snps_rand, metaAnal, by = "SNP")

    ## 6.2 there may be duplicates because of multiallelic SNPs in the random set -- check this
    dups_rand = assoc_rand[duplicated(assoc_rand$SNP),]
    all_dups_rand = assoc_rand[which(assoc_rand$SNP %in% dups_rand$SNP),]
    all_dups_rand = unique(all_dups_rand)
    assoc_rand = assoc_rand[which(!(assoc_rand$SNP %in% all_dups_rand$SNP)),]
    all_dups_rand$keep = NA
    for (i in 1:nrow(all_dups_rand)){
        alleles_gwas = toupper(c(all_dups_rand$A1[i], all_dups_rand$A2[i]))
        alleles_gwas = alleles_gwas[order(alleles_gwas)]
        alleles_gwas = paste(alleles_gwas, collapse = ":")
        alleles_ref = toupper(c(all_dups_rand$ref[i], all_dups_rand$alt[i]))
        alleles_ref = alleles_ref[order(alleles_ref)]
        alleles_ref = paste(alleles_ref, collapse = ":")
        if (alleles_gwas == alleles_ref){
            all_dups_rand$keep[i] = "YES"
        } else {
            all_dups_rand$keep[i] = "NO"
        }
    }
    kp = all_dups_rand[which(all_dups_rand$keep == "YES"),]
    kp$keep = NULL
    assoc_rand = rbind(assoc_rand, kp)

    ## 6.3 finally output the association of the random snps
    write.table(assoc_rand, paste0(MAIN, "RESULTS_with_random/snps_association_random.txt"), quote=F, row.names =F, sep = "\t")

    ## 7. run function to get dosages for these variants -- for AD snps and random snps
    message("** Getting dosages")
    chr_list <- unique(clean_snps$chr)
    all_dosages <- mclapply(chr_list, GetDosages, clean_snps=clean_snps, path_pvar_files=path_pvar_files, mc.cores=n_cores)
    chr_list_rand <- unique(clean_snps_rand$chr)
    all_dosages_rand <- mclapply(chr_list_rand, GetDosages, clean_snps=clean_snps_rand, path_pvar_files=path_pvar_files, mc.cores=n_cores)

    ## 7.1 exclude duplicated snps in the random set
    for (chrom in 1:length(all_dosages_rand)){
        tmp = all_dosages_rand[[chrom]]
        dups = colnames(tmp)[duplicated(colnames(tmp))]
        for (i in dups){
            indexes = grep(i, colnames(tmp))
            tmp[, indexes] = NULL
        }
        all_dosages_rand[[chrom]] = tmp
    }

    ## 8. merge results from the different chromosomes -- for AD snps and random snps
    message("** Merging dosages")
    dosages <- Reduce(function(x, y) merge(x, y, by="IID", all.x=T, all.y=T), all_dosages)
    dosages_rand <- Reduce(function(x, y) merge(x, y, by="IID", all.x=T, all.y=T), all_dosages_rand)

    ## 9. read covariates and phenotypes
    pheno <- fread(pheno, h=T)
    covar <- fread(covar, h=T)

    ## 10. run function to perform bootstrap -- for AD snps and random snps
    message("** Bootstrapping associations")
    boot_allSNPs <- mclapply(2:ncol(dosages), Bootstrap, dosages=dosages, pheno=pheno, covar=covar, n_boot=1000, mc.cores=n_cores)
    boot_allSNPs_rand <- mclapply(2:ncol(dosages_rand), Bootstrap, dosages=dosages_rand, pheno=pheno, covar=covar, n_boot=1000, mc.cores=n_cores)

    ## 11. merge results together -- for AD snps and random snps
    message("** Merging results")
    boot_clean <- as.data.frame(list.cbind(boot_allSNPs))
    names(boot_clean) <- names(dosages)[2:ncol(dosages)]
    boot_clean_rand <- as.data.frame(list.cbind(boot_allSNPs_rand))
    names(boot_clean_rand) <- names(dosages_rand)[2:ncol(dosages_rand)]

    ## 12. save image and bootstrap matrix
    message("** Saving outputs")
    save(boot_clean, file=paste0(MAIN, "RESULTS_with_random/Bootstrap_survival_10k_ADsnps.RData"))
    write.table(boot_clean, paste0(MAIN, "RESULTS_with_random/Bootstrap_survival_10k_ADsnps.txt"), quote=F, row.names=F, sep="\t")
    save(boot_clean_rand, file=paste0(MAIN, "RESULTS_with_random/Bootstrap_survival_10k_ADsnps_random.RData"))
    write.table(boot_clean_rand, paste0(MAIN, "RESULTS_with_random/Bootstrap_survival_10k_ADsnps_random.txt"), quote=F, row.names=F, sep="\t")