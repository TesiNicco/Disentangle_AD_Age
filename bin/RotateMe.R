##########################################
#   THIS IS THE FINAL ROTATE ME SCRIPT   #  
#                                        #
#  #  #  #  ###  ###  ###   #    ###  #  #
#  # ##  #  #    #    # #   #    # #     #
#  #  #  #  ###  ###  ###   ###  ###     #
#                                        #
#  ### ### ### ### ### ###   # # ###     #
#  ##  # #  #  ###  #  ##    ### ##      #
#  # # ###  #  # #  #  ###   # # ###     #
#                                        #
##########################################

##########################################
## CLEAN WORKSPACE
rm(list=ls())
##########################################

##########################################
## LIBRARIES
suppressPackageStartupMessages({
  library(pheatmap)
  library(dendextend)
  library(wordcloud2)
  library(RColorBrewer)
  library(tibble)
  library(dplyr)
  library(stringr)
  library(data.table)
  library(biomaRt)
  library(ggplot2)
  library(tidytext)
  library(caret)
  library(GOSim)
  library(WGCNA)
  library(GOSemSim)
  # library("RPMG")
  # library("stringr")
  # library("plotrix")
  # library("data.table")
  # library("grid")
  # library("ggplot2")
  # library("RColorBrewer")
  # library("wesanderson")
  # library("xlsx")
  # library("pheatmap")
  # library("tibble")
  # library("dplyr")
  # library(data.table)
  # library(stringr)
  # library(parallel)
  # library(lme4)
  # library(rlist)
  # library(pvclust)
  # library(ggridges)
  # library(RPMG)
  # library(plotrix)
  # library(GOSemSim)
  # library(org.Hs.eg.db)
  # library(GO.db)
  # library(pheatmap)
  # library(tidytext)
  # library(biomaRt)
  # library(dplyr)
  # library(tidyverse)  # data manipulation
  # library(cluster)    # clustering algorithms
  # library(factoextra) # clustering algorithms & visualization
  # library(grid)
  # library(reshape2)
  # library(ggplot2)
  # library(RColorBrewer)
  # library(wesanderson)
  # library(xlsx)
  # library(ngram)
  # library(viridis)
})
######################################

######################################
## FUNCTIONSline.x
# Function to sample distribution of effect sizes from beta and standard error from literature
function.distrLit <- function(lit, boot.n){
  # Define output
  lit.dist <- matrix(data = NA, nrow = nrow(lit), ncol = boot.n)
  rownames(lit.dist) <- lit$rsid
  
  # Main loop
  for (snp in 1:nrow(lit)){
    # Make normal distribution
    dist <- rnorm(n = boot.n, mean = lit$beta[snp], sd = lit$se[snp])
    
    # Put distribution in matrix
    lit.dist[snp, ] <- dist
  }
  # Assign row names
  rownames(lit.dist) <- lit$rsid
  
  return(lit.dist)
}

## Function to order two data frames in the same way
matchBy <- function(df1, by.1, by.2){ 
  df1 <- df1[match(by.1, by.2),]
  return(df1)
}

## Function to manage names of the snps and alleles so that they are the same
makeItTheSame <- function(boot_clean, lit, lit_dist, p_filter){
  ## Get the names of both datasets
  names_boot <- as.data.frame(str_split_fixed(rownames(boot_clean), "_", 2))
  
  ## Main loop
  for (i in 1:nrow(boot_clean)){
    ## Get variant name and allele
    snp_name <- as.character(names_boot$V1[i])
    snp_allele <- as.character(names_boot$V2[i])
    
    ## Check if the name is in literature ids
    lit_info <- lit[which(lit$rsid == snp_name), ]
    
    ## Check if there was a match -- cause can be different id for my variants
    if (nrow(lit_info) >0){         # if there was a match, check alleles and in case flip all my betas
      if (snp_allele != as.character(toupper(lit_info$a1))) { boot_clean[i, ] <- -boot_clean[i, ] }
      if (lit_info$beta <0){ boot_clean[i, ] <- -boot_clean[i, ]; lit_dist[which(rownames(lit_dist) == snp_name), ] <- -lit_dist[which(rownames(lit_dist) == snp_name), ] }
      rownames(boot_clean)[i] <- lit_info$rsid
    } else {                # if there was no match, need to grep by chr:pos
      lit_info <- lit[which(paste(lit$chr, lit$pos, sep=":") == snp_name),]
      if (snp_allele != as.character(toupper(lit_info$a1))) { boot_clean[i, ] <- -boot_clean[i, ] }
      if (lit_info$beta <0){ boot_clean[i, ] <- -boot_clean[i, ]; lit_dist[which(paste(lit$chr, lit$pos, sep=":") == snp_name), ] <- -lit_dist[which(paste(lit$chr, lit$pos, sep=":") == snp_name), ] }
      rownames(boot_clean)[i] <- lit_info$rsid
    }
  }
  
  ## Add a filter on significance for the literature -- keep replicated with p<5e-4
  if (p_filter == TRUE){ 
    toKeep <- lit$rsid[which(lit$p <= 5e-4 & lit$rsid != "rs75932628")] 
    ## Next, select same variants from literature and bootstrap
    lit_dist <- lit_dist[which(rownames(lit_dist) %in% toKeep), ]
    boot_clean <- boot_clean[which(rownames(boot_clean) %in% toKeep), ]
    lit_dist <- lit_dist[which(rownames(lit_dist) %in% rownames(boot_clean)),]
  }
  
  ## Match order of the rows
  lit_dist <- matchBy(lit_dist, rownames(boot_clean), rownames(lit_dist))
  
  return(list(lit_dist, boot_clean))    
}

## Function to transform bootstrap points
Transform <- function(boot.x, boot.y){  
  ## Transform to matrix
  boot.x <- as.matrix(boot.x)
  boot.y <- as.matrix(boot.y)
  
  ## Compute normalized angle
  angle <- atan2(boot.y, boot.x)/(-pi/2)+1

  ## Take care of the y<0
  angle[which(angle >1)] = angle[which(angle >1)] - 2
  
  return(angle)
}
# the below part is just to make sure everything is going ok
#tmp_x = rnorm(n = 100, mean = 0.5, sd = 0.1)
#tmp_y = rnorm(n = 100, mean = -0.01, sd = 0.1)
snp = sample(x = rownames(boot_dist_rand), size = 1, replace = F)
snp = "rs17736861"
tmp_x = boot_dist_rand[which(rownames(boot_dist_rand) == snp),]
mean(tmp_x)
tmp_y = lit_rand_dist[which(rownames(lit_rand_dist) == snp),]
mean(tmp_y)
angle = Transform(tmp_x, tmp_y)
par(mfrow=c(1,1), mar = c(4,5,4,5))
plot(tmp_x, tmp_y, xlim=c(-0.5, 0.5), ylim=c(-0.25, 0.25), pch = 16, col = viridis(n = 3, alpha = 0.5)[2],
     xlab = "Effect on Longevity", ylab = "Effect on AD", main = snp)
abline(h=0)
abline(v=0)
abline(a = 0, b = -1, lty = 2)
abline(a = 0, b = 1, lty = 2)
lines(density(angle, bw = "sj"), lwd = 2, col = "red")
median(angle)
mean(angle)
#### MEAN oR MEDIAN???
####

## Semantic similarity of GO terms and clustering of GO terms
function.GlobalProps <- function(angle){
  # Read ouput file with similarity distances -- Lin
  data <- fread("RESULTS/semSim_out_LIN.txt", h=T, sep=",")
  
  # Plot heatmap and save it
  png("RESULTS/pheatmap_lin_distance.png", height=7, width=7, res=300, units="in")
  pheatmap(data, show_rownames = F, show_colnames = F)
  dev.off()
  
  hr <- hclust(as.dist(1-data), method="ward.D2", members = NULL)
  minModuleSize = 15
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = hr, distM = (1-data),
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  tb <- as.data.frame(table(dynamicMods))
  n_clust <- max(as.numeric(as.character(tb$dynamicMods)))
  
  # Hierarchical clustering on the distance matrix
  hr <- as.dendrogram(hclust(as.dist(1-data), method="ward.D2", members = NULL))
  
  # Plot dendrogram and rectangles of the clusters
  png("RESULTS/dendrogram_GOterms.png", height=7, width=17, res=300, units="in")
  d1=color_branches(hr, k=n_clust, groupLabels = T)
  plot(d1)
  dev.off()

  # Merge with terms description
  data_clus <- as.data.frame(cutree(hr, k=n_clust))
  data_clus$term <- rownames(data_clus)
  colnames(data_clus) <- c("cluster", "term")
  desc <- fread("RESULTS/enrichent_results_sampling.txt", h=F)
  colnames(desc) <- c("term_id", "term_name", "avgP", "avgQ", "genes")
  clusters <- merge(data_clus, desc, by.x="term", by.y="term_id")
  
  # look at the labels
  for (i in 1:n_clust){ clusters$term_name[which(clusters$cluster == i)] }

  # return
  return(list(clusters, data, n_clust))
}

## Function to count the number of occurrences of each word given a vector of sentences
CountFrequency_words <- function(functional_clusters, n_clust){
  # Read the description of all GO:BP -- this will be the background to remove redundant words
  all_go_bp <- fread("INPUTS/go_terms_BP_all.txt", h=F, stringsAsFactors = F, sep="\t")
  
  # Put all different words datasets in a list
  if (n_clust == 6){
    all_dset <- list(all_go_bp$V1, functional_clusters$term_name[which(functional_clusters$cluster == 1)], 
                   functional_clusters$term_name[which(functional_clusters$cluster == 2)],
                   functional_clusters$term_name[which(functional_clusters$cluster == 3)],
                   functional_clusters$term_name[which(functional_clusters$cluster == 4)],
                   functional_clusters$term_name[which(functional_clusters$cluster == 5)],
                   functional_clusters$term_name[which(functional_clusters$cluster == 6)])
  } else if (n_clust == 4){
    all_dset <- list(all_go_bp$V1, functional_clusters$term_name[which(functional_clusters$cluster == 1)], 
                     functional_clusters$term_name[which(functional_clusters$cluster == 2)],
                     functional_clusters$term_name[which(functional_clusters$cluster == 3)],
                     functional_clusters$term_name[which(functional_clusters$cluster == 4)])
  } else if (n_clust == 5){
    all_dset <- list(all_go_bp$V1, functional_clusters$term_name[which(functional_clusters$cluster == 1)], 
                     functional_clusters$term_name[which(functional_clusters$cluster == 2)],
                     functional_clusters$term_name[which(functional_clusters$cluster == 3)],
                     functional_clusters$term_name[which(functional_clusters$cluster == 4)],
                     functional_clusters$term_name[which(functional_clusters$cluster == 5)])
  }
  
  # Define output
  word_frequency_dset <- list()
  
  # Main loop to calculate word frequency
  for (i in 1:length(all_dset)){
    # Convert to dataframe and make sure there are characters only
    text <- data.frame(text=all_dset[[i]])
    text$text <- as.character(text$text)
    
    # Covert dataframe to tibble
    text_df <- tibble(line = 1:nrow(text), text = text)
    text_df <- mutate(text_df, text = text$text)
    
    # Count word occurrence
    x <- text_df %>% unnest_tokens(word, text) %>%    # split words
      anti_join(stop_words) %>%    # take out "a", "an", "the", etc.
      count(word, sort = TRUE)    # count occurrences
    
    # Convert back to dataframe
    x <- as.data.frame(x)
    colnames(x) <- c("word", "freq")
    
    # Assign to output
    word_frequency_dset[[i]] <- x
  }
 
  # Extract background and calculate % and take top 5%
  bkg <- word_frequency_dset[[1]]
  bkg$perc <- bkg$freq/nrow(all_go_bp)
  top_5_pc <- bkg[which(bkg$perc >= 0.025),]
  
  # Now exclude these frequent words from the clusters names
  for (i in 2:length(word_frequency_dset)){
    sb <- word_frequency_dset[[i]]
    sb <- sb[which(!(sb$word %in% top_5_pc$word)),]
    word_frequency_dset[[i]] <- sb
  }
  
  # Try some plot -- seems like it is not possible to plot multiple plots in one --> make 4 plots and save them from "Export"
  for (i in 2:length(word_frequency_dset)){ wordcloud2(data = word_frequency_dset[[i]]) }
  
  return(word_frequency_dset)
}

## Function to plot densities only
plotDensities <- function(line.x, ordered){
  # color palette
  color.palette <- brewer.pal(9, "Set1")
  
  # Global graphical parameters
  mat <- matrix(c(1, 2), 1)
  layout(mat, widths = c(1, 4), heights = c(1,1))
  par(mar=c(1, 1, 8, 3))
  # ymax for the plot
  ymax <- nrow(line.x) + 1
  # define square height
  h <- 1
  # define counter for the plot
  c <- 1
  
  # PLOT 1 -- clustering
  d <- dist(line.x, method = "manhattan")
  fit <- hclust(d, method = "complete") 
  nodePar <- list(lab.cex = 1, pch = c(NA, 19), cex = 1.20, col = "coral")
  edgePar = list(col = 1, lwd = 2)
  dendr <- as.dendrogram(fit)
  names.order <- fit$labels[order.dendrogram(dendr)]
  plot(dendr, yaxt="none", horiz = T, nodePar = nodePar, edgePar = edgePar, leaflab="none",
       xpd=T, ylim = c(0.5, ymax+0.5))
  text(x = 2000, y = (ymax - c + h*1.5-0.5), labels = 'A', font=2, xpd=T, cex = 2)
  text(x = 2000, y = (ymax + 1.75), labels = 'Clustering', font=1, xpd=T, cex = 1.50)
  rev.names <- rev(names.order)
  #rev.names <- ordered$snp
  
  # PLOT 2 -- big plot
  # set graphical parameters for second plot
  par(mar=c(1, 9, 8, 2))
  
  # to automatize, need to use parameters
  xmx <- 25
  st.labels <- 0.90
  step.labels <- 1.50
  l.dens <- 9.5
  polygon.hei <- 3
  
  # prepare window to plot
  plot(1, xlim = c(0, xmx),  ylim=c(0, ymax), xlab='', ylab='', xaxt='none', 
       yaxt='none', col='white', bty='n')
  
  text(x = st.labels+step.labels, y = (ymax + st.labels), labels = "Effect on AD", adj = 0, xpd=T, srt=90, cex=1.50, font=1)
  text(x = st.labels+step.labels*2, y = (ymax + st.labels), labels = "Effect on Age", adj = 0, xpd=T, srt=90, cex=1.50, font=1)
  
  text(x = -3.5, y = (ymax + st.labels*1.5), labels = "Variant annotation", xpd=T, cex=1.50, font=1)
  
  # try to put only two polygons for the classes here
  polygon(x = c(st.labels+step.labels*3, st.labels+step.labels*3, st.labels+step.labels*3+l.dens), y = c((ymax + st.labels), (ymax + polygon.hei), (ymax + st.labels)), col = ggplot2::alpha(color.palette[5], 0.40), xpd=T)
  polygon(x = c(st.labels+step.labels*3+l.dens, st.labels+step.labels*3+l.dens*2, st.labels+step.labels*3+l.dens*2), y = c((ymax + st.labels), (ymax + polygon.hei), (ymax + st.labels)), col = ggplot2::alpha(color.palette[2], 0.40), xpd=T)
  polygon(x = c(st.labels+step.labels*3, st.labels+step.labels*3+l.dens, st.labels+step.labels*3+l.dens*2), y = c((ymax + polygon.hei+0.15), (ymax + st.labels+0.15), (ymax + polygon.hei+0.15)), col = ggplot2::alpha(color.palette[6], 0.40), xpd=T)
  text(x = st.labels+step.labels*3+0.5, y = (ymax + 1.25), labels = "Age effect", adj=0, xpd=T, font=2, cex = 1.40)
  text(x = st.labels+step.labels*3+l.dens*2-0.5, y = (ymax + 1.25), labels = "! Age effect", adj=1, xpd=T, font=2, cex = 1.40)
  text(x = st.labels+step.labels*3+l.dens, y = (ymax + 2.75), labels = "AD effect", adj=0.5, xpd=T, font=2, cex = 1.40)
  
  # add labels for expected and unexpected
  colors.unexp <- brewer.pal(n = 3, name = "Pastel1")
  rect(xleft = st.labels+step.labels*3, ybottom = ymax + 3.30, xright = st.labels+step.labels*3+l.dens, ytop = ymax + 4.30, col = ggplot2::alpha(colors.unexp[1], 0.6), xpd=T)
  rect(xleft = st.labels+step.labels*3+l.dens, ybottom = ymax + 3.30, xright = st.labels+step.labels*3+l.dens*2, ytop = ymax + 4.30, col = ggplot2::alpha(colors.unexp[2], 0.6), xpd=T)
  text(x = st.labels+step.labels*3, y = ymax + 3.8, labels = "Expected direction", adj = 0, xpd=T, font=2, cex=1.40)
  text(x = st.labels+step.labels*3+l.dens*2, y = ymax + 3.8, labels = "Unexpected direction", adj = 1, xpd=T, font=2, cex=1.40)
  
  c <- c + h
  
  # main loop across loci
  for (loc in 1:length(rev.names)){
    # select locus
    rsid_tmp <- rev.names[loc]
    # take variant density and put x and y in dataframe
    locus <- line.x[which(rownames(line.x) == rsid_tmp), ]
    #locus[is.na(locus)] <- median(na.omit(locus))
    locus <- locus[!is.na(locus)]
    locus_dens <- density(locus)
    dens <- data.frame(x = locus_dens$x, y = locus_dens$y)
    # normalize x-densities for the place where to plot [a, b]
    # in case of APOE as the distribution are cutted, I need to extend them
    if (rsid_tmp %in% c("rs429358", "rs7412")){
      tmp.up <- data.frame(x=0, y=0)
      tmp.down <- data.frame(x=2, y=0)
      dens <- rbind(tmp.up, dens, tmp.down)
    }
    b <- st.labels+step.labels*3+l.dens*2
    a <- st.labels+step.labels*3
    dens$x <- (b-a) * ((dens$x - min(dens$x))/(max(dens$x) - min(dens$x))) + a
    # also check the y and in case normalize between 0 and 0.4
    if (max(dens$y) > 2){ dens$y <- 1.5*(dens$y - min(dens$y))/(max(dens$y) - min(dens$y)) }
    lim.dens <- a
    lim.ad <- a + l.dens/2
    lim.unex <- lim.ad + l.dens/2
    lim.ext <- b
    dens$y <- dens$y + (ymax - c)
    dens <- subset(dens, ((dens$x > lim.dens) & (dens$x <= lim.ext)))
    
    # select square size
    sz <- 1.5
    
    # plot effect on AD
    #sb.lit <- lit_igap[which(lit_igap$rsid == rsid_tmp),]
    #rect(xleft = sz, ybottom = (ymax - c), xright = sz*2, ytop = (ymax - c + h), col=alpha(sb.lit$color.lit, min(sb.lit$norm.beta.lit+0.5, 1)), border = 'black')
    
    # plot effect on survival
    #sb.my <- meta[which(meta$rsid == rsid_tmp),]
    #rect(xleft = sz*2, ybottom = (ymax - c), xright = sz*3, ytop = (ymax - c + h), col=alpha(sb.my$color.mine, min(abs(sb.my$norm.beta.mine)+0.5, 1)), border = 'black')
    
    # plot locus name
    text(x = 1, y = (ymax - c + (h/2)), labels = rsid_tmp, font = 3, adj = 1, xpd=T, cex = 1)
    
    # plot chromosome:position
    #text(x = -6.5, y = (ymax - c + (h/2)), labels = sb.my$LOCUS, font = 3, adj = 0, xpd=T, cex = 1)
    
    # plot unexpected
    polygon(dens, col=color.palette[5])
    
    # plot larger effect on AD
    ad <- subset(dens, dens$x > lim.ad)
    supp <- data.frame(x = lim.ad, y = (ymax - c))
    ad.supp <- rbind(supp, ad)
    polygon(ad.supp, col=color.palette[6])
    
    # plot larger effect on survival
    unex <- subset(dens, dens$x > lim.unex)
    supp <- data.frame(x = lim.unex, y = (ymax - c))
    unex.supp <- rbind(supp, unex)
    polygon(unex.supp, col=color.palette[2])
    
    #increment c
    c <- c + h
    
    #draw line to divide loci
    segments(x0 = -12, y0 = (ymax - c + h), x1 = st.labels+step.labels*3+l.dens*2, y1 = (ymax - c + h), xpd=T)
  }
  
  #lines
  b <- st.labels+step.labels*3+l.dens*2
  a <- st.labels+step.labels*3
  segments(x0 = a+l.dens/2, y0 = 0, x1 = a+l.dens/2, y1 = ymax-1, lty = 2, lwd=2)
  segments(x0 = a+l.dens, y0 = 0, x1 = a+l.dens, y1 = ymax-1, lty = 1, lwd=3.0)
  segments(x0 = a+l.dens*2, y0 = 0, x1 = a+l.dens*2, y1 = ymax-1, lty = 1, lwd=3.0)
  
  #draw line to "divide" figure 
  segments(x0 = sz, y0 = 0, x1 = sz, y1 = ymax, xpd=T, lwd=3)
  segments(x0 = sz*3, y0 = 0, x1 = sz*3, y1 = ymax, xpd=T, lwd=3)
}

## Function to plot densities only
plotDensities_noDendro <- function(angle, lit, rand){
  # read my association file
  if (rand == FALSE){
    myass <- fread("INPUTS/myAssoc_ADsnps_GRACE.txt", h=T)
    # match same as the densities
    myass <- matchBy(myass, rownames(angle), myass$ID)
    lit <- matchBy(lit, rownames(angle), lit$rsid)
    
    # check alleles and in case change my beta
    lit$beta_adj <- lit$beta
    lit$beta_adj[which(lit$beta <0)] <- lit$beta_adj[which(lit$beta <0)]*(-1)
    lit$a1_adj <- toupper(lit$a1)
    lit$a1_adj[which(lit$beta <0)] <- toupper(lit$a2[which(lit$beta <0)])
    myass$BETA[which(myass$A1 != lit$a1_adj)] <- myass$BETA[which(myass$A1 != lit$a1_adj)]*(-1)
    
    # color palette
    color.palette <- brewer.pal(9, "Set1")
    
    # manage colors for beta -- normalize between 0 and 1
    min.beta <- min(na.omit(abs(lit$beta_adj), abs(myass$BETA)))
    max.beta <- max(na.omit(abs(lit$beta_adj), abs(myass$BETA)))
    lit$norm.beta.lit <- (abs(lit$beta_adj) - min.beta) / (max.beta - min.beta)
    myass$norm.beta.mine <- (abs(myass$BETA) - min.beta) / (max.beta - min.beta)
    lit$color.lit <- color.palette[1]
    lit$color.lit[which(lit$beta_adj > 0)] <- color.palette[4]
    myass$color.mine <- color.palette[1]
    myass$color.mine[which(myass$BETA < 0)] <- color.palette[4]
  }
  
  # PLOT
  # Global graphical parameters
  par(mar=c(1, 17, 16, 0), mfrow=c(1,1))
  # ymax for the plot
  ymax <- nrow(angle) + 1
  # define square height
  h <- 1
  # define counter for the plot
  c <- 1
  # to automatize, need to use parameters
  xmx <- 37
  st.labels <- 0.90
  step.labels <- 1.50
  l.dens <- 9.5
  polygon.hei <- 3
  
  # prepare window to plot
  plot(1, xlim = c(0, xmx),  ylim=c(0, ymax), xlab='', ylab='', xaxt='none', 
       yaxt='none', col='white', bty='n')
  
  # put text on columns
  text(x = st.labels, y = (ymax + st.labels), labels = "Beta-Amyloid", adj = 0, xpd=T, srt=90, cex=2, font=2)
  text(x = st.labels+step.labels, y = (ymax + st.labels), labels = "Lipid/Cholesterol", adj = 0, xpd=T, srt=90, cex=2, font=2)
  text(x = st.labels+step.labels*2, y = (ymax + st.labels), labels = "Endocytosis/Immune", adj = 0, xpd=T, srt=90, cex=2, font=2)
  text(x = st.labels+step.labels*3, y = (ymax + st.labels), labels = "Synaptic plasticity", adj = 0, xpd=T, srt=90, cex=2, font=2)
  
  text(x = st.labels+step.labels*4, y = (ymax + st.labels), labels = "Effect on AD", adj = 0, xpd=T, srt=90, cex=2, font=2)
  text(x = st.labels+step.labels*5, y = (ymax + st.labels), labels = "Effect on LGV", adj = 0, xpd=T, srt=90, cex=2, font=2)
  
  text(x = -4.75, y = (ymax + st.labels*1.5), labels = "Variant annotation", xpd=T, cex=2, font=2)
  
  # try to put only two polygons for the classes here
  # polygon(x = c(st.labels+step.labels*6, st.labels+step.labels*6, st.labels+step.labels*6+l.dens), y = c((ymax + st.labels), (ymax + polygon.hei), (ymax + st.labels)), col = alpha(color.palette[5], 0.40), xpd=T)
  # polygon(x = c(st.labels+step.labels*6+l.dens, st.labels+step.labels*6+l.dens*2, st.labels+step.labels*6+l.dens*2), y = c((ymax + st.labels), (ymax + polygon.hei), (ymax + st.labels)), col = alpha(color.palette[2], 0.40), xpd=T)
  # polygon(x = c(st.labels+step.labels*6, st.labels+step.labels*6+l.dens, st.labels+step.labels*6+l.dens*2), y = c((ymax + polygon.hei+0.15), (ymax + st.labels+0.15), (ymax + polygon.hei+0.15)), col = alpha(color.palette[6], 0.40), xpd=T)
  # text(x = st.labels+step.labels*6+0.5, y = (ymax + 1.50), labels = "CHA effect", adj=0, xpd=T, cex = 1.40, font=2)
  # text(x = st.labels+step.labels*6+l.dens*2-0.5, y = (ymax + 1.50), labels = "! CHA effect", adj=1, xpd=T, font=2, cex = 1.40)
  # text(x = st.labels+step.labels*6+l.dens, y = (ymax + 2.50), labels = "AD effect", adj=0.5, xpd=T, font=2, cex = 1.40)
  
  # add labels for expected and unexpected
  colors.unexp <- brewer.pal(n = 3, name = "Pastel1")
  # rect(xleft = st.labels+step.labels*6, ybottom = ymax + 3.30, xright = st.labels+step.labels*6+l.dens, ytop = ymax + 4.30, col = alpha(colors.unexp[1], 0.6), xpd=T)
  # rect(xleft = st.labels+step.labels*6+l.dens, ybottom = ymax + 3.30, xright = st.labels+step.labels*6+l.dens*2, ytop = ymax + 4.30, col = alpha(colors.unexp[2], 0.6), xpd=T)
  # text(x = st.labels+step.labels*6+l.dens/2, y = ymax + 3.8, labels = "Expected", adj = 0.5, xpd=T, font=2, cex=1.40)
  # text(x = st.labels+step.labels*6+(l.dens*3/2), y = ymax + 3.8, labels = "Unexpected", adj = 0.5, xpd=T, font=2, cex=1.40)
  # rect(xleft = st.labels+step.labels*6, ybottom = ymax + st.labels, xright = st.labels+step.labels*6+l.dens, ytop = ymax + st.labels*2, col = alpha(colors.unexp[1], 0.6), xpd=T)
  # rect(xleft = st.labels+step.labels*6+l.dens, ybottom = ymax + st.labels, xright = st.labels+step.labels*6+l.dens*2, ytop = ymax + st.labels*2, col = alpha(colors.unexp[2], 0.6), xpd=T)
  text(x = st.labels+step.labels*6+l.dens/2, y = ymax + st.labels*1.5, labels = "Expected", adj = 0.5, xpd=T, font=2, cex=2)
  text(x = st.labels+step.labels*6+(l.dens*3/2), y = ymax + st.labels*1.5, labels = "Unexpected", adj = 0.5, xpd=T, font=2, cex=2)
  
  # text about cell type
  text(x = st.labels+step.labels*7+l.dens*2, y = (ymax + st.labels), labels = "Astrocytes", adj = 0, xpd=T, srt=90, cex=2, font=2)
  text(x = st.labels+step.labels*8+l.dens*2, y = (ymax + st.labels), labels = "Oligodendrocytes", adj = 0, xpd=T, srt=90, cex=2, font=2)
  text(x = st.labels+step.labels*9+l.dens*2, y = (ymax + st.labels), labels = "Microglia", adj = 0, xpd=T, srt=90, cex=2, font=2)
  text(x = st.labels+step.labels*10+l.dens*2, y = (ymax + st.labels), labels = "Endotelial", adj = 0, xpd=T, srt=90, cex=2, font=2)
  text(x = st.labels+step.labels*11+l.dens*2, y = (ymax + st.labels), labels = "Neurons", adj = 0, xpd=T, srt=90, cex=2, font=2)
  
  # define positions
  p1 <- (st.labels+step.labels*3)/2 + st.labels/2
  p2 <- st.labels+step.labels*4.5
  p3 <- st.labels+step.labels*7+l.dens*2 + ((st.labels+step.labels*11+l.dens*2)-(st.labels+step.labels*7+l.dens*2))/2
  text(x = -4.75, y = (ymax - c + h), labels = 'A', font=2, xpd=T, cex = 2)
  text(x = p1, y = (ymax - c + h), labels = 'B', font=2, cex = 2)
  text(x = p2, y = (ymax - c + h), labels = 'C', font=2, cex = 2)
  text(x = st.labels+step.labels*6+l.dens, y = (ymax - c + h), labels = 'D', font=2, cex = 2)
  text(x = p3, y = (ymax - c + h), labels = 'E', font=2, cex = 2)
  segments(x0 = -8, y0 = (ymax - c), x1 = st.labels+step.labels*11+l.dens*2, y1 = (ymax - c), xpd=T, cex = 2)
  c <- c + h
  
  # select square size
  sz <- 1.5
  
  #draw line to "divide" figure 
  segments(x0 = 0, y0 = 0, x1 = 0, y1 = ymax, xpd=T, lwd=3)
  segments(x0 = sz*4, y0 = 0, x1 = sz*4, y1 = ymax, xpd=T, lwd=3)
  segments(x0 = sz*6, y0 = 0, x1 = sz*6, y1 = ymax, xpd=T, lwd=3)
  segments(x0 = st.labels+step.labels*6.5+l.dens*2, y0 = 0, x1 = st.labels+step.labels*6.5+l.dens*2, y1 = ymax, xpd=T, lwd=3)
  
  # last thing: big rectangles for clustering
  # rect(xleft = -8.5, ybottom = 10, xright = st.labels+step.labels*11.5+l.dens*2, ytop = 0, border="navy", lwd=6, xpd=T)
  # rect(xleft = -8.5, ybottom = 27, xright = st.labels+step.labels*11.5+l.dens*2, ytop = 10, border="orange", lwd=6, xpd=T)
  # rect(xleft = -8.5, ybottom = 27, xright = st.labels+step.labels*11.5+l.dens*2, ytop = 38, border="red", lwd=6, xpd=T)
  
  # also the labels of the groups
  # text(x = st.labels+step.labels*11.5+l.dens*2+0.5, y = 32.5, labels = "CHA-group", col="red", xpd=T, cex=1.40, font=4, srt=-90)
  # text(x = st.labels+step.labels*11.5+l.dens*2+0.5, y = 18.5, labels = "AD-group", col="orange", xpd=T, cex=1.40, font=4, srt=-90)
  # text(x = st.labels+step.labels*11.5+l.dens*2+0.5, y = 5, labels = "Unex-group", col="navy", xpd=T, cex=1.40, font=4, srt=-90)
  if (rand == FALSE){
    text(x = -10, y = 32.5, labels = "Longevity-group", col="red", xpd=T, cex=1.80, font=4, srt=-90)
    text(x = -11, y = 32.5, labels = expression(bolditalic(E[LGV]*" > "*E[AD])), col="red", xpd=T, cex=1.40, font=4, srt=-90)
    text(x = -10, y = 18.5, labels = "AD-group", col="orange", xpd=T, cex=1.80, font=4, srt=-90)
    text(x = -11, y = 18.5, labels = expression(bolditalic(E[LGV]*" < "*E[AD])), col="orange", xpd=T, cex=1.40, font=4, srt=-90)
    text(x = -10, y = 5, labels = "Unex-group", col="navy", xpd=T, cex=1.80, font=4, srt=-90)
    text(x = -11, y = 5, labels = "Unexpected direction", col="navy", xpd=T, cex=1.40, font=4, srt=-90)
  }
  
  # main loop across loci
  for (loc in 1:nrow(angle)){
    # select locus
    rsid_tmp <- rownames(angle)[loc]
    # take variant density and put x and y in dataframe
    locus <- density(angle[which(rownames(angle) == rsid_tmp), ])
    # also find median
    med_dens <- median(angle[which(rownames(angle) == rsid_tmp), ])
    dens <- data.frame(x = locus$x, y = locus$y)
    # normalize x-densities for the place where to plot [a, b]
    # in case of APOE as the distribution are cutted, I need to extend them
    if (rsid_tmp %in% c("rs429358", "rs7412")){
      tmp.up <- data.frame(x=-1, y=0)
      tmp.down <- data.frame(x=1, y=0)
      dens <- rbind(tmp.up, dens, tmp.down)
    }
    b <- st.labels+step.labels*6+l.dens*2
    a <- st.labels+step.labels*6
    dens$x <- (b-a) * ((dens$x - min(dens$x))/(max(dens$x) - min(dens$x))) + a
    # also normalize median in the same interval
    med_dens <- (b-a) * ((med_dens - min(angle))/(max(angle) - min(angle))) + a
    # also check the y and in case normalize between 0 and 0.4
    if (max(dens$y) >2){ dens$y <- 1.5*(dens$y - min(dens$y))/(max(dens$y) - min(dens$y)) }
    lim.dens <- a
    lim.ad <- a + l.dens/2
    lim.unex <- lim.ad + l.dens/2
    lim.ext <- b
    dens$y <- dens$y + (ymax - c)

    if (rand == FALSE){
      # plot effect on AD
      sb.lit <- lit[which(lit$rsid == rsid_tmp),]
      rect(xleft = sz*4, ybottom = (ymax - c), xright = sz*5, ytop = (ymax - c + h), col=ggplot2::alpha(sb.lit$color.lit, min(sb.lit$norm.beta.lit+0.5, 1)), border = 'black')
      
      # plot effect on survival
      sb.my <- myass[which(myass$ID == rsid_tmp),]
      rect(xleft = sz*5, ybottom = (ymax - c), xright = sz*6, ytop = (ymax - c + h), col=ggplot2::alpha(sb.my$color.mine, min(abs(sb.my$norm.beta.mine)+0.5, 1)), border = 'black')
      # also add * in case association was significant
      if (sb.my$P <= 0.05){ text(x =sz*5.5, y = (ymax - c + h/2), labels="*", cex=1.5, font=2) }
    }    
    # plot annotation
    # text(x = -11.25, y = (ymax - c + (h/2)), labels = path$locus, font=2, adj=0, xpd=T, cex=1)
    # text(x = -5.25, y = (ymax - c + (h/2)), labels = rsid_tmp, font=2, adj=0.5, xpd=T, cex=1)
    # text(x = -0.25, y = (ymax - c + (h/2)), labels = sb.lit$gene, font=4, adj=1, xpd=T, cex=1)
    #text(x = -9.25, y = (ymax - c + (h/2)), labels = path$locus, font=2, adj=0, xpd=T, cex=1)
    text(x = -4.5, y = (ymax - c + (h/2)), labels = rsid_tmp, font=2, adj=0.5, xpd=T, cex=1)
    text(x = -0.25, y = (ymax - c + (h/2)), labels = sb.lit$gene, font=4, adj=1, xpd=T, cex=1)
    
    # cell types
    cell.col <- "navy"
    
    if (rand == FALSE){
      if (loc == 1){
        rect(xleft = -9.5, ybottom = 38, xright = st.labels+step.labels*11.5+l.dens*2, ytop = 27.1, border="red", lwd=5, xpd=T)
      } else if (loc == 12){
        # last thing: big rectangles for clustering
        rect(xleft = -9.5, ybottom = 10.1, xright = st.labels+step.labels*11.5+l.dens*2, ytop = 26.9, border="orange", lwd=5, xpd=T)
        segments(x0 = -9.5, y0 = 27.1, x1 = st.labels+step.labels*11.5+l.dens*2, y1 = 27.1, col="red", lwd=5, xpd=T)
        segments(x0 = st.labels+step.labels*11.5+l.dens*2, y0 = 38, x1 = st.labels+step.labels*11.5+l.dens*2, y1 = 27.1, col="red", lwd=5, xpd=T)
      } else if (loc == nrow(angle)){
        # last thing: big rectangles for clustering
        segments(x0 = -9.5, y0 = 10.1, x1 = st.labels+step.labels*11.5+l.dens*2, y1 = 10.1, col="orange", lwd=5, xpd=T)
        segments(x0 = st.labels+step.labels*11.5+l.dens*2, y0 = 26.9, x1 = st.labels+step.labels*11.5+l.dens*2, y1 = 10.1, col="orange", lwd=5, xpd=T)
        rect(xleft = -9.5, ybottom = 9.9, xright = st.labels+step.labels*11.5+l.dens*2, ytop = 0, border="navy", lwd=5, xpd=T)
      }
    }    
    # plot unexpected
    polygon(dens, col=color.palette[5], border = NA)
    # add median -- need to find the y where the median intesect the density
    # calculate distance from the median and take the closest value
    tmp <- dens
    tmp$median_dist <- abs(med_dens - tmp$x)
    y_inters <- tmp$y[which(tmp$median_dist == min(tmp$median_dist))]
    
    # plot larger effect on AD
    ad <- subset(dens, dens$x > lim.ad)
    supp <- data.frame(x = lim.ad, y = (ymax - c))
    ad.supp <- rbind(supp, ad)
    polygon(ad.supp, col=color.palette[6], border = NA)
    
    # plot larger effect on survival
    unex <- subset(dens, dens$x > lim.unex)
    supp <- data.frame(x = lim.unex, y = (ymax - c))
    unex.supp <- rbind(supp, unex)
    polygon(unex.supp, col=color.palette[2], border = NA)
    
    # finally add medians
    segments(x0 = med_dens, y0 = min(dens$y), x1 = med_dens, y1 = y_inters, lwd=1, col="navy")
    
    # re-add border of the polygon
    polygon(dens, col = NULL, border = "black")
    
    #increment c
    c <- c + h
    
    #draw line to divide loci
    segments(x0 = -9.25, y0 = (ymax - c + h), x1 = st.labels+step.labels*11+l.dens*2, y1 = (ymax - c + h), xpd=T)
  }
  
  #lines
  b <- st.labels+step.labels*6+l.dens*2
  a <- st.labels+step.labels*6
  #segments(x0 = a+l.dens/2, y0 = 0, x1 = a+l.dens/2, y1 = ymax-1, lty = 2)
  segments(x0 = a+l.dens, y0 = 0, x1 = a+l.dens, y1 = ymax-1, lty = 1, lwd=1.5)
  
  # also legend at the very bottom
  text(x = -8.5, y = -0.5, labels = "LGV: longevity        *: p-value<0.05", font=2, cex=1, xpd=T, adj=0)
  segments(x0 = 3.4, y0 = -0.5, x1 = 4, y1 = -0.5, lwd = 2, col="navy", xpd=T)
  text(x = 5, y = -0.5, labels = "Median", cex=1, xpd=T, font=2)
}

## Function to do clustering on the normalized angle
ClusterJointDistribution <- function(angle){
  # Calculate distance -- manhattan as measure
  d <- dist(angle, method = "manhattan")
  
  # Do hierarchical clustering
  fit <- hclust(d) 
  plot(fit)
  
  # get names
  names.order <- fit$labels[order.dendrogram(dendr)]
  
}

## Function to compute bins intervals given bins number and x limits
function_bins <- function(bin.number, limits){
  #explicitize left and right limits
  left.l <- limits[1]
  right.l <- limits[2]
  
  #find bin size as total interval divided by number of bins
  bin.size <- (abs(left.l) + right.l) / bin.number
  
  #define bins matrix
  bins = matrix(data=0, nrow = bin.number, ncol = 3)
  
  #define starting point
  start = left.l
  
  #main loop
  for (i in 1:nrow(bins)){
    bins[i, 1] <- i
    bins[i, 2] <- start
    bins[i, 3] <- start + bin.size
    start = bins[i, 3]
  }
  
  return(bins)
}

## Function to compute the number of points per bin, per variant
PointsPerBin <- function(i, bins, angle){
  # Select variant
  tmp <- angle[i, ]
  
  # Create empty vector to fill with number of points per bin
  n_points <- rep(0, nrow(bins))
  
  # Main loop on bins
  for (j in 1:nrow(bins)){ n_points[j] <- length(tmp[which(tmp > bins[j, 2] & tmp < bins[j, 3])]) }
  
  return(n_points)
}

## Function to compare the number of points per variant: exp. vs. unexp // ad vs. aging
CompareNumber_ofPoints <- function(points_bin){
  # Define output
  out <- as.data.frame(matrix(data=NA, nrow=nrow(points_bin), ncol=4))
  colnames(out) <- c("p_exp_unexp", "direction_exp_unexp", "p_ad_cha", "direction_ad_cha")
  rownames(out) <- rownames(points_bin)
  
  # Main loop on variants
  for (i in 1:nrow(points_bin)){
    # Points in expected and unexpected direction
    exp <- points_bin[i, 1:50]
    unexp <- points_bin[i, 51:100]
    # Wilcoxon sum rank test
    wil <- wilcox.test(x = exp, y = unexp, alternative = "greater")
    # Assign to output -- this depends on the direction
    if (wil$p.value > 0.5){ 
      out[i, "p_exp_unexp"] <- 1-wil$p.value
      out[i, "direction_exp_unexp"] <- "larger_unexp"
    } else {
      out[i, "p_exp_unexp"] <- wil$p.value
      out[i, "direction_exp_unexp"] <- "larger_exp"
    }
    ####################################
    
    # Points in expected and unexpected direction
    ad <- points_bin[i, 26:50]
    cha <- points_bin[i, 1:25]
    # Wilcoxon sum rank test
    wil <- wilcox.test(x = ad, y = cha, alternative = "greater")
    # Assign to output -- this also depends on the direction
    if (wil$p.value > 0.5){
      out[i, "p_ad_cha"] <- 1-wil$p.value
      out[i, "direction_ad_cha"] <- "larger_cha"
    } else {
      out[i, "p_ad_cha"] <- wil$p.value
      out[i, "direction_ad_cha"] <- "larger_ad"
    }
    #####################################
  }
  # Correct p-values independently (exp vs unexp) and (ad vs. cha)  
  out$p_exp_unexp <- p.adjust(out$p_exp_unexp, method = "bonferroni", n = 2*nrow(out))
  out$p_ad_cha <- p.adjust(out$p_ad_cha, method = "bonferroni", n = 2*nrow(out))
  
  # We need to define groups:
  # 1: non significantly different between expected and unexpected --> noise
  # 2: significant more points in unexpcted direction --> pleiotropy
  # 3: significant more points in expected, no difference between ad and cha --> more ad
  # 4: significant more points in expected, cha side --> effect on multiple phenotypes other than ad
  non_sign_exp_unex <- out[which(out$p_exp_unexp > 0.05),]
  sign_unex <- out[which(out$p_exp_unexp <= 0.05 & out$direction_exp_unexp == "larger_unexp"),]
  sign_exp_ad <- out[which(out$p_exp_unexp <= 0.05 & out$direction_exp_unexp == "larger_exp" & (out$p_ad_cha > 0.05 | out$direction_ad_cha == "larger_ad" )),]
  sign_exp_cha <- out[which(out$p_exp_unexp <= 0.05 & out$direction_exp_unexp == "larger_exp" & out$p_ad_cha <= 0.05 & out$direction_ad_cha == "larger_cha"),]
  
  # also look at the global significance in expected vs. unexpected direction
  exp <- colSums(points_bin[, 1:50])
  unexp <- colSums(points_bin[, 51:100])
  wilcox.test(x = exp, y = unexp, alternative = "greater")
  
  # Put results together
  l = list(non_sign_exp_unex, sign_unex, sign_exp_ad, sign_exp_cha)
  return(l)
}

## Function to reorder densitites based on median value of the angle
Reorder_median <- function(angle){
  # Define output first
  tmp_df <- as.data.frame(matrix(data=NA, nrow=nrow(angle), ncol=3))
  colnames(tmp_df) <- c("snp", "median", "mean")
  
  # Main loop
  for (i in 1:nrow(angle)){ tmp_df[i, ] <- c(rownames(angle)[i], median(angle[i, ]), mean(angle[i, ])) }
  
  # Order by median
  tmp_df$median = as.numeric(tmp_df$median)
  tmp_df <- tmp_df[order(tmp_df$median),]
  
  return(tmp_df)
}

## Function to plot distribution of variant groups 
Plot_variantGroups <- function(res_groups){
  # Divide the groups from the main list
  non_sign_exp_unex <- res_groups[[1]]
  sign_unex <- res_groups[[2]]
  sign_exp_ad <- res_groups[[3]]
  sign_exp_cha <- res_groups[[4]]
  
  # colors and graphical parameters
  cols <- brewer.pal(n = 6, name = "Set1")
  par(mar=c(6, 4, 4, 6), mfrow=c(1, 1))
  
  # plot
  pdf("RESULTS/distribution_variant_groups.pdf", height = 7, width = 10)
  barplot(c(nrow(non_sign_exp_unex), nrow(sign_unex), nrow(sign_exp_ad), nrow(sign_exp_cha)),
          names=c("Non sign (Exp vs. Unex)", "Sign Unex", "Sign Exp AD>CHA", "Sign Exp CHA>AD"),
          col=cols, main="Distribution of variant groups", ylab="# of variants")
  dev.off()
}

## Function to merge sampling procedure for the functional annotation -- clusterProfiler
mergeSampling_clusProf <- function(){
  # load results
  load("RESULTS/Results_enrichment_sampling_clusPro.RData")
  tr <- enrich.res
  # extract terms common to all datasets
  all <- rbindlist(tr)
  counts <- as.data.frame(table(all$ID))
  common <- as.character(counts[which(counts$Freq == length(tr)), "Var1"])
  
  # function to do mean
  tmp_f <- function(i, all){
    # subset to term of interest
    sb <- all[which(all$ID == i),]
    
    # order by number of characters of the gene IDs
    sb <- sb[order(-nchar(sb$geneID)),]
    
    # put in separate dataframe
    df <- data.frame(ID=i, Description=unique(sb$Description), p_adj = mean(sb$p.adjust), q_value = mean(sb$qvalue), geneID=sb$geneID[1])
    
    return(df)
  }
  
  # run in multiprocessing
  terms.avg <- mclapply(common, tmp_f, all=all, mc.cores=1)
  all_avg <- rbindlist(terms.avg)
  
  # order by p
  all_avg <- all_avg[order(all_avg$p_adj),]
  dim(all_avg[which(all_avg$q_value <= 0.01),])
  dim(all_avg[which(all_avg$p_adj <= 0.05),])
  
  # write output for further analyses
  write.table(all_avg[, c("ID", "p_adj")], "RESULTS/revigo_inp_clusterProfiler.txt", quote=F, row.names=F, col.names=F, sep="\t")
  
  # also write the full list of enrichment results
  write.table(all_avg, "RESULTS/enrichent_results_sampling.txt", quote=F, row.names=F, col.names=F, sep="\t")
  
  return(all_avg)
}

## Function to link variant to pathways through gene-overlap
FindGenesPerCluster_clusProf <- function(functional_clusters, all_avg, lin_matrix){
  ## read the variant-gene mapping also
  var_gene_mapping <- fread("RESULTS/AnnotateMe/snp_annotation.txt", h=T, stringsAsFactors=F)
  geneList <- read.table("RESULTS/AnnotateMe/snp_annotation_geneList.txt", h=F, stringsAsFactors=F)
  enrich_sign <- all_avg[which(all_avg$p_adj <= 0.05),]
  var_gene_mapping <- data.frame(var_gene_mapping)
  
  # identify number of clusters
  n_clust <- length(unique(functional_clusters$cluster))
  if (n_clust == 6){
    # first identify cluster names
    clust_names <- c("Beta-Amyloid", "Lipid metabolism", "Immune activation", "Synaptic plasticity", "Protein activity", "Lipid transport")
    
    # output will be in the var_gene_mapping -- need to generate the respective columns
    for (n in clust_names){ var_gene_mapping[, n] <- 0 }    
  } else if (n_clust == 4){
    # first identify cluster names
    clust_names <- c("Beta-Amyloid", "Lipid metabolism", "Immune activation", "Synaptic plasticity")
    
    # output will be in the var_gene_mapping -- need to generate the respective columns
    for (n in clust_names){ var_gene_mapping[, n] <- 0 }
  }

  # main loop on variants
  for (i in 1:nrow(var_gene_mapping)){
    # get genes
    sb_genes <- unlist(strsplit(var_gene_mapping$geneList[i], ","))
    
    # grep in the enrichment results
    tmp_f <- function(g, all_avg){ return(all_avg[grep(g, all_avg$geneID),]) }
    res_grep <- lapply(sb_genes, tmp_f, all_avg=enrich_sign)
    res_grep <- rbindlist(res_grep)
    
    # look at which functional cluster these terms are involved
    sb_clusters <- as.data.frame(table(functional_clusters[which(functional_clusters$term %in% res_grep$ID), "cluster"]))
    sb_clusters$Prop <- sb_clusters$Freq/sum(sb_clusters$Freq)
    sb_clusters$Var1 <- as.numeric(as.character(sb_clusters$Var1))
    
    # assign in var_gene_mapping dataset
    for (cl in sb_clusters$Var1){
      var_gene_mapping[i, (15+cl)] <- as.numeric(sb_clusters$Prop[which(sb_clusters$Var1 == cl)])
    }
  }
  
  # this is nice but there are missings: the missings are those that are associated with GO that are not significantly enriched
  # now try for the missing to find the cluster that best fit to the annotation
  miss <- var_gene_mapping[which(rowSums(var_gene_mapping[, 16:(16+n_clust-1)]) ==0),]
  
  # define output
  out <- as.data.frame(matrix(data=NA, nrow=nrow(miss), ncol=n_clust))
  colnames(out) <- seq(1, n_clust)
  rownames(out) <- miss$ID
  
  # prepare training dset for the knn model
  train <- as.data.frame(lin_matrix)
  names_go <- colnames(train)
  rownames(train) <- names_go
  # prepare groups from hierarchical clustering
  train <- matchBy(train, functional_clusters$term, rownames(train))
  train <- cbind(train, cluster=functional_clusters$cluster)
  train$cluster <- as.factor(train$cluster)
  # replace : with _ cause it is causing issues for the model
  names_go <- str_replace_all(string = names_go, pattern = ":", replacement = "_")
  rownames(train) <- names_go
  colnames(train) <- c(names_go, "cluster")
  # knn model on the training
  knn_mod = knn3(cluster ~ . , data = train, k = n_clust)
  
  # set up bioMart
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  # loop on the missing
  for (i in 1:nrow(miss)){
    # get genes
    sb_genes <- unlist(strsplit(miss$geneList[i], ","))
    print(paste("   Annotating genes: ", sb_genes, sep=""))
    
    # look for the biological processes the gene is invovled into
    results <- getBM(attributes = c("external_gene_name", "go_id", "namespace_1003"), 
                     filters = "external_gene_name", values = sb_genes, mart = mart)
    
    # restrict to biological processes only
    results <- results[which(results$namespace_1003 == "biological_process"),]
    
    # grep in the enrichment results
    #res_grep <- lapply(sb_genes, tmp_f, all_avg = all_avg)
    #res_grep <- rbindlist(res_grep)
    
    # grep genes in the full enrichment profile
    go <- results$go_id
    do_annot <- TRUE
    if (length(go) == 1){
      goSim_mt <- data.frame(term=all_avg$ID, Sim=goSim(all_avg$ID, go, semData=hsGO, measure="Lin"))
    } else if (length(go) >1){
      goSim_mt <- as.data.frame(matrix(data=NA, nrow = nrow(enrich_sign), ncol = length(go)))
      rownames(goSim_mt) <- colnames(lin_matrix)
      for (x in 1:length(go)){
        goSim_tmp <- data.frame(term=enrich_sign$ID, Sim=goSim(enrich_sign$ID, go[x], semData=hsGO, measure="Lin"))
        goSim_mt[, x] <- goSim_tmp$Sim
      } 
    } else {
      do_annot <- FALSE
    }
    goSim_mt[is.na(goSim_mt)] <- 0
    
    # if the matrix is empty (colSum==0), don't do the annotation but report this as missing (NA)
    if (do_annot == FALSE){
      # assign NA to out object
      out[i, ] <- rep(NA, n_clust)
    } else {
      # use knn to classify the new data
      # test using knn model
      test <- data.frame(t(goSim_mt[, 1:ncol(goSim_mt)]))
      colnames(test) <- names_go
      # return predicted probabilities
      probs <- predict(knn_mod, test, type = "prob")
      
      # if multiple rows, do the mean
      if (nrow(probs) >1){
        probs <- t(data.frame(colMeans(probs)))
      }
      
      # finally assign to out object
      out[i, ] <- probs[1, ]
    }    
  }
  
  # assign to miss dataset
  if (n_clust == 6){
    miss$"Beta-Amyloid" <- out[, 1]
    miss$"Lipid metabolism" <- out[, 2]
    miss$"Immune activation" <- out[, 3]
    miss$"Synaptic plasticity" <- out[, 4]
    miss$"Protein activity" <- out[, 5]
    miss$"Lipid transport" <- out[, 6]
  } else if (n_clust == 4){
    miss$"Beta-Amyloid" <- out[, 1]
    miss$"Lipid metabolism" <- out[, 2]
    miss$"Immune activation" <- out[, 3]
    miss$"Synaptic plasticity" <- out[, 4]
  }  
  
  # re-assign to main dataset
  var_gene_mapping <- var_gene_mapping[which(!(var_gene_mapping$ID %in% miss$ID)),]
  var_gene_mapping <- rbind(var_gene_mapping, miss)
  
  return(var_gene_mapping)
}

# Compare the annotation of the different functional annotations across density groups
CompareClusterAnnot_mod_3grps <- function(var_gene_mapping, angle, all_gr, n_clust){
  # Define main output -- this depends on how many clusters we have
  if (n_clust == 4){
    out_means <- as.data.frame(matrix(data=NA, nrow=3, ncol=n_clust))
    out_sd <- as.data.frame(matrix(data=NA, nrow=3, ncol=n_clust))
    rownames(out_means) <- c("cha", "ad", "unexp")
    colnames(out_means) <- names(var_gene_mapping[16:19])
    rownames(out_sd) <- c("cha", "ad", "unexp")
    colnames(out_sd) <- names(var_gene_mapping[16:19])
  } else if (n_clust == 6){
    out_means <- as.data.frame(matrix(data=NA, nrow=3, ncol=n_clust))
    out_sd <- as.data.frame(matrix(data=NA, nrow=3, ncol=n_clust))
    rownames(out_means) <- c("cha", "ad", "unexp")
    colnames(out_means) <- names(var_gene_mapping[16:21])
    rownames(out_sd) <- c("cha", "ad", "unexp")
    colnames(out_sd) <- names(var_gene_mapping[16:21])
  }
    
  # Main loop on the groups
  counter <- 1
  for (x in all_gr){
    # get variants in the cluster
    snp_info <- var_gene_mapping[which(var_gene_mapping$ID %in% x),]
    
    # convert NAs to 0 and exclude empty rows
    snp_info[is.na(snp_info)] <- 0
    snp_info <- snp_info[which(rowSums(snp_info[, 16:(16+n_clust-1)]) != 0),]
    snp_info <- snp_info[!duplicated(snp_info$geneList),]
    # store mean in output dataframe
    if (n_clust == 4){
      out_means[counter, ] <- c(mean(snp_info$"Beta-Amyloid"), mean(snp_info$`Lipid metabolism`), mean(snp_info$`Immune activation`),
                              mean(snp_info$`Synaptic plasticity`))
      out_sd[counter, ] <- c(sd(snp_info$"Beta-Amyloid"), sd(snp_info$`Lipid metabolism`), sd(snp_info$`Immune activation`),
                           sd(snp_info$`Synaptic plasticity`))
      counter <- counter + 1
    } else if (n_clust == 6){
      out_means[counter, ] <- c(mean(snp_info$"Beta-Amyloid"), mean(snp_info$`Lipid metabolism`), mean(snp_info$`Immune activation`),
                                mean(snp_info$`Synaptic plasticity`), mean(snp_info$`Protein activity`), mean(snp_info$`Lipid transport`))
      out_sd[counter, ] <- c(sd(snp_info$"Beta-Amyloid"), sd(snp_info$`Lipid metabolism`), sd(snp_info$`Immune activation`),
                             sd(snp_info$`Synaptic plasticity`), sd(snp_info$`Protein activity`), sd(snp_info$`Lipid transport`))
      counter <- counter + 1
    }
  }
  
  # Now compare within groups
  within_groups <- list()
  counter <- 1
  for (x in all_gr){
    # take variants and annontation from that group
    sb <- var_gene_mapping[which(var_gene_mapping$ID %in% x),]
    # exclude NAs
    sb[is.na(sb)] <- 0
    sb <- sb[which(rowSums(sb[, 16:(16+n_clust-1)]) != 0),]
    sb <- sb[!duplicated(sb$geneList),]
    if (n_clust == 4){
      # Do the tests -- it is 6 tests per group
      t1 <- wilcox.test(x = sb$`Beta-Amyloid`, y = sb$`Lipid metabolism`, alternative = "two.sided")
      t2 <- wilcox.test(x = sb$`Beta-Amyloid`, y = sb$`Immune activation`, alternative = "two.sided")
      t3 <- wilcox.test(x = sb$`Beta-Amyloid`, y = sb$`Synaptic plasticity`, alternative = "two.sided")
      t4 <- wilcox.test(x = sb$`Lipid metabolism`, y = sb$`Immune activation`, alternative = "two.sided")
      t5 <- wilcox.test(x = sb$`Lipid metabolism`, y = sb$`Synaptic plasticity`, alternative = "two.sided")
      t6 <- wilcox.test(x = sb$`Immune activation`, y = sb$`Synaptic plasticity`, alternative = "two.sided")
      tmp_df <- data.frame(Amy_Lip=t1$p.value, Amy_Imm=t2$p.value, Amy_Syn=t3$p.value, Lip_Imm=t4$p.value, Lip_Syn=t5$p.value, Imm_Syn=t6$p.value)
      within_groups[[counter]] <- tmp_df
      counter <- counter + 1
    } else if (n_clust == 6){
      # Do the tests -- it is 6 tests per group
      t1 <- wilcox.test(x = sb$`Beta-Amyloid`, y = sb$`Lipid metabolism`, alternative = "two.sided")
      t2 <- wilcox.test(x = sb$`Beta-Amyloid`, y = sb$`Immune activation`, alternative = "two.sided")
      t3 <- wilcox.test(x = sb$`Beta-Amyloid`, y = sb$`Synaptic plasticity`, alternative = "two.sided")
      t4 <- wilcox.test(x = sb$`Beta-Amyloid`, y = sb$`Protein activity`, alternative = "two.sided")
      t5 <- wilcox.test(x = sb$`Beta-Amyloid`, y = sb$`Lipid transport`, alternative = "two.sided")
      t6 <- wilcox.test(x = sb$`Lipid metabolism`, y = sb$`Immune activation`, alternative = "two.sided")
      t7 <- wilcox.test(x = sb$`Lipid metabolism`, y = sb$`Synaptic plasticity`, alternative = "two.sided")
      t8 <- wilcox.test(x = sb$`Lipid metabolism`, y = sb$`Protein activity`, alternative = "two.sided")
      t9 <- wilcox.test(x = sb$`Lipid metabolism`, y = sb$`Lipid transport`, alternative = "two.sided")
      t10 <- wilcox.test(x = sb$`Immune activation`, y = sb$`Synaptic plasticity`, alternative = "two.sided")
      t11 <- wilcox.test(x = sb$`Immune activation`, y = sb$`Protein activity`, alternative = "two.sided")
      t12 <- wilcox.test(x = sb$`Immune activation`, y = sb$`Lipid transport`, alternative = "two.sided")
      t13 <- wilcox.test(x = sb$`Protein activity`, y = sb$`Lipid transport`, alternative = "two.sided")
      tmp_df <- data.frame(Amy_Lip=t1$p.value, Amy_Imm=t2$p.value, Amy_Syn=t3$p.value, Amy_Pro=t4$p.value, Amy_Tra=t5$p.value,
                           Lip_Imm=t6$p.value, Lip_Syn=t7$p.value, Lip_Pro=t8$p.value, Lip_Tra=t9$p.value,
                           Imm_Syn=t10$p.value, Imm_Pro=t11$p.value, Imm_Tra=t12$p.value, Pro_Tra=t13$p.value)
      within_groups[[counter]] <- tmp_df
      counter <- counter + 1
    }
  }
  within_groups <- rbindlist(within_groups)
  rownames(within_groups) <- c("cha", "ad", "unexp")
  
  # Now between groups
  between_groups <- list()
  counter <- 1
  for (x in 2:length(all_gr)-1){
    for (j in (x+1):length(all_gr)){
      # Take the two datasets to compare
      df1 <- var_gene_mapping[which(var_gene_mapping$ID %in% all_gr[[x]]),]
      df2 <- var_gene_mapping[which(var_gene_mapping$ID %in% all_gr[[j]]),]
      # Clean
      df1[is.na(df1)] <- 0
      df2[is.na(df2)] <- 0
      df1 <- df1[which(rowSums(df1[, 16:(16+n_clust-1)]) != 0),]
      df2 <- df2[which(rowSums(df2[, 16:(16+n_clust-1)]) != 0),]
      df1 <- df1[!duplicated(df1$geneList),]
      df2 <- df2[!duplicated(df2$geneList),]
      if (n_clust == 4){
        # Do the test for the 4 functional clusters
        t1 <- wilcox.test(x = df1$`Beta-Amyloid`, y = df2$`Beta-Amyloid`, alternative = "two.sided")
        t2 <- wilcox.test(x = df1$`Lipid metabolism`, y = df2$`Lipid metabolism`, alternative = "two.sided")
        t3 <- wilcox.test(x = df1$`Immune activation`, y = df2$`Immune activation`, alternative = "two.sided")
        t4 <- wilcox.test(x = df1$`Synaptic plasticity`, y = df2$`Synaptic plasticity`, alternative = "two.sided")
        tmp_df <- data.frame(comparison=paste(x, j, sep="_"), Amy=t1$p.value, Lip=t2$p.value, Imm=t3$p.value, Syn=t4$p.value)
        between_groups[[counter]] <- tmp_df
        counter <- counter + 1
      } else if (n_clust == 6){
        # Do the test for the 4 functional clusters
        t1 <- wilcox.test(x = df1$`Beta-Amyloid`, y = df2$`Beta-Amyloid`, alternative = "two.sided")
        t2 <- wilcox.test(x = df1$`Lipid metabolism`, y = df2$`Lipid metabolism`, alternative = "two.sided")
        t3 <- wilcox.test(x = df1$`Immune activation`, y = df2$`Immune activation`, alternative = "two.sided")
        t4 <- wilcox.test(x = df1$`Synaptic plasticity`, y = df2$`Synaptic plasticity`, alternative = "two.sided")
        t5 <- wilcox.test(x = df1$`Protein activity`, y = df2$`Protein activity`, alternative = "two.sided")
        t6 <- wilcox.test(x = df1$`Lipid transport`, y = df2$`Lipid transport`, alternative = "two.sided")
        tmp_df <- data.frame(comparison=paste(x, j, sep="_"), Amy=t1$p.value, Lip=t2$p.value, Imm=t3$p.value, Syn=t4$p.value, Pro=t5$p.value, Tra=t6$p.value)
        between_groups[[counter]] <- tmp_df
        counter <- counter + 1
      }
    }
  }
  between_groups <- rbindlist(between_groups)

  # For the comparison between groups, also do 1 groups vs. the others
  one_vs_all <- list()
  counter <- 1
  for (x in 1:length(all_gr)){
    # Take the two datasets to compare
    df1 <- var_gene_mapping[which(var_gene_mapping$ID %in% all_gr[[x]]),]
    df2 <- var_gene_mapping[which(!(var_gene_mapping$ID %in% all_gr[[x]])),]
    # Clean
    df1[is.na(df1)] <- 0
    df2[is.na(df2)] <- 0
    df1 <- df1[which(rowSums(df1[, 16:(16+n_clust-1)]) != 0),]
    df2 <- df2[which(rowSums(df2[, 16:(16+n_clust-1)]) != 0),]
    df1 <- df1[!duplicated(df1$geneList),]
    df2 <- df2[!duplicated(df2$geneList),]
    if (n_clust == 4){
      # Do the test for the 4 functional clusters
      t1 <- wilcox.test(x = df1$`Beta-Amyloid`, y = df2$`Beta-Amyloid`, alternative = "two.sided")
      t2 <- wilcox.test(x = df1$`Lipid metabolism`, y = df2$`Lipid metabolism`, alternative = "two.sided")
      t3 <- wilcox.test(x = df1$`Immune activation`, y = df2$`Immune activation`, alternative = "two.sided")
      t4 <- wilcox.test(x = df1$`Synaptic plasticity`, y = df2$`Synaptic plasticity`, alternative = "two.sided")
      tmp_df <- data.frame(comparison=paste(x, "_vs_all", sep=""), Amy=t1$p.value, Lip=t2$p.value, Imm=t3$p.value, Syn=t4$p.value)
      one_vs_all[[counter]] <- tmp_df
      counter <- counter + 1
    } else if (n_clust == 6){
      # Do the test for the 4 functional clusters
      t1 <- wilcox.test(x = df1$`Beta-Amyloid`, y = df2$`Beta-Amyloid`, alternative = "two.sided")
      t2 <- wilcox.test(x = df1$`Lipid metabolism`, y = df2$`Lipid metabolism`, alternative = "two.sided")
      t3 <- wilcox.test(x = df1$`Immune activation`, y = df2$`Immune activation`, alternative = "two.sided")
      t4 <- wilcox.test(x = df1$`Synaptic plasticity`, y = df2$`Synaptic plasticity`, alternative = "two.sided")
      t5 <- wilcox.test(x = df1$`Protein activity`, y = df2$`Protein activity`, alternative = "two.sided")
      t6 <- wilcox.test(x = df1$`Lipid transport`, y = df2$`Lipid transport`, alternative = "two.sided")
      tmp_df <- data.frame(comparison=paste(x, "_vs_all", sep=""), Amy=t1$p.value, Lip=t2$p.value, Imm=t3$p.value, Syn=t4$p.value, Pro=t5$p.value, Tra=t6$p.value)
      one_vs_all[[counter]] <- tmp_df
      counter <- counter + 1
    }
  }
  one_vs_all <- rbindlist(one_vs_all)
  
  l=list(out_means, out_sd, between_groups, within_groups, one_vs_all)
  return(l)
}

## function to create n sampling dsets
samplingGset <- function(i, mapping){
  all.genes <- strsplit(mapping$geneList, ",")
  tmp_f <- function(x, all.genes){ tmp <- all.genes[[x]]; g <- sample(x=tmp, size=1); return(g)}
  gset <- lapply(1:length(all.genes), tmp_f, all.genes=all.genes)
  gset.clean <- unlist(gset)
  return(gset.clean)
}

## function to perform overlap analysis over the sampling datasets
overlapAnalysis <- function(i, gene_list, source_gset){
  g <- gene_list[[i]]
  
  res <- suppressMessages(gost(g, organism = "hsapiens", ordered_query = FALSE, multi_query = FALSE, significant = FALSE, exclude_iea = TRUE,
                               measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 1, correction_method = "fdr",
                               domain_scope = "annotated", custom_bg = NULL, numeric_ns = "", sources = source_gset))
  return(res$result)
}

## function to merge sampling procedure for the functional annotation
mergeSampling <- function(enrich.res){
  n.sampl = length(enrich.res)
  all.enrich <- do.call("rbind", enrich.res)      #put all results together
  sign.only <- all.enrich[which(all.enrich$p_value <= 0.05),]         #isolate only significant findings
  if (nrow(sign.only) >0){
    sign.table <- as.data.frame(table(sign.only$term_name))
    info <- sign.only[, c("term_name", "term_id")]
    sign.table <- sign.table[order(-sign.table$Freq),]
    sign.table <- merge(sign.table, info, by.x="Var1", by.y="term_name")
  } else { sign.table <- data.frame(chr=NULL, position=NULL)}
  
  # also try with averaging pvalues
  if ("intersection" %in% colnames(all.enrich)){
    terms <- all.enrich[, c("term_name", "term_id", "intersection", "p_value")]
    slow <- TRUE
  } else {
    terms <- all.enrich[, c("term_name", "term_id", "p_value")]
    slow = FALSE
  }
  # see how many times each term was tested and resctrict to common (tested all times)
  ttt = as.data.frame(table(terms$term_id))
  ttt <- ttt[which(ttt$Freq == n.sampl),]
  terms <- terms[which(terms$term_id %in% ttt$Var1),]
  #terms <- terms[!duplicated(terms$term_name),]
  avgP <- function(i, terms, term_list, slow){
    if (isTRUE(slow)){
      all_intersections <- paste0(unique(unlist(strsplit(paste0(terms$intersection[which(terms$term_id == term_list[i])], collapse = ","), ","))), collapse = ",")
      df <- data.frame(term_name = unique(terms$term_name[which(terms$term_id == term_list[i])]), term_id = term_list[i], avgP = mean(terms$p_value[which(terms$term_id == term_list[i])]), intersection = all_intersections)
    } else {
      df <- data.frame(term_name = unique(terms$term_name[which(terms$term_id == term_list[i])]), term_id = term_list[i], avgP = mean(terms$p_value[which(terms$term_id == term_list[i])]))
    }
    return(df) 
  }
  term_list <- unique(terms$term_id)
  all.terms.avg <- rbindlist(mclapply(1:length(term_list), avgP, terms=terms, term_list = term_list, slow = slow, mc.cores=2))
  all.terms.avg <- all.terms.avg[order(all.terms.avg$avgP),]
  all.terms.avg$log10P <- -log10(all.terms.avg$avgP)
  print(head(all.terms.avg))
  #write.table(all.terms.avg, paste("RESULTS_", random_num, "/enrichent_results_sampling.txt", sep=""), quote=F, row.names=F, sep="\t")
  
  # if (nrow(all.terms.avg[which(all.terms.avg$avgP <= 0.05),]) >0){
  #   pdf(paste("RESULTS_", random_num, "/treemap_sampling_PRS_genes.pdf", sep=""), height=10, width=10)
  #   treemap(all.terms.avg[which(all.terms.avg$avgP <= 0.05),], index="term_name", vSize="log10P", type="index")
  #   invisible(dev.off())
  # }
  
  l <- list(sign.table, all.terms.avg)
  return(l)
}

# Alternative function to test between groups weights
betweenGroupsTestEnrichment = function(var_gene_mapping, all_gr){
  # add group label to var_gene_mapping
  var_gene_mapping$group = NA
  for (i in 1:nrow(var_gene_mapping)){
    if (var_gene_mapping$ID[i] %in% all_gr[[1]]){
      var_gene_mapping$group[i] = "CHA"
    } else if (var_gene_mapping$ID[i] %in% all_gr[[2]]){
      var_gene_mapping$group[i] = "AD"
    } else {
      var_gene_mapping$group[i] = "UNEX"
    }
  }
  # for each group, do gene-set enrichment with sampling
  for (gr in unique(var_gene_mapping$group)){
    tmp = var_gene_mapping[which(var_gene_mapping$group == gr),]
    n.sampl <- 500
    tmp_sam_dsets = mclapply(1:n.sampl, samplingGset, mapping = tmp, mc.cores = 4)
    enrich.res <- mclapply(1:n.sampl, overlapAnalysis, gene_list=tmp_sam_dsets, source_gset = c("GO:BP"), mc.cores=4)
    sampling.res <- mergeSampling(enrich.res)
    
  }
}

# Compare the annotation of the different functional annotations across density groups
CompareClusterAnnot_cellTypes <- function(var_gene_mapping, angle, all_gr){
  # Define main output
  out_means <- as.data.frame(matrix(data=NA, nrow=3, ncol=5))
  out_sd <- as.data.frame(matrix(data=NA, nrow=3, ncol=5))
  rownames(out_means) <- c("cha", "ad", "unexp")
  colnames(out_means) <- c("astrocytes", "endothelial", "myeloid", "neuron", "oligodendrocytes")
  rownames(out_sd) <- c("cha", "ad", "unexp")
  colnames(out_sd) <- c("astrocytes", "endothelial", "myeloid", "neuron", "oligodendrocytes")
  
  # Main loop on the groups
  counter <- 1
  for (x in all_gr){
    # get variants in the cluster
    snp_info <- var_gene_mapping[which(var_gene_mapping$ID %in% x),]
    
    # convert NAs to 0 and exclude empty rows
    snp_info[is.na(snp_info)] <- 0
    snp_info <- snp_info[which(rowSums(snp_info[, c("astrocytes", "endothelial", "myeloid", "neuron", "oligodendrocytes")]) != 0),]
    snp_info <- snp_info[!duplicated(snp_info$geneList),]
    # store mean in output dataframe
    out_means[counter, ] <- c(mean(snp_info$astrocytes), mean(snp_info$endothelial), mean(snp_info$myeloid),
                              mean(snp_info$neuron), mean(snp_info$oligodendrocytes))
    out_sd[counter, ] <- c(sd(snp_info$astrocytes), sd(snp_info$endothelial), sd(snp_info$myeloid),
                           sd(snp_info$neuron), sd(snp_info$oligodendrocytes))
    counter <- counter + 1
  }
  
  # Now compare within groups
  within_groups <- list()
  counter <- 1
  for (x in all_gr){
    # take variants and annontation from that group
    sb <- var_gene_mapping[which(var_gene_mapping$ID %in% x),]
    # exclude NAs
    sb[is.na(sb)] <- 0
    sb <- sb[which(rowSums(sb[, c("astrocytes", "endothelial", "myeloid", "neuron", "oligodendrocytes")]) != 0),]
    sb <- sb[!duplicated(sb$geneList),]
    # Do the tests -- it is 6 tests per group
    t1 <- wilcox.test(x = sb$astrocytes, y = sb$endothelial, alternative = "two.sided")
    t2 <- wilcox.test(x = sb$astrocytes, y = sb$myeloid, alternative = "two.sided")
    t3 <- wilcox.test(x = sb$astrocytes, y = sb$neuron, alternative = "two.sided")
    t4 <- wilcox.test(x = sb$astrocytes, y = sb$oligodendrocytes, alternative = "two.sided")
    t5 <- wilcox.test(x = sb$endothelial, y = sb$myeloid, alternative = "two.sided")
    t6 <- wilcox.test(x = sb$endothelial, y = sb$neuron, alternative = "two.sided")
    t7 <- wilcox.test(x = sb$endothelial, y = sb$oligodendrocytes, alternative = "two.sided")
    t8 <- wilcox.test(x = sb$myeloid, y = sb$neuron, alternative = "two.sided")
    t9 <- wilcox.test(x = sb$myeloid, y = sb$oligodendrocytes, alternative = "two.sided")
    t10 <- wilcox.test(x = sb$neuron, y = sb$oligodendrocytes, alternative = "two.sided")
    tmp_df <- data.frame(Ast_End=t1$p.value, Ast_Mye=t2$p.value, Ast_Neu=t3$p.value, Ast_Oli=t4$p.value, 
                         End_Mye=t5$p.value, End_Neu=t6$p.value, End_Oli=t7$p.value, 
                         Mye_Neu=t8$p.value, Mye_Oli=t9$p.value, Neu_Oli=t10$p.value)
    within_groups[[counter]] <- tmp_df
    counter <- counter + 1
  }
  within_groups <- rbindlist(within_groups)
  rownames(within_groups) <- rownames(out_means)
  
  # Now between groups
  between_groups <- list()
  counter <- 1
  for (x in 2:length(all_gr)-1){
    for (j in (x+1):length(all_gr)){
      # Take the two datasets to compare
      df1 <- var_gene_mapping[which(var_gene_mapping$ID %in% all_gr[[x]]),]
      df2 <- var_gene_mapping[which(var_gene_mapping$ID %in% all_gr[[j]]),]
      # Clean
      df1[is.na(df1)] <- 0
      df2[is.na(df2)] <- 0
      df1 <- df1[which(rowSums(df1[, c("astrocytes", "endothelial", "myeloid", "neuron", "oligodendrocytes")]) != 0),]
      df2 <- df2[which(rowSums(df2[, c("astrocytes", "endothelial", "myeloid", "neuron", "oligodendrocytes")]) != 0),]
      df1 <- df1[!duplicated(df1$geneList),]
      df2 <- df2[!duplicated(df2$geneList),]
      # Do the test for the 4 functional clusters
      t1 <- wilcox.test(x = df1$astrocytes, y = df2$astrocytes, alternative = "two.sided")
      t2 <- wilcox.test(x = df1$endothelial, y = df2$endothelial, alternative = "two.sided")
      t3 <- wilcox.test(x = df1$myeloid, y = df2$myeloid, alternative = "two.sided")
      t4 <- wilcox.test(x = df1$neuron, y = df2$neuron, alternative = "two.sided")
      t5 <- wilcox.test(x = df1$oligodendrocytes, y = df2$oligodendrocytes, alternative = "two.sided")
      tmp_df <- data.frame(comparison=paste(x, j, sep="_"), Ast=t1$p.value, End=t2$p.value, Mye=t3$p.value, Neu=t4$p.value, Oli=t5$p.value)
      between_groups[[counter]] <- tmp_df
      counter <- counter + 1
    }
  }
  between_groups <- rbindlist(between_groups)
  
  # For the comparison between groups, also do 1 groups vs. the others
  one_vs_all <- list()
  counter <- 1
  for (x in 1:length(all_gr)){
    # Take the two datasets to compare
    df1 <- var_gene_mapping[which(var_gene_mapping$ID %in% all_gr[[x]]),]
    df2 <- var_gene_mapping[which(!(var_gene_mapping$ID %in% all_gr[[x]])),]
    # Clean
    df1[is.na(df1)] <- 0
    df2[is.na(df2)] <- 0
    df1 <- df1[which(rowSums(df1[, c("astrocytes", "endothelial", "myeloid", "neuron", "oligodendrocytes")]) != 0),]
    df2 <- df2[which(rowSums(df2[, c("astrocytes", "endothelial", "myeloid", "neuron", "oligodendrocytes")]) != 0),]
    df1 <- df1[!duplicated(df1$geneList),]
    df2 <- df2[!duplicated(df2$geneList),]
    # Do the test for the 4 functional clusters
    t1 <- wilcox.test(x = df1$astrocytes, y = df2$astrocytes, alternative = "two.sided")
    t2 <- wilcox.test(x = df1$endothelial, y = df2$endothelial, alternative = "two.sided")
    t3 <- wilcox.test(x = df1$myeloid, y = df2$myeloid, alternative = "two.sided")
    t4 <- wilcox.test(x = df1$neuron, y = df2$neuron, alternative = "two.sided")
    t5 <- wilcox.test(x = df1$oligodendrocytes, y = df2$oligodendrocytes, alternative = "two.sided")
    tmp_df <- data.frame(comparison=paste(x, "_vs_all", sep=""), Ast=t1$p.value, End=t2$p.value, Mye=t3$p.value, Neu=t4$p.value, Oli=t5$p.value)
    one_vs_all[[counter]] <- tmp_df
    counter <- counter + 1
  }
  one_vs_all <- rbindlist(one_vs_all)
  
  l=list(out_means, out_sd, between_groups, within_groups, one_vs_all)
  return(l)
}

# Figure to draw figure 1 this needs to be finished
function.Figure1 <- function(var_gene_mapping_CT, angle, lit){
  # read my association file
  myass <- fread("INPUTS/myAssoc_ADsnps_GRACE.txt", h=T)
  
  # match same as the densities
  myass <- matchBy(myass, rownames(angle), myass$ID)
  lit <- matchBy(lit, rownames(angle), lit$rsid)
  
  # check alleles and in case change my beta
  lit$beta_adj <- lit$beta
  lit$beta_adj[which(lit$beta <0)] <- lit$beta_adj[which(lit$beta <0)]*(-1)
  lit$a1_adj <- toupper(lit$a1)
  lit$a1_adj[which(lit$beta <0)] <- toupper(lit$a2[which(lit$beta <0)])
  myass$BETA[which(myass$A1 != lit$a1_adj)] <- myass$BETA[which(myass$A1 != lit$a1_adj)]*(-1)
  
  # color palette
  color.palette <- brewer.pal(9, "Set1")
  
  # manage colors for beta -- normalize between 0 and 1
  min.beta <- min(na.omit(abs(lit$beta_adj), abs(myass$BETA)))
  max.beta <- max(na.omit(abs(lit$beta_adj), abs(myass$BETA)))
  lit$norm.beta.lit <- (abs(lit$beta_adj) - min.beta) / (max.beta - min.beta)
  myass$norm.beta.mine <- (abs(myass$BETA) - min.beta) / (max.beta - min.beta)
  lit$color.lit <- color.palette[1]
  lit$color.lit[which(lit$beta_adj > 0)] <- color.palette[4]
  myass$color.mine <- color.palette[1]
  myass$color.mine[which(myass$BETA < 0)] <- color.palette[4]
  
  # find missings to be added in proportions
  var_gene_mapping_CT[is.na(var_gene_mapping_CT)] <- 0
  
  # normalize expression data between 0 and 1
  max.expre <- max(na.omit(var_gene_mapping_CT[, 20:24]))
  min.expre <- min(na.omit(var_gene_mapping_CT[, 20:24]))
  var_gene_mapping_CT[, 20:24] <- (var_gene_mapping_CT[, 20:24] - min.expre)/(max.expre - min.expre)
  
  # PLOT
  # Global graphical parameters
  par(mar=c(0, 18, 16, 0), mfrow=c(1,1))
  # ymax for the plot
  ymax <- nrow(angle) + 1
  # define square height
  h <- 1
  # define counter for the plot
  c <- 1
  # to automatize, need to use parameters
  xmx <- 37
  st.labels <- 0.90
  step.labels <- 1.50
  l.dens <- 9.5
  polygon.hei <- 3
  
  # prepare window to plot
  plot(1, xlim = c(0, xmx),  ylim=c(0, ymax), xlab='', ylab='', xaxt='none', 
       yaxt='none', col='white', bty='n')
  
  # put text on columns
  text(x = st.labels, y = (ymax + st.labels), labels = "Beta-Amyloid", adj = 0, xpd=T, srt=90, cex=2, font=2)
  text(x = st.labels+step.labels, y = (ymax + st.labels), labels = "Lipid/Cholesterol", adj = 0, xpd=T, srt=90, cex=2, font=2)
  text(x = st.labels+step.labels*2, y = (ymax + st.labels), labels = "Endocytosis/Immune", adj = 0, xpd=T, srt=90, cex=2, font=2)
  text(x = st.labels+step.labels*3, y = (ymax + st.labels), labels = "Synaptic plasticity", adj = 0, xpd=T, srt=90, cex=2, font=2)
  
  text(x = st.labels+step.labels*4, y = (ymax + st.labels), labels = "Effect on AD", adj = 0, xpd=T, srt=90, cex=2, font=2)
  text(x = st.labels+step.labels*5, y = (ymax + st.labels), labels = "Effect on LGV", adj = 0, xpd=T, srt=90, cex=2, font=2)
  
  text(x = -4.75, y = (ymax + st.labels*1.5), labels = "Variant annotation", xpd=T, cex=2, font=2)
  
  # try to put only two polygons for the classes here
  # polygon(x = c(st.labels+step.labels*6, st.labels+step.labels*6, st.labels+step.labels*6+l.dens), y = c((ymax + st.labels), (ymax + polygon.hei), (ymax + st.labels)), col = alpha(color.palette[5], 0.40), xpd=T)
  # polygon(x = c(st.labels+step.labels*6+l.dens, st.labels+step.labels*6+l.dens*2, st.labels+step.labels*6+l.dens*2), y = c((ymax + st.labels), (ymax + polygon.hei), (ymax + st.labels)), col = alpha(color.palette[2], 0.40), xpd=T)
  # polygon(x = c(st.labels+step.labels*6, st.labels+step.labels*6+l.dens, st.labels+step.labels*6+l.dens*2), y = c((ymax + polygon.hei+0.15), (ymax + st.labels+0.15), (ymax + polygon.hei+0.15)), col = alpha(color.palette[6], 0.40), xpd=T)
  # text(x = st.labels+step.labels*6+0.5, y = (ymax + 1.50), labels = "CHA effect", adj=0, xpd=T, cex = 1.40, font=2)
  # text(x = st.labels+step.labels*6+l.dens*2-0.5, y = (ymax + 1.50), labels = "! CHA effect", adj=1, xpd=T, font=2, cex = 1.40)
  # text(x = st.labels+step.labels*6+l.dens, y = (ymax + 2.50), labels = "AD effect", adj=0.5, xpd=T, font=2, cex = 1.40)
  
  # add labels for expected and unexpected
  colors.unexp <- brewer.pal(n = 3, name = "Pastel1")
  # rect(xleft = st.labels+step.labels*6, ybottom = ymax + 3.30, xright = st.labels+step.labels*6+l.dens, ytop = ymax + 4.30, col = alpha(colors.unexp[1], 0.6), xpd=T)
  # rect(xleft = st.labels+step.labels*6+l.dens, ybottom = ymax + 3.30, xright = st.labels+step.labels*6+l.dens*2, ytop = ymax + 4.30, col = alpha(colors.unexp[2], 0.6), xpd=T)
  # text(x = st.labels+step.labels*6+l.dens/2, y = ymax + 3.8, labels = "Expected", adj = 0.5, xpd=T, font=2, cex=1.40)
  # text(x = st.labels+step.labels*6+(l.dens*3/2), y = ymax + 3.8, labels = "Unexpected", adj = 0.5, xpd=T, font=2, cex=1.40)
  # rect(xleft = st.labels+step.labels*6, ybottom = ymax + st.labels, xright = st.labels+step.labels*6+l.dens, ytop = ymax + st.labels*2, col = alpha(colors.unexp[1], 0.6), xpd=T)
  # rect(xleft = st.labels+step.labels*6+l.dens, ybottom = ymax + st.labels, xright = st.labels+step.labels*6+l.dens*2, ytop = ymax + st.labels*2, col = alpha(colors.unexp[2], 0.6), xpd=T)
  text(x = st.labels+step.labels*6+l.dens/2, y = ymax + st.labels*1.5, labels = "Expected", adj = 0.5, xpd=T, font=2, cex=2)
  text(x = st.labels+step.labels*6+(l.dens*3/2), y = ymax + st.labels*1.5, labels = "Unexpected", adj = 0.5, xpd=T, font=2, cex=2)
  
  # text about cell type
  text(x = st.labels+step.labels*7+l.dens*2, y = (ymax + st.labels), labels = "Astrocytes", adj = 0, xpd=T, srt=90, cex=2, font=2)
  text(x = st.labels+step.labels*8+l.dens*2, y = (ymax + st.labels), labels = "Oligodendrocytes", adj = 0, xpd=T, srt=90, cex=2, font=2)
  text(x = st.labels+step.labels*9+l.dens*2, y = (ymax + st.labels), labels = "Microglia", adj = 0, xpd=T, srt=90, cex=2, font=2)
  text(x = st.labels+step.labels*10+l.dens*2, y = (ymax + st.labels), labels = "Endotelial", adj = 0, xpd=T, srt=90, cex=2, font=2)
  text(x = st.labels+step.labels*11+l.dens*2, y = (ymax + st.labels), labels = "Neurons", adj = 0, xpd=T, srt=90, cex=2, font=2)
  
  # define positions
  p1 <- (st.labels+step.labels*3)/2 + st.labels/2
  p2 <- st.labels+step.labels*4.5
  p3 <- st.labels+step.labels*7+l.dens*2 + ((st.labels+step.labels*11+l.dens*2)-(st.labels+step.labels*7+l.dens*2))/2
  text(x = -4.75, y = (ymax - c + h), labels = 'A', font=2, xpd=T, cex = 2)
  text(x = p1, y = (ymax - c + h), labels = 'B', font=2, cex = 2)
  text(x = p2, y = (ymax - c + h), labels = 'C', font=2, cex = 2)
  text(x = st.labels+step.labels*6+l.dens, y = (ymax - c + h), labels = 'D', font=2, cex = 2)
  text(x = p3, y = (ymax - c + h), labels = 'E', font=2, cex = 2)
  segments(x0 = -8, y0 = (ymax - c), x1 = st.labels+step.labels*11+l.dens*2, y1 = (ymax - c), xpd=T, cex = 2)
  c <- c + h
  
  # select square size
  sz <- 1.5
  
  #draw line to "divide" figure 
  segments(x0 = 0, y0 = 0, x1 = 0, y1 = ymax, xpd=T, lwd=3)
  segments(x0 = sz*4, y0 = 0, x1 = sz*4, y1 = ymax, xpd=T, lwd=3)
  segments(x0 = sz*6, y0 = 0, x1 = sz*6, y1 = ymax, xpd=T, lwd=3)
  segments(x0 = st.labels+step.labels*6.5+l.dens*2, y0 = 0, x1 = st.labels+step.labels*6.5+l.dens*2, y1 = ymax, xpd=T, lwd=3)
  
  # last thing: big rectangles for clustering
  # rect(xleft = -8.5, ybottom = 10, xright = st.labels+step.labels*11.5+l.dens*2, ytop = 0, border="navy", lwd=6, xpd=T)
  # rect(xleft = -8.5, ybottom = 27, xright = st.labels+step.labels*11.5+l.dens*2, ytop = 10, border="orange", lwd=6, xpd=T)
  # rect(xleft = -8.5, ybottom = 27, xright = st.labels+step.labels*11.5+l.dens*2, ytop = 38, border="red", lwd=6, xpd=T)
  
  # also the labels of the groups
  # text(x = st.labels+step.labels*11.5+l.dens*2+0.5, y = 32.5, labels = "CHA-group", col="red", xpd=T, cex=1.40, font=4, srt=-90)
  # text(x = st.labels+step.labels*11.5+l.dens*2+0.5, y = 18.5, labels = "AD-group", col="orange", xpd=T, cex=1.40, font=4, srt=-90)
  # text(x = st.labels+step.labels*11.5+l.dens*2+0.5, y = 5, labels = "Unex-group", col="navy", xpd=T, cex=1.40, font=4, srt=-90)
  text(x = -10, y = 32.5, labels = "Longevity-group", col="red", xpd=T, cex=1.80, font=4, srt=-90)
  text(x = -11, y = 32.5, labels = expression(bolditalic(E[LGV]*" > "*E[AD])), col="red", xpd=T, cex=1.40, font=4, srt=-90)
  text(x = -10, y = 18.5, labels = "AD-group", col="orange", xpd=T, cex=1.80, font=4, srt=-90)
  text(x = -11, y = 18.5, labels = expression(bolditalic(E[LGV]*" < "*E[AD])), col="orange", xpd=T, cex=1.40, font=4, srt=-90)
  text(x = -10, y = 5, labels = "Unex-group", col="navy", xpd=T, cex=1.80, font=4, srt=-90)
  text(x = -11, y = 5, labels = "Unexpected direction", col="navy", xpd=T, cex=1.40, font=4, srt=-90)
  
  # main loop across loci
  for (loc in 1:nrow(var_gene_mapping_CT)){
    # select locus
    rsid_tmp <- var_gene_mapping_CT$ID[loc]
    # take information for the variant
    path <- var_gene_mapping_CT[which(var_gene_mapping_CT$ID == rsid_tmp), ]
    # take variant density and put x and y in dataframe
    locus <- density(angle[which(rownames(angle) == rsid_tmp), ])
    # also find median
    med_dens <- median(angle[which(rownames(angle) == rsid_tmp), ])
    dens <- data.frame(x = locus$x, y = locus$y)
    # normalize x-densities for the place where to plot [a, b]
    # in case of APOE as the distribution are cutted, I need to extend them
    if (rsid_tmp %in% c("rs429358", "rs7412")){
      tmp.up <- data.frame(x=0, y=0)
      tmp.down <- data.frame(x=2, y=0)
      dens <- rbind(tmp.up, dens, tmp.down)
    }
    b <- st.labels+step.labels*6+l.dens*2
    a <- st.labels+step.labels*6
    dens$x <- (b-a) * ((dens$x - min(dens$x))/(max(dens$x) - min(dens$x))) + a
    # also normalize median in the same interval
    med_dens <- (b-a) * ((med_dens - min(angle))/(max(angle) - min(angle))) + a
    # also check the y and in case normalize between 0 and 0.4
    if (max(dens$y) >2){ dens$y <- 1.5*(dens$y - min(dens$y))/(max(dens$y) - min(dens$y)) }
    lim.dens <- a
    lim.ad <- a + l.dens/2
    lim.unex <- lim.ad + l.dens/2
    lim.ext <- b
    dens$y <- dens$y + (ymax - c)
    dens <- subset(dens, ((dens$x > lim.dens) & (dens$x <= lim.ext)))
    
    # plot functional annotations
    if (sum(path[,16:19]) == 0){
      text(x = sz/2, y = (ymax - c + h/2), labels="X", font=2, col="red", cex = 1.50)
      text(x = sz/2 + sz, y = (ymax - c + h/2), labels="X", font=2, col="red", cex = 1.50)
      text(x = sz/2 + sz*2, y = (ymax - c + h/2), labels="X", font=2, col="red", cex = 1.50)
      text(x = sz/2 + sz*3, y = (ymax - c + h/2), labels="X", font=2, col="red", cex = 1.50)
      #text(x = sz/2 + sz*4, y = (ymax - c + h/2), labels="X", font=2, col="red")
      #text(x = sz/2 + sz*5, y = (ymax - c + h/2), labels="X", font=2, col="red")
      #text(x = sz/2 + sz*6, y = (ymax - c + h/2), labels="X", font=2, col="red")
    }
    rect(xleft = 0, ybottom = (ymax - c), xright = sz, ytop = (ymax - c + h), col=ggplot2::alpha("black", path$"Beta-Amyloid"), border = 'black')
    rect(xleft = sz, ybottom = (ymax - c), xright = sz*2, ytop = (ymax - c + h), col=ggplot2::alpha("black", path$"Lipid metabolism"), border = 'black')
    rect(xleft = sz*2, ybottom = (ymax - c), xright = sz*3, ytop = (ymax - c + h), col=ggplot2::alpha("black", path$"Immune activation"), border = 'black')
    rect(xleft = sz*3, ybottom = (ymax - c), xright = sz*4, ytop = (ymax - c + h), col=ggplot2::alpha("black", path$"Synaptic plasticity"), border = 'black')

    # plot effect on AD
    sb.lit <- lit[which(lit$rsid == rsid_tmp),]
    rect(xleft = sz*4, ybottom = (ymax - c), xright = sz*5, ytop = (ymax - c + h), col=ggplot2::alpha(sb.lit$color.lit, min(sb.lit$norm.beta.lit+0.5, 1)), border = 'black')
    
    # plot effect on survival
    sb.my <- myass[which(myass$ID == rsid_tmp),]
    rect(xleft = sz*5, ybottom = (ymax - c), xright = sz*6, ytop = (ymax - c + h), col=ggplot2::alpha(sb.my$color.mine, min(abs(sb.my$norm.beta.mine)+0.5, 1)), border = 'black')
    # also add * in case association was significant
    if (sb.my$P <= 0.05){ text(x =sz*5.5, y = (ymax - c + h/2), labels="*", cex=2, font=2) }
    
    # plot annotation
    # text(x = -11.25, y = (ymax - c + (h/2)), labels = path$locus, font=2, adj=0, xpd=T, cex=1)
    # text(x = -5.25, y = (ymax - c + (h/2)), labels = rsid_tmp, font=2, adj=0.5, xpd=T, cex=1)
    # text(x = -0.25, y = (ymax - c + (h/2)), labels = sb.lit$gene, font=4, adj=1, xpd=T, cex=1)
    text(x = -9.25, y = (ymax - c + (h/2)), labels = path$locus, font=2, adj=0, xpd=T, cex=0.95)
    text(x = -4.5, y = (ymax - c + (h/2)), labels = rsid_tmp, font=2, adj=0.5, xpd=T, cex=0.95)
    text(x = -0.25, y = (ymax - c + (h/2)), labels = sb.lit$gene, font=4, adj=1, xpd=T, cex=1.15)
    
    if (rowSums(path[, 20:24]) == 0){
      text(x=st.labels+step.labels*7+l.dens*2, y = (ymax-c + h/2), labels="X", col="red", font=2, cex = 1.50)
      text(x=st.labels+step.labels*8+l.dens*2, y = (ymax-c + h/2), labels="X", col="red", font=2, cex = 1.50)
      text(x=st.labels+step.labels*9+l.dens*2, y = (ymax-c + h/2), labels="X", col="red", font=2, cex = 1.50)
      text(x=st.labels+step.labels*10+l.dens*2, y = (ymax-c + h/2), labels="X", col="red", font=2, cex = 1.50)
      text(x=st.labels+step.labels*11+l.dens*2, y = (ymax-c + h/2), labels="X", col="red", font=2, cex = 1.50)
    }
    # cell types
    cell.col <- "navy"
    
    # astrocytes
    rect(xleft = st.labels+step.labels*6.5+l.dens*2, ybottom = (ymax-c), xright = st.labels+step.labels*7.5+l.dens*2, ytop = (ymax - c + h), 
         col=ggplot2::alpha(cell.col, path$astrocytes), border = 'black')
    #oligodendrocytes
    rect(xleft = st.labels+step.labels*7.5+l.dens*2, ybottom = (ymax-c), xright = st.labels+step.labels*8.5+l.dens*2, ytop = (ymax - c + h), 
         col=ggplot2::alpha(cell.col, path$oligodendrocytes), border = 'black')
    #microglia
    rect(xleft = st.labels+step.labels*8.5+l.dens*2, ybottom = (ymax-c), xright = st.labels+step.labels*9.5+l.dens*2, ytop = (ymax - c + h), 
         col=ggplot2::alpha(cell.col, path$myeloid), border = 'black')
    #endotelial cells
    rect(xleft = st.labels+step.labels*9.5+l.dens*2, ybottom = (ymax-c), xright = st.labels+step.labels*10.5+l.dens*2, ytop = (ymax - c + h), 
         col=ggplot2::alpha(cell.col, path$endothelial), border = 'black')
    #neurons
    rect(xleft = st.labels+step.labels*10.5+l.dens*2, ybottom = (ymax-c), xright = st.labels+step.labels*11.5+l.dens*2, ytop = (ymax - c + h), 
         col=ggplot2::alpha(cell.col, path$neuron), border = 'black')
    
    if (loc == 1){
      rect(xleft = -9.5, ybottom = 38, xright = st.labels+step.labels*11.5+l.dens*2, ytop = 27.1, border="red", lwd=5, xpd=T)
    } else if (loc == 12){
      # last thing: big rectangles for clustering
      rect(xleft = -9.5, ybottom = 10.1, xright = st.labels+step.labels*11.5+l.dens*2, ytop = 26.9, border="orange", lwd=5, xpd=T)
      segments(x0 = -9.5, y0 = 27.1, x1 = st.labels+step.labels*11.5+l.dens*2, y1 = 27.1, col="red", lwd=5, xpd=T)
      segments(x0 = st.labels+step.labels*11.5+l.dens*2, y0 = 38, x1 = st.labels+step.labels*11.5+l.dens*2, y1 = 27.1, col="red", lwd=5, xpd=T)
    } else if (loc == nrow(var_gene_mapping_CT)){
      # last thing: big rectangles for clustering
      segments(x0 = -9.5, y0 = 10.1, x1 = st.labels+step.labels*11.5+l.dens*2, y1 = 10.1, col="orange", lwd=5, xpd=T)
      segments(x0 = st.labels+step.labels*11.5+l.dens*2, y0 = 26.9, x1 = st.labels+step.labels*11.5+l.dens*2, y1 = 10.1, col="orange", lwd=5, xpd=T)
      rect(xleft = -9.5, ybottom = 9.9, xright = st.labels+step.labels*11.5+l.dens*2, ytop = 0, border="navy", lwd=5, xpd=T)
    }
    
    # plot unexpected
    polygon(dens, col=color.palette[5], border = NA)
    # add median -- need to find the y where the median intesect the density
    # calculate distance from the median and take the closest value
    tmp <- dens
    tmp$median_dist <- abs(med_dens - tmp$x)
    y_inters <- tmp$y[which(tmp$median_dist == min(tmp$median_dist))]
    
    # plot larger effect on AD
    ad <- subset(dens, dens$x > lim.ad)
    supp <- data.frame(x = lim.ad, y = (ymax - c))
    ad.supp <- rbind(supp, ad)
    polygon(ad.supp, col=color.palette[6], border = NA)
    
    # plot larger effect on survival
    unex <- subset(dens, dens$x > lim.unex)
    supp <- data.frame(x = lim.unex, y = (ymax - c))
    unex.supp <- rbind(supp, unex)
    polygon(unex.supp, col=color.palette[2], border = NA)
    
    # finally add medians
    segments(x0 = med_dens, y0 = min(dens$y), x1 = med_dens, y1 = y_inters, lwd=1, col="navy")
    
    # re-add border of the polygon
    polygon(dens, col = NULL, border = "black")
    
    #increment c
    c <- c + h
    
    #draw line to divide loci
    segments(x0 = -9.25, y0 = (ymax - c + h), x1 = st.labels+step.labels*11+l.dens*2, y1 = (ymax - c + h), xpd=T)
  }
  
  #lines
  b <- st.labels+step.labels*6+l.dens*2
  a <- st.labels+step.labels*6
  #segments(x0 = a+l.dens/2, y0 = 0, x1 = a+l.dens/2, y1 = ymax-1, lty = 2)
  segments(x0 = a+l.dens, y0 = 0, x1 = a+l.dens, y1 = ymax-1, lty = 1, lwd=1.5)
  
  # also legend at the very bottom
  text(x = -8.5, y = -0.5, labels = "LGV: longevity        *: p-value<0.05", font=2, cex=1, xpd=T, adj=0)
  segments(x0 = 3.4, y0 = -0.5, x1 = 4, y1 = -0.5, lwd = 2, col="navy", xpd=T)
  text(x = 5, y = -0.5, labels = "Median", cex=1, xpd=T, font=2)
}

# Function to extract expression data from publicly available single-cell dataset
cellType_expression <- function(var_gene_mapping){
  # First, read data
  expr.data <- fread("INPUTS/RNAexpr_fromBrain.tab", h=T, sep='\t', dec=',')
  colnames(expr.data) <- c("Gene.Symbol", "astrocytes", "endothelial", "myeloid", "neuron", "oligodendrocytes", "diff.analysis.adjP",
                           "diff.analysis.log2FC")
  
  # Derive list of genes associated with all variants
  geneList <- unlist(strsplit(x = var_gene_mapping$geneList, ","))
  
  # Match gene expression of genes of interest
  gene_expr_all <- expr.data[which(expr.data$Gene.Symbol %in% geneList),]
  gene_expr_all <- gene_expr_all[!duplicated(gene_expr_all$Gene.Symbol),]
  miss <- geneList[which(!(geneList %in% gene_expr_all$Gene.Symbol))]
  pie(x = c(nrow(gene_expr_all), length(miss)), labels = c("Match", "Miss"))

  # Plot heatmap of cell-type expression
  pheatmap(gene_expr_all[, c(2:6)], labels_row = gene_expr_all$Gene.Symbol)
  
  # Create empty columns regarding cell-type expression
  for (x in names(gene_expr_all)[2:6]){ var_gene_mapping[, x] <- NA }
  
  # Now need to average by the genes, per variant
  for (i in 1:nrow(var_gene_mapping)){
    # select genes of the variant of interest
    gene_sb <- unlist(strsplit(var_gene_mapping$geneList[i], ","))
    
    # take relative cell-type expression
    cell_t <- gene_expr_all[which(gene_expr_all$Gene.Symbol %in% gene_sb),]
    
    # do the mean and assign if there are entries
    if (nrow(cell_t) >=1){
      var_gene_mapping[i, c("astrocytes", "endothelial", "myeloid", "neuron", "oligodendrocytes")] <- c(mean(cell_t$astrocytes), mean(cell_t$endothelial), mean(cell_t$myeloid), mean(cell_t$neuron), mean(cell_t$oligodendrocytes))
      # also normalize between 0 and 1
      var_gene_mapping[i, c("astrocytes", "endothelial", "myeloid", "neuron", "oligodendrocytes")] <- var_gene_mapping[i, c("astrocytes", "endothelial", "myeloid", "neuron", "oligodendrocytes")] / sum(var_gene_mapping[i, c("astrocytes", "endothelial", "myeloid", "neuron", "oligodendrocytes")])
    } else {
      var_gene_mapping[i, c("astrocytes", "endothelial", "myeloid", "neuron", "oligodendrocytes")] <- NA
    }
  }
  
  # Plot heatmap of cell-type at variant level
  pheatmap(var_gene_mapping[, c("astrocytes", "endothelial", "myeloid", "neuron", "oligodendrocytes")], cluster_rows = F, labels_row = var_gene_mapping$locus)
  
  return(var_gene_mapping)
}

# Function to correct pvalues given a matrix -- return the same matrix with corrected pvalues
correctMe <- function(df, skipFirst=F){
  if (skipFirst == FALSE){
    # put all values of the matrix into a vector
    pv <- as.vector(t(df))
    # correct pvalues
    pv_adj <- p.adjust(p = pv, method = "fdr", n = length(pv))
    # restore dataframe
    newdf <- as.data.frame(matrix(data=pv_adj, nrow = nrow(df), ncol = ncol(df), byrow = T))
    rownames(newdf) <- rownames(df)
    colnames(newdf) <- colnames(df)
  } else {
    # do the same but skip first column
    toExc <- df[, 1]
    # put all values of the matrix into a vector
    pv <- as.vector(t(df[, 2:ncol(df)]))
    # correct pvalues
    pv_adj <- p.adjust(p = pv, method = "fdr", n = length(pv))
    # restore dataframe
    newdf <- as.data.frame(matrix(data=pv_adj, nrow = nrow(df), ncol = (ncol(df)-1), byrow = T))
    rownames(newdf) <- rownames(df)
    newdf <- cbind(toExc, newdf)
    colnames(newdf) <- colnames(df)
  }
  return(newdf)
}

# Function to draw forest plot with GWAS results
InDaForest <- function(var_gene_mapping, orig_lit){
  # load working framework
  load("/Users/nicco/Desktop/2k19_work/rotationPaper/FINAL_VERSION/RESULTS/analysis_GWAS_context.RData")
  # important: flip all variants to the AD-risk increasing allele (so that also direction can be seen)
  for (i in 1:nrow(all_long)){
    # grep variant association on ad
    ad_as <- orig_lit[which(orig_lit$rsid == all_long$rsid_JOR[i]),]
    if (ad_as$beta <0){
      # if beta <0 --> protective, make sure the allele tested on longevity is the other
      oth_all <- toupper(ad_as$a2)
      if (all_long$A1[i] != oth_all){
        all_long$A1[i] <- oth_all
        all_long$BETA[i] <- all_long$BETA[i]*(-1)
      }
    } else if (ad_as$beta >0) {
      all_risk <- toupper(ad_as$a1)
      if (all_long$A1[i] != all_risk){
        all_long$A1[i] <- all_risk
        all_long$BETA[i] <- all_long$BETA[i]*(-1)
      }
    }  
  }
  # use the same variants
  background <- var_gene_mapping[which(var_gene_mapping$ID %in% all_long$rsid_JOR),]
  # also need to add standard errors for my association
  se_sb <- res[, c("LOCUS", "SE", "P", "P_ADJ")]
  all_long <- merge(all_long, se_sb, by.x="locus_JOR", by.y="LOCUS")
  # define max
  max <- 34.5
  # first make sure all alleles and betas are the same
  for (i in 1:nrow(all_long)){
    # make sure to have same alleles
    my_beta <- all_long[i, "BETA"]
    my_se <- all_long[i, "SE"]
    my_a1 <- all_long[i, "A1"]
    if (my_a1 != toupper(all_long[i, "a1_JOR"])){ all_long[i, "beta_JOR"] <- all_long[i, "beta_JOR"]*(-1) }
    if (my_a1 != toupper(all_long[i, "a1_TIM"])){ all_long[i, "beta_TIM"] <- all_long[i, "beta_TIM"]*(-1) }
  }
  # first, calculate CIs
  all_long$upCI <- (all_long$BETA+(1.96*all_long$SE))
  all_long$lowCI <- (all_long$BETA-(1.96*all_long$SE))
  all_long$upCI_jor <- (all_long$beta_JOR+(1.96*all_long$se_JOR))
  all_long$lowCI_jor <- (all_long$beta_JOR-(1.96*all_long$se_JOR))
  all_long$upCI_tim <- (all_long$beta_TIM+(1.96*all_long$se_TIM))
  all_long$lowCI_tim <- (all_long$beta_TIM-(1.96*all_long$se_TIM))
  # scale betas and ci with respect to where points should be plotted
  all_long[, c("upCI_norm", "lowCI_norm", "BETA_norm")] <- 4*((all_long[, c("upCI", "lowCI", "BETA")] - min(all_long[, c("upCI", "lowCI", "BETA")]))/(max(all_long[, c("upCI", "lowCI", "BETA")])-min(all_long[, c("upCI", "lowCI", "BETA")])))+2
  all_long[, c("upCI_jor_norm", "lowCI_jor_norm", "beta_JOR_norm")] <- 4*((all_long[, c("upCI_jor", "lowCI_jor", "beta_JOR")] - min(all_long[, c("upCI_jor", "lowCI_jor", "beta_JOR")]))/(max(all_long[, c("upCI_jor", "lowCI_jor", "beta_JOR")])-min(all_long[, c("upCI_jor", "lowCI_jor", "beta_JOR")])))+7
  all_long[, c("upCI_tim_norm", "lowCI_tim_norm", "beta_TIM_norm")] <- 4*((all_long[, c("upCI_tim", "lowCI_tim", "beta_TIM")] - min(all_long[, c("upCI_tim", "lowCI_tim", "beta_TIM")]))/(max(all_long[, c("upCI_tim", "lowCI_tim", "beta_TIM")])-min(all_long[, c("upCI_tim", "lowCI_tim", "beta_TIM")])))+12
  # i also need the 0 position
  pos_0_chc <- 4*((0 - min(all_long[, c("upCI", "lowCI", "BETA")]))/(max(all_long[, c("upCI", "lowCI", "BETA")])-min(all_long[, c("upCI", "lowCI", "BETA")])))+2
  pos_0_jor <- 4*((0 - min(all_long[, c("upCI_jor", "lowCI_jor", "beta_JOR")]))/(max(all_long[, c("upCI_jor", "lowCI_jor", "beta_JOR")])-min(all_long[, c("upCI_jor", "lowCI_jor", "beta_JOR")])))+7
  pos_0_tim <- 4*((0 - min(all_long[, c("upCI_tim", "lowCI_tim", "beta_TIM")]))/(max(all_long[, c("upCI_tim", "lowCI_tim", "beta_TIM")])-min(all_long[, c("upCI_tim", "lowCI_tim", "beta_TIM")])))+12

  # graphical params
  par(mar=c(5, 4, 4, 10))
  # make background plot
  plot(0, 0, pch=16, col="white", bty='none', xaxt='none', yaxt='none', xlab="Effect size on longevity of AD risk-increasing alleles", ylab="", xlim=c(0, 16), ylim=c(0, 34), cex.lab=2)
  # add grid for each rectangle
  for (i in seq(0, 34, 3.4)){ segments(x0 = 2, y0 = i, x1 = 6, y1 = i, lwd=0.5, col="grey80"); segments(x0 = 7, y0 = i, x1 = 11, y1 = i, lwd=0.5, col="grey60"); segments(x0 = 12, y0 = i, x1 = 16, y1 = i, lwd=0.5, col="grey60") }
  for (i in seq(2, 6)){ segments(x0 = i, y0 = 0, x1 = i, y1 = 34, lwd=0.5, col="grey80") }
  for (i in seq(7, 11)){ segments(x0 = i, y0 = 0, x1 = i, y1 = 34, lwd=0.5, col="grey80") }
  for (i in seq(12, 16)){ segments(x0 = i, y0 = 0, x1 = i, y1 = 34, lwd=0.5, col="grey80") }
  # there will be 3 rectangles, each divided in 3
  rect(xleft = 2, ybottom = 0, xright = 6, ytop = 34, lwd=3, border="black")
  rect(xleft = 7, ybottom = 0, xright = 11, ytop = 34, lwd=3, border="black")
  rect(xleft = 12, ybottom = 0, xright = 16, ytop = 34, lwd=3, border="black")
  # add names on top of the rectangles
  text(x = 4, y = 35, labels = "100-plus", cex=1.50, font=2, xpd=T, col="black")
  text(x = 9, y = 35, labels = "Deelen et al.", cex=1.50, font=2, xpd=T, col="black")
  text(x = 14, y = 35, labels = "Timmers et al.", cex=1.50, font=2, xpd=T, col="black")
  # add rectangles for the groups
  rect(xleft = 2, ybottom = 0, xright = 6, ytop = 8, lwd=3, border="navy")
  rect(xleft = 7, ybottom = 0, xright = 11, ytop = 8, lwd=3, border="navy")
  rect(xleft = 12, ybottom = 0, xright = 16, ytop = 8, lwd=3, border="navy")
  rect(xleft = 2, ybottom = 8, xright = 6, ytop = 25, lwd=3, border="orange")
  rect(xleft = 7, ybottom = 8, xright = 11, ytop = 25, lwd=3, border="orange")
  rect(xleft = 12, ybottom = 8, xright = 16, ytop = 25, lwd=3, border="orange")
  rect(xleft = 2, ybottom = 34, xright = 6, ytop = 25, lwd=3, border="red")
  rect(xleft = 7, ybottom = 34, xright = 11, ytop = 25, lwd=3, border="red")
  rect(xleft = 12, ybottom = 34, xright = 16, ytop = 25, lwd=3, border="red")
  # also add main line in the middle
  pos_0_list <- list(pos_0_chc, pos_0_jor, pos_0_tim)
  for (i in 1:length(pos_0_list)){ 
    segments(x0 = pos_0_list[[i]], y0 = 0, x1 = pos_0_list[[i]], y1 = 34, lwd=1, col="black")
    text(x = pos_0_list[[i]], y = -1, labels = 0, xpd=T)
    if (i == 1){
      text(x = 2, y = -1, labels = round(min(all_long[, c("upCI", "lowCI")]), 2), xpd=T)
      text(x = 6, y = -1, labels = round(max(all_long[, c("upCI", "lowCI")]), 2), xpd=T)
    } else if (i == 2){
      text(x = 7, y = -1, labels = round(min(all_long[, c("upCI_jor", "lowCI_jor")]), 2), xpd=T)
      text(x = 11, y = -1, labels = round(max(all_long[, c("upCI_jor", "lowCI_jor")]), 2), xpd=T)
    } else {
      text(x = 12, y = -1, labels = round(min(all_long[, c("upCI_tim", "lowCI_tim")]), 2), xpd=T)
      text(x = 16, y = -1, labels = round(max(all_long[, c("upCI_tim", "lowCI_tim")]), 2), xpd=T)
    }
  }
  # main loop per variant
  for (i in 1:nrow(background)){
    # get line of interest
    sb <- all_long[which(all_long$rsid_JOR == background$ID[i]),]
    # now draw points
    # 100-plus
    if (max-i > 25){ ccol <- "red" } else if (max-i > 8) { ccol <- "orange" } else { ccol <- "navy" } 
    segments(x0 = sb$lowCI_norm, y0 = max-i, x1 = sb$upCI_norm, y1 = max-i, lwd=1.25, col=ccol)
    points(x = sb$BETA_norm, y = max-i, pch=16, col=ccol)
    if (sb$P_ADJ <= 0.05){
      text(x = 5.9, y = max-i, labels = "**", font=2, col=ccol, adj=1)
    } else if (sb$P <= 0.05){
      text(x = 5.9, y = max-i, labels = "*", font=2, col=ccol, adj=1)
    }
    # deelen
    segments(x0 = sb$lowCI_jor_norm, y0 = max-i, x1 = sb$upCI_jor_norm, y1 = max-i, lwd=1.25, col=ccol)
    points(x = sb$beta_JOR_norm, y = max-i, pch=17, col=ccol)
    if (sb$p_adj_JOR <= 0.05){
      text(x = 10.9, y = max-i, labels = "**", font=2, col=ccol, adj=1)
    } else if (sb$p_JOR <= 0.05){
      text(x = 10.9, y = max-i, labels = "*", font=2, col=ccol, adj=1)
    }
    # timmers
    segments(x0 = sb$lowCI_tim_norm, y0 = max-i, x1 = sb$upCI_tim_norm, y1 = max-i, lwd=1.25, col=ccol)
    points(x = sb$beta_TIM_norm, y = max-i, pch=15, col=ccol)
    if (sb$p_adj_TIM <= 0.05){
      text(x = 15.9, y = max-i, labels = "**", font=2, col=ccol, adj=1)
    } else if (sb$p_TIM <= 0.05){
      text(x = 15.9, y = max-i, labels = "*", font=2, col=ccol, adj=1)
    }
    # name on the side
    text(x = 0.5, y = max-i, labels = paste(sb$locus_JOR, "~"), adj = 1, xpd=T)
    text(x = 1.80, y = max-i, labels = orig_lit$gene[which(orig_lit$rsid == sb$rsid_JOR)], adj = 1, xpd=T, font=3)
  }
  # final -- groups annotation
  text(x = 16.5, y = 30.5, labels = "CHA-group", cex=1.25, col="red", srt=90, font=4, xpd=T)
  text(x = 16.5, y = 16.5, labels = "AD-group", cex=1.25, col="orange", srt=90, font=4, xpd=T)
  text(x = 16.5, y = 4, labels = "Unex-group", cex=1.25, col="navy", srt=90, font=4, xpd=T)
  # legend
  legend(x = 17, y = 34, legend = c("*  p<0.05", "**  FDR<5%"), xpd=T, col="black", pch = c("", ""), pt.cex = 1.50, bty='n')
}

# Function to make figure 2 -- comparison of weights for functional annotation and cell-type
function_figure2 <- function(var_gene_mapping_CT, final_within_groups, final_within_groups_CT, all_gr){
  # set colors
  colz_fun <- brewer.pal(n = 4, name = "Accent")
  colz_ct <- brewer.pal(n = 5, name = "Dark2")
  
  # graphical parameters
  par(mfrow=c(2, 1), mar=c(4, 5, 4, 5))
  
  # plot 1 -- functional clusters
  # empty background
  plot(0, 0, pch=16, col="white", xlim=c(0.5, 14.5), ylim=c(0, 1), ylab="Functional clusters weight", xaxt='none', cex.lab=1.50, cex.axis=1.25, bty='n', xlab="")
  
  # grid
  for (i in seq(0.5, 14.5, (14.5-0.5)/10)){ abline(v=i, lwd=0.4, col="grey60") }
  for (i in seq(0, 1, 0.1)){ abline(h=i, lwd=0.4, col="grey60") }
  
  # parameters
  k <- 1
  w <- 0.4
  stripes <- c(0, 10, 25)  
  
  # main loop
  for (gr in 1:length(all_gr)){
    # get variants and compute means across groups
    sb <- var_gene_mapping_CT[which(var_gene_mapping_CT$ID %in% all_gr[[gr]]), 16:19]
    sb <- na.omit(sb)
    mean_fun <- colMeans(sb)
    # then plot rectangles
    for (i in 1:length(mean_fun)){
      rect(xleft = k-w, ybottom = 0, xright = k+w, ytop = mean_fun[i], col = colz_fun[i], lwd=1.5)
      rect(xleft = k-w, ybottom = 0, xright = k+w, ytop = mean_fun[i], col = NA, lwd=1.5, density = stripes[gr])
      k <- k + 1
    }
    k <- k + 1
  }
  text(x = c(2.5, 7.5, 12.5), y = -0.075, labels = c("Longevity-group", "AD-group", "Unex-group"), xpd=T, cex=1.50, font=4)
  legend(x = 2.5, y = 1.25, legend = c("Beta-Amyloid", "Lipid/Cholesterol", "Endocytosis/Immune", "Synaptic plasticity"), bty='n', cex=1.25, col="black", pch=22, ncol=2, pt.lwd = 1.5, pt.bg = colz_fun, xpd=T)  
  text(x = 15.5, y = 1.2, labels = "A", cex=2, font=2, xpd=T)
  
  # add significance
  # cha
  segments(x0 = 3, y0 = 0.95, x1 = 1, y1 = 0.95, lwd=1.25); text(x = 2, y = 0.975, labels = "**", font = 2, cex = 1.25)
  segments(x0 = 3, y0 = 0.90, x1 = 2, y1 = 0.90, lwd=1.25); text(x = 2.5, y = 0.925, labels = "**", font = 2, cex = 1.25)
  segments(x0 = 3, y0 = 0.85, x1 = 4, y1 = 0.85, lwd=1.25); text(x = 3.5, y = 0.875, labels = "**", font = 2, cex = 1.25)
  # ad
  segments(x0 = 6, y0 = 0.95, x1 = 7, y1 = 0.95, lwd=1.25); text(x = 6.5, y = 0.975, labels = "**", font = 2, cex = 1.25)
  segments(x0 = 7, y0 = 0.90, x1 = 8, y1 = 0.90, lwd=1.25); text(x = 7.5, y = 0.925, labels = "**", font = 2, cex = 1.25)
  segments(x0 = 7, y0 = 0.85, x1 = 9, y1 = 0.85, lwd=1.25); text(x = 8, y = 0.875, labels = "**", font = 2, cex = 1.25)
  # unex
  segments(x0 = 11, y0 = 0.95, x1 = 12, y1 = 0.95, lwd=1.25); text(x = 11.5, y = 0.975, labels = "**", font = 2, cex = 1.25)
  segments(x0 = 11, y0 = 0.90, x1 = 14, y1 = 0.90, lwd=1.25); text(x = 12.5, y = 0.925, labels = "*", font = 2, cex = 1.25)
  segments(x0 = 12, y0 = 0.85, x1 = 13, y1 = 0.85, lwd=1.25); text(x = 12.5, y = 0.875, labels = "**", font = 2, cex = 1.25)
  segments(x0 = 13, y0 = 0.80, x1 = 14, y1 = 0.80, lwd=1.25); text(x = 13.5, y = 0.825, labels = "**", font = 2, cex = 1.25)
  
  # plot 2 -- cell-types
  # empty background
  plot(0, 0, pch=16, col="white", xlim=c(0.5, 17.5), ylim=c(0, 1), ylab="Cell-type weight", xaxt='none', cex.lab=1.50, cex.axis=1.25, bty='n', xlab="")
  
  # grid
  for (i in seq(0.5, 17.5, (17.5-0.5)/10)){ abline(v=i, lwd=0.4, col="grey60") }
  for (i in seq(0, 1, 0.1)){ abline(h=i, lwd=0.4, col="grey60") }
  
  # parameters
  k <- 1
  w <- 0.4
  stripes <- c(0, 10, 25)  
  
  # main loop
  for (gr in 1:length(all_gr)){
    # get variants and compute means across groups
    ct <- var_gene_mapping_CT[which(var_gene_mapping_CT$ID %in% all_gr[[gr]]), 20:24]
    ct <- na.omit(ct)
    mean_ct <- colMeans(ct)
    # then plot rectangles
    for (i in 1:length(mean_ct)){
      rect(xleft = k-w, ybottom = 0, xright = k+w, ytop = mean_ct[i], col = colz_ct[i], lwd=1.5)
      rect(xleft = k-w, ybottom = 0, xright = k+w, ytop = mean_ct[i], col = NA, lwd=1.5, density = stripes[gr])
      k <- k + 1
    }
    k <- k + 1
  }
  text(x = c(3, 9, 15), y = -0.075, labels = c("Longevity-group", "AD-group", "Unex-group"), xpd=T, cex=1.50, font=4)
  legend(x = 1.5, y = 1.25, legend = c("Astrocytes", "Endothelial", "Myeloid", "Neuron", "Oligodendrocytes"), bty='n', cex=1.25, col="black", pch=22, ncol=3, pt.lwd = 1.5, pt.bg = colz_ct, xpd=T)  
  text(x = 18.5, y = 1.2, labels = "B", cex=2, font=2, xpd=T)
  # add significance
  # cha
  segments(x0 = 2, y0 = 0.95, x1 = 1, y1 = 0.95, lwd=1.25); text(x = 1.5, y = 0.975, labels = "**", font = 2, cex = 1.25)
  segments(x0 = 2, y0 = 0.90, x1 = 4, y1 = 0.90, lwd=1.25); text(x = 3, y = 0.925, labels = "**", font = 2, cex = 1.25)
  segments(x0 = 2, y0 = 0.85, x1 = 5, y1 = 0.85, lwd=1.25); text(x = 3.5, y = 0.875, labels = "**", font = 2, cex = 1.25)
  segments(x0 = 3, y0 = 0.80, x1 = 1, y1 = 0.80, lwd=1.25); text(x = 2, y = 0.825, labels = "**", font = 2, cex = 1.25)
  segments(x0 = 3, y0 = 0.75, x1 = 4, y1 = 0.75, lwd=1.25); text(x = 3.5, y = 0.775, labels = "**", font = 2, cex = 1.25)
  segments(x0 = 3, y0 = 0.70, x1 = 5, y1 = 0.70, lwd=1.25); text(x = 4, y = 0.725, labels = "**", font = 2, cex = 1.25)
  # ad
  segments(x0 = 9, y0 = 0.95, x1 = 7, y1 = 0.95, lwd=1.25); text(x = 8, y = 0.975, labels = "**", font = 2, cex = 1.25)
  segments(x0 = 9, y0 = 0.90, x1 = 8, y1 = 0.90, lwd=1.25); text(x = 8.5, y = 0.925, labels = "**", font = 2, cex = 1.25)
  segments(x0 = 9, y0 = 0.85, x1 = 10, y1 = 0.85, lwd=1.25); text(x = 9.5, y = 0.875, labels = "**", font = 2, cex = 1.25)
  segments(x0 = 9, y0 = 0.80, x1 = 11, y1 = 0.80, lwd=1.25); text(x = 10, y = 0.825, labels = "**", font = 2, cex = 1.25)
  # unex
  segments(x0 = 14, y0 = 0.95, x1 = 16, y1 = 0.95, lwd=1.25); text(x = 15, y = 0.975, labels = "**", font = 2, cex = 1.25)
  segments(x0 = 14, y0 = 0.90, x1 = 17, y1 = 0.90, lwd=1.25); text(x = 15.5, y = 0.925, labels = "**", font = 2, cex = 1.25)
  segments(x0 = 15, y0 = 0.85, x1 = 17, y1 = 0.85, lwd=1.25); text(x = 16, y = 0.875, labels = "*", font = 2, cex = 1.25)
  # finally annotation for the p-value and groups
  text(x = 0, y = -0.2, labels = "**: FDR<5%    *: FDR<10%", font=2, xpd=T, adj=0, cex=1)
  #text(x = 0, y = -0.25, labels = "LGV: Longevity", font=2, xpd=T, adj=0, cex=0.60)
}

# Function to draw forest plot with GWAS results -- only 100plus and Timmers
InDaForest_only100plus_and_Timmers <- function(var_gene_mapping, orig_lit){
  # load working framework
  load("/Users/nicco/Desktop/2k19_work/rotationPaper/FINAL_VERSION/RESULTS/analysis_GWAS_context.RData")
  # important: flip all variants to the AD-risk increasing allele (so that also direction can be seen)
  for (i in 1:nrow(all_long)){
    # grep variant association on ad
    ad_as <- orig_lit[which(orig_lit$rsid == all_long$rsid_JOR[i]),]
    if (ad_as$beta <0){
      # if beta <0 --> protective, make sure the allele tested on longevity is the other
      oth_all <- toupper(ad_as$a2)
      if (all_long$A1[i] != oth_all){
        all_long$A1[i] <- oth_all
        all_long$BETA[i] <- all_long$BETA[i]*(-1)
      }
    } else if (ad_as$beta >0) {
      all_risk <- toupper(ad_as$a1)
      if (all_long$A1[i] != all_risk){
        all_long$A1[i] <- all_risk
        all_long$BETA[i] <- all_long$BETA[i]*(-1)
      }
    }  
  }
  # use the same variants
  background <- var_gene_mapping[which(var_gene_mapping$ID %in% all_long$rsid_JOR),]
  # also need to add standard errors for my association
  se_sb <- res[, c("LOCUS", "SE", "P", "P_ADJ")]
  all_long <- merge(all_long, se_sb, by.x="locus_JOR", by.y="LOCUS")
  # define max
  max <- 34.5
  # first make sure all alleles and betas are the same
  for (i in 1:nrow(all_long)){
    # make sure to have same alleles
    my_beta <- all_long[i, "BETA"]
    my_se <- all_long[i, "SE"]
    my_a1 <- all_long[i, "A1"]
    if (my_a1 != toupper(all_long[i, "a1_JOR"])){ all_long[i, "beta_JOR"] <- all_long[i, "beta_JOR"]*(-1) }
    if (my_a1 != toupper(all_long[i, "a1_TIM"])){ all_long[i, "beta_TIM"] <- all_long[i, "beta_TIM"]*(-1) }
  }
  # first, calculate CIs
  all_long$upCI <- (all_long$BETA+(1.96*all_long$SE))
  all_long$lowCI <- (all_long$BETA-(1.96*all_long$SE))
  all_long$upCI_jor <- (all_long$beta_JOR+(1.96*all_long$se_JOR))
  all_long$lowCI_jor <- (all_long$beta_JOR-(1.96*all_long$se_JOR))
  all_long$upCI_tim <- (all_long$beta_TIM+(1.96*all_long$se_TIM))
  all_long$lowCI_tim <- (all_long$beta_TIM-(1.96*all_long$se_TIM))
  # scale betas and ci with respect to where points should be plotted
  all_long[, c("upCI_norm", "lowCI_norm", "BETA_norm")] <- 4*((all_long[, c("upCI", "lowCI", "BETA")] - min(all_long[, c("upCI", "lowCI", "BETA")]))/(max(all_long[, c("upCI", "lowCI", "BETA")])-min(all_long[, c("upCI", "lowCI", "BETA")])))+2
  all_long[, c("upCI_tim_norm", "lowCI_tim_norm", "beta_TIM_norm")] <- 4*((all_long[, c("upCI_tim", "lowCI_tim", "beta_TIM")] - min(all_long[, c("upCI_tim", "lowCI_tim", "beta_TIM")]))/(max(all_long[, c("upCI_tim", "lowCI_tim", "beta_TIM")])-min(all_long[, c("upCI_tim", "lowCI_tim", "beta_TIM")])))+7
  # i also need the 0 position
  pos_0_chc <- 4*((0 - min(all_long[, c("upCI", "lowCI", "BETA")]))/(max(all_long[, c("upCI", "lowCI", "BETA")])-min(all_long[, c("upCI", "lowCI", "BETA")])))+2
#  pos_0_jor <- 4*((0 - min(all_long[, c("upCI_jor", "lowCI_jor", "beta_JOR")]))/(max(all_long[, c("upCI_jor", "lowCI_jor", "beta_JOR")])-min(all_long[, c("upCI_jor", "lowCI_jor", "beta_JOR")])))+7
  pos_0_tim <- 4*((0 - min(all_long[, c("upCI_tim", "lowCI_tim", "beta_TIM")]))/(max(all_long[, c("upCI_tim", "lowCI_tim", "beta_TIM")])-min(all_long[, c("upCI_tim", "lowCI_tim", "beta_TIM")])))+7
  
  # graphical params
  par(mar=c(5, 4, 4, 5))
  # make background plot
  plot(0, 0, pch=16, col="white", bty='none', xaxt='none', yaxt='none', xlab="Effect size on longevity of AD risk-increasing alleles", ylab="", xlim=c(0, 11), ylim=c(0, 34), cex.lab=2)
  # add grid for each rectangle
  for (i in seq(0, 34, 3.4)){ segments(x0 = 2, y0 = i, x1 = 6, y1 = i, lwd=0.5, col="grey80"); segments(x0 = 7, y0 = i, x1 = 11, y1 = i, lwd=0.5, col="grey60")}
  for (i in seq(2, 6)){ segments(x0 = i, y0 = 0, x1 = i, y1 = 34, lwd=0.5, col="grey80") }
  for (i in seq(7, 11)){ segments(x0 = i, y0 = 0, x1 = i, y1 = 34, lwd=0.5, col="grey80") }
  # there will be 3 rectangles, each divided in 3
  rect(xleft = 2, ybottom = 0, xright = 6, ytop = 34, lwd=3, border="black")
  rect(xleft = 7, ybottom = 0, xright = 11, ytop = 34, lwd=3, border="black")
  # add names on top of the rectangles
  text(x = 4, y = 35, labels = "100-plus", cex=1.50, font=4, xpd=T, col="black")
  text(x = 9, y = 35, labels = "Timmers et al.", cex=1.50, font=4, xpd=T, col="black")
  # add rectangles for the groups
  rect(xleft = 2, ybottom = 0, xright = 6, ytop = 7.9, lwd=3, border="navy")
  rect(xleft = 7, ybottom = 0, xright = 11, ytop = 7.9, lwd=3, border="navy")
  rect(xleft = 2, ybottom = 8.1, xright = 6, ytop = 24.9, lwd=3, border="orange")
  rect(xleft = 7, ybottom = 8.1, xright = 11, ytop = 24.9, lwd=3, border="orange")
  rect(xleft = 2, ybottom = 34, xright = 6, ytop = 25.1, lwd=3, border="red")
  rect(xleft = 7, ybottom = 34, xright = 11, ytop = 25.1, lwd=3, border="red")
  # also add main line in the middle
  pos_0_list <- list(pos_0_chc, pos_0_tim)
  for (i in 1:length(pos_0_list)){ 
    segments(x0 = pos_0_list[[i]], y0 = 0, x1 = pos_0_list[[i]], y1 = 34, lwd=1, col="black")
    text(x = pos_0_list[[i]], y = -1, labels = 0, xpd=T)
    if (i == 1){
      text(x = 2, y = -1, labels = round(min(all_long[, c("upCI", "lowCI")]), 2), xpd=T)
      text(x = 6, y = -1, labels = round(max(all_long[, c("upCI", "lowCI")]), 2), xpd=T)
    } else if (i == 2){
      text(x = 7, y = -1, labels = round(min(all_long[, c("upCI_tim", "lowCI_tim")]), 2), xpd=T)
      text(x = 11, y = -1, labels = round(max(all_long[, c("upCI_tim", "lowCI_tim")]), 2), xpd=T)
    }
  }
  # main loop per variant
  for (i in 1:nrow(background)){
    # get line of interest
    sb <- all_long[which(all_long$rsid_JOR == background$ID[i]),]
    # now draw points
    # 100-plus
    if (max-i > 25){ ccol <- "red" } else if (max-i > 8) { ccol <- "orange" } else { ccol <- "navy" } 
    segments(x0 = sb$lowCI_norm, y0 = max-i, x1 = sb$upCI_norm, y1 = max-i, lwd=1.50, col=ccol)
    points(x = sb$BETA_norm, y = max-i, pch=16, col=ccol, cex=1.50)
    if (sb$P_ADJ <= 0.05){
      text(x = 5.9, y = max-i, labels = "**", font=2, col=ccol, adj=1, cex=1.50)
    } else if (sb$P <= 0.05){
      text(x = 5.9, y = max-i, labels = "*", font=2, col=ccol, adj=1, cex=1.50)
    }
    # timmers
    segments(x0 = sb$lowCI_tim_norm, y0 = max-i, x1 = sb$upCI_tim_norm, y1 = max-i, lwd=1.50, col=ccol)
    points(x = sb$beta_TIM_norm, y = max-i, pch=17, col=ccol, cex=1.50)
    if (sb$p_adj_TIM <= 0.05){
      text(x = 10.9, y = max-i, labels = "**", font=2, col=ccol, adj=1, cex=1.50)
    } else if (sb$p_TIM <= 0.05){
      text(x = 10.9, y = max-i, labels = "*", font=2, col=ccol, adj=1, cex=1.50)
    }
    # name on the side
    text(x = 0.5, y = max-i, labels = paste(sb$locus_JOR, "~"), adj = 1, xpd=T)
    text(x = 1.80, y = max-i, labels = orig_lit$gene[which(orig_lit$rsid == sb$rsid_JOR)], adj = 1, xpd=T, font=3)
  }
  # final -- groups annotation
  text(x = 11.2, y = 29.5, labels = "Longevity-group", cex=1.25, col="red", srt=90, font=4, xpd=T)
  text(x = 11.2, y = 16.5, labels = "AD-group", cex=1.25, col="orange", srt=90, font=4, xpd=T)
  text(x = 11.2, y = 4, labels = "Unex-group", cex=1.25, col="navy", srt=90, font=4, xpd=T)
  # legend
  legend(x = 11.1, y = 34, legend = c("*  p<0.05", "**  FDR<5%"), xpd=T, col="black", pch = c("", ""), pt.cex = 1.50, bty='n')
}

######################################

######################################
## MAIN
# 1. First, set the working directory
setwd("/Users/nicco/Desktop/2k19_work/rotationPaper/FINAL_VERSION/")

######################################
# 2. Create distribution of effect-sizes from literature
# 2.1 Read literature AD effects
lit <- fread("INPUTS/ADsnps_GRACE.txt", h=F, dec=",")
colnames(lit) <- c("chr", "pos", "a1", "a2", "beta", "se", "p", "rsid", "gene")

# 2.2 Make distribution for literature AD betas
set.seed(1234)
lit_dist <- function.distrLit(lit, 1000)

# 2.3 Also read random SNPs
rand_assoc = fread("RESULTS_with_random/snps_association_random.txt", h=T, sep = "\t")
rand_assoc = rand_assoc[, c("chr", "pos", "A1", "A2", "beta", "SE", "P", "RS")]
rand_assoc$gene = NA
colnames(rand_assoc) <- c("chr", "pos", "a1", "a2", "beta", "se", "p", "rsid", "gene")

# 2.4 Make distribution for literature AD random betas
lit_rand_dist = function.distrLit(rand_assoc, 1000)
######################################

######################################
# 3. Now load the bootstrap and combine effects
# 3.1 Load bootstrapped dataset of effects on aging -- this includes also random SNPs
#load("INPUTS/Bootstrap_survival_10k_ADsnps.RData")
load("RESULTS_with_random/Bootstrap_survival_10k_ADsnps.RData")
load("RESULTS_with_random/Bootstrap_survival_10k_ADsnps_random.RData")

# 3.2 Do some cleaning: make the transpose and remove two rare variants cause are causing problems
boot_clean <- t(boot_clean)
boot_clean_rand <- t(boot_clean_rand)

# 3.3 Check alleles and variant names -- output is the two datasets ordered the same way
ordered_list <- makeItTheSame(boot_clean, lit, lit_dist, p_filter = TRUE)
lit_dist <- ordered_list[[1]]
boot_dist <- ordered_list[[2]]
# for the random snps make sure the ids are ok
names_boot_rand <- as.data.frame(str_split_fixed(rownames(boot_clean_rand), "_", 2))
names_boot_rand$ID = paste(names_boot_rand$V1, names_boot_rand$V2, sep = "_")
rand_assoc = rand_assoc[which(rand_assoc$rsid %in% names_boot_rand$V1),]
names_boot_rand = names_boot_rand[which(names_boot_rand$V1 %in% rand_assoc$rsid),]
boot_clean_rand = boot_clean_rand[which(rownames(boot_clean_rand) %in% names_boot_rand$ID),]
lit_rand_dist = lit_rand_dist[which(rownames(lit_rand_dist) %in% rand_assoc$rsid),]
ordered_list_rand <- makeItTheSame(boot_clean_rand, rand_assoc, lit_rand_dist, p_filter = FALSE)
lit_dist_rand <- ordered_list_rand[[1]]
boot_dist_rand <- ordered_list_rand[[2]]

# 4. Combine data and get normalized angles
angle <- Transform(boot.x = boot_dist, boot.y = lit_dist)
angle_rand <- Transform(boot.x = boot_dist_rand, boot.y = lit_dist_rand)

# 5. Calculate bins in the region of the densities
bins <- function_bins(bin.number = 100, limits = c(0, 2))

# 6. Then calculate number of points per bin, per variant
# AD SNPs
res_pointsPerBin <- lapply(1:nrow(angle), PointsPerBin, bins=bins, angle=angle)
points_bin <- matrix(data=unlist(res_pointsPerBin), nrow = nrow(angle), ncol = nrow(bins), byrow = T)
rownames(points_bin) <- rownames(angle)
# random snps
res_pointsPerBin_rand <- lapply(1:nrow(angle_rand), PointsPerBin, bins=bins, angle=angle_rand)
points_bin_rand <- matrix(data=unlist(res_pointsPerBin_rand), nrow = nrow(angle_rand), ncol = nrow(bins), byrow = T)
rownames(points_bin_rand) <- rownames(angle_rand)

# 7. Also get the ordered dataframe of angles according to median
ordered_angle_median <- Reorder_median(angle)
angle = angle[match(ordered_angle_median$snp, rownames(angle)),]
ordered_angle_median_rand <- Reorder_median(angle_rand)
angle_rand = angle_rand[match(ordered_angle_median_rand$snp, rownames(angle_rand)),]

# 8. define the three groups based on median threshold -- <-0.5 -- -0.5<x<0 -- >0
gr1_snps = ordered_angle_median[which(as.numeric(ordered_angle_median$median) <= -0.5),]
gr1_snps_mean = ordered_angle_median[which(as.numeric(ordered_angle_median$mean) <= -0.5),]
gr2_snps = ordered_angle_median[which(as.numeric(ordered_angle_median$median) <= 0 & as.numeric(ordered_angle_median$median) > -0.5),]
gr2_snps_mean = ordered_angle_median[which(as.numeric(ordered_angle_median$mean) <= 0 & as.numeric(ordered_angle_median$mean) > -0.5),]
gr3_snps = ordered_angle_median[which(as.numeric(ordered_angle_median$median) > 0),]
gr3_snps_mean = ordered_angle_median[which(as.numeric(ordered_angle_median$mean) > 0),]
# also for random snps
gr1_snps_rand = ordered_angle_median_rand[which(as.numeric(ordered_angle_median_rand$median) <= -0.5),]
gr1_snps_rand_mean = ordered_angle_median_rand[which(as.numeric(ordered_angle_median_rand$mean) <= -0.5),]
gr2_snps_rand = ordered_angle_median_rand[which(as.numeric(ordered_angle_median_rand$median) <= 0 & as.numeric(ordered_angle_median_rand$median) > -0.5),]
gr2_snps_rand_mean = ordered_angle_median_rand[which(as.numeric(ordered_angle_median_rand$mean) <= 0 & as.numeric(ordered_angle_median_rand$mean) > -0.5),]
gr3_snps_rand = ordered_angle_median_rand[which(as.numeric(ordered_angle_median_rand$median) > 0),]
gr3_snps_rand_mean = ordered_angle_median_rand[which(as.numeric(ordered_angle_median_rand$mean) > 0),]
# explore the snps groups
  tmp_x = boot_dist_rand[which(rownames(boot_dist_rand) == "rs62069508"),]
  tmp_y = lit_dist_rand[which(rownames(lit_dist_rand) == "rs62069508"),]
  plot(tmp_x, tmp_y, 
       pch = 16, xlim = c(min(tmp_x, -1), max(tmp_x, 1)), col = viridis(n = 1, alpha = 0.5, option = "plasma"), cex = 1, xlab = "Bootstrapped effect on longevity",
       ylab = "Sampled effect on AD", ylim = c(-1, 1))
  b = 1
  a = -1
  angle = (b-a) * ((angle-min(angle))/(max(angle)-min(angle))) + a 
  angle_rand = (b-a) * ((angle_rand - min(na.omit(angle_rand)))/(max(na.omit(angle_rand))-min(na.omit(angle_rand)))) + a 
  abline(h=0)
  abline(v=0)
  abline(v=-0.5, lty = 2)
  d = density(angle_rand[which(rownames(angle_rand) == "rs62069508"),])
  #d = density(angle_rand)
  b = 0.80
  a = 0
  d$y = (b-a) * ((d$y - min(d$y))/(max(d$y) - min(d$y))) + a
  points(x = d$x, y = d$y, type = "l", col = viridis(n = 3, alpha = 0.8, option = "plasma")[2], lwd = 3)
  
  
  hist(angle_rand[which(rownames(angle_rand) == "rs7188"),])
  
# 9. graphical verification of weird behaviors
# 9.1 plot random snps with weird behavior
target = "rs7306136"
ad_dist = lit_dist_rand[which(rownames(lit_dist_rand) == target),]
lg_dist = boot_dist_rand[which(rownames(boot_dist_rand) == target),]
plot(lg_dist, ad_dist)
abline(h=0)
abline(v=0)
# this is weird because AD effects are negavite also, which causes the flipping
# 9.2 barplot of percentagle of snps in the different groups
opar <- par(lwd = 4)
pos = barplot(c(nrow(gr1_snps)/nrow(angle), nrow(gr2_snps)/nrow(angle), nrow(gr3_snps)/nrow(angle), nrow(gr1_snps_rand)/nrow(angle_rand), nrow(gr2_snps_rand)/nrow(angle_rand), nrow(gr3_snps_rand)/nrow(angle_rand)), 
        names.arg = c("Longevity group", "AD group", "Unexpected group", "Longevity group", "AD group", "Unexpected group"),
        col = rep(viridis(n = 3, alpha = 0.8, option = "viridis"), 2), ylim = c(0, 1), ylab = "Percentage (%)",
        density = c(0,0,0, 10,10,10), border = rep(viridis(n = 3, alpha = 0.8, option = "viridis"), 2))
abline(v = 3.7, lwd = 3, lty = 2)
text(x = 1.9, y = 1, labels = "AD SNPs from GWAS", font = 2, pos = 1, offset = 0.5, cex = 2)
text(x = 5.5, y = 1, labels = "Random SNPs", font = 2, pos = 1, offset = 0.5, cex = 2)
######################################
# Plots
# Plot densities ordered according to median
plotDensities_noDendro(angle, lit = lit, rand=F)
plotDensities_noDendro(angle = angle_rand, lit = rand_assoc, rand = TRUE)
plotDensities(line.x = angle_rand[sample(x = 1:nrow(angle_rand), size = 20, replace = F),], ordered_angle_median_rand)

# Barplot of the variant groups
Plot_variantGroups(res_groups)

# Final figure 1
png("RESULTS/figure_1_v3.png", height = 14, width = 17, res=600, units="in")
pdf("RESULTS/figure_1_v3.pdf", height = 14, width = 17)
#png("/Users/nicco/Desktop/THESIS/Thesis_NT/chapter_5_rotation/figure_1_v4.png", height = 19, width = 17, res=600, units="in")
function.Figure1(var_gene_mapping_CT, angle, lit)
dev.off()

# Final figure 2
png("RESULTS/figure_2_v3.png", height = 10, width = 10, res=600, units="in")
pdf("RESULTS/figure_2_v3.pdf", height = 10, width = 10)
function_figure2(var_gene_mapping_CT, final_within_groups, final_within_groups_CT, all_gr)
dev.off()

######################################
# Now, with the input variants I will run AnnotateMe tool 
# The gene-set overlap analysis is done differently on the cluster, after that need to merge results
all_avg <- mergeSampling_clusProf()

# Instead of running REVIGO, I will do semantic similarity myself with a python script
# The command is very simple:
# python3 Alternative_REVIGO.py ../RESULTS/AnnotateMe/newMeth.revigo.txt [INPUT] ../RESULTS/AnnotateMe/Lin_distance_measure.txt [OUTPUT - LIN] ../RESULTS/AnnotateMe/Resnik_distance_measure.txt [OUTPUT - RESNIK]
# Here i start from the output of the semantic similarity matrix
results_functional <- function.GlobalProps(angle)
functional_clusters <- results_functional[[1]]
lin_matrix <- results_functional[[2]]
n_clust <- results_functional[[3]]

# Would be nice to implement a wordcloud visualisation for the clusters
# First, need to remove frequent words. For this, need to get all terms from GO:BP
# The following bash command can be used to do that
# wget http://current.geneontology.org/ontology/go-basic.obo .
# grep -B 3 -w "namespace: biological_process" go-basic.obo | grep "name: " | sed 's/name: //g' > go_terms_BP_all.txt
# rm go-basic.obo
# Then use the function to count the frequency of words in different scenarios
mostFreq_words <- CountFrequency_words(functional_clusters, n_clust)
#######################################

#######################################
# Now it is time to integrate densities and functional clusters
# Find all genes underlying each cluster
hsGO <- godata('org.Hs.eg.db', ont="BP")
var_gene_mapping <- FindGenesPerCluster_clusProf(functional_clusters, all_avg, lin_matrix)

# Order in the same way as the ordered densities densities
var_gene_mapping <- matchBy(var_gene_mapping, ordered_angle_median$snp, var_gene_mapping$ID)

# Run comparison between groups
results_comparison <- CompareClusterAnnot_mod_3grps(var_gene_mapping, angle, all_gr, n_clust)
functional_means_perGroup <- results_comparison[[1]]
functional_sd_perGroup <- results_comparison[[2]]
between_groups_comparison <- results_comparison[[3]]
within_groups_comparison <- results_comparison[[4]]
one_vs_all_groups_comparison <- results_comparison[[5]]

# Finally correct pvalues
final_between_groups <- correctMe(between_groups_comparison, skipFirst = T)
final_within_groups <- correctMe(within_groups_comparison, skipFirst = F)
final_oneVsAll_groups <- correctMe(one_vs_all_groups_comparison, skipFirst = T)
######################################

######################################
# Now it is time to look at the single-cell expression
# Find expression of genes in hippocampus in the different cell-types
var_gene_mapping_CT <- cellType_expression(var_gene_mapping)

# Then the final thing is to run comparisons at cell-type level
results_comparison_CT <- CompareClusterAnnot_cellTypes(var_gene_mapping_CT, angle, all_gr)
CT_means_perGroup <- results_comparison_CT[[1]]
CT_sd_perGroup <- results_comparison_CT[[2]]
between_groups_CT <- results_comparison_CT[[3]]
within_groups_CT <- results_comparison_CT[[4]]
one_vs_all_groups_CT <- results_comparison_CT[[5]]

# Finally correct pvalues
final_between_groups_CT <- correctMe(between_groups_CT, skipFirst = T)
final_within_groups_CT <- correctMe(within_groups_CT, skipFirst = F)
final_oneVsAll_groups_CT <- correctMe(one_vs_all_groups_CT, skipFirst = T)
####################################

####################################
# Also good to lokk into the association of these variants on large gwas of longevity
# Then do a forest plot based on the effect sizes
png("RESULTS/Figure_3_v3.png", height = 8, width = 10, res=300, units = "in")
#pdf("RESULTS/Figure_3_v3.pdf", height = 8, width = 10)
png("/Users/nicco/Desktop/THESIS/Thesis_NT/chapter_5_rotation/figure_3_v4.png", height = 12, width = 10, res=300, units = "in")
InDaForest_only100plus_and_Timmers(var_gene_mapping, orig_lit = lit)
dev.off()
####################################

####################################
# Last, it is good to save this workspace as next time don't need to run the whole thing again
save.image("Rotation_paper_analysis_workspace.RData")
load("/Users/nicco/Desktop/2k19_work/rotationPaper/FINAL_VERSION/Rotation_paper_analysis_workspace.RData")

####################################
####################################
# does it make sense to look at expression of genes in AD and aging?
allen_all <- fread("/Users/nicco/Desktop/2k19_work/rotationPaper/FINAL_VERSION/gene_expression_matrix_2016-03-03/fpkm_table_normalized.csv")
annot_genes <- fread("/Users/nicco/Desktop/2k19_work/rotationPaper/FINAL_VERSION/gene_expression_matrix_2016-03-03/rows-genes.csv")
annot_samples <- fread("/Users/nicco/Desktop/2k19_work/rotationPaper/FINAL_VERSION/gene_expression_matrix_2016-03-03/columns-samples.csv")
full_samples_info <- fread("/Users/nicco/Desktop/2k19_work/rotationPaper/FINAL_VERSION/gene_expression_matrix_2016-03-03/DonorInformation.csv")

# merge annotation of the samples with the full information
samples_info <- merge(full_samples_info, annot_samples, by="donor_id")

# put all ad samples somewhere
ad <- samples_info[which(samples_info$nincds_arda_diagnosis == "Probable Alzheimer'S Disease"),]
nd <- samples_info[which(samples_info$nincds_arda_diagnosis == "No Dementia"),]

# let's look at hippocampus
ad_ <- ad[which(ad$structure_name == "hippocampus (hippocampal formation)"),]
nd_hipp <- nd[which(nd$structure_name == "hippocampus (hippocampal formation)"),]

# samples are set -- let take only those
kk <- names(allen_all)
samples_index <- c()
for (i in 1:nrow(ad_hipp)){ tmp <- ad_hipp$rnaseq_profile_id[i]; samples_index <- c(samples_index, grep(ad_hipp$rnaseq_profile_id[i], kk)) }
for (i in 1:nrow(nd_hipp)){ tmp <- nd_hipp$rnaseq_profile_id[i]; samples_index <- c(samples_index, grep(nd_hipp$rnaseq_profile_id[i], kk)) }
allen_samples_interest <- allen_all[, ..samples_index]
allen_samples_interest <- cbind(allen_all[, 1], allen_samples_interest)

# then need to select the genes

genes_ids <- allen_samples_interest[, 1]
myGenes <- fread("/Users/nicco/Desktop/2k19_work/rotationPaper/FINAL_VERSION/RESULTS/AnnotateMe/snp_annotation_geneList.txt", h=F)
annot_genes_my <- annot_genes[which(annot_genes$gene_symbol %in% myGenes$V1),]
colnames(allen_samples_interest)[1] <- "gene_id"
expr_myGenes <- merge(annot_genes_my, allen_samples_interest, by="gene_id")
head(expr_myGenes)

## function to average expression of genes mapping to the same variant
avgExpr <- function(expr_myGenes, var_gene_mapping){
  # define list for the results
  res_list = list()
  for (snp in 1:nrow(var_gene_mapping)){
    # get gene list
    g_list <- unlist(strsplit(var_gene_mapping$geneList[snp], ","))
    # get the expression of those genes
    expr_g <- expr_myGenes[which(expr_myGenes$gene_symbol %in% g_list),]
    # take avg expression
    avg_exp <- as.data.frame(t(colMeans(expr_g[, 6:ncol(expr_g)])))
    # create output df
    avg_exp <- cbind(avg_exp, expr_g[1, 1:5])
    avg_exp$snp_id <- var_gene_mapping$ID[snp]
    # add to results list
    res_list[[snp]] <- avg_exp
  }
  # put things together
  res_list <- rbindlist(res_list)
  res_list <- res_list[!duplicated(res_list$gene_symbol),]
  return(res_list)
}
avg_expr <- avgExpr(expr_myGenes, var_gene_mapping)

## need to plot expression wrt age
# check age
all_samples <- rbind(ad_hipp, nd_hipp)
length(table(all_samples$age))
order_age <- c("78", "79", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90-94", "95-99", "100+")
## function to make plot
plotOverAge <- function(order_age, all_samples, avg_expr){
  # take groups of variants
  cha_gr <- avg_expr[which(avg_expr$snp_id %in% c(all_gr[[1]])),]
  ad_gr <- avg_expr[which(avg_expr$snp_id %in% c(all_gr[[2]])),]
  unex_gr <- avg_expr[which(avg_expr$snp_id %in% c(all_gr[[3]])),]
  
  # first, cha-group
  # modify expression data
  mod_exp <- as.data.frame(t(cha_gr[, 1:67]))
  colnames(mod_exp) <- cha_gr$gene_symbol
  # add ages
  ages_sb <- all_samples[, c("rnaseq_profile_id", "age", "sex", "education_years", "apo_e4_allele", "nincds_arda_diagnosis")]
  # merge all together
  all_tog <- merge(mod_exp, ages_sb, by.x="row.names", by.y="rnaseq_profile_id")
  all_tog$pch <- 21
  all_tog$pch[which(all_tog$nincds_arda_diagnosis == "No Dementia")] <- 23
  # background plot
  par(mar=c(5, 5, 2, 8))
  plot(0, 0, pch=16, col="white", xlab="Age", ylab="RNA expression", xlim=c(1, 14), ylim=c(0, max(all_tog[, 2:11])), xaxt='n')
  for (i in seq(0, max(all_tog[, 2:11]), max(all_tog[, 2:11])/10)){abline(h=i, lwd=0.4, col="grey60")}
  for (i in seq(1, 14, 1.3)){abline(v=i, lwd=0.4, col="grey60")}
  axis(side = 1, at = seq(1, 14), labels = order_age, cex.axis=1.25)
  # then loop on ages
  x_pos <- 1
  gene_cols <- brewer.pal(n = 10, name = "Set3")
  for (x in order_age){
    ad <- colMeans(all_tog[which(all_tog$age == x & all_tog$nincds_arda_diagnosis != "No Dementia"), 2:11])
    nd <- colMeans(all_tog[which(all_tog$age == x & all_tog$nincds_arda_diagnosis == "No Dementia"), 2:11])
    for (g in 2:11){
      points(x = x_pos, y = nd[g], pch=21, bg=gene_cols[g-1], cex=1.50)
      points(x = x_pos, y = ad[g], pch=23, bg=gene_cols[g-1], cex=1.50)
    }
    x_pos <- x_pos + 1
  }
  legend(x = 15, y = 30, legend = names(all_tog)[2:11], xpd=T, pch=21, pt.bg=gene_cols)
  legend(x = 15, y = 15, legend = c("AD", "Control"), pch=c(21, 23), xpd=T)
  
  #######################################################################
  # alternatively, do regression over time and show only regression lines
  # make the age a numeric value
  all_tog$age[which(all_tog$age == "78")] <- 78
  all_tog$age[which(all_tog$age == "79")] <- 79
  all_tog$age[which(all_tog$age == "81")] <- 81
  all_tog$age[which(all_tog$age == "82")] <- 82
  all_tog$age[which(all_tog$age == "83")] <- 83
  all_tog$age[which(all_tog$age == "84")] <- 84
  all_tog$age[which(all_tog$age == "85")] <- 85
  all_tog$age[which(all_tog$age == "86")] <- 86
  all_tog$age[which(all_tog$age == "87")] <- 87
  all_tog$age[which(all_tog$age == "88")] <- 88
  all_tog$age[which(all_tog$age == "89")] <- 89
  all_tog$age[which(all_tog$age == "90-94")] <- 92
  all_tog$age[which(all_tog$age == "95-99")] <- 97
  all_tog$age[which(all_tog$age == "100+")] <- 101
  all_tog$age <- as.numeric(all_tog$age)
  par(mfrow=c(5, 2), mar=c(1, 3, 3, 1))
  # main loop across genes
  for (g in 2:11){
    plot(0, 0, pch=16, col="white", xlab="Age", ylab="RNA expression", xlim=c(75, 105), ylim=c(0, max(all_tog[, 2:11])), bty='n', main=names(all_tog)[g])
    for (i in seq(0, max(all_tog[, 2:11]), max(all_tog[, 2:11])/10)){abline(h=i, lwd=0.4, col="grey60")}
    for (i in seq(75, 105, (105-75)/10)){abline(v=i, lwd=0.4, col="grey60")}
    axis(side = 1, at = seq(1, 14), labels = order_age, cex.axis=1.25)
    # get data
    sb <- all_tog[, c(g, 12:16)]
    colnames(sb) <- c("gene", "age", "sex", "education", "apoe", "pheno")
    ad <- sb[which(sb$pheno != "No Dementia"),]
    nd <- sb[which(sb$pheno == "No Dementia"),]
    # plot points
    model_ad <- lm(gene~age+sex+education, data=ad)
    points(ad$age, ad$gene, pch=16, col="red")
    abline(lm(gene~age+sex+education, data=ad), col="red", lty=1, lwd=2)
    if (summary(model_ad)$coefficient[2, 4] <0.05){text(x = 100, y = 30, labels = "sign.assoc in AD")}
    points(nd$age, nd$gene, pch=16, col="navy")
    model_nd <- lm(gene~age+sex+education, data=nd)
    abline(lm(gene~age+sex+education, data=nd), lwd=2, col="navy")
    if (summary(model_nd)$coefficient[2, 4] <0.05){text(x = 100, y = 30, labels = "sign.assoc in ND")}
  }
  legend(x = 15, y = 30, legend = names(all_tog)[2:11])
}

########################################
########################################


