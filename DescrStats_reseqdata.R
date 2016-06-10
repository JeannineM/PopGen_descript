if (F) {setwd("C:/Users/Alex/Desktop/bluebell/script_vAug11-15")} # choose wd for alex
if (F) {setwd("/Users/Jeannine/Dropbox/Parkplatz blue bells/6. Data analyses/nuclear_genotypes/script_vAug15")} # choose wd for jeannine

# load packages
if (F){
  library("fields")
  library("geosphere")
  library("ade4")
  library("RColorBrewer")
  library("gplots")
  library("ggplot2")
  library('ggbiplot')
  library('MASS')
  library('ape')
  library('MultinomialCI')
} # new 25/08

#import function
source('function.R')

#import data

if (F){
  raw <- read.table('/Users/Jeannine/Documents/NGS_transcriptomeSWA/12_Re_seq_data/04_analyses/organelle/org_GTpopfreq.dat3.txt', sep='\t', header = T, na.strings="NA")
  
  nc_2228Loc_dat <- read.table('/Users/Jeannine/Dropbox/Parkplatz blue bells/1. Data/final_Snpcalls/nc_hatkgHC.012',sep='\t', header = F, na.strings="-1")

  pos <- read.table('/Users/Jeannine/Dropbox/Parkplatz blue bells/1. Data/final_Snpcalls/nc_hatkgHC.012.MOD2.txt',sep='\t', header = T, stringsAsFactors = FALSE)

  ind <- read.table('/Users/Jeannine/Dropbox/Parkplatz blue bells/1. Data/final_Snpcalls/nc_hatkgHC.012.indv',sep='\t', header = F, stringsAsFactors = FALSE)
  
  # add data columns to the data frame and save as table file
  
  popinfo <- read.table('/Users/Jeannine/Dropbox/Parkplatz blue bells/1. Data/final_Snpcalls/nc_hatkgHC.012.popinfo' ,sep='\t', header = T, stringsAsFactors = FALSE)
}# jeannine

if (T){
  nc_2228Loc_dat <- read.table('nc_hatkgHC.012',sep='\t', header = F, na.strings="-1")

  ind <- read.table('nc_hatkgHC.012.indv',sep='\t', header = F, stringsAsFactors = FALSE)
  
  pos <- read.table('nc_hatkgHC.012.MOD2.txt',sep='\t', header = T, stringsAsFactors = FALSE)
  
  popinfo <- read.table('nc_hatkgHC.012.popinfo' ,sep='\t', header = T, stringsAsFactors = FALSE)
}# alex

# add proper column and row names to data
if (T){
  nc_dat.df <- nc_2228Loc_dat[,2:dim(nc_2228Loc_dat)[2][1]]
  dimnames(nc_dat.df)[[1]] <- ind$V1 
  dimnames(nc_dat.df)[[2]] <- pos$CHROM_POS
  
  dimnames(popinfo)[[1]] <- ind$V1 
  
  rownames(pos) <- pos$CHROM_POS
  posinfo <- t(pos)

  targetSNPs <- nc_dat.df[,pos$set == "Target"]
  
  nc_dat.dfv2 <- data.frame(popinfo, nc_dat.df)
  
  addinfo_col_data=dim(popinfo)[2]
  
  # first gene: nc_dat.dfv2[, 11]
  # last gene: dim(nc_dat.dfv2)[[2]] -> 2238
  # write.table(nc_dat.dfv2, "/Users/Jeannine/Documents/NGS_transcriptomeSWA/12_Re_seq_data/03_snpcalls/snpcall_11-06/GT_1069SNPs_mod.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
  

 if (F){hzdat <- subset(raw, subset = raw$Order == "hz")}
}

# explore the H-ns frequency in the organelles
if (F){
  hist(hzdat$freq_SWA2, main="Histogram of allele proportions per population", xlab = "proportion ns genotype p.pop")
  hist(hzdat$freq_SWA2[hzdat$freq_SWA2 != 1 & hzdat$freq_SWA2 != 0], main = "Histogram of allele proportions in hybrid pops", xlab = "proportion ns genotype p.pop")
  abline(v=mean(hzdat$freq_SWA2[hzdat$freq_SWA2 != 1 & hzdat$freq_SWA2 != 0]), col="red")

# divived the H-ns genotype frequency
  color <- hzdat$freq_SWA2
  (color <- ifelse(color == 0, "blue", # Hisp GT
                 ifelse(color > 0 & color <= 0.5, "green", # HyS-GT
                        ifelse(color > 0.5 & color <= 0.9, "orange", # HyN-GT
                               ifelse(color > 0.9 , "yellow", NA))))) # NonScr GT

# Map H-ns-GT frequency per population
  plot(hzdat$lat ~ hzdat$lon, ylab="latitude", xlab="longitude", main="Geographic distribution of the Hns genotype frequency in each location", col=color, pch=19, cex =0.8)

# H-ns-GT frequency per population along latitude
# plot(hzdat$lat ~ hzdat$freq_SWA2, xlab = "proportion ns genotype", ylab = "S --> N", main="H-ns-GT frequency per population along latitude", col=color, pch=19)
# H-ns-GT frequency per population along longitude
# plot(hzdat$freq_SWA2 ~ hzdat$lon, xlab = "W --> E", ylab = "proportion ns genotype", main="H-ns-GT frequency per population along longitude",col=color, pch=19)

# Cline of H-ns-GT frequency per population along latitude close to the hybrid zone
  plot(hzdat$freq_SWA2 ~ hzdat$lat, xlab="latitude", ylab = "Frequency of Hns genotype per population", main="Cline of Hns genotype along latitude", col=color, pch=19)

# Cline inlcuding distant samples
# color2 <- raw$freq_SWA2
# (color2 <- ifelse(color2 == 0, "blue", # Hisp GT
#                  ifelse(color2 > 0 & color2 <= 0.5, "green", # HyS-GT
#                         ifelse(color2 > 0.5 & color2 <= 0.9, "orange", # HyN-GT
#                                ifelse(color2 > 0.9 , "yellow", NA))))) # NonScr GT
# plot(raw$freq_SWA2 ~ raw$lat, ylab = "proportion ns genotype per population", col=color2, pch=19)
}

# Transform 012 to Minor allele frequency (MAF), name data
if (T){
  data=popinfo
  counter=0
  for ( i in 1:dim(nc_dat.df)[2]){
    if (mean(nc_dat.df[,i],na.rm=TRUE) <1){
      data=cbind(data,nc_dat.df[,i])
    }
    else{
      data=cbind(data,abs(nc_dat.df[,i]-2))
      counter=counter+1
    }
    if (mean(nc_dat.df[,i],na.rm=TRUE)==1){print(1)}
  }
  
  colnames(data) <- colnames(nc_dat.dfv2)
  
  # represent data as matrix
  testMAF <- as.matrix(data[,(addinfo_col_data+1):dim(data)[2]])
  
  # write.table(data,file="/Users/Jeannine/Dropbox/Parkplatz blue bells/6. Data analyses/Popgen/nc_hatkgHC.012_MAF.tab.txt",row.names = TRUE, col.names = TRUE, sep = "\t")
  
} # create a table, name data, with 2 always corresponding to the Minor Allele

# number of gene, snps and origin and subpop
if (T){
  gene_info=array(0,c(1,3))
  for (i in 1:dim(pos)[[1]]){
    nqme=unlist(strsplit(pos[i,1],'|',fixed=TRUE))
    if (gene_info[dim(gene_info)[1],1]==nqme[1]){
      gene_info[dim(gene_info)[1],2]=as.numeric(gene_info[dim(gene_info)[1],2])+1
    } else{
      gene_info=rbind(gene_info,c(nqme[1], 1 ,substr(unlist(strsplit(nqme[2],'_'))[1],4,4)))
    }
  }
  # write.table(we, "/Users/Jeannine/Dropbox/Parkplatz blue bells/6. Data analyses/Popgen/genlist_info_vector.tab.txt",row.names = TRUE, col.names = TRUE, sep = "\t")
  
  list_ns=seq(312)[data$Sp2=="ns"]
  list_hisp=seq(312)[data$Sp2=="hisp"]
  list_hyN=seq(312)[data$Sp2=="hyN"]
  list_hyS=seq(312)[data$Sp2=="hyS"]
  
} # variable: gene_info, list_hisp, list_ns, list_hyN, list_hyS

# Estimate allele frequencies per sampling population
if (T){
  f=c()
  g=c()
  CPF=c()
  N=c()
  for (j in levels(as.factor(data$pop.ID))){
    test=data[data$pop.ID==j,]
    values=sapply(test[,(addinfo_col_data+1):dim(data)[2]],mean,na.rm=TRUE)/2
    freqcp=mean(test$cp.GT,na.rm=TRUE)
    Nsamp=dim(test)[1]
    f=rbind(f,values)
    CPF=rbind(CPF,freqcp)
    N=rbind(N, Nsamp)
    g=rbind(g,test[1,2:dim(popinfo)[2]])
  }
  
  data_pop2=data.frame(g,f,row.names=g$pop.ID)
  data_pop3=data.frame(g$pop.ID,N, CPF, g[,3:9],f,row.names=g$pop.ID)

} # variable: data_pop2 + data_pop3 * jeannine 07-08

# Estimate mean MAF frequency * jeannine 10-08

if (F) {
  
  MAFmean51 <- rowMeans(data_pop3[,10:dim(data_pop3)[2]], na.rm = TRUE)
  # colMAF=gsub("hyS","orange",gsub("hyN","green",gsub("ns","blue",gsub("hisp","red",data_pop2$Sp2[order(MAFmean51)]))))
  # plot(MAFmean51[order(MAFmean51)], col=colMAF)

  popsumtab <- data.frame(data_pop3[,1:9], MAFmean51)

} # ... per sampling location # popsumtab # 51 subsets

if (F){
  k=c()
  CPF=c()
  g=c()
  for (j in levels(as.factor(data$Sp2))){
    test=data[data$Sp2==j,]
    values=sapply(test[,(addinfo_col_data+1):dim(data)[2]],mean,na.rm=TRUE)/2
    freqcp=mean(test$cp.GT,na.rm=TRUE)
    k=rbind(k,values)
    CPF=rbind(CPF,freqcp)
    g=rbind(g,test[1,2:dim(popinfo)[2]])
  }
  
  
  MAFmean4 <- rowMeans(k, na.rm = TRUE)
  data_subset4=data.frame(CPF,MAFmean4,row.names=g$Sp2)
  
} # ...per subpopulation # 4 subsets

# Estimate observed heterozygosity
if (F){
  if (T){
save_Ho=c()
  for (j in 1:2228) {
    example=data[,addinfo_col_data+j]
    save_Ho=c(save_Ho, length(subset(example, example==1)))
  } # total 
  save_Ho=(save_Ho/312)/2
  par(mfrow=c(2,2))
  plot(save_Ho,main="Total observed heterozygosity", ylim=c(0,0.5))
  hist(save_Ho, breaks=100, xlim=c(0,0.5))
  abline(v = mean(save_Ho), col = "blue")
  plot(save_Ho[pos$set=="Target"],main="Total observed Heterozygosity for the target SNPs", ylim=c(0,0.5))
  hist(save_Ho[pos$set=="Target"], breaks=100, xlim=c(0,0.5))
  abline(v = mean(save_Ho[pos$set=="Target"]), col = "blue")
  
  # mean(save_Ho)

} # represent total observed heterozygosity jeannine 10-08
  
  if (F){
  g=c()
  f=c()
  for (m in levels(as.factor(data$pop.ID))){
      example=data[data$pop.ID==m,]
      het=colSums(example[,(addinfo_col_data+1):2238]==1,na.rm=T)
      f=rbind(f, het)
  }
  
  rownames(f) <- levels(as.factor(data$pop.ID))
  Ho_mean <- rowMeans(f)
  
  # popsumtab["HoMean"] <- Ho_mean
  
  hist(f, main="Total observed heterozygosity \n per 51 sampling location")
  abline(v = mean(f), col = "blue")

  } # count observed heterozygosity per 51 sampling location, add column to popsumtab
}

# A heatmap of genotypes 
if (F){

  my_palette <- colorRampPalette(c("darkgreen", "gold", "darkorchid1"))(n = 5)
  
  col_breaks = c(seq(-0.5,0.5, length=2),             # for major allele
                 seq(0.5001,1.5,length=2),            # for heterozygote
                 seq(1.5001,2.5,length=2))            # for minor allele
  
  pop_display=gsub("hisp","red",data$Sp2)
  pop_display=gsub("ns","blue",pop_display)
  pop_display=gsub("hyN","green",pop_display)
  pop_display=gsub("hyS","orange",pop_display)
  
  pop_chlr_display=gsub("0","blue",data$cp.GT)
  pop_chlr_display=gsub("1","red",pop_chlr_display)
  
  # with euclidian distance (not good)
  if (F){
    heatmap_MAF=heatmap.2(as.matrix(data[,(addinfo_col_data+1):dim(data)[2]]), main="Clustering of samples by MAF genotype", trace="none",  Colv=FALSE, margins =c(1.09,1.09), breaks=col_breaks, col=my_palette,density.info = "density",ylab = "" , xlab = "", distfun=dist, hclustfun=hclust, dendrogram="row", labCol = NA, labRow=NA,RowSideColors=pop_display)
    heatmap_MAF_chlr=heatmap.2(as.matrix(data[,(addinfo_col_data+1):dim(data)[2]]), main="Clustering of samples by MAF genotype", trace="none",  Colv=FALSE, margins =c(1.09,1.09), breaks=col_breaks, col=my_palette,density.info = "density",ylab = NA , xlab = NA, hclustfun=hclust, dendrogram="row", labCol = NA, labRow=NA,RowSideColors=pop_chlr_display)
  }
  
  # with manhattan distance ( much better)
  if (T){
    d <- dist(testMAF, method="manhattan")
    hr <- hclust(d, method="ward.D2")
    # manhattan distance = sum of pairwise absolute distance between each locus
    # Ward's minimum variance hier. clustering ward.D2 implements Ward's 1963 clustering criterion
    heatmap.2(testMAF, Rowv=as.dendrogram(hr), Colv=FALSE, margins =c(1.09,1.09), breaks=col_breaks, col=my_palette, density.info = "density",ylab = NA , xlab = NA, dendrogram="row", labCol = NA, labRow=NA, RowSideColors=pop_display, scale="none", trace="none")
    
    heatmap.2(testMAF, Rowv=as.dendrogram(hr), Colv=FALSE, margins =c(1.09,1.09), breaks=col_breaks, col=my_palette, density.info = "density",ylab = NA , xlab = NA, dendrogram="row", labCol = NA, labRow=NA, RowSideColors=pop_chlr_display, scale="none", trace="none")
    
    heatmap.2(as.matrix(d), Rowv=as.dendrogram(hr), Colv=rev(as.dendrogram(hr)), trace="none",density.info = "density", ColSideColors=pop_chlr_display, RowSideColors = pop_display, dendrogram = "row")
    
    # Target MAF SNPs
    
    targetSNPs
    d_T <- dist(as.matrix(targetSNPs), method="manhattan")
    hr_T <- hclust(d_T, method="ward.D2") 
    
    heatmap.2(as.matrix(targetSNPs), Rowv=as.dendrogram(hr_T), Colv=FALSE, margins =c(1.09,1.09), breaks=col_breaks, col=my_palette, density.info = "density",ylab = NA , xlab = NA, dendrogram="row", labCol = NA, labRow=NA, RowSideColors=pop_display, scale="none", trace="none")
    
    heatmap.2(as.matrix(d_T), Rowv=as.dendrogram(hr_T), Colv=rev(as.dendrogram(hr_T)), trace="none",density.info = "density", ColSideColors=pop_chlr_display, RowSideColors = pop_display, dendrogram = "row")
    
    ## Citation
    # Gao X, Starmer J (2007) Human population structure detection via multilocus genotype clustering. BMC Genet 8, 34.
    
    
  }
  
}

# Isolation by distance
if (F){
  
  # remove the British BB-506 sample 
  # genetic distance matrix
  display_IBD(testMAF,TRUE,popinfo, "total") # the whole data set
  display_IBD(testMAF,list_ns, popinfo, "non-scripta")
  display_IBD(testMAF,list_hyN, popinfo, "hyN")
  display_IBD(testMAF,list_hyS,popinfo, "hyS")
  display_IBD(testMAF,list_hisp,popinfo, "hispanica")
  
  display_IBD(testMAF,c(list_hyN, list_hyS), popinfo, "hybrids")
  display_IBD(testMAF,c(list_ns, list_hisp), popinfo, "parents") 
  
  # Mantel test to check for auto-correlation of spatial distance and genetic distance
  # object of class 'dist' (done in the display_IBD function)

}

if (F) {

  # genetic manhattan distance matrix
  
  gDist.ns <- as.matrix(dist(testMAF[list_ns,], method="manhattan"))
  gDist.hisp <- as.matrix(dist(testMAF[list_hisp,], method="manhattan"))
  gDist.hyN <- as.matrix(dist(testMAF[list_hyN,], method="manhattan"))
  gDist.hyS <- as.matrix(dist(testMAF[list_hyS,], method="manhattan"))
  
  gDist.hyb <- as.matrix(dist(testMAF[c(list_hyN,list_hyS),], method="manhattan"))
  gDist.hyShisp <- as.matrix(dist(testMAF[c(list_hyS, list_hisp),], method="manhattan"))
  gDist.hyNns <- as.matrix(dist(testMAF[c(list_ns, list_hyN),], method="manhattan"))
  
  gDist.T <- as.matrix(dist(testMAF, method="manhattan"))
  
  # geographic Haversine distance matrix in km
 
  ds.ns <- data.frame(popinfo$lat, popinfo$lon)[list_ns,]
  ds.hisp <- data.frame(popinfo$lat, popinfo$lon)[list_hisp,]
  ds.hyN <- data.frame(popinfo$lat, popinfo$lon)[list_hyN,]
  ds.hyS <- data.frame(popinfo$lat, popinfo$lon)[list_hyS,]
  
  ds.hyb <- data.frame(popinfo$lat, popinfo$lon)[c(list_hyN,list_hyS),]
  ds.hyShisp <- data.frame(popinfo$lat, popinfo$lon)[c(list_hyS, list_hisp),]
  ds.hyNns <- data.frame(popinfo$lat, popinfo$lon)[c(list_ns, list_hyN),]
  
  ds.T <- data.frame(popinfo$lat, popinfo$lon)
  
  sDist.ns <-as.matrix(apply(ds.ns, 1, FUN=function(X) round(distHaversine(X, ds.ns)/1000, digits=3)))
  sDist.hyN <-as.matrix(apply(ds.hyN, 1, FUN=function(X) round(distHaversine(X, ds.hyN)/1000, digits=3)))
  sDist.hyS <-as.matrix(apply(ds.hyS, 1, FUN=function(X) round(distHaversine(X, ds.hyS)/1000, digits=3)))
  sDist.hisp <-as.matrix(apply(ds.hisp, 1, FUN=function(X) round(distHaversine(X, ds.hisp)/1000, digits=3)))
  
  sDist.hyb <-as.matrix(apply(ds.hyb, 1, FUN=function(X) round(distHaversine(X, ds.hyb)/1000, digits=3)))
  sDist.hyShisp <-as.matrix(apply(ds.hyShisp, 1, FUN=function(X) round(distHaversine(X, ds.hyShisp)/1000, digits=3)))
  sDist.hyNns <-as.matrix(apply(ds.hyNns, 1, FUN=function(X) round(distHaversine(X, ds.hyNns)/1000, digits=3)))
  
  sDist.T <-as.matrix(apply(ds.T, 1, FUN=function(X) round(distHaversine(X, ds.T)/1000, digits=3)))
  
  plot(gDist.T~sDist.T, xlim=c(0,180), ylim=c(0,700), pch='.', xlab="sDist in km", ylab="gDist", main="Isolation by distance", col="darkgrey")
  
  points(gDist.hyb~sDist.hyb, pch='.', col = "black" )
  points(gDist.ns~sDist.ns, pch='.', col = "darkgrey" )
  points(gDist.hyN~sDist.hyN, pch='.', col = "green" )
  points(gDist.hyS~sDist.hyS, pch='.', col = "orange" )
  points(gDist.hisp~sDist.hisp, pch='.', col = "darkgrey" )
  
  
  points(gDist.hyNns~sDist.hyNns, pch='.', col="gold")
  points(gDist.hyShisp~sDist.hyShisp, pch='.', col="black")
  
} # plot the IBD * jeannine 07-08

# pca * jeannine 12-Aug
if (F){
  # pca on pop
  list_def=seq(dim(data_pop2)[2]-addinfo_col_data)[!is.na(sapply(data_pop2[,(addinfo_col_data+1):dim(data_pop2)[2]], mean))]
  
  # number of included sites without NA calls:
  length(list_def)
  
  pca_2158loci=prcomp(data_pop2[1:51,addinfo_col_data+list_def],scaled=TRUE)
  summary(pca_2158loci)
  scores2 <- data.frame(data_pop2$Sp2[1:51], pca_2158loci$x[,1:3])
  
  data_pop2$PCA1=(scores2$PC1-min(scores2$PC1))/max(scores2$PC1-min(scores2$PC1))
  
  pc1.1 <- signif(summary(pca_2158loci)$importance[2,1]*100, digits=4)
  pc1.2 <- signif(summary(pca_2158loci)$importance[2,2]*100, digits=3)
  pc1.3 <- signif(summary(pca_2158loci)$importance[2,3]*100, digits=3)
  
  
  colpca=gsub("hyS","orange",gsub("hyN","green",gsub("ns","blue",gsub("hisp","red",data_pop2$Sp2[1:51]))))
  
  pcaplotPC12_pop <- ggplot(scores2,aes(x=scores2$PC1,y=scores2$PC2)) +
    geom_text(aes(x=scores2$PC1,y=scores2$PC2,label=data_pop2[,1], colour=data_pop2$Sp2[1:51])) + 
    scale_colour_manual(values=c("red", "green", "orange","blue")) + 
    labs(x= bquote("PC1 (" ~ .(pc1.1) ~ ") %"), y= bquote("PC2 (" ~ .(pc1.2) ~ ") %"))
  
  pcaplotPC13_pop <- ggplot(scores2,aes(x=scores2$PC1,y=scores2$PC3)) +
    geom_text(aes(x=scores2$PC1,y=scores2$PC3,label=data_pop2[,1], colour=data_pop2$Sp2[1:51])) + 
    scale_colour_manual(values=c("red", "green", "orange","blue")) + 
    labs(x= bquote("PC1 (" ~ .(pc1.1) ~ ") %"), y= bquote("PC2 (" ~ .(pc1.3) ~ ") %"))
  
  ## pca on the individuals 
  list_def2=seq(dim(data)[2]-addinfo_col_data)[!is.na(sapply(data[,(addinfo_col_data+1):dim(data)[2]], mean))]
  length(list_def2)
  
  pca_1125SNPs=prcomp(data[1:312,addinfo_col_data+list_def2],scaled=TRUE)
  
  scores3=data.frame(data$Sp2, pca_1125SNPs$x[,1:3])
  
  # plot the variances
  
  variances2 <- data.frame(variances=pca_1125SNPs$sdev**2, pcomp=1:length(pca_1125SNPs$sdev))
  
  varPlot2 <- ggplot(variances2, aes(pcomp, variances)) 
  varPlot2 <- varPlot2 + geom_bar(stat="identity", fill="gray") 
  varPlot2 <- varPlot2 + geom_line()
  # varPlot2
  
  pc2.1 <- signif(summary(pca_1125SNPs)$importance[2,1]*100, digits=4)
  pc2.2 <- signif(summary(pca_1125SNPs)$importance[2,2]*100, digits=3)
  pc2.3 <- signif(summary(pca_1125SNPs)$importance[2,3]*100, digits=3)
  
  pca1plot <- ggplot(scores3,aes(x=scores3$PC1,y=scores3$PC2)) +
    geom_text(aes(x=scores3$PC1,y=scores3$PC2,label=data$pop.ID,colour=data$Sp2)) + 
    scale_colour_manual(values=c("red", "green", "orange","blue")) +
    labs(x= bquote("PC1 (" ~ .(pc2.1) ~ ") %"), y= bquote("PC2 (" ~ .(pc2.2) ~ ") %"))
  
  pca2plot <- ggplot(scores3,aes(x=scores3$PC1,y=scores3$PC3)) +
    geom_text(aes(x=scores3$PC1,y=scores3$PC3,label=data$pop.ID,colour=data$Sp2)) + 
    scale_colour_manual(values=c("red", "green", "orange","blue")) +
    labs(x= bquote("PC1 (" ~ .(pc2.1) ~ ") %"), y= bquote("PC3 (" ~ .(pc2.3) ~ ") %"))
  
  # pca3plot <- ggplot(scores3,aes(x=scores3$PC1,y=scores3$PC2)) +
  # geom_point(aes(x=scores3$PC1,y=scores3$PC2, label=data[,1], col=data$Sp2)) +
  # scale_colour_manual(values=c("red", "green", "orange","blue")) +
  # labs(x= bquote("PC1 (" ~ .(pc2.1) ~ ") %"), y= bquote("PC2 (" ~ .(pc2.2) ~ ") %"))
  
  # pca4plot <- ggplot(scores3,aes(x=scores3$PC1,y=scores3$PC3)) +
  # geom_point(aes(x=scores3$PC1,y=scores3$PC3,label=data$pop.ID,colour=data$Sp2)) + 
  # scale_colour_manual(values=c("red", "green", "orange","blue")) +
  # labs(x= bquote("PC1 (" ~ .(pc2.1) ~ ") %"), y= bquote("PC3 (" ~ .(pc2.3) ~ ") %"))
  
  multiplot(pca1plot, pca2plot, pcaplotPC12_pop, pcaplotPC13_pop, cols=2)
}


# calculate Wright stats

if(T){
    save_fst=c()
    for (j in 1:2228){ 
      save_fst=c(save_fst,calculate_fst(data[,addinfo_col_data+j],data$Sp2))
    }
    
} # calculate Fst for total pop (1)  ** updated 11-Aug

if(T){
  save_fst_ns=c()
  for (j in 1:2228){ 
    save_fst_ns=c(save_fst_ns,calculate_fst(data[list_ns,addinfo_col_data+j],data$pop.ID[list_ns]))
  }
  
  save_fst_hyN=c()
  for (j in 1:2228){ 
    save_fst_hyN=c(save_fst_hyN,calculate_fst(data[list_hyN,addinfo_col_data+j],data$pop.ID[list_hyN]))
  }
  
  save_fst_hyS=c()
  for (j in 1:2228){ 
    save_fst_hyS=c(save_fst_hyS,calculate_fst(data[list_hyS,addinfo_col_data+j],data$pop.ID[list_hyS]))
  }
  
  save_fst_hisp=c()
  for (j in 1:2228){ 
    save_fst_hisp=c(save_fst_hisp,calculate_fst(data[list_hisp,addinfo_col_data+j],data$pop.ID[list_hisp]))
  }
  
  table_local_fst=rbind(save_fst_ns,save_fst_hyN,save_fst_hyS,save_fst_hisp)
  
} # calculate all Fst for the 4 pop ** updated 11-Aug output: table_local_fst

if(T){} # 

if(T){
  save_fst_target=save_fst[pos$set=="Target"]
  save_fst_ns_target=save_fst_ns[pos$set=="Target"]
  save_fst_hyN_target=save_fst_hyN[pos$set=="Target"]
  save_fst_hyS_target=save_fst_hyS[pos$set=="Target"]
  save_fst_hisp_target=save_fst_hisp[pos$set=="Target"]
} # quick subset of target snps for all fst, keep F because it's done the lazy way **** new 12/08 alex non lazy****

if (T){
  split.screen(c(3,2))
  screen(1)
  plot(save_fst)
  screen(2)
  hist(save_fst,breaks=100,freq=T)
  abline(v=mean(save_fst,na.rm=T),col="blue")
  abline(v=quantile(save_fst,p=.95,na.rm=T),col="red")
  split.screen(c(1,2),screen=3)
  screen(7)
  plot(save_fst_ns)
  screen(8)
  hist(save_fst_ns,breaks=100,freq=T)
  abline(v=mean(save_fst_ns,na.rm=T),col="blue")
  abline(v=quantile(save_fst_ns,p=.95,na.rm=T),col="red")
  split.screen(c(1,2),screen=4)
  screen(9)
  plot(save_fst_hisp)
  screen(10)
  hist(save_fst_hisp,breaks=100,freq=T)
  abline(v=mean(save_fst_hisp,na.rm=T),col="blue")
  abline(v=quantile(save_fst_hisp,p=.95,na.rm=T),col="red")
  split.screen(c(1,2),screen=5)
  screen(11)
  plot(save_fst_hyN)
  screen(12)
  hist(save_fst_hyN,breaks=100,freq=T)
  abline(v=mean(save_fst_hyN,na.rm=T),col="blue")
  abline(v=quantile(save_fst_hyN,p=.95,na.rm=T),col="red")
  split.screen(c(1,2),screen=6)
  screen(13)
  plot(save_fst_hyS)
  screen(14)
  hist(save_fst_hyS,breaks=100,freq=T)
  abline(v=mean(save_fst_hyS,na.rm=T),col="blue")
  abline(v=quantile(save_fst_hyS,p=.95,na.rm=T),col="red")
  readline()
  dev.off()
} # allow to plot Fst, its distribution in the total pop and the subdivide one

if (T){
  save_fis=c()
  save_fis_sign=c()
  for (i in 1:2228){
    example=data[,addinfo_col_data+i]
    new_vect=c(sum(example[list_ns]/2,na.rm=TRUE),sum(example[list_hyN]/2,na.rm=TRUE),sum(example[list_hyS]/2,na.rm=TRUE),sum(example[list_hisp]/2,na.rm=TRUE))
    het_obs=c(sum(example[list_ns]%%2,na.rm=TRUE),sum(example[list_hyN]%%2,na.rm=TRUE),sum(example[list_hyS]%%2,na.rm=TRUE),sum(example[list_hisp]%%2,na.rm=TRUE))
    pop_size=c(sum(!is.na(example[list_ns])),sum(!is.na(example[list_hyN])),sum(!is.na(example[list_hyS])),sum(!is.na(example[list_hisp])))
    Fis=1-het_obs/pop_size/(2*new_vect/pop_size*(1-new_vect/pop_size)) # Alexandre is an idiot
    Fis[is.nan(Fis)]<-NA
    save_fis=cbind(save_fis,Fis)
    sign_fis=c()
    for (k in 1:4){
      sign_fis=c(sign_fis,dist_het(pop_size[k],2*new_vect[k]/pop_size[k]*(1-new_vect[k]/pop_size[k]),het_obs[k],c(1422,1344,1370,1369)[k]))
    }
    save_fis_sign=cbind(save_fis_sign,sign_fis)
    #sign_fis[is.na(sign_fis)]<-FALSE
  } #calculate fis
  
  par(mfrow=c(4,2))
  for(i in 1:4){
    plot(save_fis[i,],main=paste("Fis",c("ns","hyN","hyS","hisp")[i]))
    hist(save_fis[i,],main=paste("Fis",c("ns","hyN","hyS","hisp")[i]))
    abline(v=mean(save_fis[i,]),col="blue")} # represent Fis
  
  par(mfrow=c(4,2))
  for(i in 1:4){
    plot(save_fis[i,pos$set=="Target"],main=paste("Fis",c("ns","hyN","hyS","hisp")[i]))
    hist(save_fis[i,pos$set=="Target"],main=paste("Fis",c("ns","hyN","hyS","hisp")[i]))
    abline(v=mean(save_fis[i,]),col="blue")} # represent Fis for the target SNPs
  
  sign_loc_fis=c(sum(save_fis_sign[1,],na.rm=T), sum(save_fis_sign[2,],na.rm=T), sum(save_fis_sign[3,],na.rm=T), sum(save_fis_sign[4,],na.rm=T))
  par(mfrow=c(4,2))
  for(i in 1:4){
    plot(save_fis[i,],main=paste("Fis",c("ns","hyN","hyS","hisp")[i]),col=as.factor(save_fis_sign[i,]))
    hist(save_fis[i,],main=paste("Fis",c("ns","hyN","hyS","hisp")[i],"\n number of sign. loci ",sign_loc_fis[i]),breaks=seq(-10,10)/10,ylim=c(0,150))
    par(new=T)
    hist(save_fis[i,save_fis_sign[i,]],col="red",main=paste("Fis",c("ns","hyN","hyS","hisp")[i],"\n number of sign. loci ",sign_loc_fis[i]),breaks=seq(-10,10)/10,ylim=c(0,150))
  } # represent Fis
  
} # calculate and represent Fis for all SNPs and the target snps ****COPIED 10-Aug*****

if (T){
  table_sign_fis_ns=cbind(colnames(data)[(addinfo_col_data+1):dim(data)[2]][save_fis_sign[1,]],save_fis[1,save_fis_sign[1,]],save_fst_ns[save_fis_sign[1,]])
  colnames(table_sign_fis_ns)=c('SNP_name','Fis','Fst')
  table_sign_fis_ns=table_sign_fis_ns[!is.na(table_sign_fis_ns[,1]),]
  fst_threshold_val=c()
  for (i in table_sign_fis_ns[,1]){
    y=data[,seq(dim(data)[2])[colnames(data)==i]]
    y=y[list_ns]
    fst_threshold_val=c(fst_threshold_val,calculate_random_fst(y,data$pop.ID[list_ns],length(table_sign_fis_ns[,1])))
  }
  table_sign_fis_ns=cbind(table_sign_fis_ns,fst_threshold_val)
  table_sign_fis_ns=table_sign_fis_ns[!((table_sign_fis_ns[,3]>table_sign_fis_ns[,4] )& table_sign_fis_ns[,2]>0),]
  
  
  table_sign_fis_hyN=cbind(colnames(data)[(addinfo_col_data+1):dim(data)[2]][save_fis_sign[2,]],save_fis[2,save_fis_sign[2,]],save_fst_hyN[save_fis_sign[2,]])
  colnames(table_sign_fis_hyN)=c('SNP_name','Fis','Fst')
  table_sign_fis_hyN=table_sign_fis_hyN[!is.na(table_sign_fis_hyN[,1]),]
  fst_threshold_val=c()
  for (i in table_sign_fis_hyN[,1]){
    y=data[,seq(dim(data)[2])[colnames(data)==i]]
    y=y[list_hyN]
    fst_threshold_val=c(fst_threshold_val,calculate_random_fst(y,data$pop.ID[list_hyN],length(table_sign_fis_hyN[,1])))
  }
  table_sign_fis_hyN=cbind(table_sign_fis_hyN,fst_threshold_val)
  table_sign_fis_hyN=table_sign_fis_hyN[!((table_sign_fis_hyN[,3]>table_sign_fis_hyN[,4] )& table_sign_fis_hyN[,2]>0),]
  
  
  table_sign_fis_hyS=cbind(colnames(data)[(addinfo_col_data+1):dim(data)[2]][save_fis_sign[3,]],save_fis[3,save_fis_sign[3,]],save_fst_hyS[save_fis_sign[3,]])
  colnames(table_sign_fis_hyS)=c('SNP_name','Fis','Fst')
  table_sign_fis_hyS=table_sign_fis_hyS[!is.na(table_sign_fis_hyS[,1]),]
  fst_threshold_val=c()
  for (i in table_sign_fis_hyS[,1]){
    y=data[,seq(dim(data)[2])[colnames(data)==i]]
    y=y[list_hyS]
    fst_threshold_val=c(fst_threshold_val,calculate_random_fst(y,data$pop.ID[list_hyS],length(table_sign_fis_hyS[,1])))
  }
  table_sign_fis_hyS=cbind(table_sign_fis_hyS,fst_threshold_val)
  table_sign_fis_hyS=table_sign_fis_hyS[!((table_sign_fis_hyS[,3]>table_sign_fis_hyS[,4] )& table_sign_fis_hyS[,2]>0),]
  
  
  table_sign_fis_hisp=cbind(colnames(data)[(addinfo_col_data+1):dim(data)[2]][save_fis_sign[4,]],save_fis[4,save_fis_sign[4,]],save_fst_hisp[save_fis_sign[4,]])
  colnames(table_sign_fis_hisp)=c('SNP_name','Fis','Fst')
  table_sign_fis_hisp=table_sign_fis_hisp[!is.na(table_sign_fis_hisp[,1]),]
  fst_threshold_val=c()
  for (i in table_sign_fis_hisp[,1]){
    y=data[,seq(dim(data)[2])[colnames(data)==i]]
    y=y[list_hisp]
    fst_threshold_val=c(fst_threshold_val,calculate_random_fst(y,data$pop.ID[list_hisp],length(table_sign_fis_hisp[,1])))
  }
  table_sign_fis_hisp=cbind(table_sign_fis_hisp,fst_threshold_val)
  table_sign_fis_hisp=table_sign_fis_hisp[!((table_sign_fis_hisp[,3]>table_sign_fis_hisp[,4] )& table_sign_fis_hisp[,2]>0),]
  
} # check the fis outliers for substrucuture ** updated 11-Aug

if (T){
  save_He=c()
  for (i in 1:2228){
    example=data[,addinfo_col_data+i]
    save_He=c(save_He,2*mean(example/2,na.rm=T)*(1-mean(example/2,na.rm=T)))
  }
  
  save_He_subpop=array(0,dim=c(4,2228))
  for (i in 1:2228){
    for (j in 1:4){
      example=data[get(c("list_ns","list_hyN","list_hyS","list_hisp")[j]),addinfo_col_data+i]
      save_He_subpop[j,i]=2*mean(example/2,na.rm=T)*(1-mean(example/2,na.rm=T))
    }
  } # new 13/08
  
  save_Ho_subpop=array(0,dim=c(5,2228))
  for (i in 1:2228){
    for (j in 1:5){
      example=data[get(c("list_ns","list_hyN","list_hyS","list_hisp","T")[j]),addinfo_col_data+i]
      save_Ho_subpop[j,i]=sum(example==1,na.rm=T)/sum(example*0+1,na.rm=T)
    }
  }  # new 13/08
  
  par(mfrow=c(2,2))
  plot(save_He,main="Expected Heterozygosity")
  hist(save_He, breaks=100)
  plot(save_He[pos$set=="Target"],main="Expected Heterozygosity for the target SNPs")
  hist(save_He[pos$set=="Target"], breaks=100)
  
  par(mfrow=c(1,2))
  plot(density(save_He_subpop[1,]),col="blue",ylim=c(0,16),xlim=c(0,.5), lwd=3)
  par(new=T)
  plot(density(save_He_subpop[2,]),col="green",ylim=c(0,16),xlim=c(0,.5), lwd=3)
  par(new=T)
  plot(density(save_He_subpop[3,]),col="orange",ylim=c(0,16),xlim=c(0,.5), lwd=3)
  par(new=T)
  plot(density(save_He_subpop[4,]),col="red",ylim=c(0,16),xlim=c(0,.5), lwd=3)
  
  plot(density(save_Ho_subpop[1,]),col="blue",ylim=c(0,16),xlim=c(0,1), lwd=3)
  par(new=T)
  plot(density(save_Ho_subpop[2,]),col="green",ylim=c(0,16),xlim=c(0,1), lwd=3)
  par(new=T)
  plot(density(save_Ho_subpop[3,]),col="orange",ylim=c(0,16),xlim=c(0,1), lwd=3)
  par(new=T)
  plot(density(save_Ho_subpop[4,]),col="red",ylim=c(0,16),xlim=c(0,1), lwd=3)
  
  
} # calculate and represent HET (expected heterozygosity for the total pop) **** new 13/08

if (F){
  # ________NA = 0__________
  #   pairwise_fst=c()
  #   counter=1
  #   for (pop1 in data_pop2$pop.ID[1:50]){
  #     val=rep(NA,counter)
  #     for (pop2 in data_pop2$pop.ID[(counter+1):51]){
  #       fst_gen=c()
  #       for (k in 11:2238){
  #         list_ind <- c(seq(312)[data$pop.ID == pop1 ],seq(312)[data$pop.ID==pop2])
  #         fst_gen=c(fst_gen, calculate_fst(data[list_ind,k],data$pop.ID[list_ind]))
  #       }
  #     fst_gen[is.na(fst_gen)]=0
  #     val= c(val, mean(fst_gen, na.rm=TRUE)) 
  #     }
  #   counter=counter+1
  #   pairwise_fst=rbind(pairwise_fst,val)
  #   }
  # 
  #   pairwise_fst2=rbind(pairwise_fst,rep(NA,counter))
  #   pairwise_fst2[is.na(pairwise_fst2)]=0
  #   t(pairwise_fst2)+pairwise_fst2
  # 
  #   colnames(pairwise_fst2)=data_pop2$pop.ID
  #   rownames(pairwise_fst2)=data_pop2$pop.ID
  
  # _______NA -> NA ___________
  pairwise_fst3=c()
  counter=1
  for (pop1 in data_pop2$pop.ID[1:50]){
    val=rep(NA,counter)
    for (pop2 in data_pop2$pop.ID[(counter+1):51]){
      fst_gen=c()
      for (k in 11:2238){
        list_ind <- c(seq(312)[data$pop.ID == pop1 ],seq(312)[data$pop.ID==pop2])
        fst_gen=c(fst_gen, calculate_fst(data[list_ind,k],data$pop.ID[list_ind]))
      }
      # fst_gen[is.na(fst_gen)]=0 # new
      val= c(val, mean(fst_gen, na.rm=TRUE)) 
    }
    counter=counter+1
    pairwise_fst3=rbind(pairwise_fst3,val)
  }
  
  pairwise_fst5=rbind(pairwise_fst3,rep(NA,counter))
  pairwise_fst5[is.na(pairwise_fst5)]=0
  t(pairwise_fst5)+pairwise_fst5
  
  colnames(pairwise_fst5)=data_pop2$pop.ID
  rownames(pairwise_fst5)=data_pop2$pop.ID
  # __________________
  
  
  # write.table(pairwise_fst2, file = "../pairwiseFst2.txt", sep = "\t", row.names = TRUE, col.names = TRUE )
  
  
  par(mfrow=c(1,2))
  plot(pairwise_fst2[1,], main="pairwise Fst estimates \n for sampling locations")
  for (i in 2:51){points(pairwise_fst2[i,])}
  abline(h=mean(pairwise_fst2), col = "blue")
  hist(pairwise_fst2)
  abline(v=mean(pairwise_fst2), col = "blue")
  
  
} # mean pairwise Fst per sampling location & plot * Aug 11 eve Jeannine

if (F){
  pairwise_fst_4pop=c()
  counter=1
  for (pop1 in unique(data_pop2$Sp2)[1:3]){
    val=rep(NA,counter)
    for (pop2 in unique(data_pop2$Sp2)[(counter+1):4]){
      fst_gen=c()
      for (k in 11:2238){
        list_ind <- c(seq(312)[data$Sp2 == pop1],seq(312)[data$Sp2 == pop2])
        fst_gen=c(fst_gen, calculate_fst(data[list_ind,k],data$Sp2[list_ind]))
      }
      val= c(val, mean(fst_gen, na.rm=TRUE))
    }
    counter=counter+1
    pairwise_fst_4pop=rbind(pairwise_fst_4pop,val)
  }
  
  pairwise_fst_4pop2=rbind(pairwise_fst_4pop,rep(NA,counter))
  pairwise_fst_4pop2[is.na(pairwise_fst_4pop2)]=0
  t(pairwise_fst_4pop2)+pairwise_fst_4pop2
  
  colnames(pairwise_fst_4pop2)=unique(data_pop2$Sp2)
  rownames(pairwise_fst_4pop2)=unique(data_pop2$Sp2)
  
} # mean pairwise Fst per population * Aug 10 Jeannine

if (F){
  pairwise_fst_4popT=c()
  counter=1
  for (pop1 in unique(data_pop2$Sp2)[1:3]){
    val=rep(NA,counter)
    for (pop2 in unique(data_pop2$Sp2)[(counter+1):4]){
      fst_gen=c()
      for (k in (11:2238)[pos$set=="Target"]){
        list_ind <- c(seq(312)[data$Sp2 == pop1],seq(312)[data$Sp2 == pop2])
        fst_gen=c(fst_gen, calculate_fst(data[list_ind,k],data$Sp2[list_ind]))
      }
      val= c(val, mean(fst_gen, na.rm=TRUE))
    }
    counter=counter+1
    pairwise_fst_4popT=rbind(pairwise_fst_4popT,val)
  }
  
  pairwise_fst_4pop2=rbind(pairwise_fst_4pop,rep(NA,counter))
  pairwise_fst_4pop2[is.na(pairwise_fst_4pop2)]=0
  t(pairwise_fst_4pop2)+pairwise_fst_4pop2
  
  colnames(pairwise_fst_4pop2)=unique(data_pop2$Sp2)
  rownames(pairwise_fst_4pop2)=unique(data_pop2$Sp2)
  
} # mean pairwise Fst per population * Aug 13 ALex for target SNPs only

if (F){
  counter=1
  fst_gen=c()
  for (pop1 in unique(data_pop2$Sp2)[1:3]){
    for (pop2 in unique(data_pop2$Sp2)[(counter+1):4]){
      val=c()
      for (k in 11:2238){
        list_ind <- c(seq(312)[data$Sp2 == pop1],seq(312)[data$Sp2 == pop2])
        val=c(val, calculate_fst(data[list_ind,k],data$Sp2[list_ind]))
      }
      fst_gen= rbind(fst_gen, val)
    }
    counter=counter+1
  }
  
  colnames(fst_gen)=colnames(data)[11:2238]
  rownames(fst_gen)=c("ns_hisp", "ns_hyN", "ns_hyS", "hisp_hyN", "hisp_hyS", "hyN_hyS")
  # signif(rowMeans(val, na.rm = TRUE), 2)
  
} # pairwise Fst per population for all loci ** Aug 11 Jeannine

if (F){
  split.screen(c(3,2))
  screen(1)
  plot(save_fst)
  screen(2)
  hist(save_fst,breaks=100,freq=T)
  abline(v=mean(save_fst,na.rm=T),col="blue")
  abline(v=quantile(save_fst,p=.95,na.rm=T),col="red")
  
  split.screen(c(1,2),screen=3)
  screen(7)
  plot(fst_gen[1,], main=paste(rownames(fst_gen)[1]))
  screen(8)
  hist(fst_gen[1,],breaks=100,freq=T)
  abline(v=mean(fst_gen[1,],na.rm=T),col="blue")
  abline(v=quantile(fst_gen[1,],p=.95,na.rm=T),col="red")
  
  split.screen(c(1,2),screen=4)
  screen(9)
  plot(fst_gen[2,], main=paste(rownames(fst_gen)[2]))
  screen(10)
  hist(fst_gen[2,],breaks=100,freq=T)
  abline(v=mean(fst_gen[2,],na.rm=T),col="blue")
  abline(v=quantile(fst_gen[2,],p=.95,na.rm=T),col="red")
  
  split.screen(c(1,2),screen=5)
  screen(11)
  plot(fst_gen[3,], main=paste(rownames(fst_gen)[3]))
  screen(12)
  hist(fst_gen[3,],breaks=100,freq=T)
  abline(v=mean(fst_gen[3,],na.rm=T),col="blue")
  abline(v=quantile(fst_gen[3,],p=.95,na.rm=T),col="red")
  
  split.screen(c(1,2),screen=6)
  screen(13)
  plot(fst_gen[4,], main=paste(rownames(fst_gen)[4]))
  screen(14)
  hist(fst_gen[4,],breaks=100,freq=T)
  abline(v=mean(fst_gen[4,],na.rm=T),col="blue")
  abline(v=quantile(fst_gen[4,],p=.95,na.rm=T),col="red")
  
  # part 2 
  split.screen(c(3,2))
  screen(1)
  plot(save_fst)
  screen(2)
  hist(save_fst,breaks=100,freq=T)
  abline(v=mean(save_fst,na.rm=T),col="blue")
  abline(v=quantile(save_fst,p=.95,na.rm=T),col="red")
  
  split.screen(c(1,2),screen=3)
  screen(7)
  plot(fst_gen[5,], main=paste(rownames(fst_gen)[5]))
  screen(8)
  hist(fst_gen[5,],breaks=100,freq=T)
  abline(v=mean(fst_gen[5,],na.rm=T),col="blue")
  abline(v=quantile(fst_gen[5,],p=.95,na.rm=T),col="red")
  
  split.screen(c(1,2),screen=4)
  screen(9)
  plot(fst_gen[6,], main=paste(rownames(fst_gen)[6]))
  screen(10)
  hist(fst_gen[6,],breaks=100,freq=T)
  abline(v=mean(fst_gen[6,],na.rm=T),col="blue")
  abline(v=quantile(fst_gen[6,],p=.95,na.rm=T),col="red")
  readline()
  dev.off()
} # plot Fst, its distribution in the total pop and the pairwise comparisons ** Aug 11 Jeannine

if (F){
  
  # distance matrix estimation (Fst is already the distance ...)
  # ManDist = dist(x = data[11:2238], method = "manhattan")
  # genetDist = dist.gene(x = data[11:2238], method = "pairwise")
  
  #   #NA -> 0
  #   pairwise_fst2.m=t(pairwise_fst2)
  #   colnames(pairwise_fst2.m)=data_pop2$pop.ID
  #   rownames(pairwise_fst2.m)=data_pop2$pop.ID
  #   
  #   FstDist = as.dist(pairwise_fst2.m)
  #    
  #NA -> NA 
  pairwise_fst5.m=t(pairwise_fst5)
  colnames(pairwise_fst5.m)=data_pop2$pop.ID
  rownames(pairwise_fst5.m)=data_pop2$pop.ID
  
  FstDistNA = as.dist(pairwise_fst5.m)
  
  #hierarchical clustering
  color <- data_pop2$Sp2
  (color <- ifelse(color == "hisp", "red",
                   ifelse(color == "hyS", "orange", 
                          ifelse(color == "hyN", "green", 
                                 ifelse(color == "ns" , "blue", NA)))))
  
  par(mfrow=c(1,2))
  method = "ward.D2"
  plot(as.phylo(hclust(FstDist, method)), tip.color= c(color), main="Fst NA -> 0")
  
  plot(as.phylo(hclust(FstDistNA, method)), tip.color= c(color), main="Fst NA -> NA")
  
  # neighbor-joining tree
  par(mfrow=c(1,2))
  
  plot(as.phylo(nj(pairwise_fst2.m)), tip.color= c(color), main="Fst NA -> 0")
  plot(as.phylo(nj(pairwise_fst5.m)), tip.color= c(color), main="Fst NA -> NA")
  
  # pcoa of the fst
  ## ape package
  col2r <- data$Sp2
  (col2r <- ifelse(col2r == "hisp", "red",
                   ifelse(col2r == "hyS", "orange", 
                          ifelse(col2r == "hyN", "green", 
                                 ifelse(col2r == "ns" , "blue", NA)))))
  
  # Fst distance pcoa
  pcoaFst = pcoa(FstDist)
  scores5 = data.frame(data_pop2$Sp2, pcoaFst$vectors[,1:3])
  pcoaA1 = pcoaFst$values$Cum_corr_eig[1]*100
  pcoaA2 = (pcoaFst$values$Cum_corr_eig[2]-pcoaFst$values$Cum_corr_eig[1])*100
  pcoaA3 = (pcoaFst$values$Cum_corr_eig[3]-pcoaFst$values$Cum_corr_eig[2])*100
  
  pcaplot_pop2 <- ggplot(scores5,aes(x=scores5$Axis.1,y=scores5$Axis.2)) + 
    geom_text(aes(x=scores5$Axis.1, y=scores5$Axis.2,label=data_pop2[,1], colour=color)) + 
    scale_colour_manual(values=c("blue", "green", "orange","red")) + 
    labs(x= bquote("Axis 1 (" ~ .(pcoaA1) ~ ") %"), y= bquote("Axis 2 (" ~ .(pcoaA2) ~ ") %")) + 
    ggtitle("Fst -> 0")
  plot(pcaplot_pop2)
  
  pcoaFstNA = pcoa(FstDistNA)
  scores6 = data.frame(data_pop2$Sp2, pcoaFstNA$vectors[,1:3])
  
  pcoaA1 = pcoaFstNA$values$Cum_corr_eig[1]*100
  pcoaA2 = (pcoaFstNA$values$Cum_corr_eig[2]-pcoaFstNA$values$Cum_corr_eig[1])*100
  pcoaA3 = (pcoaFstNA$values$Cum_corr_eig[3]-pcoaFstNA$values$Cum_corr_eig[2])*100
  
  pcaplot_pop3 <- ggplot(scores6,aes(x=scores6$Axis.1,y=scores6$Axis.2)) + 
    geom_text(aes(x=scores6$Axis.1, y=scores6$Axis.2,label=data_pop2[,1], colour=color)) + 
    scale_colour_manual(values=c("blue", "green", "orange","red")) + 
    labs(x= bquote("Axis 1 (" ~ .(pcoaA1) ~ ") %"), y= bquote("Axis 2 (" ~ .(pcoaA2) ~ ") %")) + 
    ggtitle("Fst -> NA")
  plot(pcaplot_pop3)
  
  multiplot(pcaplot_pop2, pcaplot_pop3)
  
  # heatmap of pairwise fst
  
  pop2_display=gsub("hisp","red",data_pop2$Sp2)
  pop2_display=gsub("ns","blue",pop2_display)
  pop2_display=gsub("hyN","green",pop2_display)
  pop2_display=gsub("hyS","orange",pop2_display)
  
  
  pop2_chlr_display <- data_pop3$CPF
  (pop2_chlr_display <- ifelse(pop2_chlr_display  == 1, "red",
                               ifelse(pop2_chlr_display > 0.5 , "orange", 
                                      ifelse(pop2_chlr_display  == 0 , "blue", 
                                             ifelse(pop2_chlr_display < 0.5, "green", "purple")))))
  
  heatmap.2(as.matrix(as.dist(pairwise_fst5.m)), trace="none", Rowv=as.dendrogram(hclust(as.dist(pairwise_fst5.m), method="ward.D2")), Colv=rev(as.dendrogram(hclust(as.dist(pairwise_fst5.m), method="ward.D2"))), density.info = "density", ColSideColors = pop2_display, RowSideColors=pop2_chlr_display, dendrogram = "row", main="Fst matrix using Ward's clustering")
  
  
  
} # playing with the Fst pairwise sampl.loc matrix * jeannine 12-Aug

#-------------------------------------------------------------------------------
# display potential cline and additional info per gene
displayallelefrequency(gene_info,data_pop2,pos$set, save_fst,save_fis ,save_He,table_local_fst,table_sign_fis_ns,table_sign_fis_hyN,table_sign_fis_hyS,table_sign_fis_hisp,sign_fst) # ** updated 21 Aug
while(T){dev.off()}

# calculate significant Fst (per permutations) # old
if (F){
  # We randomise the sample assignments into the taxon cluster 
  list_nsv2=seq(length(list_ns))
  list_hispv2=seq(length(list_ns)+1,length(list_ns)+length(list_hisp))
  list_hyNv2=seq(length(list_ns)+length(list_hisp)+1,length(list_ns)+length(list_hisp)+length(list_hyN))
  list_hySv2=seq(length(list_ns)+length(list_hisp)+length(list_hyN)+1,312)
  
  save_Array_fst=c()
  for (r in 1:1000){
    random_order=order(rnorm(312))
    data_shuffled=data[random_order,]
    list_Fst=c()
    for (i in 1:2228){
      example=data_shuffled[,addinfo_col_data+i]
      new_vect=c(mean(example[list_nsv2]/2,na.rm=TRUE),mean(example[list_hyNv2]/2,na.rm=TRUE),mean(example[list_hySv2]/2,na.rm=TRUE),mean(example[list_hispv2]/2,na.rm=TRUE),mean(example/2,na.rm=TRUE))
    vect_Size=c(sum(!is.na(example[list_nsv2])),sum(!is.na(example[list_hyNv2])),sum(!is.na(example[list_hySv2])),sum(!is.na(example[list_hispv2])))
    Fst=1-sum((2*new_vect*(1-new_vect))[1:4]*vect_Size/sum(vect_Size))/(2*new_vect*(1-new_vect))[5]
    list_Fst=c(list_Fst,Fst)
    }
    save_Array_fst=rbind(save_Array_fst,list_Fst)
  }
  # write.table(save_Array_fst, )
  # plot(save_Array_fst)
}

# calculate significant Fst (per permutations) 
if (T){
  # We randomise the sample assignments into the taxon cluster 
  threshold_fst=c()
  for (i in (addinfo_col_data+1):dim(data)[2]){
    threshold_fst=c(threshold_fst,calculate_random_fst(data[,i],data$Sp2,2228))
  }
  
  sign_fst=save_fst>threshold_fst
  # write.table(save_Array_fst, )
  # plot(save_Array_fst)
} # new 10/08   ***** alex*****


# plot Fst significant SNPs over each other
if (F){
  a=read.table("/Users/Jeannine/Dropbox/Parkplatz blue bells/1. Data/final_Snpcalls/save_Array_fst.tab.txt",header=TRUE,sep='\t',row.names=1)
  n=dim(a)
  a= a[,2:n[2]]
  n=dim(a)
  
  # Perform bonferoni correction for multiple testing
  threshold=c()
  for (i in 1:n[2]){
    test= quantile(a[,i],p=1-.05/n[2])
    threshold=c(threshold,test)
  }
  
  print(paste("Bonferroni correction, p=",.05/n[2]))
  plot(threshold, main=paste("p=",.05/n[2]))
  
  #fst=threshold+rnorm(n[2],0,10^-5)
  
  plot(save_fst,col=as.numeric(save_fst>threshold)+1, ylab="Fst", xlab="locus",main=paste("Fst outliers per locus \n p=",.05/n[2]," based on 1000 permutations"))
  
  
  data_popSignf <- data.frame(data_pop2[1:9],data_pop2[10:2237][save_fst>threshold])
  total = 9
  plot(data_popSignf$lat, data_popSignf[,total],xlim=c(41.8,43),ylim=c(0,1),col=color, ylab = "allele frequency", xlab = "latitude")
  for (i in 7:dim(data_popSignf)[2]-7){
    points(data_popSignf$lat, data_popSignf[,total],xlim=c(41.8,43),ylim=c(0,1),col=color)
    i = i + 1
    total = total + 1
  }
} # untested


# extract Fst for the flowering genes and incompatible? genes
if(F){
  new_save_fst=matrix(data=save_fst,nrow=1,ncol=2228,dimnames=list(c("Fst"),genelist))

  list_flower_gene=c("GSMUA_Achr1G01980","GSMUA_Achr5G18070","GSMUA_Achr6G33840","GSMUA_Achr7G15560")
  list_flowering_snp=c()
  for (i in 1:4){
    print(new_save_fst[grep(list_flower_gene[i], colnames(new_save_fst))])
    print(new_save_fst[grep(list_flower_gene[i], colnames(new_save_fst))][new_save_fst[grep(list_flower_gene[i], colnames(new_save_fst))]>threshold[grep(list_flower_gene[i], colnames(new_save_fst))]])
    list_flowering_snp=c(list_flowering_snp,grep(list_flower_gene[i], colnames(new_save_fst)))
  }
  
  
  x11()  
  plot(save_fst*as.numeric(save_fst>threshold),col=as.numeric(save_fst>threshold)+1, ylab="Fst", xlab="locus",main=paste("Fst outliers \n p=",.05/n[2]," based on 1000 permutations"))
  points(list_flowering_snp,save_fst[list_flowering_snp],pch=15,col="green")
  ###################
  list_organel_inter_gene=c("GSMUA_Achr10G00370", "GSMUA_Achr10G26390", "GSMUA_Achr11G15710", "GSMUA_Achr11G22050", "GSMUA_Achr2G08260", "GSMUA_Achr2G09300", "GSMUA_Achr2G13970", "GSMUA_Achr3G03330", "GSMUA_Achr3G04190", "GSMUA_Achr3G23800", "GSMUA_Achr4G18340", "GSMUA_Achr4G28160", "GSMUA_Achr6G01600", "GSMUA_Achr6G09110", "GSMUA_Achr6G25590", "GSMUA_Achr6G33840", "GSMUA_Achr6G35510", "GSMUA_Achr7G10900", "GSMUA_Achr7G16120", "GSMUA_Achr7G24470", "GSMUA_Achr8G28550", "GSMUA_AchrUn_randomG12680", "GSMUA_AchrUn_randomG18870")
  
  list_organel_snp=c()
  for (i in 1:23){
    print(new_save_fst[grep(list_organel_inter_gene[i], colnames(new_save_fst))])
    print(new_save_fst[grep(list_organel_inter_gene[i], colnames(new_save_fst))][new_save_fst[grep(list_organel_inter_gene[i], colnames(new_save_fst))]>threshold[grep(list_organel_inter_gene[i], colnames(new_save_fst))]])
    list_organel_snp=c(list_organel_snp,grep(list_organel_inter_gene[i], colnames(new_save_fst)))
  }
  
  
  x11()  
  plot(save_fst*as.numeric(save_fst>threshold),col=as.numeric(save_fst>threshold)+1, ylab="Fst", xlab="locus",main=paste("Fst outliers \n p=",.05/n[2]," based on 1000 permutations"))
  points(list_organel_snp,save_fst[list_organel_snp],pch=15,col="blue")
} # untested again

# removing singleton
if(F){
  #a=save_Array_fst
  
  # plot(save_fst*as.numeric(save_fst>threshold),col=as.numeric(save_fst>threshold)+1, ylab="Fst", xlab="locus",main=paste("Fst outliers \n p=",.05/n[2]," based on 1000 permutations"))
  
  # plot(save_fst*as.numeric(save_fst>threshold), main=paste("p=",.05/n[2],"based on 1000 permutations"))
  
  # ## to remove singleton -----
  # filter=c()
  # for (k in 7+1:7+n[2]){
  #   if (sum(data[,k],na.rm=TRUE)==1||sum(data[,k]-2,na.rm=TRUE)==-1){filter=c(filter,FALSE)} else {filter=c(filter,TRUE)}  
  # }
  # b=a[,filter]
  # n2=dim(b)
  # save_fst_without_singletons=save_fst[filter]
  # 
  # threshold2=c()
  # for (i in 1:n2[2]){
  #   test= quantile(b[,i],p=1-.05/n[2])
  #   threshold2=c(threshold2,test)
  # }
  # 
  # plot(save_fst_without_singletons,col=as.numeric(save_fst_without_singletons>threshold2)+1, main=paste("p=",.05/n[2],"based on 1000 permutations"))
  # plot(save_fst_without_singletons*as.numeric(save_fst_without_singletons>threshold2), main=paste("p=",.05/n[2],"based on 1000 permutations"))
}# untested

# test singleton on fis
if (T){
  filter=c()
  for (k in (addinfo_col_data+1):dim(data)[2]){
    if (sum(data[,k],na.rm=TRUE)==1){filter=c(filter,FALSE)} else {filter=c(filter,TRUE)}  
  }
  par(mfrow=c(2,2))
  for(i in 1:4){
    plot(save_fis[i,],main=paste("Fis",c("ns","hyN","hyS","hisp")[i]),col=as.factor(filter))
  } # represent Fis for the target SNPs
} # ****COPIED 10-AUG*****

if (F){
  genelistHighFst <- colnames(data)[11:2238][save_fst>threshold]
  
  uniq_genelistHighFst=array(0,c(1,2))
  for (i in 1:length(genelistHighFst)){
    nqme=unlist(strsplit(genelistHighFst[i],'.',fixed=TRUE))[1]
    if (uniq_genelistHighFst[dim(uniq_genelistHighFst)[1],1]==nqme){
      uniq_genelistHighFst[dim(uniq_genelistHighFst)[1],2]=as.numeric(uniq_genelistHighFst[dim(uniq_genelistHighFst)[1],2])+1
    } else{
      uniq_genelistHighFst=rbind(uniq_genelistHighFst,c(nqme, 1))
    }
  } 
  uniq_genelistHighFst=uniq_genelistHighFst[2:dim(uniq_genelistHighFst)[1],]
  
  genelist <- colnames(data)[8:2235]
  
  genelist_names=c()
  for (i in 1:length(genelist)){
    nqme=unlist(strsplit(genelist[i],'.',fixed=TRUE))[1]
    if (!sum(nqme==genelist_names)){
      genelist_names=c(genelist_names,nqme)
    }
  } 
  
  length(uniq_genelistHighFst[,1])
} # untested again

if (T){
  o1=unique(unlist(strsplit(table_sign_fis_ns[,1],'.',fixed=TRUE))[seq(1,2*dim(table_sign_fis_ns)[1],2)])
  o2=unique(unlist(strsplit(table_sign_fis_hyN[,1],'.',fixed=TRUE))[seq(1,2*dim(table_sign_fis_hyN)[1],2)])
  o3=unique(unlist(strsplit(table_sign_fis_hyS[,1],'.',fixed=TRUE))[seq(1,2*dim(table_sign_fis_hyS)[1],2)])
  o4=unique(unlist(strsplit(table_sign_fis_hisp[,1],'.',fixed=TRUE))[seq(1,2*dim(table_sign_fis_hisp)[1],2)])
  oo=unique(c(o1,o2,o3,o4))
  all_sign_gene=Reduce(intersect,list(o1,o2,o3,o4))
  for (i in oo)
  {
    display_specific_gene_frequency(i,gene_info,data_pop2,pos$set, save_fst,save_fis,save_He, table_local_fst,table_sign_fis_ns,table_sign_fis_hyN,table_sign_fis_hyS,table_sign_fis_hisp,sign_fst)
  }
  
  for (i in all_sign_gene)
  {
    display_specific_gene_frequency(i,gene_info,data_pop2,pos$set, save_fst,save_fis,save_He, table_local_fst,table_sign_fis_ns,table_sign_fis_hyN,table_sign_fis_hyS,table_sign_fis_hisp,sign_fst)
  }
  
  for (i in list_flower_gene)
  {
    display_specific_gene_frequency(i,gene_info,data_pop2,pos$set, save_fst,save_fis,save_He, table_local_fst,table_sign_fis_ns,table_sign_fis_hyN,table_sign_fis_hyS,table_sign_fis_hisp,sign_fst)
  }
  
  for (i in list_organel_inter_gene)
  {
    display_specific_gene_frequency(i,gene_info,data_pop2,pos$set, save_fst,save_fis,save_He, table_local_fst,table_sign_fis_ns,table_sign_fis_hyN,table_sign_fis_hyS,table_sign_fis_hisp,sign_fst)
  }
  
  list_gene_fst=unique(unlist(strsplit( pos$CHROM_POS[sign_fst],'|',fixed=TRUE))[seq(1,2*length(pos$CHROM_POS[sign_fst]),2)])
  
  for (i in list_gene_fst)
  {
    display_specific_gene_frequency(i,gene_info,data_pop2,pos$set, save_fst,save_fis,save_He, table_local_fst,table_sign_fis_ns,table_sign_fis_hyN,table_sign_fis_hyS,table_sign_fis_hisp,sign_fst)
  }
} # display gene with significant Fis in one pop, flower gene & org *** new 10/08 alex ***


if(T){
  # AMstr2= data.frame(data$Sp2, data$Sp1)
  # colnames(AMstr2)<-c("sep hy","merge hy")
  
  AMsamples2 = as.data.frame(na.omit(t(data[11:2238]))) # 1125 SNPs remaining
  
  for (i in 1:1125){
    empty_array=array(0,dim=c(2,51))
    for (j in 1:51){
      empty_array[2,j]=sum(AMsamples2[i,data$pop.ID==unique(data$pop.ID)[j]],na.rm=T)
      empty_array[1,j]=sum(AMsamples2[i,data$pop.ID==unique(data$pop.ID)[j]]*0+2,na.rm=T)-sum(AMsamples2[i,data$pop.ID==unique(data$pop.ID)[j]],na.rm=T)
    }
    rownames(empty_array)=c(0,1)
    assign(paste("amova",i,sep='_'),amova(as.data.frame(empty_array),dist=NULL,data.frame(data_pop2$Sp2,data_pop2$Sp1)))
  }
  
  variations=c()
  for (i in 1:1125){
    empty_stuff=c()  
    empty_stuff = t(eval(parse(text = paste('amova_',i,'$componentsofcovariance[,2]',sep = ''))))
    variations=rbind(variations, empty_stuff)
  }
  
  colnames(variations) <- rownames(amova_1$componentsofcovariance)
  
  x11()
  par(mfrow=c(5,1))
  plot(variations[variations[,4] >= 80,1], main="Between data_pop2.Sp1")
  plot(variations[variations[,4] >= 80,2], main="Between data_pop2.Sp2 Within data_pop2.Sp1")
  plot(variations[variations[,4] >= 80,3], main="Between samples Within data_pop2.Sp2")
  plot(variations[variations[,4] >= 80,4], main="Within samples")
  plot(variations[variations[,4] >= 80,5], main="Total")
  
  x11()
  par(mfrow=c(5,1))
  plot(variations[variations[,4] <= 80,1], main="Between data_pop2.Sp1")
  plot(variations[variations[,4] <= 80,2], main="Between data_pop2.Sp2 Within data_pop2.Sp1")
  plot(variations[variations[,4] <= 80,3], main="Between samples Within data_pop2.Sp2")
  plot(variations[variations[,4] <= 80,4], main="Within samples")
  plot(variations[variations[,4] <= 80,5], main="Total")
  
  x11()
  par(mfrow=c(5,1))
  hist(variations[,1], main="Between data_pop2.Sp1")
  hist(variations[,2], main="Between data_pop2.Sp2 Within data_pop2.Sp1")
  hist(variations[,3], main="Between samples Within data_pop2.Sp2")
  hist(variations[,4], main="Within samples")
  hist(variations[,5], main="Total")
  
  phis=c()
  for (i in 1:1125){
    empty_stuff=c()  
    empty_stuff = t(eval(parse(text = paste('amova_',i,'$statphi',sep = ''))))
    phis=rbind(phis, empty_stuff)
  }
  
  colnames(phis) <- rownames(amova_1$statphi)
  
  
  
}# AMOVA * 15-Aug


while(T){dev.off()}

if (T){
#generated from mathematica
  cumul_neut=read.table('NeutralCumulDist.csv',sep='\n')[,1]
  cumul_selected_snp=read.table('SelectedSNPDist.csv',sep='\n')[,1]

  save_result1=array(-1,dim=c(1,3))
  while(dim(save_result1)[1]<1001){
    temp=drawing_into_pdist(cumul_neut,312,312,1000)
    save_result1=rbind(save_result1,temp[(temp[,3]!=0 & temp[,3]!=1),])
  }
  save_result1=save_result1[2:1001,]
  arb_Fst1=1-(save_result1[,1]*(1-save_result1[,1])+save_result1[,2]*(1-save_result1[,2]))/(2*save_result1[,3]*(1-save_result1[,3]))
  
  save_result2=drawing_into_pdist(cumul_selected_snp,312,312,1000)
  sel_SNP_Fst=1-(save_result2[,1]*(1-save_result2[,1])+save_result2[,2]*(1-save_result2[,2]))/(2*save_result2[,3]*(1-save_result2[,3]))
  
  
  par(mfrow=c(3,2))
  hist(save_result1[,1],main="dist of allele freq in pop1")
  hist(save_result1[,2],main="dist of allele freq in pop2")
  hist(save_result1[,3],main="dist of allele freq in tot pop")
  hist(save_result2[,1],main="dist of allele freq in pop1")
  hist(save_result2[,2],main="dist of allele freq in pop2")
  hist(save_result2[,3],main="dist of allele freq in tot pop")

  
  par(mfrow=c(3,1))
  hist(arb_Fst1,main="dist. of Fst for random neutral polymorphic SNPs",col="blue",breaks=20,freq=FALSE,xlim=c(0,.02),ylim=c(0,600))
  hist(sel_SNP_Fst,main="dist of Fst for SNPs following the choice process",col="green",breaks=20,freq=FALSE,xlim=c(0,.02),ylim=c(0,600))
  hist(save_fst_target,main="dist Fst of the actual target SNPs",col="red",breaks=20,freq=FALSE)

 # ks test p-values <2.2e-16



 ks.test(sel_SNP_Fst,arb_Fst1)
 ks.test(save_fst_target,arb_Fst1)
 ks.test(save_fst_target,sel_SNP_Fst)
} #calculate dist of Fst due to the SNP design***** new 31/08

if (T){
  list_pure_ns=read.table('ns.list')
  list_pure_hisp=read.table('hisp.S.list')
  pure_ns=rep(F,312)
  pure_hisp=rep(F,312)
  for (i in 1:dim(data)[1]){
    pure_ns[i]=any(data[i,1]==list_pure_ns)
    pure_hisp[i]=any(data[i,1]==list_pure_hisp)
  }
    
  table_amp=read.table('amplicon_pos.txt')
  
  mean(table_amp$V3-table_amp$V2)
  
  new_col=c()
  for (j in 1:dim(table_amp)[1]){
    new_col=c(new_col,unlist(strsplit(paste(table_amp[j,1]),'|',fixed=T))[1])
  }
  table_amp=cbind(table_amp,new_col)
  
  exon_vec=c()
  for (j in seq(2228)){
    i=pos$CHROM_POS[j]
    extract=unlist(strsplit(i,'|',fixed=T))
    if (sum(extract[1]==table_amp[,6])==1){
      exon_vec=c(exon_vec,j)
    } else{
      temp_list=table_amp[extract[1]==table_amp[,6],]
      fgt=unlist(strsplit(extract[2],'_',fixed=T))
      if (fgt[length(fgt)]>= temp_list[1,2] & fgt[length(fgt)]<=temp_list[1,3]){
        exon_vec=c(exon_vec,j)
      }
    }
  }
  exon_vec=c(seq(10),(exon_vec+10))
  
  uniq_table_amp=table_amp[1,]
  for (k in 2:279){
    if (table_amp[k,6]!=table_amp[(k-1),6]){
      uniq_table_amp=rbind(uniq_table_amp,table_amp[k,])
    }
  } 
  mean(uniq_table_amp$V3-uniq_table_amp$V2)
  
  calculate_stat_2pop_adap_data(data[pure_ns,exon_vec],data[pure_hisp,exon_vec])

  c(calculate_stat_2pop_adap_data(data[data$Sp1=='ns',exon_vec],data[data$Sp1=='hy',exon_vec]),calculate_stat_2pop_adap_data(data[data$Sp1=='ns',exon_vec],data[data$Sp1=='hisp',exon_vec]),calculate_stat_2pop_adap_data(data[data$Sp1=='hy',exon_vec],data[data$Sp1=='hisp',exon_vec]))
 
  c(calculate_stat_2pop_adap_data(data[data$Sp2=='ns',exon_vec],data[data$Sp2=='hyN',exon_vec]),calculate_stat_2pop_adap_data(data[data$Sp2=='ns',exon_vec],data[data$Sp2=='hyS',exon_vec]),calculate_stat_2pop_adap_data(data[data$Sp2=='ns',exon_vec],data[data$Sp2=='hisp',exon_vec]),calculate_stat_2pop_adap_data(data[data$Sp2=='hyN',exon_vec],data[data$Sp2=='hyS',exon_vec]),calculate_stat_2pop_adap_data(data[data$Sp2=='hyN',exon_vec],data[data$Sp2=='hisp',exon_vec]),calculate_stat_2pop_adap_data(data[data$Sp2=='hyS',exon_vec],data[data$Sp2=='hisp',exon_vec]))
  
}# folded jsfs of the data*** 18/08 alex modified end august 2 pure species add toehr comb****




data=popinfo
counter=0
for ( i in 1:dim(nc_dat.df)[2]){
  if (mean(nc_dat.df[,i],na.rm=TRUE) <1){
    data=cbind(data,nc_dat.df[,i])
  }
  else{
    data=cbind(data,abs(nc_dat.df[,i]-2))
    counter=counter+1
  }
  if (mean(nc_dat.df[,i],na.rm=TRUE)==1){print(1)}
}

colnames(data) <- colnames(nc_dat.dfv2)

# represent data as matrix
testMAF <- as.matrix(data[,(addinfo_col_data+1):dim(data)[2]])

if(T){
  all_sign_SNP=Reduce(intersect,list(table_sign_fis_ns[,1],table_sign_fis_hyN[,1],table_sign_fis_hyS[,1],table_sign_fis_hisp[,1]))
  for (i in all_sign_SNP){
    print(paste(i,paste(hist(data[,seq(2238)[i==colnames(data)]],breaks=c(-.5,.5,1.5,2.5),plot=F)$counts,sep=" ",collapse=' '),sep='  ...  '))
  }
  
  multinomialCI(c(39,273,0),.05)
  
  genotypes=array(-1,dim=c(312,8))
  counter=1
  for (i in all_sign_SNP){
    genotypes[,counter]=data[,seq(2238)[i==colnames(data)]]
    counter=counter+1
  }
  cor(genotypes,use="complete.obs")
  heatmap(cor(genotypes,use="complete.obs"))
  cor(genotypes,use="pairwise.complete.obs")
  pseudo_gen=c()
  for (i in 1:312){
    pseudo_gen= c(pseudo_gen,paste(genotypes[i,],collapse='') )
  }
  barplot(table(pseudo_gen),las=2)
  
  chisquaretest=array(NA,dim=c(8,8))
  for (i in 1:7){
    for (j in (i+1):8){
      t1=genotypes[,i]+10*genotypes[,j]
      t1bis=c(t1,0,1,2,10,11,12,20,21,22)
      t2=table(as.factor(t1bis))-1
      t3=(table(c(genotypes[!is.na(t1),i],0,1,2))-1)/sum(!is.na(t1))
      t4=(table(c(genotypes[!is.na(t1),j],0,1,2))-1)/sum(!is.na(t1))
      t5=as.vector(as.matrix(t3) %*%t(as.matrix(t4)))
      t6=chisq.test(t2[t5>0],p=t5[t5>0])$p.value
      chisquaretest[i,j]=t6
    }
  }
  chisquaretest
  colnames(chisquaretest)=c("G1.1","G1.2","G1.3","G2","G3","G4","G5.1","G5.2")
  rownames(chisquaretest)=colnames(chisquaretest)
  par(cex.main=1,oma=c(0,1,1,0))
  heatmap(chisquaretest,Rowv = NA, Colv = NA, breaks = c(0, 10^-4 ,0.001785714, .05, 1), col = c("magenta","red", "orange", "yellow"),scale="none",main="pairwise associations of genotypes (chi square)",cex.main=.25)
  legend("topleft", legend = c("p<10^-4","p<.05/28 (bonf)", "p<.05","p>=.05"), col = c("magenta","red", "orange","yellow"), pch = 15)
}# SNP sign for all FIs new 24/08

if(T){
  potential_cline=c()
  for (i in 11:2238){
    if((mean(data[data$Sp2=='ns',i],na.rm=T)/2<.2 & mean(data[data$Sp2=='hisp',i],na.rm=T)/2>.8)|(mean(data[data$Sp2=='ns',i],na.rm=T)/2>.8 & mean(data[data$Sp2=='hisp',i],na.rm=T)/2<.2)){
      potential_cline=c(potential_cline,colnames(data)[i])
    }  
  }
  gene_list_pc=unique(unlist(strsplit(potential_cline,'.',fixed=T))[c(T,F)])
  for (i in gene_list_pc)
  {
    display_specific_gene_frequency(i,gene_info,data_pop2,pos$set, save_fst,save_fis,save_He, table_local_fst,table_sign_fis_ns,table_sign_fis_hyN,table_sign_fis_hyS,table_sign_fis_hisp,sign_fst)
  }
  save_cline=list()
  par(mfrow=c(6,5))
  for ( j in potential_cline){
    k=seq(2238)[colnames(data_pop2)==j]
    x=data_pop2$lat
    y=data_pop2[,k]
    mod2=nls(y ~ 1 / (1 + exp(c * (x-d))), start = list(c =10,d=42))
    h=summary(mod2)
    save_cline=c(save_cline,list(j,h$coefficients))
  }
  
  for (j in 1:length(potential_cline)){
    k=seq(2238)[colnames(data_pop2)==potential_cline[j]]
    display_gene_freq(1,k-10,gene_info,data_pop2,pos$set, save_fst,save_fis,save_He, table_local_fst,table_sign_fis_ns,table_sign_fis_hyN,table_sign_fis_hyS,table_sign_fis_hisp,sign_fst)
    par(new=T)
    plot( seq(41,43,.001),1/(1 + exp(save_cline[[2*j]][1] * ( seq(41,43,.001)-save_cline[[2*j]][2]))),xlim=c(41.8,43),ylim=c(0,1),type='l',col='black',xaxt='n',yaxt='n')
    legend("bottomleft",c(paste('c= ',signif(save_cline[[2*j]],4)[1],', d= ',signif(save_cline[[2*j]],4)[2])),col='white')
  }
  

  for (j in 1:length(potential_cline)){
    if (save_cline[[2*j]][1]>0){
      plot( seq(41,43,.001),1/(1 + exp(-save_cline[[2*j]][1] * ( seq(41,43,.001)-save_cline[[2*j]][2]))),xlim=c(41.8,43),ylim=c(0,1),type='l',col='red',xaxt='n',yaxt='n',ylab='',xlb='')
      par(new=T)    
    } else{
      plot( seq(41,43,.001),1/(1 + exp(save_cline[[2*j]][1] * ( seq(41,43,.001)-save_cline[[2*j]][2]))),xlim=c(41.8,43),ylim=c(0,1),type='l',col='blue',xaxt='n',yaxt='n',ylab='',xlab='')
      par(new=T)    
    }
  }
  
  width_cline=c()
  center_cline=c()
  for (j in 1:length(potential_cline)){
    width_cline=c(width_cline,save_cline[[2*j]][1])
    center_cline=c(center_cline,save_cline[[2*j]][2])
  }
  par(mfrow=c(1,2))
  hist(center_cline,breaks=10)
  hist(abs(width_cline),breaks=10)
  # 1° roughly 110 km, ie width of the clines goes from 4km to 25 km
  save_cline_tan=list()
  for ( j in potential_cline){
    k=seq(2238)[colnames(data_pop2)==j]
    x=data_pop2$lat
    y=data_pop2[,k]
    mod3=nls(y ~ .5*  (1 + tanh(2 * (x-d)/c)), start = list(c =0.2,d=42))
    h=summary(mod2)
    save_cline_tan=c(save_cline_tan,list(j,h$coefficients))
  }
  
} #new 2/09 modified 22/009 potential cline

if (T){
  #impossible data are unphased
  test_ld=cov(data[,11:2238],use="pairwise.complete.obs")
  test2_ld=cor(data[,11:2238],use="pairwise.complete.obs")
  image(seq(2228), seq(2228), test_ld, zlim = c(0, max(test_ld)), col = gray(seq(1, 0, -0.01)))
  image(seq(2228), seq(2228), test2_ld, zlim = c(-1, 1), col = gray(seq(1, 0, -0.01)))
}# LD per loci

if(T){
  snp_fis_all=Reduce(intersect,list(as.vector(table_sign_fis_ns[,1]),as.vector(table_sign_fis_hisp[,1]),as.vector(table_sign_fis_hyN[,1]),as.vector(table_sign_fis_hyS[,1])))
  par(mfrow=c(4,2))
  for (i in snp_fis_all)
  {
    k=seq(2238)[colnames(data_pop2)==i]
    display_gene_freq(1,k-10,gene_info,data_pop2,pos$set, save_fst,save_fis,save_He, table_local_fst,table_sign_fis_ns,table_sign_fis_hyN,table_sign_fis_hyS,table_sign_fis_hisp,sign_fst)
  }
}


if(T){
  t0=seq(216)[gene_info[,1]=="GSMUA_Achr7G07160"]
  t1=sum(as.integer(gene_info[1:(t0-1),2]))+11
  data[1,t1:(t1+as.integer(gene_info[t0,2])-1)]
}# new 28 sep, for haplotype