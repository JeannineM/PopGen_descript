
# load packages
if (T){
  library("fields")
  library("geosphere")
  library("ade4")
  library("RColorBrewer")
  library("gplots")
  library("ggplot2")
  library('ggbiplot')
  library('MASS')
  library('ape')
}

#import function
source('function-1.R')
x11()

# Map GT frequency per population as pie charts
  library("plotrix") 
  library("mapplots")
  library(rgdal)
  library(raster)
  srtm_35_04 <- getData(name='SRTM', download = TRUE,lon=-7, lat=42)
  srtm_36_04 <- getData(name='SRTM', download = TRUE,lon=-4, lat=42)
  srtm <- mosaic(srtm_35_04, srtm_36_04, fun=mean)
  
  plot(srtm, xlim=c(-7.6, -5.9), ylim=c(41.8,43.0), main="genotype proportions", col=rgb(colorRamp(c("chartreuse3","darkorange4", "darkgrey","whitesmoke"))((1:200)/200)/255)[20:200])
  for (i in 1:52){
    add.pie(z=c(as.numeric(hzdat$freq_SWA2[i]),as.numeric(hzdat$freq_ALT[i])),x=hzdat$lon[i], y=hzdat$lat[i], labels="", radius=0.015, col=c(alpha("blue", 1.0), alpha("red", 0.6)))
    }

 
 if(F){
   x11()
   par(mfrow=c(1,2)) 
   # chloroplast
   plot(srtm, xlim=c(-7.6, -5.9), ylim=c(41.8,43.0), col=rgb(colorRamp(c("chartreuse3","darkorange4", "darkgrey","whitesmoke"))((1:200)/200)/255)[20:200], legend=FALSE)
   
   for (i in 1:52){
     add.pie(z=c(as.numeric(hzdat$freq_SWA2[i]),as.numeric(hzdat$freq_ALT[i])),x=hzdat$lon[i], y=hzdat$lat[i], labels="", radius=0.015, col=c("blue","red"))
   }
   
   # nuclear
   plot(srtm, xlim=c(-7.6, -5.9), ylim=c(41.8,43.0), col=rgb(colorRamp(c("chartreuse3","darkorange4", "darkgrey","whitesmoke"))((1:200)/200)/255)[20:200], legend=FALSE)
   
   for (i in 1:51){
     add.pie(z=c(as.numeric(admixDat$V1[i]),as.numeric(admixDat$V2[i])),x=admixDat$lon[i], y=admixDat$lat[i], labels="", radius=0.015, col=c(col=c("blue","red")))
   }
   
   par(mfrow=c(1,1))
   par(old.par)
   
   plot(srtm, xlim=c(-7.6, -5.9), ylim=c(41.8,43.0), col=rgb(colorRamp(c("chartreuse3","darkorange4", "darkgrey","whitesmoke"))((1:200)/200)/255)[20:200])
   
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
  # write.table(we, "~/genlist_info_vector.tab.txt",row.names = TRUE, col.names = TRUE, sep = "\t")
  
  list_ns=seq(nrow(data))[data$Sp2=="ns"]
  list_hisp=seq(nrow(data))[data$Sp2=="hisp"]
  list_hyN=seq(nrow(data))[data$Sp2=="hyN"]
  list_hyS=seq(nrow(data))[data$Sp2=="hyS"]
  
} # variable: gene_info, list_hisp, list_ns, list_hyN, list_hyS

# represent data as matrix
testMAF <- as.matrix(data[,(addinfo_col_data+1):dim(data)[2]])

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
  data_pop3=data.frame(pop.ID=g$pop.ID,N, CPF, g[,3:9],f,row.names=g$pop.ID)

} # variable: data_pop2 + data_pop3

# Estimate mean MAF frequency

if (T) {
  
  MAFmean51 <- rowMeans(data_pop3[,10:dim(data_pop3)[2]], na.rm = TRUE)

  popsumtab <- data.frame(data_pop3[,1:9], MAFmean51)

} # ... per sampling location # popsumtab # 51 subsets
if (T) {
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

# Estimate observed heterozygosity: Ho_mean
if (T) {
  if (T){
save_Ho=c()
  for (j in ((addinfo_col_data+1):ncol(data))) {
    example=data[,j]
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

} # represent total observed heterozygosity
  
  if (T){
  g=c()
  f=c()
  for (m in levels(as.factor(data$pop.ID))){
      example=data[data$pop.ID==m,]
      het=colSums(example[,(addinfo_col_data+1):ncol(data)]==1,na.rm=T)/2
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
if (F) {

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
  
  # with manhattan distance
  if (T){
    d <- dist(testMAF, method="manhattan")
    hr <- hclust(d, method="ward.D2")
    heatmap.2(testMAF, Rowv=as.dendrogram(hr), Colv=FALSE, margins =c(1.09,1.09), breaks=col_breaks, col=my_palette, density.info = "density",ylab = NA , xlab = NA, dendrogram="row", labCol = NA, labRow=NA, RowSideColors=pop_display, scale="none", trace="none")

    # Target MAF SNPs
    x11()
    dim(targetSNPs)
    d_T <- dist(as.matrix(targetSNPs), method="manhattan")
    hr_T <- hclust(d_T, method="ward.D2") 
    
    heatmap.2(as.matrix(targetSNPs), Rowv=as.dendrogram(hr_T), Colv=FALSE, margins =c(1.09,1.09), breaks=col_breaks, col=my_palette, density.info = "density",ylab = NA , xlab = NA, dendrogram="row", labCol = NA, labRow=NA, RowSideColors=pop_display, scale="none", trace="none")
    
    # plotting trees of the hierarchical clustering
    library("ggdendro")
    dendr_hz <- dendro_data(hr, type = "rectangle")
    tree_hz <- ggplot() + 
      geom_segment(data=segment(dendr_hz), aes(x=x, y=y, xend=xend, yend=yend)) + 
      geom_text(data=label(dendr_hz), aes(x=x, y=y, label=label, hjust=0), size=3) +
      coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
      theme(axis.line.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank())
    x11(); tree_hz
    ## Citation
    # Gao X, Starmer J (2007) Human population structure detection via multilocus genotype clustering. BMC Genet 8, 34.
    
    
  }
  
}

# Isolation by distance
if (F) {
  # genetic distance matrix
  par(mfrow=c(3,2))
  display_IBD(testMAF,TRUE,popinfo, "total", "black") # the whole data set
  display_IBD(testMAF,list_ns, popinfo, "non-scripta", "blue")
  display_IBD(testMAF,list_hyN, popinfo, "hyN", "green4")
  display_IBD(testMAF,list_hyS,popinfo, "hyS", "orange")
  display_IBD(testMAF,list_hisp,popinfo, "hispanica", "red")
  
  display_IBD(testMAF,c(list_hyN, list_hyS), popinfo, "hybrids", "deeppink4")
  display_IBD(testMAF,c(list_ns, list_hisp), popinfo, "parents", "blueviolet") 
  
  # Mantel test to check for auto-correlation of spatial distance and genetic distance
  # object of class 'dist' (done in the display_IBD function)

par(mfrow=c(3,2))
  display_IBD(data_pop3[,((addinfo_col_data+1):ncol(data_pop3))],TRUE,data_pop3[,1:addinfo_col_data], "total", "black")
  display_IBD(data_pop3[,((addinfo_col_data+1):ncol(data_pop3))],data_pop3$Sp2 == "ns",data_pop3[,1:addinfo_col_data], "non-scripta", "blue")
  display_IBD(data_pop3[,((addinfo_col_data+1):ncol(data_pop3))],data_pop3$Sp2 == "hyN",data_pop3[,1:addinfo_col_data], "hyN", "green")
  display_IBD(data_pop3[,((addinfo_col_data+1):ncol(data_pop3))],data_pop3$Sp2 == "hyS",data_pop3[,1:addinfo_col_data], "hyS", "orange")
  
}

# pca
if (T) {
  # pca on pop
  list_def=seq(dim(data_pop2)[2]-addinfo_col_data)[!is.na(sapply(data_pop2[,(addinfo_col_data+1):dim(data_pop2)[2]], mean))]
  
  # number of included sites without NA calls:
  length(list_def)
  
  pca_2158loci=prcomp(data_pop2[1:nrow(data_pop2),addinfo_col_data+list_def],scaled=TRUE)
  summary(pca_2158loci)
  scores2 <- data.frame(data_pop2$Sp2[1:nrow(data_pop2)], pca_2158loci$x[,1:3])
  
  pc1.1 <- signif(summary(pca_2158loci)$importance[2,1]*100, digits=4)
  pc1.2 <- signif(summary(pca_2158loci)$importance[2,2]*100, digits=3)
  pc1.3 <- signif(summary(pca_2158loci)$importance[2,3]*100, digits=3)
  
  
  colpca=gsub("hyS","orange",gsub("hyN","green",gsub("ns","blue",gsub("hisp","red",data_pop2$Sp2))))
  
  pcaplotPC12_pop <- ggplot(scores2,aes(x=scores2$PC1,y=scores2$PC2)) +
    geom_text(aes(x=scores2$PC1,y=scores2$PC2,label=data_pop2[,1], colour=data_pop2$Sp2)) + 
    scale_colour_manual(values=c("red", "green", "orange","blue")) + 
    labs(x= bquote("PC1 (" ~ .(pc1.1) ~ ") %"), y= bquote("PC2 (" ~ .(pc1.2) ~ ") %"))
  
  pcaplotPC13_pop <- ggplot(scores2,aes(x=scores2$PC1,y=scores2$PC3)) +
    geom_text(aes(x=scores2$PC1,y=scores2$PC3,label=data_pop2[,1], colour=data_pop2$Sp2)) + 
    scale_colour_manual(values=c("red", "green", "orange","blue")) + 
    labs(x= bquote("PC1 (" ~ .(pc1.1) ~ ") %"), y= bquote("PC2 (" ~ .(pc1.3) ~ ") %"))
  
  ## pca on the individuals 
  list_def2=seq(dim(data)[2]-addinfo_col_data)[!is.na(sapply(data[,(addinfo_col_data+1):dim(data)[2]], mean))]
  length(list_def2) # 1120
  
  pca_1125SNPs=prcomp(data[,addinfo_col_data+list_def2],scaled=TRUE)
  
  scores3=data.frame(data$Sp2, pca_1125SNPs$x[,1:10])
  
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
}

# AMOVA (locus-by-locus and phi-stats) * June-2016

if(F){
  # hierarchical clustering using amova of non-NA data
  AM304samples_subSp3 <- amova(samples = as.data.frame(na.omit(t(data[,11:ncol(data)]))),
                               distances = dist(as.data.frame(na.omit(t(data[,11:ncol(data)])))), 
                               structures = data.frame(Pop=as.factor(data$pop.ID), Sp3 = data$Sp1))
  
  AM304samples_subSp4 <- amova(samples =  as.data.frame(na.omit(t(data[,11:ncol(data)]))),
                               distances =  dist(as.data.frame(na.omit(t(data[,11:ncol(data)])))), 
                               structures = data.frame(Pop=as.factor(data$pop.ID), Sp4 = data$Sp2))
  
  AM304samples_subTargSp4 <- amova(samples =  as.data.frame(na.omit(t(data[,(which(pos$set == "Target") + 10)]))),
                                   distances =  dist(as.data.frame(na.omit(t(data[,(which(pos$set == "Target") + 10)])))), 
                                   structures = data.frame(Pop=as.factor(data$pop.ID), Sp4 = data$Sp2))
  
  AM304samples_subTargSp3  <- amova(samples =  as.data.frame(na.omit(t(data[,(which(pos$set == "Target") + 10)]))),
                                    distances =  dist(as.data.frame(na.omit(t(data[,(which(pos$set == "Target") + 10)])))), 
                                    structures = data.frame(Pop=as.factor(data$pop.ID), Sp3 = data$Sp1))
  
  set.seed(1999)
  sign_AM304samples_subSp4 <- randtest(AM304samples_subSp4, nrepet = 100) # very slow
 
  # locus-by-locus AMOVA accounting for either 3 or 4 demes
  
  AMarraySp3=c()
  for (i in 11:(ncol(data) - addinfo_col_data)){
    Locus = na.omit(data.frame(Loc=data[,i], Pop=as.factor(data$pop.ID), Sp3 = data$Sp1, Sp4 = data$Sp2))
    AMovaLocus <- amova(samples = as.data.frame(t(cbind(Locus$Loc, 2-Locus$Loc))), 
                        distances = NULL,
                        structures = as.data.frame(Locus[,2:3]))
    AMarraySp3=rbind(AMarraySp3, c(t(AMovaLocus$statphi), t(AMovaLocus$componentsofcovariance$Sigma), t(AMovaLocus$componentsofcovariance$`%`)))
  }
  
  AMarraySp3_phi <- as.data.frame(AMarraySp3, row.names = colnames(data[11:(ncol(data) - addinfo_col_data)]))
  colnames(AMarraySp3_phi) <- c("Phi-samples-total", "Phi-samples-Pop", "Phi-Pop-Sp3", "Phi-Sp3-total", "Sig-btw-Sp3", "Sig-btw-Pop-win-Sp3", "Sig-btw-Ind-win-Pop", "Sig-win-Ind", "Sig-total", "Pct-btw-Sp3", "Pct-btw-Pop-win-Sp3", "Pct-btw-Ind-win-Pop", "Pct-win-Ind", "Pct-Total")
  
  par(mfrow=c(2,2))
  for (i in 1:4) {
    plot(AMarraySp3_phi[,i], cex=0.5, main=colnames(AMarraySp3_phi)[i], xlab="Snps", ylab="phi-stat")
    abline(h=mean(AMarraySp3_phi[,i]), col="red")
    abline(h=quantile(AMarraySp3_phi[,i])[4], col="blue")
  }
  
  AMarraySp4=c()
  for (i in 11:(ncol(data) - addinfo_col_data)){
    Locus = na.omit(data.frame(Loc=data[,i], Pop=as.factor(data$pop.ID), Sp4 = data$Sp2))
    AMovaLocus <- amova(samples = as.data.frame(t(cbind(Locus$Loc, 2-Locus$Loc))), 
                        distances = NULL,
                        structures = as.data.frame(Locus[,2:3]))
    AMarraySp4=rbind(AMarraySp4, c(t(AMovaLocus$statphi), t(AMovaLocus$componentsofcovariance$Sigma), t(AMovaLocus$componentsofcovariance$`%`)))
  }
  
  AMarraySp4_phi <- as.data.frame(AMarraySp4, row.names = colnames(data[11:(ncol(data) - addinfo_col_data)]))
  colnames(AMarraySp4_phi) <- c("Phi-samples-total", "Phi-samples-Pop", "Phi-Pop-Sp4", "Phi-Sp4-total", "Sig-btw-Sp4", "Sig-btw-Pop-win-Sp4", "Sig-btw-Ind-win-Pop", "Sig-win-Ind", "Sig-total", "Pct-btw-Sp4", "Pct-btw-Pop-win-Sp4", "Pct-btw-Ind-win-Pop", "Pct-win-Ind", "Pct-Total")
  
  par(mfrow=c(2,2))
  for (i in 1:4) {
    plot(AMarraySp4_phi[,i], cex=0.5, main=colnames(AMarraySp4_phi)[i], xlab="Snps", ylab="phi-stat")
    abline(h=mean(AMarraySp4_phi[,i]), col="red")
    abline(h=quantile(AMarraySp4_phi[,i])[4], col="blue")
  }
  
  plot(density(AMarraySp3_phi$`Pct-btw-Sp3`), col="green", type = "l", xlim=c(-3, 75), 
       main="% of variance between Sp3 (green) and Sp4 (blue)")
  
  points(density(AMarraySp4_phi$`Pct-btw-Sp4`), col="blue", type = "l")
  
}# AMOVA

# calculate Wright stats

if (T) {
    save_fst=c()
    for (j in 1:(ncol(data) - addinfo_col_data)){ 
      save_fst=c(save_fst,calculate_fst(data[,addinfo_col_data+j],data$Sp2))
    }
    
} # calculate Fst for total pop (1)  ** updated 11-Aug

if (T) {
  save_fst_ns=c()
  for (j in 1:(ncol(data) - addinfo_col_data)){ 
    save_fst_ns=c(save_fst_ns,calculate_fst(data[list_ns,addinfo_col_data+j],data$pop.ID[list_ns]))
  }
  
  save_fst_hyN=c()
  for (j in 1:(ncol(data) - addinfo_col_data)){ 
    save_fst_hyN=c(save_fst_hyN,calculate_fst(data[list_hyN,addinfo_col_data+j],data$pop.ID[list_hyN]))
  }
  
  save_fst_hyS=c()
  for (j in 1:(ncol(data) - addinfo_col_data)){ 
    save_fst_hyS=c(save_fst_hyS,calculate_fst(data[list_hyS,addinfo_col_data+j],data$pop.ID[list_hyS]))
  }
  
  save_fst_hisp=c()
  for (j in 1:(ncol(data) - addinfo_col_data)){ 
    save_fst_hisp=c(save_fst_hisp,calculate_fst(data[list_hisp,addinfo_col_data+j],data$pop.ID[list_hisp]))
  }
  
  table_local_fst=rbind(save_fst_ns,save_fst_hyN,save_fst_hyS,save_fst_hisp)
  
} # calculate all Fst for the 4 pop; output: table_local_fst

if (T) {
  save_fst_target=save_fst[pos$set=="Target"]
  save_fst_ns_target=save_fst_ns[pos$set=="Target"]
  save_fst_hyN_target=save_fst_hyN[pos$set=="Target"]
  save_fst_hyS_target=save_fst_hyS[pos$set=="Target"]
  save_fst_hisp_target=save_fst_hisp[pos$set=="Target"]
} # quick subset of target snps for all fst

if (T) {
  x11()
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

  loci_ns   <- sum(colSums(data[list_ns,10:ncol(data)], na.rm = TRUE) != 0)
  loci_hyN  <- sum(colSums(data[list_hyN,10:ncol(data)], na.rm = TRUE) != 0)
  loci_hyS  <- sum(colSums(data[list_hyS,10:ncol(data)], na.rm = TRUE) != 0)
  loci_hisp <- sum(colSums(data[list_hisp,10:ncol(data)], na.rm = TRUE) != 0)
  
  for (i in 1:(ncol(data) - addinfo_col_data)){
    example=data[,addinfo_col_data+i]
    new_vect=c(sum(example[list_ns]/2,na.rm=TRUE),sum(example[list_hyN]/2,na.rm=TRUE),sum(example[list_hyS]/2,na.rm=TRUE),sum(example[list_hisp]/2,na.rm=TRUE))
    het_obs=c(sum(example[list_ns]%%2,na.rm=TRUE),sum(example[list_hyN]%%2,na.rm=TRUE),sum(example[list_hyS]%%2,na.rm=TRUE),sum(example[list_hisp]%%2,na.rm=TRUE))
    pop_size=c(sum(!is.na(example[list_ns])),sum(!is.na(example[list_hyN])),sum(!is.na(example[list_hyS])),sum(!is.na(example[list_hisp])))
    Fis=1-het_obs/pop_size/(2*new_vect/pop_size*(1-new_vect/pop_size)) # Alexandre is an idiot
    Fis[is.nan(Fis)]<-NA
    save_fis=cbind(save_fis,Fis)
    sign_fis=c()
    for (k in 1:4){
      sign_fis=c(sign_fis,dist_het(pop_size[k],2*new_vect[k]/pop_size[k]*(1-new_vect[k]/pop_size[k]),het_obs[k],c(loci_ns,loci_hyN,loci_hyS,loci_hisp)[k]))
    }
    save_fis_sign=cbind(save_fis_sign,sign_fis)
    rownames(save_fis_sign) <- c("ns", "hyN", "hyS", "hisp")
    rownames(save_fis) <- c("ns", "hyN", "hyS", "hisp")
    #sign_fis[is.na(sign_fis)]<-FALSE
  } #calculate fis
  
  x11()
  par(mfrow=c(4,2))
  for(i in 1:4){
    plot(save_fis[i,],main=paste("Fis",c("ns","hyN","hyS","hisp")[i]))
    hist(save_fis[i,],main=paste("Fis",c("ns","hyN","hyS","hisp")[i]))
    abline(v=mean(save_fis[i,]),col="blue")} # represent Fis
  
  x11()
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
  
} # calculate and represent Fis for all SNPs and the target snps

if (F){
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
  
} # check the fis outliers for substructure * very slow

if (T){
  save_He=c()
  for (i in 1:(ncol(data) - addinfo_col_data)){
    example=data[,addinfo_col_data+i]
    save_He=c(save_He,2*mean(example/2,na.rm=T)*(1-mean(example/2,na.rm=T)))
  }
  
  save_He_subpop=array(0,dim=c(4,(ncol(data) - addinfo_col_data)))
  for (i in 1:(ncol(data) - addinfo_col_data)){
    for (j in 1:4){
      example=data[get(c("list_ns","list_hyN","list_hyS","list_hisp")[j]),addinfo_col_data+i]
      save_He_subpop[j,i]=2*mean(example/2,na.rm=T)*(1-mean(example/2,na.rm=T))
    }
  } # estimate expected heterozygosity
  
  save_Ho_subpop=array(0,dim=c(5,(ncol(data) - addinfo_col_data)))
  for (i in 1:(ncol(data) - addinfo_col_data)){
    for (j in 1:5){
      example=data[get(c("list_ns","list_hyN","list_hyS","list_hisp","T")[j]),addinfo_col_data+i]
      save_Ho_subpop[j,i]=sum(example==1,na.rm=T)/sum(example*0+1,na.rm=T)
    }
  }  # estimate observed heterozygosity
  
  x11()
  par(mfrow=c(2,2))
  plot(save_He,main="Expected Heterozygosity")
  hist(save_He, breaks=100)
  plot(save_He[pos$set=="Target"],main="Expected Heterozygosity for the target SNPs")
  hist(save_He[pos$set=="Target"], breaks=100)
  
  x11()
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
  
  
} # calculate and represent HET (expected heterozygosity for the total pop)

if (F){

  # _______NA -> NA ___________
  pairwise_fst3=c()
  counter=1
  for (pop1 in data_pop2$pop.ID[1:(nrow(data_pop2)-1)]){
    val=rep(NA,counter)
    for (pop2 in data_pop2$pop.ID[(counter+1):nrow(data_pop2)]){
      fst_gen=c()
      for (k in ((1+addinfo_col_data):ncol(data))){
        list_ind <- c(seq(nrow(data))[data$pop.ID == pop1],seq(nrow(data))[data$pop.ID==pop2])
        fst_gen=c(fst_gen, calculate_fst(data[list_ind,k],data$pop.ID[list_ind]))
      }
  
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

  # plot
  par(mfrow=c(1,1))
  plot(pairwise_fst5[1,], main="pairwise Fst estimates \n for sampling locations", las=2, col="white", axes=FALSE, ylab="distance", xlab="sampling location")
    axis(side = 1, labels=rownames(pairwise_fst5), at = c(seq(1:nrow(data_pop2))), las=2 )
  
  for (i in 1:nrow(data_pop2)){points(pairwise_fst5[i,][pairwise_fst5[i,] >0])}
  abline(h=mean(pairwise_fst5[pairwise_fst5 > 0]), col = "blue")
  
  hist(pairwise_fst5[pairwise_fst5 > 0], breaks = nrow(data_pop2)) 
  abline(v=mean(pairwise_fst5), col = "blue")
  abline(v=mean(pairwise_fst5[pairwise_fst5 > 0]), col = "green")
  
} # mean pairwise Fst per sampling location & plot * very slow

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
  
} # mean pairwise Fst per 4 populations

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
  
} # mean pairwise Fst per population for target SNPs only

if (T){
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
  
} # pairwise Fst per population for all loci

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
} # plot Fst, its distribution in the total pop and the pairwise comparisons

if (F){

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
  
  # isolation by distance using Fst
  par(mfrow=c(3,2))
  display_IBD_Fst(pairwise_fst5, data_pop2$Sp2 == "hyS" | data_pop2$Sp2 == "hyN", data_pop2, "hybrids", "deeppink2")
  display_IBD_Fst(pairwise_fst5, data_pop2$Sp2 == "hisp" | data_pop2$Sp2 == "ns", data_pop2, "parents", "gray")
  display_IBD_Fst(pairwise_fst5, data_pop2$Sp2 == "hyN", data_pop2, "hyN", "green")
  display_IBD_Fst(pairwise_fst5, data_pop2$Sp2 == "ns", data_pop2, "non-scripta", "blue")
  display_IBD_Fst(pairwise_fst5, data_pop2$Sp2 == "hyS", data_pop2, "hyS", "orange")
  display_IBD_Fst(pairwise_fst5, data_pop2$Sp2 == "hisp", data_pop2, "hispanica", "red")
  
} # playing with the Fst pairwise sampl.loc matrix

# calculate significant Fst (per permutations) 
if (F){
  # We randomise the sample assignments into the taxon cluster 
  threshold_fst=c()
  for (i in (addinfo_col_data+1):dim(data)[2]){
    threshold_fst=c(threshold_fst,calculate_random_fst(data[,i],data$Sp2,(ncol(data) - addinfo_col_data)))
  }
  
  sign_fst=save_fst>threshold_fst
  # write.table(save_Array_fst, )
  # plot(save_Array_fst)
} # calculate significant Fst (per permutations) 


