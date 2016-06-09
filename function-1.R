# isolation by D

display_IBD=function(dat,subset,popinf, main, colour){
  # UPDATE: use Fst as genetic divergence between sampling sites vs position
  d2 <- dist(dat[subset,], method="manhattan")
  gDist <- as.matrix(d2)
  #  rownames(gDist) <- row.names(popinfo[subset])
  #  colnames(gDist) <- row.names(popinfo[subset])
  
  # geographic Haversine distance matrix in km
  df.x2 <- data.frame(popinf$lat, popinf$lon)[subset,]
  ma.hav<-apply(df.x2, 1, FUN=function(X) round(distHaversine(X, df.x2)/1000, digits=3))
  #ma.hav[ma.hav == 0] <- Inf
  sDist <- as.matrix(ma.hav)
  plot(gDist~sDist, pch=19, main=paste(main), col=colour)
  test=mantel.rtest(as.dist(d2), as.dist(ma.hav), nrepet=1000)
  print(test)
}

# xlim=c(0,180), ylim=c(0,700),
# plot function
if (F){
  # displayallelefrequency=function(gene_info,data_pop2,targetSNP,fst_values,fst1,fst2,fst3,fst4,He){
  #   color <- data_pop2$Sp2
  #   (color <- ifelse(color == "hisp", "red",
  #                    ifelse(color == "hyS", "orange", 
  #                           ifelse(color == "hyN", "green", 
  #                                  ifelse(color == "ns" , "blue", NA)))))
  #   total=0
  #   test=expression(plot(data_pop2$lat, data_pop2[,addinfo_col_data-1+total+k],xlim=c(41.8,43),ylim=c(0,1),col=color, pch=16,main=paste(pos$CHROM_POS[total+k]),col.main=colortitle, ylab = "allele frequency",xlab = paste("latitude\n Fst=" ,signif(fst_values[total+k],3),"he=", signif(He[total+k],3), " Fis in pop hisp",signif(fst4[total+k],3),"hyS",signif(fst3[total+k],3),"hyN",signif(fst2[total+k],3),"ns",signif(fst1[total+k],3))))
  #   for (i in 2:dim(gene_info)[1]){
  #     x11()
  #     n=as.numeric(gene_info[i,2])
  #     if (n<=12){
  #       par(mfrow=switch(n,c(1,1),c(2,1),c(2,2),c(2,2),c(3,2),c(3,2),c(3,3),c(3,3),c(3,3),c(4,3),c(4,3),c(4,3)))
  #       for (k in 1:n) {
  #         if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'} 
  #         eval(test)
  #       }
  #     }
  #     if (n>12 & n<=24){
  #       par(mfrow=c(4,3))
  #       for (k in 1:12) {
  #         if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'} 
  #         eval(test)
  #       }
  #       x11()
  #       par(mfrow=switch(n-12,c(1,1),c(2,1),c(2,2),c(2,2),c(3,2),c(3,2),c(3,3),c(3,3),c(3,3),c(4,3),c(4,3),c(4,3)))
  #       for (k in 13:n) {
  #         if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'} 
  #         eval(test)
  #       }
  #     }
  #     if (n>24 & n<=36){
  #       par(mfrow=c(4,3))
  #       for (k in 1:12) {
  #         if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'} 
  #         eval(test)
  #       }
  #       x11()
  #       par(mfrow=c(4,3))
  #       for (k in 13:24) {
  #         if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'} 
  #         eval(test)
  #       }
  #       x11()
  #       par(mfrow=switch(n-24,c(1,1),c(2,1),c(2,2),c(2,2),c(3,2),c(3,2),c(3,3),c(3,3),c(3,3),c(4,3),c(4,3),c(4,3)))
  #       for (k in 25:n) {
  #         if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'} 
  #         eval(test)
  #       }
  #     }
  #     
  #     total=total+n
  #     u=readline()
  #     if (u==1){break}
  #     dev.off()
  #   }
  # } 
}# old version


###################################################################
# Multiple plot function
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
###################################################################

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


dist_het=function(n,h,i,cutoff){
  if (h==0){
    bool=NA
  }else{
    j=0:n
    t=choose(n,j)*h^j*(1-h)^(n-j)
    t=cumsum(t)
    bool=((i<=max((0:n)[t<.025/cutoff])) || (i>= min((0:n)[t>1-.025/cutoff])) )
  }
  return(bool)
} #****COPIED*****

calculate_fst=function(gene_seq,pop_vector){
  if ( sum(is.na(gene_seq))==length(gene_seq)) {return(NA)} 
  else {
    name_list=c() # list of subpop
    for (i in unique(pop_vector)){
      assign(paste("list",i,sep="_"),seq(length(pop_vector))[pop_vector==i])
      name_list=c(name_list,paste("list",i,sep="_"))
    }
    
    new_vect=c() # frequency of the minor allele in eahc subpop
    vect_size=c() # number of indivudual per subpop
    for (i in name_list){
      new_vect=c(new_vect,mean(gene_seq[get(i)]/2,na.rm=TRUE))
      vect_size=c(vect_size,sum(!is.na(gene_seq[get(i)])))
    }
    if(2*mean(gene_seq/2,na.rm=TRUE)*(1-mean(gene_seq/2,na.rm=TRUE))==0){
      Fst=NA
    } else{
      #(2*new_vect*(1-new_vect)) heterozigosity within each subpop
      Fst=1-sum((2*new_vect*(1-new_vect))*vect_size/sum(vect_size))/(2*mean(gene_seq/2,na.rm=TRUE)*(1-mean(gene_seq/2,na.rm=TRUE))) # FST if the initial pop is split in n subpop. 
    }
    return(Fst)
    }
}

calculate_random_fst=function(gene_seq, pop_vector,bf_cor){
  save_random_fst=c()
  for (i in 1:1000){
    pop_vector2=pop_vector[order(rnorm(length(pop_vector)))]
    save_random_fst=c(save_random_fst,calculate_fst(gene_seq,pop_vector2))
  }
  return(quantile(save_random_fst,p=1-.05/bf_cor,na.rm=T))
} # 10/08 ***alex****

if (F){
  # display_specific_gene_allele_frequency=function(allele_name,gene_info,data_pop2,targetSNP,fst_values,fis,He,local_fst){
  #   color <- data_pop2$Sp2
  #   (color <- ifelse(color == "hisp", "red",
  #                    ifelse(color == "hyS", "orange", 
  #                           ifelse(color == "hyN", "green", 
  #                                  ifelse(color == "ns" , "blue", NA)))))
  # 
  #   test=expression(plot(data_pop2$lat, data_pop2[,addinfo_col_data-1+total+k],xlim=c(41.8,43),ylim=c(0,1), col=color, pch=16,main=paste(pos$CHROM_POS[total+k]),col.main=colortitle, sub = paste(add_info),ylab = "allele frequency", xlab = paste("latitude\n Fst=" ,signif(fst_values[total+k],3),"he=", signif(He[total+k],3)," Fis in pop hisp",signif(fis[4,total+k],3),"hyS",signif(fis[3,total+k],3),"hyN",signif(fis[2,total+k],3),"ns",signif(fis[1,total+k],3),"\n Fst in pop hisp ",signif(local_fst[4,total+k],3),"hyS",signif(local_fst[3,total+k],3),"hyN",signif(local_fst[2,total+k],3),"ns",signif(local_fst[1,total+k],3))))
  #   
  #   gene_name=unlist(strsplit(allele_name,".",fixed=TRUE))[1]
  #   n1=seq(dim(gene_info)[1])[gene_info[,1]==gene_name]
  #   
  #   total=sum(as.numeric(gene_info[1:(n1-1),2]))
  #   x11()
  #   n=as.numeric(gene_info[n1,2])
  #   
  #   if (n<=12){
  #     par(mfrow=switch(n,c(1,1),c(2,1),c(2,2),c(2,2),c(3,2),c(3,2),c(3,3),c(3,3),c(3,3),c(4,3),c(4,3),c(4,3)))
  #     for (k in 1:n) {
  #       if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'} 
  #       if (){add_info=}else {addinfo=''}
  #       eval(test)
  #     }
  #   }
  #   if (n>12 & n<=24){
  #     par(mfrow=c(4,3))
  #     for (k in 1:12) {
  #       if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'} 
  #       eval(test)
  #     }
  #     x11()
  #     par(mfrow=switch(n-12,c(1,1),c(2,1),c(2,2),c(2,2),c(3,2),c(3,2),c(3,3),c(3,3),c(3,3),c(4,3),c(4,3),c(4,3)))
  #     for (k in 13:n) {
  #       if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'} 
  #       eval(test)
  #     }
  #   }
  #   if (n>24 & n<=36){
  #     par(mfrow=c(4,3))
  #     for (k in 1:12) {
  #       if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'} 
  #       eval(test)
  #     }
  #     x11()
  #     par(mfrow=c(4,3))
  #     for (k in 13:24) {
  #       if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'} 
  #       eval(test)
  #     }
  #     x11()
  #     par(mfrow=switch(n-24,c(1,1),c(2,1),c(2,2),c(2,2),c(3,2),c(3,2),c(3,3),c(3,3),c(3,3),c(4,3),c(4,3),c(4,3)))
  #     for (k in 25:n) {
  #       if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'} 
  #       eval(test)
  #     }
  #   }
  # }
} # old version specific gene

# ---- new update ---

add_significant_test=function(name,sign_fis_ns,sign_fis_hyN,sign_fis_hyS,sign_fis_hisp,sign_fst){
  add_info="SIGN:"
  if (sign_fst){add_info=paste(add_info," FST; ", sep='')}
  if (sum(name==sign_fis_ns[,1])){add_info=paste(add_info," fis in ns; ", sep='')}
  if (sum(name==sign_fis_hyN[,1])){add_info=paste(add_info," fis in hyN; ", sep='')}
  if (sum(name==sign_fis_hyS[,1])){add_info=paste(add_info," fis in hyS; ", sep='')}
  if (sum(name==sign_fis_hisp[,1])){add_info=paste(add_info," fis in hisp; ", sep='')}
  return(add_info)
}

display_gene_freq=function(n,total,gene_info,data_pop2,targetSNP,fst_values,fis,He,local_fst,sign_fis_ns,sign_fis_hyN,sign_fis_hyS,sign_fis_hisp,sign_fst){
  color <- data_pop2$Sp2
  (color <- ifelse(color == "hisp", "red",
                   ifelse(color == "hyS", "orange", 
                          ifelse(color == "hyN", "green", 
                                 ifelse(color == "ns" , "blue", NA)))))
  
  test=expression(plot(data_pop2$lat, data_pop2[,addinfo_col_data-1+total+k],xlim=c(41.8,43),ylim=c(0,1), col=color, pch=16,main=paste(pos$CHROM_POS[total+k]),col.main=colortitle, sub = paste(add_info),ylab = "allele frequency", xlab = paste("latitude\n Fst=" ,signif(fst_values[total+k],3),"he=", signif(He[total+k],3)," Fis in pop hisp",signif(fis[4,total+k],3),"hyS",signif(fis[3,total+k],3),"hyN",signif(fis[2,total+k],3),"ns",signif(fis[1,total+k],3),"\n Fst in pop hisp ",signif(local_fst[4,total+k],3),"hyS",signif(local_fst[3,total+k],3),"hyN",signif(local_fst[2,total+k],3),"ns",signif(local_fst[1,total+k],3))))
  x11()
  if (n<=12){
    par(mfrow=switch(n,c(1,1),c(2,1),c(2,2),c(2,2),c(3,2),c(3,2),c(3,3),c(3,3),c(3,3),c(4,3),c(4,3),c(4,3)))
    for (k in 1:n) {
      if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'}
      add_info=add_significant_test(colnames(data_pop2)[addinfo_col_data-1+total+k],sign_fis_ns,sign_fis_hyN,sign_fis_hyS,sign_fis_hyS,sign_fst[total+k])
      eval(test)
    }
  }
  if (n>12 & n<=24){
    par(mfrow=c(4,3))
    for (k in 1:12) {
      if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'}
      add_info=add_significant_test(colnames(data_pop2)[addinfo_col_data-1+total+k],sign_fis_ns,sign_fis_hyN,sign_fis_hyS,sign_fis_hyS,sign_fst[total+k])
      eval(test)
    }
    x11()
    par(mfrow=switch(n-12,c(1,1),c(2,1),c(2,2),c(2,2),c(3,2),c(3,2),c(3,3),c(3,3),c(3,3),c(4,3),c(4,3),c(4,3)))
    for (k in 13:n) {
      if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'}
      add_info=add_significant_test(colnames(data_pop2)[addinfo_col_data-1+total+k],sign_fis_ns,sign_fis_hyN,sign_fis_hyS,sign_fis_hyS,sign_fst[total+k])
      eval(test)
    }
  }
  if (n>24 & n<=36){
    par(mfrow=c(4,3))
    for (k in 1:12) {
      if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'}
      add_info=add_significant_test(colnames(data_pop2)[addinfo_col_data-1+total+k],sign_fis_ns,sign_fis_hyN,sign_fis_hyS,sign_fis_hyS,sign_fst[total+k])
      eval(test)
    }
    x11()
    par(mfrow=c(4,3))
    for (k in 13:24) {
      if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'}
      add_info=add_significant_test(colnames(data_pop2)[addinfo_col_data-1+total+k],sign_fis_ns,sign_fis_hyN,sign_fis_hyS,sign_fis_hyS,sign_fst[total+k])
      eval(test)
    }
    x11()
    par(mfrow=switch(n-24,c(1,1),c(2,1),c(2,2),c(2,2),c(3,2),c(3,2),c(3,3),c(3,3),c(3,3),c(4,3),c(4,3),c(4,3)))
    for (k in 25:n) {
      if (targetSNP[total+k]=="Target"){colortitle='red'}else{colortitle='black'}
      add_info=add_significant_test(colnames(data_pop2)[addinfo_col_data-1+total+k],sign_fis_ns,sign_fis_hyN,sign_fis_hyS,sign_fis_hyS,sign_fst[total+k])
      eval(test)
    }
  }
  
} # ******* 10/08 alex

# plot function

displayallelefrequency=function(gene_info,data_pop2,targetSNP,fst_values,fis,He,local_fst,sign_fis_ns,sign_fis_hyN,sign_fis_hyS,sign_fis_hisp,sign_fst){
  total=0
  for (i in 2:dim(gene_info)[1]){
    n=as.numeric(gene_info[i,2])
    display_gene_freq(n,total,gene_info,data_pop2,targetSNP,fst_values,fis,He,local_fst,sign_fis_ns,sign_fis_hyN,sign_fis_hyS,sign_fis_hisp,sign_fst)
    total=total+n
    u=readline()
    if (u==1){break}
    dev.off()
  }
} # ******** modified alex 10/08

display_specific_gene_frequency=function(gene_name,gene_info,data_pop2,targetSNP,fst_values,fis,He,local_fst,sign_fis_ns,sign_fis_hyN,sign_fis_hyS,sign_fis_hisp,sign_fst){
  n1=seq(dim(gene_info)[1])[gene_info[,1]==gene_name]
  
  total=sum(as.numeric(gene_info[1:(n1-1),2]))
  n=as.numeric(gene_info[n1,2])
  
  display_gene_freq(n,total,gene_info,data_pop2,targetSNP,fst_values,fis,He,local_fst,sign_fis_ns,sign_fis_hyN,sign_fis_hyS,sign_fis_hisp,sign_fst)
  u=readline()
  if (u==1){break}
  dev.off()
  
} # specific gene *******10/08 alex

#--------------------------------------------

drawing_into_pdist=function(cumul_dist,n1,n2){
  random=runif(1000,0,1)
  p_tab=array(0,dim=c(1000,3))
  for (i in 1:1000){
    p=min(seq(10001)[random[i]<=cumul_dist])/10000
    seq=c(rep(1,floor(p*(n1+n2))),rep(0,n1+n2-floor(p*(n1+n2))))[order(runif(n1+n2,0,1))]
    p_tab[i,1]=mean(seq[1:n1])
    p_tab[i,2]=mean(seq[(n1+1):(n1+n2)])
    p_tab[i,3]=mean(seq)
  }
  
  return(p_tab)
}

calculate_stat_2pop_adap_data<-function(data1,data2){
  jsfs2D <- array(0,dim=c(101,101))
  for(l in (addinfo_col_data+1):dim(data1)[2]) {
    if (mean(c(data1[,l],data2[,l]),na.rm=T)/2<.5){
    j <- round(100*mean(data1[,l],na.rm=T)/2)
    h <- round(100*mean(data2[,l],na.rm=T)/2)
    }
    else{
      j <- 100-round(100*mean(data1[,l],na.rm=T)/2)
      h <- 100-round(100*mean(data2[,l],na.rm=T)/2)
    }
    jsfs2D[j+1,h+1] <- jsfs2D[j+1,h+1]+1
  }
  
  
  i <- 0
  statistic <- numeric()
  for(j in list(1,2:5,6:94,95:99,100)) {
    for(h in list(1,2:5,6:94,95:99,100)) {
      i <- i+1
      statistic[i] <- sum(jsfs2D[j,h])
    }
  }
  statistic <- statistic[2:24]
  
  return(statistic)
} # *** 18 /08 alex ***
