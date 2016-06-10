library(ape)
merge_iden_uptona=function(cprime){
  if (dim(cprime)[1]>1){
    for (i in 1:(dim(cprime)[1]-1)){
      hap1=unlist(strsplit(rownames(cprime)[i],split=''))
      for (j in (i+1):dim(cprime)[1]){
	hap2=unlist(strsplit(rownames(cprime)[j],split=''))
	filter=hap1!=hap2
	if (sum(filter)>0 && (hap1[filter]=="-" || hap2[filter]=="-") && sum(cprime[i,])>0){
	  cprime[i,]=cprime[i,]+cprime[j,]
	  cprime[j,]=0*cprime[j,]
	}
      }
    }
  }
  return(cprime)
}

display_relevant=function(data){
  step=data[rowSums(data)>0,]
  step=step[order(rowSums(step),decreasing=T),]
  return(step)
}

display_bar=function(data,do_plot=F){
  info=t(t(data)/rowSums(t(data)))
  if (do_plot){
    x11()
    list_color5=rep('gray',dim(info)[1])
    list_color5[1:8]=c('purple','blue','cyan','green','yellow','orange','red','brown')
    barplot(info,col=list_color5)}
  nu_hap=c()
  for (i in colnames(data)){
    loc_read=data[data[,colnames(data)==i]>0,colnames(data)==i]
    loc_hap=rownames(data)[data[,colnames(data)==i]>0]
    nu_hap=c(nu_hap,length(loc_hap))
  }
  return(nu_hap)
}

export_to_fasta=function(data,name){
  output=file(name,open="w")
  for (i in 1:dim(data)[2]){
    loc_read=data[data[,i]>0,i]
    loc_hap=rownames(data)[data[,i]>0][order(loc_read,decreasing=T)]
    loc_read=loc_read[order(loc_read,decreasing=T)]
    for (j in 1:length(loc_hap)){
      line=paste(c('>',colnames(data)[i],'_hap_',j,'_freq_',signif(loc_read[j],2)),sep='',collapse='')
      writeLines(line,con=output)
      writeLines(loc_hap[j],con=output)
    }
  }
  close(output)
}

export_haplotype=function(data,name){
  output=file(name,open="w")
  for (i in 1:length(rownames(data))){
    line=paste(c('> hap_',i,'_freq_',signif(rowSums(data)[i]/sum(data),2),'\n',rownames(data)[i]),sep='',collapse='')
    writeLines(line,con=output)
  }
  close(output)
}

highlight_snp=function(file_name,list_hap){
  if (length(list_snp_fis[corresponding_amp==file_name])>0) {
    for (i in 1:length(list_hap)){
      for (k in list_snp_fis[corresponding_amp==file_name]){
	test=unlist(strsplit(list_hap[i],split=''))
	test[k]=paste(c('*',test[k],'*'),collapse='')
	test=paste(test,collapse='')
      }
    list_hap[i]=test  
    }
  }
  return(list_hap)
}

make_tree=function(data,title='',tree_for_ind=F,ind=ind,activate_bar=F,file_name=file_data){
  lim=0
  hap=rownames(data)
  ref_hap=rownames(data)[data[,1]>0]
  mat_hap_dist=adist(hap)
  colnames(mat_hap_dist)=mat_hap_dist[hap==ref_hap,]
  rownames(mat_hap_dist)=hap
  
  hc=hclust(as.dist(mat_hap_dist),method="ward.D2")

  hc2=as.phylo(hc)

  empty_count=rep(0,max(hc2$edge))
  for (i in 1:sum(sort(colSums(data))>lim)){
    empty_count[hc2$edge[hc2$edge[,2]==i,2]]=sort(colSums(data),decreasing=T)[i]
    empty_count[hc2$edge[hc2$edge[,2]==i,1]]=empty_count[hc2$edge[hc2$edge[,2]==i,1]]+empty_count[hc2$edge[hc2$edge[,2]==i,2]]
  }
  for (i in max(hc2$edge):(sum(sort(colSums(data))>lim)+1)){
    empty_count[hc2$edge[hc2$edge[,2]==i,1]]=empty_count[hc2$edge[hc2$edge[,2]==i,1]]+empty_count[hc2$edge[hc2$edge[,2]==i,2]]
  }


  value_node=c()
  for (i in 1:dim(hc2$edge)[1]){
    value_node=c(value_node,empty_count[hc2$edge[i,2]])}
  
  hc2$tip.label[hc2$tip.label==ref_hap]=paste(c(hc2$tip.label[hc2$tip.label==ref_hap],'*ref'),collapse='')
  
  hc2$tip.label=highlight_snp=highlight_snp(file_name,hc2$tip.label)
  
  list_color6=rep('black',dim(data)[1])
  list_color6[1:8]=c('purple','blue','cyan','green','orange','red','magenta','brown')
  par(mfrow=c(1,2))
  plot(hc2,edge.color=gray(1-value_node/max(value_node)),tip.color=list_color6,type = "cladogram",main=title,cex=.8,no.margin=T)
  hc2=as.phylo(hc)
  hc2$tip.label[hc2$tip.label==ref_hap]=paste(c(hc2$tip.label[hc2$tip.label==ref_hap],'*ref'),collapse='')
  plot(hc2,edge.color='black',tip.color=list_color6,type = "unrooted",lab4ut="axial")
  
  x11()
  
  heatmap(mat_hap_dist,Rowv=as.dendrogram(hc),Colv=as.dendrogram(hc))
  
  if (activate_bar){
    print('number of group')
    nb_grp=readline()
    print('number of element in each group: c(n1,n2,...)')
    com_grp=readline()
    
    plot_group(data,as.phylo(hc)$tip.label,as.integer(nb_grp),eval(parse(text=com_grp)))
  }
  
  if(tree_for_ind){
    x11()
    hc3=as.phylo(hc)
    hap_tree=hc3$tip.label
    pop='1'
    list_pop=rapply(strsplit(as.character(ind),split='-'), function(x) x[1], how = "unlist")
    while (pop!=''){
      print('choose a pop')
      pop=readline()
      if (sum(list_pop==pop)){ 
	hc3$tip.label=rep('',length(hap_tree))
	hc3$tip.label[hap_tree==ref_hap]=paste(c(hc3$tip.label[hap_tree==ref_hap],'*ref*'),collapse='')
	sub_data=display_relevant(data[,list_pop==pop])

	
	for (j in 1:dim(sub_data)[2]){
	  for (k in 1:dim(sub_data)[1]){
	    if (sub_data[k,j]>0){
	      hc3$tip.label[hap_tree==rownames(sub_data)[k]]=paste(c(hc3$tip.label[hap_tree==rownames(sub_data)[k]],colnames(sub_data)[j]),collapse=' ')
	    }
	  }
	}
	plot(hc3,edge.color='black',type = "unrooted",lab4ut="axial",tip.color='red')
      }
      else{ print("doesn't exist") }
    }
  }
}

plot_group=function(data,order_hap,nb_grp,com_grp){
  new_mat=matrix(0,nb_grp,dim(data)[2])
  k=1
  for (i in 1:length(order_hap)){
    new_mat[k,]=new_mat[k,]+data[rownames(data)==order_hap[i],]
    if (i==sum(com_grp[1:k])){
      k=k+1
    }
  }
  vect_color=rainbow(nb_grp)
  for (l in nb_grp:1){
    new_mat[l,]=as.integer(new_mat[l,]>0)
  }
  for (l in nb_grp:1){
    barplot(new_mat[nb_grp-l+1,]*0+l,ylim=c(0,nb_grp),col='white',border='white')
    par(new=T)
    barplot(new_mat[nb_grp-l+1,]+(l-1),ylim=c(0,nb_grp),col=vect_color[l])
    par(new=T)
  }
  barplot(new_mat[1,]*0,ylim=c(0,nb_grp),names.arg=colnames(data),las=2,cex.names=0.4)
}

rare_hap=function(data,loc_hap,loc_read_cor,i,filter=10,method='delete'){
  output=data
  if (method=='delete'){
    for (j in 1:length(loc_hap)){
      if (loc_read_cor[j]<filter){
	output[rownames(data)==loc_hap[j],i]=0
      } else {
	output[rownames(data)==loc_hap[j],i]=loc_read_cor[j]
      }
    }
  } else if(method=='merge_or_delete'){
      for (j in 1:length(loc_hap)){
	if (loc_read_cor[j]<filter){
	  if(sum(adist(loc_hap)[1:j,j]==1)>0){
	    k=seq(j)[adist(loc_hap)[1:j,j]==1][1]
	    if (loc_read_cor[k]>filter){
	      output[rownames(data)==loc_hap[k],i]=output[rownames(data)==loc_hap[k],i]+loc_read_cor[j]
	    }
	  }
	  output[rownames(data)==loc_hap[j],i]=0
	}
	else {
	  output[rownames(data)==loc_hap[j],i]=loc_read_cor[j]
	}
      }
  } else if (method=='merge_or_keep'){
      for (j in 1:length(loc_hap)){
	if (loc_read_cor[j]<filter){
	  if(sum(adist(loc_hap)[1:j,j]==1)>0){
	    k=seq(j)[adist(loc_hap)[1:j,j]==1][1]
	    if (loc_read_cor[k]>filter){
	      output[rownames(data)==loc_hap[k],i]=output[rownames(data)==loc_hap[k],i]+loc_read_cor[j]
	      output[rownames(data)==loc_hap[j],i]=0
	    }
	    else{
	      output[rownames(data)==loc_hap[j],i]=loc_read_cor[j]
	    }
	  }
	}
	else {
	  output[rownames(data)==loc_hap[j],i]=loc_read_cor[j]
	}
    }
  } else {
    stop("method unvalid",call.=F)} 
  return(output) 
}

filter_data=function(b,filter,method){
  output=b
  for (i in 1:dim(b)[2]){
    if (colnames(b)[i]!="REF"){
      loc_read=b[b[,i]>0,i]
      loc_hap=rownames(b)[b[,i]>0][order(loc_read,decreasing=T)]
      loc_read=loc_read[order(loc_read,decreasing=T)]
      loc_read_cor=merge_iden_uptona(as.matrix(loc_read,nrow=1))
      loc_hap=loc_hap[order(loc_read_cor,decreasing=T)]
      loc_read_cor=loc_read_cor[order(loc_read_cor,decreasing=T)]
      output=rare_hap(output,loc_hap,loc_read_cor,i,filter,method)
    }
  }
  output=display_relevant(output)
  return(output)
}

fast_analysis=function(file_ind,file_data,freq=F,method='delete',lim=10,tree_for_ind=F,activate_bar=F){
  ind=read.table(file_ind,sep='\t',header=F,colClasses="character")
  a=read.table(file_data,header=T,check.names=FALSE)
  a=a[rowSums(a)>0,]
  rownames(a)=ind
  b=t(a[,order(colSums(a),decreasing=T)])

  if (!freq){
    b_del=filter_data(b,lim,"delete")
    b_mer_or_del=filter_data(b,lim,"merge_or_delete")
    b_mer_or_keep=filter_data(b,lim,"merge_or_keep")
    b_stat=display_bar(b)
  } else {
    if(lim>=1){
      lim=.1
      print("filter=10%")}
    bprime=t(t(b)/rowSums(t(b)))
    b_del=filter_data(bprime,lim,"delete")
    b_mer_or_del=filter_data(bprime,lim,"merge_or_delete")
    b_mer_or_keep=filter_data(bprime,lim,"merge_or_keep")
    b_stat=display_bar(bprime)
  }
  
  b_del_stat=display_bar(b_del)
  b_mer_or_del_stat=display_bar(b_mer_or_del)
  b_mer_or_keep_stat=display_bar(b_mer_or_keep)

  x11()
  par(mfrow=c(2,2))
  hist(b_stat,col='black',main=file_data)
  hist(b_mer_or_keep_stat,col='blue')
  hist(b_mer_or_del_stat,col='magenta')
  hist(b_del_stat,col='red')
  
  
  x11()
  if (method=='delete'){
    make_tree(b_del,file_data,tree_for_ind,ind,activate_bar,file_data)
    return(b_del)
  } else if (method=='merge_or_delete'){
    make_tree(b_mer_or_del,file_data,tree_for_ind,ind,activate_bar,file_data)
    return(b_mer_or_del)
  } else if (method=='merge_or_keep'){
    make_tree(b_mer_or_keep,file_data,tree_for_ind,ind,activate_bar,file_data)
    return(b_mer_or_keep)
  }
  
}
#  export_to_fasta(b_mer_or_del,"GSMUA_Achr7G07160_SWA2_45268_uniq_hap_mer_or_del.fasta")



#  b_del_5=display_relevant(b_del[,b_del_stat==5])

  # with freq



