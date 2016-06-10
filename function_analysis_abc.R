library(abc)

filter_table=function(table){
  v=c()
  for (i in 1:dim(table)[2]){
    if (max(table[,i])>0){
      v=c(v,i)}
      }
      return(v)
      }

filter_table_v2=function(table,n1,n2,func=mean){
  v=c(1:n1)
  for (i in (n1+1):dim(table)[2]){
    if (eval(func)(table[,i])>n2){
      v=c(v,i)}
    }
    return(v)
  }
  
  
display_true_data=function(data,simul,n,param_name){
  max_table=c()
  mean_table=c()
  variance_table=c()
  for (i in 1:(length(data)+n)){
    mean_table=c(mean_table,mean(simul[,i]))
    variance_table=c(variance_table,var(simul[,i]))
    max_table=c(max_table,max(simul[,i]))
  }

  vect=seq(length(data))[max_table[1:length(data)+n]>0 | data>0]
  
  x11(width=20,height=12)
  if (length(vect)<100){
    par(mfrow=c(4,(length(vect)+n+3)%/%4))  
    } else {
    par(mfrow=c(6,(length(vect)+n+5)%/%6))  
    }
  for (i in c(1:n,vect+n)){
    if (i<=n){
      boxplot(simul[,i],xlab=param_name[i])
    } else {
      boxplot(simul[,i],ylim=c(min(simul[,i],data[i-n]),max(simul[,i],data[i-n])))
      par(new=T)
      plot(data[i-n],col="red",ylim=c(min(simul[,i],data[i-n]),max(simul[,i],data[i-n])),pch=16,ylab='',xlab=paste('S',i-n),xaxt='n')
    }
  }
  return(vect)
}
  
analysis_any_m=function(m){
  filter_m=filter_table_v2(table_3pop,6,m)
  table_3pop_m1=table_3pop[,filter_m]

  test2=cv4abc(param=table_3pop_m1[,1:6],sumstat=table_3pop_m1[,7:dim(table_3pop_m1)[2]],tols=.05, method="loclinear",nval=20)
  print(summary(test2))
  par(mfrow=c(2,3))
  plot(test2)

  u=readline('')

  true_data_3pop_m1=true_data_3pop[(filter_m[7:length(filter_m)]-6)]
  a=abc(true_data_3pop_m1,param=table_3pop_m1[,1:6],table_3pop_m1[,7:dim(table_3pop_m1)[2]],tol=.05,method="loclinear")
  summary(a)
  plot(a,param=table_3pop[,1:6])

  par(mfrow=c(2,3))
  for (i in 1:6){
    plot(density(table_3pop_m1[,i]),xlim=c(min(table_3pop_m1[,i]),max(table_3pop_m1[,i])),ylim=c(0,max(density(a$adj.values[,i])$y,density(table_3pop_m1[,i])$y)),xlab=c('theta','m12','m23','mt2','t1','t2')[i],main='')
    par(new=T)
    plot(density(a$adj.values[,i]),xlim=c(min(table_3pop_m1[,i]),max(table_3pop_m1[,i])),ylim=c(0,max(density(a$adj.values[,i])$y,density(table_3pop_m1[,i])$y)),col='red',xaxt='n',xlab='',main='')
  }
}

analysis_any_m_4pop=function(m){
  filter_m=filter_table_v2(table_4pop,7,m)
  table_3pop_m1=table_4pop[,filter_m]

  test2=cv4abc(param=table_3pop_m1[,1:7],sumstat=table_3pop_m1[,8:dim(table_3pop_m1)[2]],tols=.05, method="loclinear",nval=20)
  print(summary(test2))
  par(mfrow=c(3,3))
  plot(test2)

  u=readline('')

  true_data_3pop_m1=true_data_4pop[(filter_m[8:length(filter_m)]-7)]
  a=abc(true_data_3pop_m1,param=table_3pop_m1[,1:7],table_3pop_m1[,8:dim(table_3pop_m1)[2]],tol=.05,method="loclinear")
  summary(a)
  plot(a,param=table_4pop[,1:7])
  
    par(mfrow=c(3,3))
  for (i in 1:7){
    plot(density(table_3pop_m1[,i]),xlim=c(min(table_3pop_m1[,i]),max(table_3pop_m1[,i])),ylim=c(0,max(density(a$adj.values[,i])$y,density(table_3pop_m1[,i])$y)),xlab=c('theta','m12','m23','m34','t1','t2','t3')[i],main='')
    par(new=T)
    plot(density(a$adj.values[,i]),xlim=c(min(table_3pop_m1[,i]),max(table_3pop_m1[,i])),ylim=c(0,max(density(a$adj.values[,i])$y,density(table_3pop_m1[,i])$y)),col='red',xaxt='n',xlab='',main='')
  }
}

analysis_any_m_4pop_v2=function(m){
  filter_m=filter_table_v2(table_4pop_v2,9,m)
  table_3pop_m1=table_4pop_v2[,filter_m]

  test2=cv4abc(param=table_3pop_m1[,1:9],sumstat=table_3pop_m1[,10:dim(table_3pop_m1)[2]],tols=.05, method="loclinear",nval=20)
  print(summary(test2))
  par(mfrow=c(3,3))
  plot(test2)

  u=readline('')

  true_data_3pop_m1=true_data_4pop[(filter_m[10:length(filter_m)]-9)]
  a=abc(true_data_3pop_m1,param=table_3pop_m1[,1:9],table_3pop_m1[,10:dim(table_3pop_m1)[2]],tol=.05,method="loclinear")
  summary(a)
  plot(a,param=table_4pop_v2[,1:9])
  
  par(mfrow=c(3,3))
  for (i in 1:9){
    plot(density(table_3pop_m1[,i]),xlim=c(min(table_3pop_m1[,i]),max(table_3pop_m1[,i])),ylim=c(0,max(density(a$adj.values[,i])$y,density(table_3pop_m1[,i])$y)),xlab=c('theta','m12','m23','m34','mt1','mt2','t1','t2','t3')[i],main='')
    par(new=T)
    plot(density(a$adj.values[,i]),xlim=c(min(table_3pop_m1[,i]),max(table_3pop_m1[,i])),ylim=c(0,max(density(a$adj.values[,i])$y,density(table_3pop_m1[,i])$y)),col='red',xaxt='n',xlab='',main='')
  }
}


analysis_4pop=function(data,true_data,tolerance=.05,sample_cv=50,n_param=7,methode="loclinear",estimate=T){
  test2=cv4abc(param=data[,1:n_param],sumstat=data[,(n_param+1):dim(data)[2]],tols=tolerance, method=methode,nval=sample_cv)
  print(summary(test2))
  par(mfrow=c(3,3))
  plot(test2)
  u=readline('')
  a=abc(true_data,param=data[,1:n_param],data[,(n_param+1):dim(data)[2]],tol=tolerance,method=methode)
  summary(a)
  if (methode!='rejection'){plot(a,param=data[,1:n_param])}
  par(mfrow=c(3,3))
  for (i in 1:n_param){
    plot(density(data[,i]),xlim=c(min(data[,i]),max(data[,i])),ylim=c(0,max(density(a$adj.values[,i])$y,density(a$unadj.values[,i])$y,density(data[,i])$y)),xlab=c('theta','m12','m23','m34','t1','t2','t3')[i],main='')
    par(new=T)
    plot(density(a$adj.values[,i]),xlim=c(min(data[,i]),max(data[,i])),ylim=c(0,max(density(a$adj.values[,i])$y,density(a$unadj.values[,i])$y,density(data[,i])$y)),col='red',xaxt='n',xlab='',main='')
    par(new=T)
    plot(density(a$unadj.values[,i]),xlim=c(min(data[,i]),max(data[,i])),ylim=c(0,max(density(a$adj.values[,i])$y,density(a$unadj.values[,i])$y,density(data[,i])$y)),col='green',xaxt='n',xlab='',main='')
  }
#   hist(a$dist)
}

random_abc=function(simul,choice=0,tolerance=.05,n_param=7,methode="loclinear"){
  if (choice==0){
    h=ceiling(runif(1, 0, dim(simul)[1]))
  } else {
    h=choice
  }
  data=simul[-h,]
  a=abc(simul[h,(n_param+1):dim(simul)[2]],param=data[,1:n_param],data[,(n_param+1):dim(data)[2]],tol=tolerance,method=methode)
  summary(a)
#  plot(a,param=data[,1:n_param])
  par(mfrow=c(3,3))
  for (i in 1:n_param){
    plot(density(data[,i]),xlim=c(min(data[,i]),max(data[,i])),ylim=c(0,max(density(a$adj.values[,i])$y,density(a$unadj.values[,i])$y,density(data[,i])$y)),xlab=c('theta','m12','m23','m34','t1','t2','t3')[i],main='')
    par(new=T)
    plot(density(a$adj.values[,i]),xlim=c(min(data[,i]),max(data[,i])),ylim=c(0,max(density(a$adj.values[,i])$y,density(a$unadj.values[,i])$y,density(data[,i])$y)),col='red',xaxt='n',xlab='',main='')
    par(new=T)
    plot(density(a$unadj.values[,i]),xlim=c(min(data[,i]),max(data[,i])),ylim=c(0,max(density(a$adj.values[,i])$y,density(a$unadj.values[,i])$y,density(data[,i])$y)),col='green',xaxt='n',xlab='',main='')
    abline(v=simul[h,i],col='blue')
  }
}



real_data_4D=c(192,5,0,116,61,8,0,0,1,3,0,0,0,0,0,82,29,1,0,14,36,26,0,0,9,24,0,0,0,0,0,0,3,0,0,0,4,2,0,0,1,21,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,233,23,1,0,35,19,0,0,0,2,3,0,0,0,0,0,58,15,4,0,39,51,22,0,1,2,37,0,0,0,0,0,2,2,0,0,1,9,2,0,1,6,82,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,14,4,0,0,23,16,0,0,4,2,4,0,0,0,0,0,5,1,0,0,29,31,11,0,26,63,236,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
)

fst_real_data_4D=c(0.04327951,0.01908287,0.04098857,0.05863092,0.01863359,0.03411949,0.01416637,0.04327951,0.01491895,0.03375938,0.05160449,0.01383036,0.02669383,0.01100316)

diagonal_jsfs_4D=c(1006,1045,930,896,1029,917,1122)

  potential_filter=c(T,T,F,F,T,T,T,F,F,T,T,T,F,F,F,F,F,F,F,F,F,F,F)
  true_data_4pop=c(275,61,0,0,96,140,88,0,0,3,62,439,0,0,0,0,0,0,0,0,0,0,0,313,107,0,0,124,101,99,0,0,31,95,378,0,0,0,0,0,0,0,0,0,0,0,378,145,0,0,156,59,105,2,0,109,115,271,7,1,0,0,0,0,0,0,0,0,0,237,11,0,0,196,194,87,0,0,24,78,486,0,0,0,0,0,0,0,0,0,0,0,292,27,0,0,215,148,114,0,0,86,121,381,0,0,0,0,0,0,0,0,0,0,0,266,15,0,0,250,191,68,0,0,41,104,439,0,0,0,0,0,0,0,0,0,0,0)
  
convert_to_pairwise_jsfs=function(stat2,include_0=T){
  stat2=as.double(c(0,stat2))
  jsfs4D2 <- array(0,dim=c(4,4,4,4))
  for (j in 1:4) {
    for (h in 1:4) {
      jsfs4D2[j,h,,]=t(array(stat2[1:16],dim=c(4,4)))
      stat2=stat2[17:length(stat2)]
    }
  }
  stat45=c()
  for (i in 1:4){
    for (j in 1:4){
      stat45=c(stat45,sum(jsfs4D2[i,j,,]))
      }}
  for (i in 1:4){
    for (j in 1:4){
      stat45=c(stat45,sum(jsfs4D2[i,,j,]))
      }}
  for (i in 1:4){
    for (j in 1:4){
      stat45=c(stat45,sum(jsfs4D2[i,,,j]))
      }}
  for (i in 1:4){
    for (j in 1:4){
      stat45=c(stat45,sum(jsfs4D2[,i,j,]))
      }}
  for (i in 1:4){
    for (j in 1:4){
      stat45=c(stat45,sum(jsfs4D2[,i,,j]))
      }}
  for (i in 1:4){
    for (j in 1:4){
      stat45=c(stat45,sum(jsfs4D2[,,i,j]))
      }}
      
  return(stat45[c(include_0,T,T,F,T,T,T,F,T,T,T,F,F,F,F,F)])
}

convert_to_joachim_stat=function(stat2){
  stat2=as.double(c(0,stat2))
  jsfs4D2 <- array(0,dim=c(4,4,4,4))
  for (j in 1:4) {
    for (h in 1:4) {
      jsfs4D2[j,h,,]=t(array(stat2[1:16],dim=c(4,4)))
      stat2=stat2[17:length(stat2)]
    }
  }
  stat45=c()
  stat45=c(stat45,sum(jsfs4D2[2:4,1,1,1]))
  stat45=c(stat45,sum(jsfs4D2[1,2:4,1,1]))
  stat45=c(stat45,sum(jsfs4D2[1,1,2:4,1]))
  stat45=c(stat45,sum(jsfs4D2[1,1,1,2:4]))
  stat45=c(stat45,sum(jsfs4D2[2:4,2:4,2:4,2:4]))
  stat45=c(stat45,sum(jsfs4D2[2:4,2:4,,]))
  stat45=c(stat45,sum(jsfs4D2[2:4,,2:4,]))
  stat45=c(stat45,sum(jsfs4D2[2:4,,,2:4]))
  stat45=c(stat45,sum(jsfs4D2[,2:4,2:4,]))
  stat45=c(stat45,sum(jsfs4D2[,2:4,,2:4]))
  stat45=c(stat45,sum(jsfs4D2[,,2:4,2:4]))
  stat45=c(stat45,sum(jsfs4D2[2,2,,],jsfs4D2[3,3,,],jsfs4D2[4,4,,]))
  stat45=c(stat45,sum(jsfs4D2[2,,2,],jsfs4D2[3,,3,],jsfs4D2[4,,4,]))
  stat45=c(stat45,sum(jsfs4D2[2,,,2],jsfs4D2[3,,,3],jsfs4D2[4,,,4]))
  stat45=c(stat45,sum(jsfs4D2[,2,2,],jsfs4D2[,3,3,],jsfs4D2[,4,4,]))
  stat45=c(stat45,sum(jsfs4D2[,2,,2],jsfs4D2[,3,,3],jsfs4D2[,4,,4]))
  stat45=c(stat45,sum(jsfs4D2[,,2,2],jsfs4D2[,,3,3],jsfs4D2[,,4,4]))
  return(stat45)
}

convert_to_statv3=function(stat2){
  stat2=as.double(c(0,stat2))
  jsfs4D2 <- array(0,dim=c(4,4,4,4))
  for (j in 1:4) {
    for (h in 1:4) {
      jsfs4D2[j,h,,]=t(array(stat2[1:16],dim=c(4,4)))
      stat2=stat2[17:length(stat2)]
    }
  }
  stat45=c()
  stat45=c(stat45,sum(jsfs4D2[,1,1,1]))
  stat45=c(stat45,sum(jsfs4D2[1,,1,1]))
  stat45=c(stat45,sum(jsfs4D2[1,1,,1]))
  stat45=c(stat45,sum(jsfs4D2[1,1,1,]))
  stat45=c(stat45,sum(jsfs4D2[2:4,2:4,2:4,2:4]))
  stat45=c(stat45,sum(jsfs4D2[2:4,2:4,,]))
  stat45=c(stat45,sum(jsfs4D2[2:4,,2:4,]))
  stat45=c(stat45,sum(jsfs4D2[2:4,,,2:4]))
  stat45=c(stat45,sum(jsfs4D2[,2:4,2:4,]))
  stat45=c(stat45,sum(jsfs4D2[,2:4,,2:4]))
  stat45=c(stat45,sum(jsfs4D2[,,2:4,2:4])) 
  stat45=c(stat45,sum(jsfs4D2[2,2,,],jsfs4D2[3,3,,],jsfs4D2[4,4,,]))
  stat45=c(stat45,sum(jsfs4D2[2,,2,],jsfs4D2[3,,3,],jsfs4D2[4,,4,]))
  stat45=c(stat45,sum(jsfs4D2[2,,,2],jsfs4D2[3,,,3],jsfs4D2[4,,,4]))
  stat45=c(stat45,sum(jsfs4D2[,2,2,],jsfs4D2[,3,3,],jsfs4D2[,4,4,]))
  stat45=c(stat45,sum(jsfs4D2[,2,,2],jsfs4D2[,3,,3],jsfs4D2[,4,,4]))
  stat45=c(stat45,sum(jsfs4D2[,,2,2],jsfs4D2[,,3,3],jsfs4D2[,,4,4]))
  stat45=c(stat45,sum(jsfs4D2[2:4,2:4,1,1]))
  stat45=c(stat45,sum(jsfs4D2[2:4,1,2:4,1]))
  stat45=c(stat45,sum(jsfs4D2[2:4,1,1,2:4]))
  stat45=c(stat45,sum(jsfs4D2[1,2:4,2:4,1]))
  stat45=c(stat45,sum(jsfs4D2[1,2:4,1,2:4]))
  stat45=c(stat45,sum(jsfs4D2[1,1,2:4,2:4]))
  return(stat45)
}

convert_to_statv4=function(stat2){
  stat2=as.double(c(0,stat2))
  jsfs4D2 <- array(0,dim=c(4,4,4,4))
  for (j in 1:4) {
    for (h in 1:4) {
      jsfs4D2[j,h,,]=t(array(stat2[1:16],dim=c(4,4)))
      stat2=stat2[17:length(stat2)]
    }
  }
  stat45=c()
  stat45=c(stat45,sum(jsfs4D2[,1,1,1]))
  stat45=c(stat45,sum(jsfs4D2[1,,1,1]))
  stat45=c(stat45,sum(jsfs4D2[1,1,,1]))
  stat45=c(stat45,sum(jsfs4D2[1,1,1,]))
  stat45=c(stat45,sum(jsfs4D2[2:4,2:4,2:4,2:4]))
  stat45=c(stat45,sum(jsfs4D2[2:4,2:4,,]))
  stat45=c(stat45,sum(jsfs4D2[2:4,,2:4,]))
  stat45=c(stat45,sum(jsfs4D2[2:4,,,2:4]))
  stat45=c(stat45,sum(jsfs4D2[,2:4,2:4,]))
  stat45=c(stat45,sum(jsfs4D2[,2:4,,2:4]))
  stat45=c(stat45,sum(jsfs4D2[,,2:4,2:4])) 
  stat45=c(stat45,sum(jsfs4D2[2:4,2:4,1,1]))
  stat45=c(stat45,sum(jsfs4D2[2:4,1,2:4,1]))
  stat45=c(stat45,sum(jsfs4D2[2:4,1,1,2:4]))
  stat45=c(stat45,sum(jsfs4D2[1,2:4,2:4,1]))
  stat45=c(stat45,sum(jsfs4D2[1,2:4,1,2:4]))
  stat45=c(stat45,sum(jsfs4D2[1,1,2:4,2:4]))
  return(stat45)
}


plot_identity=function(n1,n2,lev=.9,color='red'){
  plot(-1,xlim=c(0,1),ylim=c(0,1))
  for (i in seq(199)){
    u=rbinom(1000,n1,i/200)/n1
    v=rbinom(1000,n2,i/200)/n2
    dataEllipse(u,v,levels=lev,xlim=c(0,1),ylim=c(0,1),center.pch=F,plot.points=F,col=color)
  }
}

plot_identity_v2=function(n1,n2,lev=.9,color='red',plot_ellipe=T,plot_all=T){
  if (plot_all){plot(-1,xlim=c(0,1),ylim=c(0,1))}
  store_val=c()
  for (i in seq(199)){
    u=rbinom(10000,n1,i/200)/n1
    v=rbinom(10000,n2,i/200)/n2
    if (plot_ellipe &plot_all){dataEllipse(u,v,levels=lev,xlim=c(0,1),ylim=c(0,1),center.pch=F,plot.points=F,col=color)}
    val=dataEllipse(u,v,levels=lev,draw=F)
    store_val=rbind(store_val,val)
  }
  store_val=t(round(t(store_val)*c(n1,n2)))
  result=array(0,dim=c(n1+1,2))
  for (i in seq(0,n1)){
    temp=store_val[store_val[,1]==i,2]
    result[i+1,]=c(max(min(temp),0),min(max(temp),n2))
  }
  if (plot_all){
    par(new=T)
    plot(seq(0,n1)/n1,result[,1]/n2,xlim=c(0,1),ylim=c(0,1),type='l',lwd=2, col=if(plot_ellipe){'black'}else{color})
    par(new=T)
    plot(seq(0,n1)/n1,result[,2]/n2,xlim=c(0,1),ylim=c(0,1),type='l',lwd=2, col=if(plot_ellipe){'black'}else{color})
    return(result)
  }
}
  
  
calculate_prediction_error=function(table,do_plot,cv){
  cv.reg5=cv4abc(param=table[,1:3],sumstat=table[,4:dim(table)[2]],tols=c(.005,.01,.02,.05), method="loclinear",nval=cv)

  error_1=summary(cv.reg5)
  error_1=signif(error_1,3)
  if(do_plot){
    legend=c()
    for (i in seq(3)){
      legend=c(legend,paste(param_name[i],': tol=.005 ',error_1[[4*i-3]],'; tol=.01 ',error_1[[4*i-2]],'; tol=.02 ',error_1[[4*i-1]],'; tol=.05 ',error_1[[4*i]]))
    }
    legend[1]=paste(c('num_simul=',dim(table)[1],' cv n=',cv,' loclinear\n',legend[1]),sep='')
    x11()
    par(mfrow=c(1,3),mar = c(4, 4,2, 0.5))
    plot(cv.reg5,caption=legend)
  }
  return(error_1)
}
  

display_true_data=function(data,simul,n,param_name){
  max_table=c()
  mean_table=c()
  variance_table=c()
  for (i in 1:(23+n)){
    mean_table=c(mean_table,mean(simul[,i]))
    variance_table=c(variance_table,var(simul[,i]))
    max_table=c(max_table,max(simul[,i]))
  }

  vect=seq(23)[max_table[1:23+n]>0 | data>0]
  
  x11(width=20,height=12)  
  par(mfrow=c(2,(length(vect)+n+1)%/%2))  
  for (i in c(1:n,vect+n)){
    if (i<=n){
      boxplot(simul[,i],xlab=param_name[i])
    } else {
      boxplot(simul[,i],ylim=c(min(simul[,i],data[i-n]),max(simul[,i],data[i-n])))
      par(new=T)
      plot(data[i-n],col="red",ylim=c(min(simul[,i],data[i-n]),max(simul[,i],data[i-n])),pch=16,ylab='',xlab=paste('S',i-n),xaxt='n')
    }
  }
  
  return(vect)
}

plot_result=function(data,simul,tol,test,true){
  result=abc(data,param=simul[,1:3],sumstat=simul[,4:dim(simul)[2]],tol=tol, method="loclinear")
  summary(result)
  if (test) {
    plot(result,param=simul[,1:3],true=true)
  } else {
   plot(result,param=simul[,1:3])
  }
  
  u=input()
  
  par(mfrow=c(1,4))
  for (i in 1:3){

    plot(density(simul[,i]),xlim=c(min(simul[,i]),max(simul[,i])),ylim=c(0,max(density(result$adj.values[,i])$y,density(simul[,i])$y)),xlab=c('Theta','m','T')[i],main='')
    par(new=T)
    plot(density(result$adj.values[,i]),xlim=c(min(simul[,i]),max(simul[,i])),ylim=c(0,max(density(result$adj.values[,i])$y,density(simul[,i])$y)),col='red',xaxt='n',xlab='',main='')
  }

  plot(result$adj.values[,2],result$adj.values[,3],xlab='m',ylab='T',main=paste('correlation between\n m and T:',signif(cor(result$adj.values[,2],result$adj.values[,3]),3)))

  
  return(result)
}  