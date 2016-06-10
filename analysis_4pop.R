source('function_analysis_abc.R')


if (T){
param_name=c('theta','m12','m23','m34','t1','t2','t3')

  table1=read.table('result/model1_fullstat_set1.txt')
  table1bis=read.table('result/model1_fullstat_set2.txt')
  table1=rbind(table1,table1bis)
  table1bis=read.table('result/model1_fullstat_set3.txt')
  table1=rbind(table1,table1bis)
  table1[table1[,6]<table1[,5],7]=table1[table1[,6]<table1[,5],7]-table1[table1[,6]<table1[,5],5]
  table1[table1[,6]>table1[,5],7]=table1[table1[,6]>=table1[,5],7]-table1[table1[,6]>=table1[,5],6]
  
  mask=filter_table_v2(table1[1:262],7,0)
  mask=c(mask,seq(263,269),seq(277,283))
  test=table1[,mask]
  z=cv4abc(param=test[,1:7], sumstat=test[,8:dim(test)[2]],tols=.01,method="loclinear",nval=50)
  best=sum(summary(z)[1:7])
  ind=1
  while (ind>0){
    ind=0
    for (i in 8:dim(test)[2]){
      print(i)
      z=cv4abc(param=test[,1:7], sumstat=test[,8:dim(test)[2]][,-i],tols=.01,method="loclinear",nval=50)
      if (sum(summary(z)[1:7])<best){
	ind=i
	best=sum(summary(z)[1:7])
      }
    }
    if(ind>0){test=test[,-ind]}
    print(dim(test))
  } # not working because cv4abc call summary()
  
  

  w=abc(c(real_data_4D,fst_real_data_4D,diagonal_jsfs_4D)[mask[8:length(mask)]-7],param=test[,1:7],test[,8:dim(test)[2]],tol=.01,method="loclinear")  
  
  pca=prcomp(test[,8:dim(test)[2]])
  z=cv4abc(param=test[,1:7], sumstat=pca$x[,1:3],tols=.01,method="loclinear",nval=200)
  
  
  if (F){
    new_table6=data.frame()
    new_table6=rbind(new_table6,seq(38))
    for (i in 1:dim(table1)[1]){
      new_table6=rbind(new_table6,c(as.double(table1[i,1:7]),convert_to_statv4(c(table1[i,8:262])),c(as.double(table1[i,263:269])),c(as.double(table1[i,277:283]))))
      }
    new_table6=new_table6[-1,]
    
    write.table(new_table6,col.names=F,file="result/joachim_stat_model1_fullstat_set1_2.txt",sep='\t',row.names=F,append=T)
    }

  new_table6=read.table("result/joachim_stat_model1_fullstat_set1_2.txt",sep='\t')
  
  mod_real_data=c(c(convert_to_statv4(real_data_4D),fst_real_data_4D[1:7],diagonal_jsfs_4D))
  
  j6=abc(mod_real_data,param=new_table6[,1:7],new_table6[,8:dim(new_table6)[2]],tol=.01,method="loclinear")
  
  
  filter1=c(rep(T,11),rep(F,6),rep(T,13))
  j7=abc(mod_real_data[filter1],param=new_table6[,1:7],new_table6[,8:dim(new_table6)[2]][,filter1],tol=.01,method="loclinear")
  plot(j7,param=new_table6[,1:7])
  analysis_4pop(new_table6[,c(rep(T,7),filter1)],mod_real_data[filter1],tolerance=.002,sample_cv=200,estimate=F)
  
  
  
  
  filter2=c(rep(T,11),rep(F,6),rep(T,7),rep(F,6))
  j8=abc(mod_real_data[filter2],param=new_table6[,1:7],new_table6[,8:dim(new_table6)[2]][,filter2],tol=.01,method="loclinear") # sound promising
  analysis_4pop(new_table6[,c(rep(T,7),filter2)],mod_real_data[filter2],tolerance=.01,sample_cv=200,estimate=F)
  
#                 V1          V2          V3          V4          V5          V6
# 0.01 0.008721517 0.043150319 0.026496582 0.051695751 0.663519836 0.518266358
#               V7
# 0.01 1.028621832


#               V1          V2          V3          V4          V5          V6
# 0.01 0.007733636 0.046840559 0.026255645 0.060703296 0.576852405 0.560013739
#               V7
# 0.01 1.001109778

  check_tol=cv4abc(param=new_table6[,1:7], sumstat=new_table6[,8:dim(new_table6)[2]][,filter2],tols=c(.0005,.001,.005,.01),method="loclinear",nval=200)
#   summary(check_tol)
  
  signif(summary(check_tol),3)
Prediction error based on a cross-validation sample of 200

           V1      V2      V3      V4      V5      V6      V7
5e-04 0.01310 0.03290 0.03270 0.03180 0.87800 0.87400 1.65000
0.001 0.00949 0.02820 0.02760 0.02400 0.67800 0.71800 1.32000
0.005 0.0129  0.0305  0.0233  0.0405  0.5050  0.5180  1.0200
0.01  0.00825 0.02920 0.02960 0.02290 0.53400 0.56900 1.03000
0.05  0.00884 0.04350 0.04460 0.03670 0.54900 0.56600 1.02000


  
  
  random_abc(new_table6,tolerance=.01)
  
  filter3=c(rep(T,11),rep(F,6),rep(F,7),rep(T,6))
  j9=abc(mod_real_data[filter3],param=new_table6[,1:7],new_table6[,8:dim(new_table6)[2]][,filter3],tol=.01,method="loclinear")
  analysis_4pop(new_table6[,c(rep(T,7),filter3)],mod_real_data[filter3],tolerance=.002,sample_cv=200,estimate=F)
  
  
}# last version

if (T){
  table2=read.table('result/model2_fullstat_set1.txt')
  table2[table2[,8]<table2[,7],9]=table2[table2[,8]<table2[,7],9]-table2[table2[,8]<table2[,7],7]
  table2[table2[,8]>table2[,7],9]=table2[table2[,8]>=table2[,7],9]-table2[table2[,8]>=table2[,7],8]
  
  
    if (F){
    new_table7=data.frame()
    new_table7=rbind(new_table7,seq(40))
    for (i in 1:dim(table2)[1]){
      new_table7=rbind(new_table7,c(as.double(table2[i,1:9]),convert_to_statv4(c(table2[i,10:264])),c(as.double(table2[i,265:271])),c(as.double(table2[i,279:285]))))
      }
    new_table7=new_table7[-1,]
    
    write.table(new_table7,col.names=F,file="result/joachim_stat_model2_fullstat_set1.txt",sep='\t',row.names=F,append=T)
    }
    
#     filterv2=c(rep(T,1),rep(F,6),rep(T,7),rep(F,6))
    check_tol=cv4abc(param=new_table7[,1:9], sumstat=new_table7[,10:dim(new_table7)[2]][,filter2],tols=c(.005,.01),method="loclinear",nval=100)
Prediction error based on a cross-validation sample of 100

          X1L     X2L     X3L     X4L     X5L     X6L     X7L     X8L     X9L
0.005 0.01040 0.02380 0.03900 0.05850 1.25000 1.39000 0.71500 0.58800 1.42000
0.01  0.00816 0.02750 0.03980 0.04890 1.09000 1.20000 0.60400 0.48400 1.15000


 j12=abc(mod_real_data[filter2],param=new_table7[,1:9],new_table7[,10:dim(new_table7)[2]][,filter2],tol=.01,method="loclinear")
 plot(j12,param=new_table7[,1:9])
}

if (T){
  table1_cd=read.table('result/model1_comp_dist.txt')
  table1_cd[table1_cd[,6]<table1_cd[,5],7]=table1_cd[table1_cd[,6]<table1_cd[,5],7]-table1_cd[table1_cd[,6]<table1_cd[,5],5]
  table1_cd[table1_cd[,6]>table1_cd[,5],7]=table1_cd[table1_cd[,6]>=table1_cd[,5],7]-table1_cd[table1_cd[,6]>=table1_cd[,5],6]

  if (F){
    new_table6=data.frame()
    new_table6=rbind(new_table6,seq(38))
    for (i in 1:dim(table1_cd)[1]){
      new_table6=rbind(new_table6,c(as.double(table1_cd[i,1:7]),convert_to_statv4(c(table1_cd[i,8:262])),c(as.double(table1_cd[i,263:269])),c(as.double(table1_cd[i,277:283]))))
      }
    new_table6=new_table6[-1,]
    
    write.table(new_table6,col.names=F,file="result/joachim_stat_model1_comp_dist.txt",sep='\t',row.names=F,append=T)
  }
  
  
  filter2=c(rep(T,11),rep(F,6),rep(T,7),rep(F,6))
  j8=abc(mod_real_data[filter2],param=new_table6[,1:7],new_table6[,8:dim(new_table6)[2]][,filter2],tol=.01,method="loclinear") # sound promising
  analysis_4pop(new_table6[,c(rep(T,7),filter2)],mod_real_data[filter2],tolerance=.01,sample_cv=200,estimate=F)
 Prediction error based on a cross-validation sample of 200

            X1L        X2L        X3L        X4L        X5L        X6L
0.01 0.01525828 0.20233104 0.07177256 0.17264175 0.63449829 0.47591102
            X7L
0.01 1.09172351

  
}