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
  

display_true_data=function(data,simul,n){
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

  
library(abc)

param_name=c('theta','m', 't')

table0a=read.table('simul_1pop_SNP.txt')

table0b=read.table('simul_1pop_SNP_v2.txt')

table1=read.table('simul_2pop_SNP_vA1.txt')
table2=read.table('simul_2pop_SNP_vA2.txt')

table=rbind(table1,table2)

table3=read.table('simul_2pop_SNP_vB1.txt')

table4=read.table('simul_2pop_SNP_vC1.txt')

table5=read.table('simul_2pop_SNP_vE1.txt')

v0a= display_true_data(almost_real_data,table0a,1)

v0b= display_true_data(real_data,table0b,1)

v1= display_true_data(almost_real_data,table,3)

v3= display_true_data(almost_real_data,table3,3)

v4= display_true_data(almost_real_data,table4,3)

v5= display_true_data(almost_real_data,table5,3)


v1f= display_true_data(real_data,table,3)
# 1000 cv for 60000 simuls, tol=c(.005,.01,.02,05)

table1bis=table[table$V3<=1,]
 cor(table1bis[,1:3])
 


#--------------------test---------------
almost_real_data=c(328 ,165,   2  , 0, 429 ,101 ,172 ,  3  , 0, 195 , 85, 309 ,  3 ,  3  , 0 ,  0 , 0,0 ,  0  , 0  , 0 ,  0  , 0)
real_data=c(257 ,132,   2  , 0, 355 ,71 ,140 ,  3  , 0, 147 , 70, 253 ,  3 ,  3  , 0 ,  0 , 0,0 ,  0  , 0  , 0 ,  0  , 0)
 v1bis= display_true_data(almost_real_data,table1bis,3)
#almost_real_data=c(328 ,165,   2  , 0, 429 ,101 ,172 ,  3  , 0, 195 , 85, 309 ,  3 ,  3  , 0 ,  0 , 0,0 ,  0  , 0  , 0 ,  0  , 0)[max_table[4:26]>0]

 # test2<-simulator_fixed_differences_2pop(.8,.8,.15)
 

# test2=c(383, 216,   0,   0, 414,  30, 106,   0,   0, 221 , 69, 327,   2 ,  0 ,  0 ,  0  , 0 ,  0,   0, 0 ,  0,   0 ,  0)[max_table[4:26]>0]

# test 3 <- simulator_fixed_differences_2pop(.3,1.76,1.42)
# test3=c(139 ,78, 0,  0,162 ,20 ,85, 2,  0, 82, 68,271, 3,  2,  0,  0, 0 , 0 , 0,0 , 0 , 0,  0)[max_table[4:26]>0]



use_table=table[,c(1:3,v1+3)]

use_table2=table[, c(1,2,3,4,5,8,9,10,13,14,15)] # keep S only fi the observed values is within the .1 -- .9 of the dist
use_table3=table[, c(1,2,3,4,8,9,10,14,15)] # keep S pnly if not weird right tail
# use_table=use_table2

colnames(use_table2)<-c('theta','m','T','S1','S2','S5','S6','S7','S10','S11','S12')
colnames(use_table)<-c('theta','m','T','S1','S2','S3','S5','S6','S7','S8','S10','S11','S12','S13','S14')
colnames(use_table3)<-c('theta','m','T','S1','S5','S6','S7','S11','S12')  

  
summary(cv4abc(param=use_table[,1:3],sumstat=use_table[,4:dim(use_table)[2]],tols=c(.001,.005,.01,.05), method="loclinear",nval=500))
  
summary(cv4abc(param=use_table[,1:3],sumstat=use_table[,4:dim(use_table)[2]],tols=c(.001,.005,.01,.05), method="rejection",nval=500))

summary(cv4abc(param=use_table2[,1:3],sumstat=use_table2[,4:dim(use_table2)[2]],tols=c(.001,.005,.01,.05), method="loclinear",nval=500))
  
summary(cv4abc(param=use_table2[,1:3],sumstat=use_table2[,4:dim(use_table2)[2]],tols=c(.001,.005,.01,.05), method="rejection",nval=500))

summary(cv4abc(param=use_table3[,1:3],sumstat=use_table3[,4:dim(use_table3)[2]],tols=c(.001,.005,.01,.05), method="loclinear",nval=500))# best way to go
# 
# Prediction error based on a cross-validation sample of 500
# 
#                V1          V2          V3
# 0.001 0.005619178 0.080826536 0.389001508
# 0.005 0.005830721 0.076793209 0.397060528
# 0.01  0.005844993 0.077876314 0.399413137
# 0.05  0.006389909 0.085256143 0.420262147


summary(cv4abc(param=use_table3[,1:3],sumstat=use_table3[,4:dim(use_table3)[2]],tols=c(.001,.005,.01,.05), method="rejection",nval=500))

summary(cv4abc(param=table0a[,1],sumstat=table0a[,c( 2 , 5 , 6 , 7 ,10 ,11, 12)],tols=c(,.005,.01,.05,.1), method="loclinear",nval=500))#decent prediction of theta


summary(cv4abc(param=table0b[,1],sumstat=table0b[,c(4,5,8,9,10,13,14,15)-2],tols=c(,.005,.01,.05,.1), method="loclinear",nval=500))#decent prediction of theta with more sim
#a=abc(almost_real_data[c(4,8,9,10,14,15)-3],param=use_table2[,1:3],sumstat=use_table3[,4:11],tol=.01,method="loclinear")


a=abc(almost_real_data[c( 2 , 5 , 6 , 7 ,10 ,11, 12)-1],param=table0b[,1],sumstat=table0b[,c( 2 , 5 , 6 , 7 ,10 ,11, 12)],tol=.05,method="loclinear")



a=plot_result(almost_real_data[c(4,5,8,9,10,13,14,15)-3],use_table2,.05,F)
a=plot_result(almost_real_data[c(4,5,8,9,10,13,14,15)-3],use_table2,.01,F)
a=plot_result(almost_real_data[c(4,5,8,9,10,13,14,15)-3],use_table2,.005,F)
a=plot_result(almost_real_data[c(4,5,8,9,10,13,14,15)-3],use_table2,.001,F)

a=plot_result(real_data[c(4,5,8,9,10,13,14,15)-3],use_table2,.05,F)
a=plot_result(real_data[c(4,5,8,9,10,13,14,15)-3],use_table2,.01,F)
a=plot_result(real_data[c(4,5,8,9,10,13,14,15)-3],use_table2,.005,F)
a=plot_result(real_data[c(4,5,8,9,10,13,14,15)-3],use_table2,.001,F)

a=plot_result(real_data[c(4,5,8,9,10,13,14,15)-3],use_table2bis,.005,F)

a=plot_result(real_data[c(4,5,8,9,10,13,14,15)-3],table0b[,c(3,4,5,8,9,10,13,14,15)-2],.01,F)

a=abc(real_data[c(4,5,8,9,10,13,14,15)-3],param=table0b[,1],table0b[,c(4,5,8,9,10,13,14,15)-2],tol=.05,method="loclinear")

plot(density(table0b[,1]),xlim=c(min(table0b[,1]),max(table0b[,1])),ylim=c(0,max(density(a$adj.values[,1])$y,density(table0b[,1])$y)),xlab=c('Theta','m','T')[1],main='')
par(new=T)
plot(density(a$adj.values[,1]),xlim=c(min(table0b[,1]),max(table0b[,1])),ylim=c(0,max(density(a$adj.values[,1])$y,density(table0b[,1])$y)),col='red',xaxt='n',xlab='',main='')
  
  
par(mfrow=c(1,4))
for (i in 1:3){

  plot(density(use_table[,i]),xlim=c(min(use_table[,i]),max(use_table[,i])),ylim=c(0,max(density(a$adj.values[,i])$y,density(use_table[,i])$y)),xlab=c('Theta','m','T')[i],main='')
  par(new=T)
  plot(density(a$adj.values[,i]),xlim=c(min(use_table[,i]),max(use_table[,i])),ylim=c(0,max(density(a$adj.values[,i])$y,density(use_table[,i])$y)),col='red',xaxt='n',xlab='',main='')
}

plot(a$adj.values[,2],a$adj.values[,3],xlab='m',ylab='T',main=paste('correlation between\n m and T:',signif(cor(a$adj.values[,2],a$adj.values[,3]),3)))


par(mfrow=c(1,4))
for (i in 1:3){

  plot(density(use_table[,i]),xlim=c(min(use_table[,i]),max(use_table[,i])),ylim=c(0,max(density(a$unadj.values[,i])$y,density(use_table[,i])$y)),xlab=c('Theta','m','T')[i],main='')
  par(new=T)
  plot(density(a$unadj.values[,i]),xlim=c(min(use_table[,i]),max(use_table[,i])),ylim=c(0,max(density(a$unadj.values[,i])$y,density(use_table[,i])$y)),col='red',xaxt='n',xlab='',main='')
}
plot(a$unadj.values[,2],a$unadj.values[,3],xlab='m',ylab='T',main=paste('correlation between\n m and T:',signif(cor(a$unadj.values[,2],a$unadj.values[,3]),3)))



#summary(cv4abc(param=use_table[,1:3],sumstat=use_table[,4:dim(use_table)[2]],tols=c(.01,.03,.05,.1), method="loclinear",nval=500))

use_table1bis=table1bis[,c(1:3,v1bis+3)]

use_table2bis=table1bis[, c(1,2,3,4,5,8,9,10,13,14,15)] # keep S only fi the observed values is within the .1 -- .9 of the dist
use_table3bis=table1bis[, c(1,2,3,4,8,9,10,14,15)]
 
summary(cv4abc(param=use_table2bis[,1:3],sumstat=use_table2bis[,4:dim(use_table2bis)[2]],tols=c(.001,.005,.01,.05), method="loclinear",nval=500))
  
summary(cv4abc(param=use_table2bis[,1:3],sumstat=use_table2bis[,4:dim(use_table2bis)[2]],tols=c(.001,.005,.01,.05), method="rejection",nval=500))
                                                                                                                                                                           

a=plot_result(almost_real_data[c(4,5,8,9,10,13,14,15)-3],use_table2bis,.01,F)
-------------------------------------------------------------------------------------
# asym mig

asym_sim=read.table('output_sym_10e6.txt')


sub_select=asym_sim[, c(1,2,3,4,5,6,9,10,11,14,15,16)]

summary(cv4abc(param=sub_select[,1:4],sumstat=sub_select[,5:dim(sub_select)[2]],tols=c(.001,.005,.01,.05), method="loclinear",nval=50))
# Prediction error based on a cross-validation sample of 500
# 
#                V1          V2          V3          V4
# 0.001 0.004215712 0.074236494 0.087316696 1.149918178
# 0.005 0.003969473 0.068891809 0.083783161 1.035321352
# 0.01  0.004019862 0.070713885 0.083986720 1.023483403
# 0.05  0.004352191 0.074931268 0.088882151 1.015606420

a=abc(real_data[c(5,6,9,10,11,14,15,16)-4],param=sub_select[,1:4],asym_sim[, c(5,6,9,10,11,14,15,16)],tol=.001,method="loclinear")
plot(a,param=sub_select[,1:4])
summary(a)
# Call: 
# abc(target = real_data[c(5, 6, 9, 10, 11, 14, 15, 16) - 4], param = sub_select[, 
#     1:4], sumstat = asym_sim[, c(5, 6, 9, 10, 11, 14, 15, 16)], 
#     tol = 0.001, method = "loclinear")
# Data:
#  abc.out$adj.values (100 posterior samples)
# Weights:
#  abc.out$weights
# 
#                            V1     V2     V3     V4
# Min.:                  0.0436 1.0168 1.1592 3.0767
# Weighted 2.5 % Perc.:  0.2948 1.0904 1.2960 3.3109
# Weighted Median:       0.5132 1.2695 1.5053 3.5835
# Weighted Mean:         0.5163 1.2883 1.5196 3.5814
# Weighted Mode:         0.5089 1.2141 1.4990 3.5870
# Weighted 97.5 % Perc.: 0.7206 1.5824 1.7621 3.8154
# Max.:                  0.7327 1.7470 2.2136 3.9758

model=c(rep(1,76125),rep(2,100000))
colnames(use_table2)<-c('theta','m','T','S1','S2','S5','S6','S7','S10','S11','S12')
colnames(sub_select)<-c('theta','m12','m21','T','S1','S2','S5','S6','S7','S10','S11','S12')
rbind(sub_select[, 5:dim(sub_select)[2]],use_table2[,4:dim(use_table2)[2]])

cv.modsel =cv4postpr(model,rbind(use_table2[,4:dim(use_table2)[2]],sub_select[, 5:dim(sub_select)[2]]),nval=1000,tol=.001,method='mnlogistic')

summary(cv.modsel)
Confusion matrix based on 1000 samples for each model.

$tol0.001
    1   2
1 811 189
2 221 779


Mean model posterior probabilities (mnlogistic)

$tol0.001
       1      2
1 0.7349 0.2651
2 0.2319 0.7681


bb=postpr(target = real_data[c(5,6,9,10,11,14,15,16)-4], index = model,sumstat=rbind(use_table2[,4:dim(use_table2)[2]],sub_select[, 5:dim(sub_select)[2]]),tol = 0.001, method = "mnlogistic")

summary(bb)
# Call: 
# postpr(target = real_data[c(5, 6, 9, 10, 11, 14, 15, 16) - 4], 
#     index = model, sumstat = rbind(use_table2[, 4:dim(use_table2)[2]], 
#         sub_select[, 5:dim(sub_select)[2]]), tol = 0.001, method = "mnlogistic")
# Data:
#  postpr.out$values (177 posterior samples)
# Models a priori:
#  1, 2
# Models a posteriori:
#  1, 2
# Warning: Posterior model probabilities are corrected for unequal number of simulations per models.
# 
# 
# Proportion of accepted simulations (rejection):
#      1      2 
# 0.8748 0.1252 
# 
# Bayes factors:
#        1      2
# 1 1.0000 5.3214
# 2 0.1879 1.0000
# 
# 
# Posterior model probabilities (mnlogistic):
#      1      2 
# 0.9847 0.0153 
# 
# Bayes factors:
#         1       2
# 1  1.0000 64.4228
# 2  0.0155  1.0000


a=abc(real_data[c(5,6,9,10,11,14,15,16)-4],param=use_table2[,1:3],use_table2[, 4:dim(use_table2)[2]],tol=.001,method="loclinear")

sub_select2=asym_sim[,c(1,2,3,4,5,9,10,11,15,16)]

summary(cv4abc(param=sub_select2[,1:4],sumstat=sub_select2[,5:dim(sub_select2)[2]],tols=c(.001,.005,.01,.05), method="loclinear",nval=50))
# Prediction error based on a cross-validation sample of 500
# 
#                V1          V2          V3          V4
# 0.001 0.004423523 0.168071770 0.201207627 1.065368196
# 0.005 0.004317429 0.164054264 0.191483444 1.020318077
# 0.01  0.004311560 0.166376091 0.193397901 1.014682320
# 0.05  0.004457035 0.177667514 0.212946113 1.003919103

a=abc(real_data[c(5,9,10,11,15,16)-4],param=sub_select2[,1:4],sub_select2[,5:dim(sub_select2)[2]],tol=.001,method="loclinear")
summary(a)
Call: 
abc(target = real_data[c(5, 6, 9, 10, 11, 14, 15, 16) - 4], param = sub_select2[, 
    1:4], sumstat = asym_sim[, c(5, 6, 9, 10, 11, 14, 15, 16)], 
    tol = 0.005, method = "loclinear")
Data:
 abc.out$adj.values (500 posterior samples)
Weights:
 abc.out$weights

                            V1      V2      V3      V4
Min.:                   0.4455  0.7082 -1.3033  1.4378
Weighted 2.5 % Perc.:   0.4735  0.8671  0.0517  1.8627
Weighted Median:        0.5240  1.0487  1.8599  3.1695
Weighted Mean:          0.5230  1.0520  1.9584  3.2024
Weighted Mode:          0.5322  1.0487  1.7890  3.8749
Weighted 97.5 % Perc.:  0.5712  1.2675  4.4820  4.4465
Max.:                   0.6091  1.7161  7.4592  4.9106


test=gfit(target=real_data[c(5,6,9,10,11,14,15,16)-4], sumstat=asym_sim[, c(5,6,9,10,11,14,15,16)],tol=.001, nb.replicate=100)
plot(test)


