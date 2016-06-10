#source('/Users/Alex/Dropbox/Bluebells/simulator2.R')
# model 1
args=(commandArgs(TRUE))


if(length(args)==0){
    stop("No arguments supplied.")
} else{
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
    print(seed)
    print(output)
}

set.seed(seed)
source('simulator_cluster.R')
N=2500
param_1=cbind(runif(N,0.2,.8),runif(N,0,20),runif(N,0,50),runif(N,0,20),runif(N,0,2),runif(N,0,2),runif(N,0,2))
for (i in seq(N)){
  param_1[i,7]=max(param_1[i,5],param_1[i,6])+param_1[i,7]
  test<-simulator_fixed_differences_4pop(param_1[i,1],param_1[i,2],param_1[i,3],param_1[i,4],param_1[i,5],param_1[i,6],param_1[i,7])
  test_mat=c(param_1[i,],test)
  write(test_mat,output,sep='\t', append =TRUE, ncolumns = 283)}

#test_mat=cbind(param_1,stat_ex)
#write.table(test_mat,'simul_fixed_differences.txt',sep='\t', append =TRUE, quote = FALSE, col.names = FALSE, row.names=FALSE)
#write.table(test_mat,'/Users/Alex/Dropbox/Bluebells/simul_fixed_differences.txt',sep='\t', append =TRUE, quote = FALSE, col.names = FALSE, row.names=FALSE)

