# .libPaths(new='/home/lv70700/blanckaert/libR')
library(scrm)

calculate_stat_4pop_4Djsfs_v2<-function(data)
{
  jsfs4D <- array(0,dim=c(4,4,4,4))
  for(M in data) {
    if(length(M)>0) {
      for(i in 1:(dim(M)[2])) {
	if (mean(M[3:626,i])<=.5){
	  j <- sum(M[3:222,i])
	  h <- sum(M[223:344,i])
	  k <- sum(M[345:464,i])
	  l <- sum(M[465:626,i])
	} else {
	  j <-220- sum(M[3:222,i])
	  h <-122- sum(M[223:344,i])
	  k <-120- sum(M[345:464,i])
	  l <-162- sum(M[465:626,i])
	}
	j=c(1,rep(2,11),rep(3,198),rep(4,11))[seq(221)==(j+1)]
	h=c(1,rep(2,6),rep(3,109),rep(4,7))[seq(123)==(h+1)]
	k=c(1,rep(2,6),rep(3,107),rep(4,7))[seq(121)==(k+1)]
	l=c(1,rep(2,8),rep(3,145),rep(4,9))[seq(163)==(l+1)]
	jsfs4D[j,h,k,l] <- jsfs4D[j,h,k,l]+1
      }
    }
  }
  statistic <- c()
  for (j in 1:4) {
    for (h in 1:4) {
      statistic =c(statistic,c(t(jsfs4D[j,h,,])))
    }
  }
  statistic <- statistic[2:256] # remove first one because no polymorphism (even if not 0, if singleton in one of the 2 individuals)

  return(statistic)
}



calculate_stat_4pop_4Djsfs_v5<-function(data)
{
  deltaF=read.table('jsfs_identity.csv',sep='\t',fill=NA)
  jsfs4D <- array(0,dim=c(4,4,4,4))
  fst_table=c(0,0,0,0,0,0,0)
  similar_freq=c(0,0,0,0,0,0,0)
  for(M in data) {
    if(length(M)>0) {
      for(i in 1:(dim(M)[2])) {
	if (mean(M[3:626,i])<=.5){
	  j <- sum(M[3:222,i])
	  h <- sum(M[223:344,i])
	  k <- sum(M[345:464,i])
	  l <- sum(M[465:626,i])
	} else {
	  j <-220- sum(M[3:222,i])
	  h <-122- sum(M[223:344,i])
	  k <-120- sum(M[345:464,i])
	  l <-162- sum(M[465:626,i])
	}
	
	# order fst : total, 12, 13 ,14, 23, 24, 34
	fst_marker=c(1-(j*(1-j/220)+h*(1-h/122)+k*(1-k/120)+l*(1-l/162))/((j+h+k+l)*(1-(j+h+k+l)/624)), 1-(j*(1-j/220)+h*(1-h/122))/((j+h)*(1-(j+h)/342)), 1-(j*(1-j/220)+k*(1-k/120))/((j+k)*(1-(j+k)/340)), 1-(j*(1-j/220)+l*(1-l/162))/((j+l)*(1-(j+l)/382)), 1-(h*(1-h/122)+k*(1-k/120))/((h+k)*(1-(h+k)/242)), 1-(h*(1-h/122)+l*(1-l/162))/((h+l)*(1-(h+l)/284)), 1-(k*(1-k/120)+l*(1-l/162))/((k+l)*(1-(k+l)/282)))
	fst_table=rbind(fst_table,fst_marker)
		
	similar_freq=similar_freq+c(deltaF[1,j+1]<=h & h<=deltaF[2,j+1]&deltaF[3,j+1]<=k & k<=deltaF[4,j+1]&deltaF[5,j+1]<=l & l<=deltaF[6,j+1]&deltaF[7,h+1]<=k & k<=deltaF[8,h+1]&deltaF[9,h+1]<=l & l<=deltaF[10,h+1]&deltaF[11,k+1]<=l & l<=deltaF[12,k+1],(h+j!=0)& deltaF[1,j+1]<=h & h<=deltaF[2,j+1],(k+j!=0) & deltaF[3,j+1]<=k & k<=deltaF[4,j+1], (l+j!=0)& deltaF[5,j+1]<=l & l<=deltaF[6,j+1],(h+k!=0)& deltaF[7,h+1]<=k & k<=deltaF[8,h+1],(h+l!=0)& deltaF[9,h+1]<=l & l<=deltaF[10,h+1],(k+l!=0)&deltaF[11,k+1]<=l & l<=deltaF[12,k+1])
	
	
	
	j=c(1,rep(2,11),rep(3,198),rep(4,11))[seq(221)==(j+1)]
	h=c(1,rep(2,6),rep(3,109),rep(4,7))[seq(123)==(h+1)]
	k=c(1,rep(2,6),rep(3,107),rep(4,7))[seq(121)==(k+1)]
	l=c(1,rep(2,8),rep(3,145),rep(4,9))[seq(163)==(l+1)]
	jsfs4D[j,h,k,l] <- jsfs4D[j,h,k,l]+1
      }
    }
  }
  fst_table=fst_table[-1,]
  table2=fst_table
  table2[is.nan(table2)]<-0

  
  statistic <- c()
  for (j in 1:4) {
    for (h in 1:4) {
      statistic =c(statistic,c(t(jsfs4D[j,h,,])))
    }
  }
  statistic=c(statistic,colMeans(fst_table,na.rm=T),colMeans(table2),similar_freq)
  statistic <- statistic[2:277] # remove first one because no polymorphism (even if not 0, if singleton in one of the 2 individuals)

  return(statistic)
}


calculate_stat_4pop_4Djsfs_v6<-function(data)
{
  deltaF=read.table('jsfs_identity.csv',sep='\t',fill=NA)
  jsfs4D <- array(0,dim=c(4,4,4,4))
  fst_table=c(0,0,0,0,0,0,0)
  similar_freq=c(0,0,0,0,0,0,0)
  for(M in data) {
    if(length(M)>0) {
      for(i in 1:(dim(M)[2])) {
	if (mean(M[3:626,i])<=.5){
	  j <- sum(M[3:222,i])
	  h <- sum(M[223:344,i])
	  k <- sum(M[345:464,i])
	  l <- sum(M[465:626,i])
	} else {
	  j <-220- sum(M[3:222,i])
	  h <-122- sum(M[223:344,i])
	  k <-120- sum(M[345:464,i])
	  l <-162- sum(M[465:626,i])
	}
	
	# order fst : total, 12, 13 ,14, 23, 24, 34
	fst_marker=c(1-(j*(1-j/220)+h*(1-h/122)+k*(1-k/120)+l*(1-l/162))/((j+h+k+l)*(1-(j+h+k+l)/624)), 1-(j*(1-j/220)+h*(1-h/122))/((j+h)*(1-(j+h)/342)), 1-(j*(1-j/220)+k*(1-k/120))/((j+k)*(1-(j+k)/340)), 1-(j*(1-j/220)+l*(1-l/162))/((j+l)*(1-(j+l)/382)), 1-(h*(1-h/122)+k*(1-k/120))/((h+k)*(1-(h+k)/242)), 1-(h*(1-h/122)+l*(1-l/162))/((h+l)*(1-(h+l)/284)), 1-(k*(1-k/120)+l*(1-l/162))/((k+l)*(1-(k+l)/282)))
	fst_table=rbind(fst_table,fst_marker)
		
	similar_freq=similar_freq+c(deltaF[1,j+1]<=h & h<=deltaF[2,j+1]&deltaF[3,j+1]<=k & k<=deltaF[4,j+1]&deltaF[5,j+1]<=l & l<=deltaF[6,j+1]&deltaF[7,h+1]<=k & k<=deltaF[8,h+1]&deltaF[9,h+1]<=l & l<=deltaF[10,h+1]&deltaF[11,k+1]<=l & l<=deltaF[12,k+1],(h+j!=0)& deltaF[1,j+1]<=h & h<=deltaF[2,j+1],(k+j!=0) & deltaF[3,j+1]<=k & k<=deltaF[4,j+1], (l+j!=0)& deltaF[5,j+1]<=l & l<=deltaF[6,j+1],(h+k!=0)& deltaF[7,h+1]<=k & k<=deltaF[8,h+1],(h+l!=0)& deltaF[9,h+1]<=l & l<=deltaF[10,h+1],(k+l!=0)&deltaF[11,k+1]<=l & l<=deltaF[12,k+1])
	
	
	
	j=c(1,rep(2,11),rep(3,198),rep(4,11))[seq(221)==(j+1)]
	h=c(1,rep(2,6),rep(3,109),rep(4,7))[seq(123)==(h+1)]
	k=c(1,rep(2,6),rep(3,107),rep(4,7))[seq(121)==(k+1)]
	l=c(1,rep(2,8),rep(3,145),rep(4,9))[seq(163)==(l+1)]
	jsfs4D[j,h,k,l] <- jsfs4D[j,h,k,l]+1
      }
    }
  }
  fst_table=fst_table[-1,]

  stat45=c()
  stat45=c(stat45,sum(jsfs4D[,1,1,1]))
  stat45=c(stat45,sum(jsfs4D[1,,1,1]))
  stat45=c(stat45,sum(jsfs4D[1,1,,1]))
  stat45=c(stat45,sum(jsfs4D[1,1,1,]))
  stat45=c(stat45,sum(jsfs4D[2:4,2:4,2:4,2:4]))
  stat45=c(stat45,sum(jsfs4D[2:4,2:4,,]))
  stat45=c(stat45,sum(jsfs4D[2:4,,2:4,]))
  stat45=c(stat45,sum(jsfs4D[2:4,,,2:4]))
  stat45=c(stat45,sum(jsfs4D[,2:4,2:4,]))
  stat45=c(stat45,sum(jsfs4D[,2:4,,2:4]))
  stat45=c(stat45,sum(jsfs4D[,,2:4,2:4])) 
  stat45=c(stat45,sum(jsfs4D[2:4,2:4,1,1]))
  stat45=c(stat45,sum(jsfs4D[2:4,1,2:4,1]))
  stat45=c(stat45,sum(jsfs4D[2:4,1,1,2:4]))
  stat45=c(stat45,sum(jsfs4D[1,2:4,2:4,1]))
  stat45=c(stat45,sum(jsfs4D[1,2:4,1,2:4]))
  stat45=c(stat45,sum(jsfs4D[1,1,2:4,2:4]))
  
  statistic=c(stat45,colMeans(fst_table,na.rm=T),similar_freq)
 
  return(statistic)
}

simulator_fixed_differences_4pop <- function(theta,m12,m23,m34,t1,t2,t3) {
  filteredloci <- list()
  counter.filtered <- 215
  while(counter.filtered){
     
    scrm.out <- scrm(paste("628 ", counter.filtered," -t ", theta," -I 4 222 122 120 164 0 -m 1 2 ", m12,
              " -m 2 1 ",m12, "-m 3 2 ",m23," -m 2 3 ",m23,
                           " -m 3 4 ",m34,
                           " -m 4 3 ",m34,
              " -ej ", t1, " 1 2 "," -ej ", t2, " 4 3 "," -ej ", t3, " 1 4 ",collapse=""))
    for (i in seq(counter.filtered)){
      if(length(scrm.out$seg_sites[[i]])) {
        M <- as.matrix(scrm.out$seg_sites[[i]][c(1,2,627,628),])
        criterion <- any(M[1,]==M[2,] & M[1,]!=M[3,] &  M[3,]==M[4,])
        if(criterion) {
          filteredloci[[counter.filtered]] <- as.matrix(scrm.out$seg_sites[[i]])
          counter.filtered <- counter.filtered - 1
        }
      }
    }
  }
  
   statistic<-calculate_stat_4pop_4Djsfs_v5(filteredloci)
   return(statistic)
#   return(filteredloci)
}

simulator_fixed_differences_4pop_jaatha <- function(theta,m12,m23,m34,t1,t2,t3) {
  filteredloci <- list()
  counter.filtered <- 215
  while(counter.filtered){
     
    scrm.out <- scrm(paste("628 ", counter.filtered," -t ", theta," -I 4 222 122 120 164 0 -m 1 2 ", m12,
              " -m 2 1 ",m12, "-m 3 2 ",m23," -m 2 3 ",m23,
                           " -m 3 4 ",m34,
                           " -m 4 3 ",m34,
              " -ej ", t1, " 1 2 "," -ej ", t2, " 4 3 "," -ej ", t3, " 1 4 ",collapse=""))
    for (i in seq(counter.filtered)){
      if(length(scrm.out$seg_sites[[i]])) {
        M <- as.matrix(scrm.out$seg_sites[[i]][c(1,2,627,628),])
        criterion <- any(M[1,]==M[2,] & M[1,]!=M[3,] &  M[3,]==M[4,])
        if(criterion) {
          filteredloci[[counter.filtered]] <- as.matrix(scrm.out$seg_sites[[i]])
          counter.filtered <- counter.filtered - 1
        }
      }
    }
  }
   statistic<-calculate_stat_4pop_4Djsfs_v6(filteredloci)
   print("done")
   return(statistic)
}

simulator_fixed_differences_4pop_verbose <- function(theta,m12,m23,m34,t1,t2,t3) {
  filteredloci <- list()
  counter.filtered <- 215
  while(counter.filtered){
    print("call_scrm")
    print(counter.filtered)
    print(paste("628 ", counter.filtered," -t ", theta," -I 4 222 122 120 164 0 -m 1 2 ", m12,
              " -m 2 1 ",m12, "-m 3 2 ",m23," -m 2 3 ",m23,
                           " -m 3 4 ",m34,
                           " -m 4 3 ",m34,
              " -ej ", t1, " 1 2 "," -ej ", t2, " 4 3 "," -ej ", t3, " 1 4 ",collapse=""))
    scrm.out <- scrm(paste("628 ", counter.filtered," -t ", theta," -I 4 222 122 120 164 0 -m 1 2 ", m12,
              " -m 2 1 ",m12, "-m 3 2 ",m23," -m 2 3 ",m23,
                           " -m 3 4 ",m34,
                           " -m 4 3 ",m34,
              " -ej ", t1, " 1 2 "," -ej ", t2, " 4 3 "," -ej ", t3, " 1 4 ",collapse=""))
    print("successful_call")
    
    for (i in seq(counter.filtered)){
      if(length(scrm.out$seg_sites[[i]])) {
        M <- as.matrix(scrm.out$seg_sites[[i]][c(1,2,627,628),])
        criterion <- any(M[1,]==M[2,] & M[1,]!=M[3,] &  M[3,]==M[4,])
        if(criterion) {
          filteredloci[[counter.filtered]] <- as.matrix(scrm.out$seg_sites[[i]])
          counter.filtered <- counter.filtered - 1
        }
      }
    }
  }
  
   statistic<-calculate_stat_4pop_4Djsfs_v5(filteredloci)
   return(statistic)
#   return(filteredloci)
}

simulator_fixed_differences_4pop_v2 <- function(theta,m12,m23,m34,mt1,mt2,t1,t2,t3) {
  filteredloci <- list()
  counter.filtered <- 215
  while(counter.filtered){
     
    scrm.out <- scrm(paste("628 ", counter.filtered," -t ", theta," -I 4 222 122 120 164 0 -m 1 2 ", m12,
              " -m 2 1 ",m12, "-m 3 2 ",m23," -m 2 3 ",m23,
                           " -m 3 4 ",m34,
                           " -m 4 3 ",m34,
              " -ej ", t1, " 1 2 ", " -em ", t1," 2 3 ", mt1 , " -em ", t1, "3 2 ", mt1 , " -ej ", t2, " 4 3 ", " -em", t2, " 2 3 ", mt2, " -em ", t2, " 3 2 ", mt2,  " -ej ", t3, " 1 4 ",collapse=""))
    for (i in seq(counter.filtered)){
      if(length(scrm.out$seg_sites[[i]])) {
        M <- as.matrix(scrm.out$seg_sites[[i]][c(1,2,627,628),])
        criterion <- any(M[1,]==M[2,] & M[1,]!=M[3,] &  M[3,]==M[4,])
        if(criterion) {
          filteredloci[[counter.filtered]] <- as.matrix(scrm.out$seg_sites[[i]])
          counter.filtered <- counter.filtered - 1
        }
      }
    }
  }
  
  statistic<-calculate_stat_4pop_4Djsfs_v5(filteredloci)
  return(statistic)
#   return(filteredloci)
}
