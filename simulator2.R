 .libPaths(new='/home/lv70700/blanckaert/libR')
library(scrm)
#library(ape)

########################################################

### HERE COMES THE STATISTICS ##########################

########################################################

calculate_stat<-function(filteredloci,otherloci)
{
jsfs3D <- array(0,dim=c(101,201,101))
for(M in filteredloci) {
  if(length(M)>0) {
    for(i in 1:(dim(M)[2])) {
      j <- sum(M[3:102,i])
      h <- sum(M[103:302,i])
      k <- sum(M[303:402,i])
      jsfs3D[j+1,h+1,k+1] <- jsfs3D[j+1,h+1,k+1]+1
    }
  }
}
i <- 0
statistic <- numeric()
for(j in list(1,2:4,5:100)) {
  for(h in list(1,2:4,5:200)) {
    for(k in list(1,2:4,5:100)) {
      i <- i+1
      statistic[i] <- sum(jsfs3D[j,h,k])
    }
  }
}
statistic <- statistic[2:27]
statistic[27] <- sum(jsfs3D[101,201,1:100])
statistic[28] <- sum(jsfs3D[101,1:200,1:100])
statistic[29] <- sum(jsfs3D[101,1:200,101])
statistic[30] <- sum(jsfs3D[1:100,1:200,101])
statistic[31] <- sum(jsfs3D[1:100,201,101])
statistic[32] <- sum(jsfs3D[1:100,201,1:100])

jsfs3D <- array(0,dim=c(3,401,3))
for(M in filteredloci) {
  if(length(M)>0) {
    for(i in 1:(dim(M)[2])) {
      j <- sum(M[1:2,i])
      h <- sum(M[3:402,i])
      k <- sum(M[403:404,i])
      jsfs3D[j+1,h+1,k+1] <- jsfs3D[j+1,h+1,k+1]+1
    }
  }
}
i <- 32
for(j in 1:3) {
  for(h in list(1,2:4,5:400,401)) {
    for(k in 1:3) {
      i <- i+1
      statistic[i] <- sum(jsfs3D[j,h,k])
    }
  }
}
jsfs <- array(0,dim=c(3,3))
for(M in otherloci) {
  if(length(M)>0) {
    for(i in 1:(dim(M)[2])) {
      j <- sum(M[1:2,i])
      h <- sum(M[3:4,i])
      jsfs[j+1,h+1] <- jsfs[j+1,h+1]+1
    }
  }
}
statistic[69] <- jsfs[3,1]+jsfs[1,3]
statistic[70] <- jsfs[2,3]+jsfs[2,1]
statistic[71] <- jsfs[3,2]+jsfs[1,2]
statistic[72] <- jsfs[2,2]
return(statistic)
}

calculate_stat_2pop<-function(filteredloci)
{
  jsfs2D <- array(0,dim=c(191,107))
  for(M in filteredloci) {
    if(length(M)>0) {
      for(i in 1:(dim(M)[2])) {
	if (mean(M[3:298,i])<=.5){
	  j <- sum(M[3:192,i])
	  h <- sum(M[193:298,i])
	} else {
	  j <- 190-sum(M[3:192,i])
	  h <- 106-sum(M[193:298,i])
	}
	jsfs2D[j+1,h+1] <- jsfs2D[j+1,h+1]+1
      }
    }
  }

  i <- 0
  statistic <- numeric()
  for(j in list(1,2:10,11:181,182:190,191)) {
    for(h in list(1,2:6,7:101,102:106,107)) {
	i <- i+1
	statistic[i] <- sum(jsfs2D[j,h])
    }
  }
  statistic <- statistic[2:24]

  return(statistic)
}


calculate_stat_4pop_4Djsfs<-function(filteredloci)
{
  jsfs4D <- array(0,dim=c(221, 123, 121, 163))
  for(M in filteredloci) {
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
	jsfs4D[j+1,h+1,k+1,l+1] <- jsfs4D[j+1,h+1,k+1,l+1]+1
      }
    }
  }
  i <- 0
  statistic <- numeric()
  for (j in list(1,2:12,13:210,211:221)) {
    for (h in list(1,2:7,8:116,117:123)) {
      for (k in list(1,2:7,8:114,115:121)) {
	for (l in list(1,2:9,10:154,155:163)) {
	  i <- i+1
	  statistic[i] <- sum(jsfs4D[j,h,k,l])
	}
      }
    }
  }

  statistic <- statistic[2:256]

  return(statistic)
}

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


calculate_stat_4pop_4Djsfs_v3<-function(data)
{
  jsfs4D <- array(0,dim=c(4,4,4,4))
  fst_table=c(0,0,0,0,0,0,0)
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
  statistic=c(statistic,colMeans(fst_table,na.rm=T),colMeans(table2))
  statistic <- statistic[2:270] # remove first one because no polymorphism (even if not 0, if singleton in one of the 2 individuals)

  return(statistic)
}

calculate_stat_4pop_4Djsfs_v4<-function(data)
{
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
	
	delta=.1
	similar_freq=similar_freq+c(abs(j/220-h/122)< delta & abs(j/220-k/120)< delta & abs(j/220-l/162)< delta & abs(h/122-k/120)< delta & abs(h/122-l/162)< delta & abs(k/120-l/162)< delta,abs(j/220-h/122)< delta,abs(j/220-k/120)< delta,abs(j/220-l/162)< delta,abs(h/122-k/120)< delta,abs(h/122-l/162)< delta,abs(k/120-l/162)< delta)
	
	
	
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
  jsfs4D <- array(0,dim=c(6,6,6,6))
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
	
	
	
	j=c(1,2,rep(3,22),rep(4,66),rep(5,43),rep(6,88))[seq(221)==(j+1)]
	h=c(1,2,rep(3,12),rep(4,36),rep(5,25),rep(6,48))[seq(123)==(h+1)]
	k=c(1,2,rep(3,12),rep(4,36),rep(5,23),rep(6,48))[seq(121)==(k+1)]
	l=c(1,2,rep(3,16),rep(4,48),rep(5,33),rep(6,64))[seq(163)==(l+1)]
	jsfs4D[j,h,k,l] <- jsfs4D[j,h,k,l]+1
      }
    }
  }
  fst_table=fst_table[-1,]
  table2=fst_table
  table2[is.nan(table2)]<-0

  
  statistic <- c()
  for (j in 1:6) {
    for (h in 1:6) {
      statistic =c(statistic,c(t(jsfs4D[j,h,,])))
    }
  }
  statistic=c(statistic,colMeans(fst_table,na.rm=T),colMeans(table2),similar_freq)
  statistic <- statistic[2:length(statistic)] # remove first one because no polymorphism (even if not 0, if singleton in one of the 2 individuals)

  return(statistic)
}


calculate_stat_4pop_fast_pairwise_jsfs<-function(filteredloci)
{
  jsfs4D <- array(0,dim=c(221, 123, 121, 163))
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
	jsfs4D[j+1,h+1,k+1,l+1] <- jsfs4D[j+1,h+1,k+1,l+1]+1
      }
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
  return(stat45)
}

calculate_stat_2pop_detailv2<-function(data,dat1,dat2)
{
  n1=length(dat1)
  n2=length(dat2)
  jsfs2D <- array(0,dim=c((n1+1),(n2+2)))
  for (i in 1:dim(data)[1]) {
    j <- sum(data[i,dat1])
    h <- sum(data[i,dat2])
    jsfs2D[j+1,h+1] <- jsfs2D[j+1,h+1]+1
  }

  i <- 0 
  
  statistic <- numeric()
  for(j in list(1,2:round(quantile(seq(n1),p=.05)),(round(quantile(seq(n1),p=.05))+1):round(quantile(seq(n1),p=.95)),(round(quantile(seq(n1),p=.95))+1):n1)) {
    for(h in list(1,2:round(quantile(seq(n2),p=.05)),(round(quantile(seq(n2),p=.05))+1):round(quantile(seq(n2),p=.95)),(round(quantile(seq(n2),p=.95))+1):n2)) {
	i <- i+1
	statistic[i] <- sum(jsfs2D[j,h])
    }
  }

  return(statistic)
}

calculate_stat_2pop_detail<-function(data,n1,n2)
{
  jsfs2D <- array(0,dim=c(101,101))
  for (i in 1:dim(data)[1]) {
    j <- round(100*mean(data[i,n1]))
    h <- round(100*mean(data[i,n2]))
    jsfs2D[j+1,h+1] <- jsfs2D[j+1,h+1]+1
  }

  i <- 0
  statistic <- numeric()
  for(j in list(1,2:5,7:94,95:99,100)) {
    for(h in list(1,2:5,7:94,95:99,100)) {
	i <- i+1
	statistic[i] <- sum(jsfs2D[j,h])
    }
  }
  statistic <- statistic[2:24]

  return(statistic)
}

calculate_stat_4pop<-function(filteredloci)
{
  local_data=c()
  for(M in filteredloci) {
    if(length(M)>0) {
      for(i in 1:(dim(M)[2])) {
	if (mean(M[3:626,i])>.5){
	  M[3:626,i]=1-M[3:626,i]
	}
      local_data=rbind(local_data,M[3:626,i])
      }
    }
  }
  pop=list((1:220),(221:342),(343:462),(463:624))
  stat=c()
  for (i in 1:3){
    for (j in (i+1):4){
      stat=c(stat,calculate_stat_2pop_detailv2(local_data,pop[[i]],pop[[j]]))
      }
    }
   return(stat)
}

calculate_stat_3pop<-function(filteredloci)
{
  local_data=c()
  for(M in filteredloci) {
    if(length(M)>0) {
      for(i in 1:(dim(M)[2])) {
	if (mean(M[3:626,i])>.5){
	  M[3:626,i]=1-M[3:626,i]
	}
      local_data=rbind(local_data,M[3:626,i])
      }
    }
  }
  pop=list((1:220),(221:462),(463:624))
  stat=c()
  for (i in 1:2){
    for (j in (i+1):3){
      stat=c(stat,calculate_stat_2pop_detail(local_data,pop[[i]],pop[[j]]))
      }
    }
   return(stat)
}

########################################################

### HERE COMES THE SIMULATOR FOR FIXED DIFFERENCES######

########################################################

simulator_fixed_differences <- function(theta,m12,m23,m34,m45,t1,t2,t3,t4) {
  filteredloci <- list()
  otherloci <- list()
  counter.filtered <- 300
  print(paste("404 ", counter.filtered," -t ", theta," -I 5 2 100 200 100 2 0 -m 1 2 ", m12,
              " -m 2 1 ",m12," -m 3 2 ",m23," -m 2 3 ",m23,
              " -m 3 4 ",m34,
              " -m 4 3 ",m34,
              " -m 4 5 ",m45,
              " -m 5 4 ",m45,
              " -es ", t1, " 3 0.5 -ej ", t1+0.00001," 3 2 -ej ", t1+0.00001,
              " 6 4 -ej ", t2, " 1 2 -ej ", t3, " 5 4 -ej ", t4 ,"4 2",collapse=""))
  while(counter.filtered){
    scrm.out <- scrm(paste("404 ", counter.filtered," -t ", theta," -I 5 2 100 200 100 2 0 -m 1 2 ", m12,
                         " -m 2 1 ",m12," -m 3 2 ",m23," -m 2 3 ",m23,
                         " -m 3 4 ",m34,
                         " -m 4 3 ",m34,
                         " -m 4 5 ",m45,
                         " -m 5 4 ",m45,
                         " -es ", t1, " 3 0.5 -ej ", t1+0.00001," 3 2 -ej ", t1+0.00001,
                         " 6 4 -ej ", t2, " 1 2 -ej ", t3, " 5 4 -ej ", t4 ,"4 2",collapse=""))
    for (i in seq(counter.filtered)){
      if(length(scrm.out$seg_sites[[i]])) {
        M <- as.matrix(scrm.out$seg_sites[[i]][c(1,2,403,404),])
        criterion <- any(M[1,]==M[2,] & M[1,]!=M[3,] &  M[3,]==M[4,])
        if(criterion) {
          filteredloci[[counter.filtered]] <- as.matrix(scrm.out$seg_sites[[i]])
          counter.filtered <- counter.filtered - 1
        }
      }
    }
  }
  
  scrm.out <- scrm(paste("4 1000 -t ", theta," -I 5 2 0 0 0 2 0 -m 1 2 ", m12,
                           " -m 2 1 ",m12," -m 3 2 ",m23," -m 2 3 ",m23,
                           " -m 3 4 ",m34,
                           " -m 4 3 ",m34,
                           " -m 4 5 ",m45,
                           " -m 5 4 ",m45,
                           " -es ", t1, " 3 0.5 -ej ", t1+0.00001," 3 2 -ej ", t1+0.00001,
                           " 6 4 -ej ", t2, " 1 2 -ej ", t3, " 5 4 -ej ", t4 ,"4 2",collapse=""))
  for (i in seq(1000)){
    if(length(scrm.out$seg_sites[[i]])) {
      otherloci[[i]] <- as.matrix(scrm.out$seg_sites[[i]])
    }
  }
  statistic<-calculate_stat(filteredloci,otherloci)
  statistic
}

simulator_fixed_differencesv2 <- function(theta,m12,m23,m32,m34,m43,m45,t1,t2,t3,t4) {
  filteredloci <- list()
  otherloci <- list()
  counter.filtered <- 300
  while(counter.filtered){
    scrm.out <- scrm(paste("404 ", counter.filtered," -t ", theta," -I 5 2 100 200 100 2 0 -m 1 2 ", m12,
                           " -m 2 1 ",m12," -m 3 2 ",m32," -m 2 3 ",m23,
                           " -m 3 4 ",m34,
                           " -m 4 3 ",m43,
                           " -m 4 5 ",m45,
                           " -m 5 4 ",m45,
                           " -es ", t1, " 3 0.5 -ej ", t1+0.00001," 3 2 -ej ", t1+0.00001,
                           " 6 4 -ej ", t2, " 1 2 -ej ", t3, " 5 4 -ej ", t4 ,"4 2",collapse=""))
    for (i in seq(counter.filtered)){
      if(length(scrm.out$seg_sites[[i]])) {
        M <- as.matrix(scrm.out$seg_sites[[i]][c(1,2,403,404),])
        criterion <- any(M[1,]==M[2,] & M[1,]!=M[3,] &  M[3,]==M[4,])
        if(criterion) {
          filteredloci[[counter.filtered]] <- as.matrix(scrm.out$seg_sites[[i]])
          counter.filtered <- counter.filtered - 1
        }
      }
    }
  }
  
  scrm.out <- scrm(paste("4 1000 -t ", theta," -I 5 2 0 0 0 2 0 -m 1 2 ", m12,
                         " -m 2 1 ",m12," -m 3 2 ",m32," -m 2 3 ",m23,
                         " -m 3 4 ",m34,
                         " -m 4 3 ",m43,
                         " -m 4 5 ",m45,
                         " -m 5 4 ",m45,
                         " -es ", t1, " 3 0.5 -ej ", t1+0.00001," 3 2 -ej ", t1+0.00001,
                         " 6 4 -ej ", t2, " 1 2 -ej ", t3, " 5 4 -ej ", t4 ,"4 2",collapse=""))
  for (i in seq(1000)){
    if(length(scrm.out$seg_sites[[i]])) {
      otherloci[[i]] <- as.matrix(scrm.out$seg_sites[[i]])
    }
  }
  statistic<-calculate_stat(filteredloci,otherloci)
  statistic
}

simulator_fixed_differences_2pop <- function(theta,m,t) {
  filteredloci <- list()
  counter.filtered <- 215
#  print(paste("148 ", counter.filtered," -t ", theta," -I 2 95 53 0 -m 1 2 ", m,
#              " -m 2 1 ",m,
#              " -es ", t, " 1 0.5 ",collapse=""))
  while(counter.filtered){
#     scrm.out <- scrm(paste("296 ", counter.filtered," -t ", theta," -I 2 190 106 0 -m 1 2 ", m,
#               " -m 2 1 ",m,
#               " -ej ", t, " 1 2 ",collapse=""))
              
    scrm.out <- scrm(paste("300 ", counter.filtered," -t ", theta," -I 2 192 108 0 -m 1 2 ", m,
              " -m 2 1 ",m,
              " -ej ", t, " 1 2 ",collapse=""))
    for (i in seq(counter.filtered)){
      if(length(scrm.out$seg_sites[[i]])) {
        M <- as.matrix(scrm.out$seg_sites[[i]][c(1,2,299,300),])
        criterion <- any(M[1,]==M[2,] & M[1,]!=M[3,] &  M[3,]==M[4,])
        if(criterion) {
          filteredloci[[counter.filtered]] <- as.matrix(scrm.out$seg_sites[[i]])
          counter.filtered <- counter.filtered - 1
        }
      }
    }
  }
  
  statistic<-calculate_stat_2pop(filteredloci)
  print('done')
  return(statistic)
#   return(filteredloci)
}


simulator_fixed_differences_2pop_asym <- function(theta,m1,m2,t) {
  filteredloci <- list()
  counter.filtered <- 215
#  print(paste("148 ", counter.filtered," -t ", theta," -I 2 95 53 0 -m 1 2 ", m,
#              " -m 2 1 ",m,
#              " -es ", t, " 1 0.5 ",collapse=""))
  while(counter.filtered){
#     scrm.out <- scrm(paste("296 ", counter.filtered," -t ", theta," -I 2 190 106 0 -m 1 2 ", m,
#               " -m 2 1 ",m,
#               " -ej ", t, " 1 2 ",collapse=""))
              
    scrm.out <- scrm(paste("300 ", counter.filtered," -t ", theta," -I 2 192 108 0 -m 1 2 ", m1,
              " -m 2 1 ",m2,
              " -ej ", t, " 1 2 ",collapse=""))
    for (i in seq(counter.filtered)){
      if(length(scrm.out$seg_sites[[i]])) {
        M <- as.matrix(scrm.out$seg_sites[[i]][c(1,2,299,300),])
        criterion <- any(M[1,]==M[2,] & M[1,]!=M[3,] &  M[3,]==M[4,])
        if(criterion) {
          filteredloci[[counter.filtered]] <- as.matrix(scrm.out$seg_sites[[i]])
          counter.filtered <- counter.filtered - 1
        }
      }
    }
  }
  
  statistic<-calculate_stat_2pop(filteredloci)
  print('done')
  return(statistic)
#   return(filteredloci)
}

simulator_fixed_differences_2pop_2contact <- function(theta,m,t1,t2) {
  filteredloci <- list()
  counter.filtered <- 215
#  print(paste("148 ", counter.filtered," -t ", theta," -I 2 95 53 0 -m 1 2 ", m,
#              " -m 2 1 ",m,
#              " -es ", t, " 1 0.5 ",collapse=""))
  while(counter.filtered){
#     scrm.out <- scrm(paste("296 ", counter.filtered," -t ", theta," -I 2 190 106 0 -m 1 2 ", m,
#               " -m 2 1 ",m,
#               " -ej ", t, " 1 2 ",collapse=""))
              
    scrm.out <- scrm(paste("300 ", counter.filtered," -t ", theta," -I 2 192 108 0 -m 1 2 ", m,
              " -m 2 1 ",m, "-em ",t1, " 1 2 0 -em ", t1, " 2 1 0 -ej ", t2, " 1 2 ",collapse=""))
    for (i in seq(counter.filtered)){
      if(length(scrm.out$seg_sites[[i]])) {
        M <- as.matrix(scrm.out$seg_sites[[i]][c(1,2,299,300),])
        criterion <- any(M[1,]==M[2,] & M[1,]!=M[3,] &  M[3,]==M[4,])
        if(criterion) {
          filteredloci[[counter.filtered]] <- as.matrix(scrm.out$seg_sites[[i]])
          counter.filtered <- counter.filtered - 1
        }
      }
    }
  }
  
  statistic<-calculate_stat_2pop(filteredloci)
  print('done')
  return(statistic)
#   return(filteredloci)
}

simulator_fixed_differences_2pop_test <- function(theta,m,t) {
  filteredloci <- list()
  counter.filtered <- 215
#  print(paste("148 ", counter.filtered," -t ", theta," -I 2 95 53 0 -m 1 2 ", m,
#              " -m 2 1 ",m,
#              " -es ", t, " 1 0.5 ",collapse=""))
  while(counter.filtered){
#     scrm.out <- scrm(paste("296 ", counter.filtered," -t ", theta," -I 2 190 106 0 -m 1 2 ", m,
#               " -m 2 1 ",m,
#               " -ej ", t, " 1 2 ",collapse=""))
              
    scrm.out <- scrm(paste("300 ", counter.filtered," -t ", theta," -I 2 192 108 0 -m 1 2 ", m,
              " -m 2 1 ",m,
              " -ej ", t, " 1 2 ",collapse=""))
    for (i in seq(counter.filtered)){
      if(length(scrm.out$seg_sites[[i]])) {
        M <- as.matrix(scrm.out$seg_sites[[i]][c(1,2,299,300),])
        criterion <- any(M[1,]==M[2,] & M[1,]!=M[3,] &  M[3,]==M[4,])
        if(criterion) {
          filteredloci[[counter.filtered]] <- as.matrix(scrm.out$seg_sites[[i]])
          counter.filtered <- counter.filtered - 1
        }
      }
    }
  }
  
#   statistic<-calculate_stat_2pop(filteredloci)
#   print('done')
#   return(statistic)
  return(filteredloci)
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
  
   statistic<-calculate_stat_4pop_4Djsfs_v6(filteredloci)
   return(statistic)
#    return(filteredloci)
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
  
  statistic<-calculate_stat_4pop(filteredloci)
  return(statistic)
#   return(filteredloci)
}

simulator_fixed_differences_3pop <- function(theta,m12,m23,mt2,t1,t2) {
  filteredloci <- list()
  counter.filtered <- 215
  while(counter.filtered){
     
    scrm.out <- scrm(paste("628 ", counter.filtered," -t ", theta," -I 3 222 242 164 0 -m 1 2 ", m12,
              " -m 2 1 ",m12, "-m 3 2 ",m23," -m 2 3 ",m23,
                           " -es ", t1, " 2 0.5 -ej ", t1+0.00001," 1 2 -ej ", t1+0.00001,
              " 4 3", " -em ", t1+0.00001," 2 3 ", mt2, " -em ", t1+0.00001," 3 2 ", mt2, " -ej ",
              t2, " 2 3 ",collapse=""))
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
  
  statistic<-calculate_stat_3pop(filteredloci)
  return(statistic)
#   return(filteredloci)
}


simulator_fixed_differences_1pop <- function(theta) {
  filteredloci <- list()
  counter.filtered <- 215
#  print(paste("148 ", counter.filtered," -t ", theta," -I 2 95 53 0 -m 1 2 ", m,
#              " -m 2 1 ",m,
#              " -es ", t, " 1 0.5 ",collapse=""))
  while(counter.filtered){
#     scrm.out <- scrm(paste("296 ", counter.filtered," -t ", theta," -I 2 190 106 0 -m 1 2 ", m,
#               " -m 2 1 ",m,
#               " -ej ", t, " 1 2 ",collapse=""))
              
    scrm.out <- scrm(paste("300 ", counter.filtered," -t ", theta," -I 1 300 0"))
    for (i in seq(counter.filtered)){
      if(length(scrm.out$seg_sites[[i]])) {
        M <- as.matrix(scrm.out$seg_sites[[i]][c(1,2,299,300),])
        criterion <- any(M[1,]==M[2,] & M[1,]!=M[3,] &  M[3,]==M[4,])
        if(criterion) {
          filteredloci[[counter.filtered]] <- as.matrix(scrm.out$seg_sites[[i]])
          counter.filtered <- counter.filtered - 1
        }
      }
    }
  }
  
  statistic<-calculate_stat_2pop(filteredloci)
  return(statistic)
#   return(filteredloci)
}





########################################################

### HERE COMES THE SIMULATOR FOR SHARED POLYMORPHISM####

########################################################

simulator_both_polymorphic <- function(theta,m12,m23,m34,m45,t1,t2,t3,t4) {
  filteredloci <- list()
  otherloci <- list()
  counter.filtered <- 300 
  while(counter.filtered){
    scrm.out <- scrm(paste("404 ", counter.filtered," -t ", theta," -I 5 2 100 200 100 2 0 -m 1 2 ", m12,
                           " -m 2 1 ",m12," -m 3 2 ",m23," -m 2 3 ",m23,
                           " -m 3 4 ",m34,
                           " -m 4 3 ",m34,
                           " -m 4 5 ",m45,
                           " -m 5 4 ",m45,
                           " -es ", t1, " 3 0.5 -ej ", t1+0.00001," 3 2 -ej ", t1+0.00001,
                           " 6 4 -ej ", t2, " 1 2 -ej ", t3, " 5 4 -ej ", t4 ,"4 2",collapse=""))
    for (i in seq(counter.filtered)){
      if(length(scrm.out$seg_sites[[i]])) {
        M <- as.matrix(scrm.out$seg_sites[[i]][c(1,2,403,404),])
        criterion <- any(M[1,]!=M[2,] & M[3,]!=M[4,])
        if(criterion) {
          filteredloci[[counter.filtered]] <- as.matrix(scrm.out$seg_sites[[i]])
          counter.filtered <- counter.filtered - 1
        }
      }
    }
  }
  
  scrm.out <- scrm(paste("4 1000 -t ", theta," -I 5 2 0 0 0 2 0 -m 1 2 ", m12,
                         " -m 2 1 ",m12," -m 3 2 ",m23," -m 2 3 ",m23,
                         " -m 3 4 ",m34,
                         " -m 4 3 ",m34,
                         " -m 4 5 ",m45,
                         " -m 5 4 ",m45,
                         " -es ", t1, " 3 0.5 -ej ", t1+0.00001," 3 2 -ej ", t1+0.00001,
                         " 6 4 -ej ", t2, " 1 2 -ej ", t3, " 5 4 -ej ", t4 ,"4 2",collapse=""))
  for (i in seq(1000)){
    if(length(scrm.out$seg_sites[[i]])) {
      otherloci[[i]] <- as.matrix(scrm.out$seg_sites[[i]])
    }
  }
  statistic<-calculate_stat(filteredloci,otherloci)
  statistic
}

########################################################

### HERE COMES THE SIMULATOR FOR PRIVATE POLYMORPHISM###

########################################################

simulator_private_polymorphic <- function(theta,m12,m23,m34,m45,t1,t2,t3,t4) {
  filteredloci <- list()
  otherloci <- list()
  counter.filtered <- 300
  
  while(counter.filtered){
    scrm.out <- scrm(paste("404 ", counter.filtered," -t ", theta," -I 5 2 100 200 100 2 0 -m 1 2 ", m12,
                           " -m 2 1 ",m12," -m 3 2 ",m23," -m 2 3 ",m23,
                           " -m 3 4 ",m34,
                           " -m 4 3 ",m34,
                           " -m 4 5 ",m45,
                           " -m 5 4 ",m45,
                           " -es ", t1, " 3 0.5 -ej ", t1+0.00001," 3 2 -ej ", t1+0.00001,
                           " 6 4 -ej ", t2, " 1 2 -ej ", t3, " 5 4 -ej ", t4 ,"4 2",collapse=""))
    for (i in seq(counter.filtered)){
      if(length(scrm.out$seg_sites[[i]])) {
        M <- as.matrix(scrm.out$seg_sites[[i]][c(1,2,403,404),])
        criterion <- any((M[1,]!=M[2,] & M[3,]==M[4,]) | (M[1,]==M[2,] & M[3,]!=M[4,]))
        if(criterion) {
          filteredloci[[counter.filtered]] <- as.matrix(scrm.out$seg_sites[[i]])
          counter.filtered <- counter.filtered - 1
        }
      }
    }
  }
  
  scrm.out <- scrm(paste("4 1000 -t ", theta," -I 5 2 0 0 0 2 0 -m 1 2 ", m12,
                         " -m 2 1 ",m12," -m 3 2 ",m23," -m 2 3 ",m23,
                         " -m 3 4 ",m34,
                         " -m 4 3 ",m34,
                         " -m 4 5 ",m45,
                         " -m 5 4 ",m45,
                         " -es ", t1, " 3 0.5 -ej ", t1+0.00001," 3 2 -ej ", t1+0.00001,
                         " 6 4 -ej ", t2, " 1 2 -ej ", t3, " 5 4 -ej ", t4 ,"4 2",collapse=""))
  for (i in seq(1000)){
    if(length(scrm.out$seg_sites[[i]])) {
      otherloci[[i]] <- as.matrix(scrm.out$seg_sites[[i]])
    }
  }
  statistic<-calculate_stat(filteredloci,otherloci)
  statistic
}

########################################################

simulator_private_polymorphic_v2 <- function(theta,m12,m23,m32,m34,m43,m45,t1,t2,t3,t4) {
  filteredloci <- list()
  otherloci <- list()
  counter.filtered <- 300
  
  while(counter.filtered){
    scrm.out <- scrm(paste("404 ", counter.filtered," -t ", theta," -I 5 2 100 200 100 2 0 -m 1 2 ", m12,
                           " -m 2 1 ",m12," -m 3 2 ",m32," -m 2 3 ",m23,
                           " -m 3 4 ",m34,
                           " -m 4 3 ",m43,
                           " -m 4 5 ",m45,
                           " -m 5 4 ",m45,
                           " -es ", t1, " 3 0.5 -ej ", t1+0.00001," 3 2 -ej ", t1+0.00001,
                           " 6 4 -ej ", t2, " 1 2 -ej ", t3, " 5 4 -ej ", t4 ,"4 2",collapse=""))
    for (i in seq(counter.filtered)){
      if(length(scrm.out$seg_sites[[i]])) {
        M <- as.matrix(scrm.out$seg_sites[[i]][c(1,2,403,404),])
        criterion <- any((M[1,]!=M[2,] & M[3,]==M[4,]) | (M[1,]==M[2,] & M[3,]!=M[4,]))
        if(criterion) {
          filteredloci[[counter.filtered]] <- as.matrix(scrm.out$seg_sites[[i]])
          counter.filtered <- counter.filtered - 1
        }
      }
    }
  }
  
  scrm.out <- scrm(paste("4 1000 -t ", theta," -I 5 2 0 0 0 2 0 -m 1 2 ", m12,
                         " -m 2 1 ",m12," -m 3 2 ",m32," -m 2 3 ",m23,
                         " -m 3 4 ",m34,
                         " -m 4 3 ",m43,
                         " -m 4 5 ",m45,
                         " -m 5 4 ",m45,
                         " -es ", t1, " 3 0.5 -ej ", t1+0.00001," 3 2 -ej ", t1+0.00001,
                         " 6 4 -ej ", t2, " 1 2 -ej ", t3, " 5 4 -ej ", t4 ,"4 2",collapse=""))
  for (i in seq(1000)){
    if(length(scrm.out$seg_sites[[i]])) {
      otherloci[[i]] <- as.matrix(scrm.out$seg_sites[[i]])
    }
  }
  statistic<-calculate_stat(filteredloci,otherloci)
  statistic
}
########################################################

### HERE COMES THE SIMULATOR ###########################

########################################################

simulator <- function(theta,m12,m23,m34,m45,t1,t2,t3,t4) {
  filteredloci <- list()
  otherloci <- list()
  counter.filtered <- 300
  while(counter.filtered){
    scrm.out <- scrm(paste("404 ", counter.filtered," -t ", theta," -I 5 2 100 200 100 2 0 -m 1 2 ", m12,
                           " -m 2 1 ",m12," -m 3 2 ",m23," -m 2 3 ",m23,
                           " -m 3 4 ",m34,
                           " -m 4 3 ",m34,
                           " -m 4 5 ",m45,
                           " -m 5 4 ",m45,
                           " -es ", t1, " 3 0.5 -ej ", t1+0.00001," 3 2 -ej ", t1+0.00001,
                           " 6 4 -ej ", t2, " 1 2 -ej ", t3, " 5 4 -ej ", t4 ,"4 2",collapse=""))
    for (i in seq(counter.filtered)){
      if(length(scrm.out$seg_sites[[i]])) {
        M <- as.matrix(scrm.out$seg_sites[[i]][c(1,2,403,404),])
        filteredloci[[counter.filtered]] <- as.matrix(scrm.out$seg_sites[[i]])
        counter.filtered <- counter.filtered - 1
      }
    }
  }
  
  scrm.out <- scrm(paste("4 1000 -t ", theta," -I 5 2 0 0 0 2 0 -m 1 2 ", m12,
                         " -m 2 1 ",m12," -m 3 2 ",m23," -m 2 3 ",m23,
                         " -m 3 4 ",m34,
                         " -m 4 3 ",m34,
                         " -m 4 5 ",m45,
                         " -m 5 4 ",m45,
                         " -es ", t1, " 3 0.5 -ej ", t1+0.00001," 3 2 -ej ", t1+0.00001,
                         " 6 4 -ej ", t2, " 1 2 -ej ", t3, " 5 4 -ej ", t4 ,"4 2",collapse=""))
  for (i in seq(1000)){
    if(length(scrm.out$seg_sites[[i]])) {
      otherloci[[i]] <- as.matrix(scrm.out$seg_sites[[i]])
    }
  }
  statistic<-calculate_stat(filteredloci,otherloci)
  statistic
}


simulator_mixed_SNP  <- function(theta,m12,m23,m34,m45,t1,t2,t3,t4) {
  filteredloci <- list()
  otherloci <- list()
  counter.filtered <- 150
  while(counter.filtered){
    scrm.out <- scrm(paste("404 ", counter.filtered," -t ", theta," -I 5 2 100 200 100 2 0 -m 1 2 ", m12,
                           " -m 2 1 ",m12," -m 3 2 ",m23," -m 2 3 ",m23,
                           " -m 3 4 ",m34,
                           " -m 4 3 ",m34,
                           " -m 4 5 ",m45,
                           " -m 5 4 ",m45,
                           " -es ", t1, " 3 0.5 -ej ", t1+0.00001," 3 2 -ej ", t1+0.00001,
                           " 6 4 -ej ", t2, " 1 2 -ej ", t3, " 5 4 -ej ", t4 ,"4 2",collapse=""))
    for (i in seq(counter.filtered)){
      if(length(scrm.out$seg_sites[[i]])) {
        M <- as.matrix(scrm.out$seg_sites[[i]][c(1,2,403,404),])
        criterion <- any(M[1,]==M[2,] & M[1,]!=M[3,] &  M[3,]==M[4,])
        if(criterion) {
          filteredloci[[counter.filtered]] <- as.matrix(scrm.out$seg_sites[[i]])
          counter.filtered <- counter.filtered - 1
        }
      }
    }
  }
  
  filteredloci2 <- list()
  counter.filtered <- 150
  while(counter.filtered){
    scrm.out <- scrm(paste("404 ", counter.filtered," -t ", theta," -I 5 2 100 200 100 2 0 -m 1 2 ", m12,
                           " -m 2 1 ",m12," -m 3 2 ",m23," -m 2 3 ",m23,
                           " -m 3 4 ",m34,
                           " -m 4 3 ",m34,
                           " -m 4 5 ",m45,
                           " -m 5 4 ",m45,
                           " -es ", t1, " 3 0.5 -ej ", t1+0.00001," 3 2 -ej ", t1+0.00001,
                           " 6 4 -ej ", t2, " 1 2 -ej ", t3, " 5 4 -ej ", t4 ,"4 2",collapse=""))
    for (i in seq(counter.filtered)){
      if(length(scrm.out$seg_sites[[i]])) {
        M <- as.matrix(scrm.out$seg_sites[[i]][c(1,2,403,404),])
        criterion <- any(M[1,]==M[2,] & M[1,]!=M[3,] &  M[3,]==M[4,])
        if(criterion) {
          filteredloci2[[counter.filtered]] <- as.matrix(scrm.out$seg_sites[[i]])
          counter.filtered <- counter.filtered - 1
        }
      }
    }
  }
  
  scrm.out <- scrm(paste("4 1000 -t ", theta," -I 5 2 0 0 0 2 0 -m 1 2 ", m12,
                         " -m 2 1 ",m12," -m 3 2 ",m23," -m 2 3 ",m23,
                         " -m 3 4 ",m34,
                         " -m 4 3 ",m34,
                         " -m 4 5 ",m45,
                         " -m 5 4 ",m45,
                         " -es ", t1, " 3 0.5 -ej ", t1+0.00001," 3 2 -ej ", t1+0.00001,
                         " 6 4 -ej ", t2, " 1 2 -ej ", t3, " 5 4 -ej ", t4 ,"4 2",collapse=""))
  for (i in seq(1000)){
    if(length(scrm.out$seg_sites[[i]])) {
      otherloci[[i]] <- as.matrix(scrm.out$seg_sites[[i]])
    }
  }
  statistic<-calculate_stat(filteredloci,otherloci)
  statistic2<-calculate_stat(filteredloci2,otherloci)
  statistic<-c(statistic[1:68],statistic2)
  return(statistic)
}
#
#
#
#
# description_stat=c('low_pol in 3', 'pol in 3', 'low_pol in 2', 'low_pol in 2 and 3','low_pol in 2 pol in 3',' pol in 2', 'pol in 2 low_pol in 3', 'pol in 2 and 3','low_pol in 1' , 'low_pol in 1 and 3',  'low_pol in 1 pol in 3', 'low_pol in 1 and 2' ,'low_pol in 1 and 2 and 3' ,'low_pol in 1 and 2 pol in 3','low_pol in 1 pol in 2' , 'pol in 2 low_pol in 1 and 3', 'low_pol in 1 pol in 2 and 3 pol in 1' , 'pol in 1 low_pol in 3',  'pol in 1 and 3', 'pol in 1 low_pol in 2', 'pol in 1 low_pol in 2 and 3' ,'low_pol in 2 pol in 1 and 3', 'pol in 1 and 2',  'pol in 1 and 2 low_pol in 3' , 'pol in 1, 2 and 3')
# description_stat=c(description_stat,'fixed in 1 and 2, pol in 3', 'fixed in 1, pol in 2 and 3',  'fixed in 1 and 3 pol in 2','pol in 1 and 2, fixed in 3','pol in 1, fixed in 2 and 3','pol in 1 and 3, fixed in 2')
# description_stat=c(description_stat,'no pol at all', 'pol in s2', 'fixed in s2', 'low_pol in all', 'low_pol in all, pol in s2','low_pol in all fixed in s2', 'pol in all','pol in all and in s2','pol in all fixed in s2','fixed in all','fixed in all pol in s2','fixed in all and in s2')
# description_stat=c(description_stat,'pol in s1','pol in s1 and s2', 'pol in s1 fixed in s2', 'pol in s1 low_pol in all', 'low_pol in all, pol in s1 and s2','pol in s1 low_pol in all fixed in s2', 'pol in s1 and all','pol in s1, all and in s2','pol in s1 and all fixed in s2','pol in s1 fixed in all','fixed in all pol in s1 and s2','pol in s1 fixed in all and in s2')
# description_stat=c(description_stat,'fixed in s1', 'fixed in s1 pol in s2', 'fixed in s1 and s2', 'fixed in s1 low_pol in all', 'fixed in s1 low_pol in all', 'fixed in s1 pol in s2','low_pol in all fixed in s1 and s2', 'fixed in s1 pol in all','fixed in s1 pol in all and in s2','pol in all fixed in s1 and s2','fixed in s1 and all','fixed in s1 and all pol in s2','fixed in s1, all and in s2')
# description_stat=c(description_stat, 'fixed diff between s1 s2', 'pol in s1, fixed or no pol in s2', 'pol in s2, fixed or no pol in s1','double het in s1 s2')
# description_stat=c(description_stat,'from 1 to 32, we consider main pop of each species + the HZ; from 1 to 26, low_pol means 1 2 or 3 copies of the derived allele, pol between 4 and n-1 copies of the derived allele; from 27 to 32 pol means from 1 to (n-1) copies of the ferived allele;  from 33 to 68 we consider the two samples individueals (s1 s2) and all the other indivius (all), low_pol means 1 2 or 3 copies of the derived alleles, pol from 4 to (n-1) for the all population; from 69 to 72, sjfs between the 2 individuals (s1 and s2) ')