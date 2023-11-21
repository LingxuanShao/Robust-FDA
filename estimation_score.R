library("foreach")
library("doParallel")

#eigen function
varphi=function(k,s){
  return(sqrt(2)*sin(k*pi*s));
}

relamse=function(cut,xi,hatxi){
  error=array(0,cut);
  for(k in 1:cut){
    error[k]=min(sum((hatxi[k,]-xi[k,])*(hatxi[k,]-xi[k,])),sum((-hatxi[k,]-xi[k,])*(-hatxi[k,]-xi[k,])))
  }
  return(sum(error)/sum(xi[1:cut,]*xi[1:cut,]));
}


###########################################################################
#mento carlo run
mento=function(it){
  write.csv(it,"process.csv");
  #generate score
  xi=array(0,c(cut,n));
  if(betascore==0){
    for(k in 1:cut){
      xi[k,]=rnorm(n,0,1)*(cut+1-k);
    } 
  }else{
  for(k in 1:cut){
    xi[k,]=(rbeta(n,betascore,betascore)-0.5)*(cut+1-k);
  }}
  #generate data
  Z=array(0,c(n,m));
  t=array(0,c(n,m));
  for(i in 1:n){
    t[i,]=sort(runif(m,0,1));
    for(j in 1:m){
      Z[i,j]=0;
      for(k in 1:cut){
        Z[i,j]=Z[i,j]+xi[k,i]*varphi(k,t[i,j]);
      }
      Z[i,j]=Z[i,j]+rnorm(1,0,sigma);
    }
  }
  ###########################################################
  error_ni=NA;
  error_pace=NA;
  library("fdapace");
  #Numerical integration method
  temp_ni=FPCA(Ly=split(Z, row(Z)),Lt=split(t, row(t)),list(methodXi="IN",kFoldMuCov=2,methodSelectK=1));
  temp_ni2=t(temp_ni$xiEst);
  error_ni=relamse(1,xi,temp_ni2);
  #PACE method
  temp_pace=FPCA(Ly=split(Z, row(Z)),Lt=split(t, row(t)),list(methodXi="CE",kFoldMuCov=2,methodSelectK=1));
  temp_pace2=t(temp_pace$xiEst);
  error_pace=relamse(1,xi,temp_pace2);
  return(list(error_ni,error_pace));
}
  
#####################################################################
###  main
#settings
for(scorenumber in c(1,2,3)){
for(n in c(100,200)){
for(m in c(5,7,9)){
set.seed(97);
cut=2;#dimension of score
sigma=0.1;
betascore=scorenumber;

  itermax=100;
  cl <- makeCluster(50);
  registerDoParallel(cl);
  result <- foreach(it=c(1:itermax), .combine='c') %dopar% mento(it);
  stopCluster(cl);
  
  result_mento=array(0,itermax);
  result_pace=array(0,itermax);
  for(it in 1:itermax){
    result_mento[it]=result[[2*it-1]];
    result_pace[it]=result[[2*it]];
  }
  
     cat(scorenumber,n,m,mean(result_mento),sd(result_mento)/sqrt(itermax),mean(result_pace),sd(result_pace)/sqrt(itermax),"\n")
output=c(scorenumber,n,m,mean(result_mento),sd(result_mento)/sqrt(itermax),mean(result_pace),sd(result_pace)/sqrt(itermax))
write.csv(output,paste0("~/score_",scorenumber,"_",n,"_",m,sep=""));
}}}





