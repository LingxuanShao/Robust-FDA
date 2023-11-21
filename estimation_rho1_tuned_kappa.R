library("foreach")
library("doParallel")
require("ggplot2")
library("expm")

###########################################################################
#mento carlo run
mento=function(it){
  write.csv(it,"process.csv");
  #prepare
  
  #eigen function
  phi=function(k,s){
    return(sqrt(2)*sin(k*pi*s));
  }
  
  #kernel
  Ker=function(s,h){
    temp=ifelse((s-h)>0,0,ifelse((s+h)<0,0,((1-(abs(s/h))**3)**3)*70/(h*81)))
    return(temp);
  }
  
  #robust function
  rho=function(x,robu){
    #un-rescaled version
    if(robu==0){
      return(x*x);
    } 
    if(robu==1){
      temp=ifelse(x>kappa,x,ifelse((x+kappa)<0,-x, (-x*x*x*x+6*kappa*kappa*x*x+3*kappa*kappa*kappa*kappa)/(8*kappa*kappa*kappa) ));
      return(temp);
    }
    if(robu==2){
      temp=ifelse(abs(x)>10,abs(x)-log(2),log(cosh(x)));
      return(temp);
    }
    if(robu==3){
      return(2*x*atan(x)/pi-log(1+x*x)/pi);
    }
    if(robu==5){ #robu=5 stands for |x|
      return(abs(x));
    }
  }
  
  psi=function(x,robu){
    if(robu==0){
      return(2*x);
    } 
    if(robu==1){
      temp=ifelse(x>kappa,1,ifelse((x+kappa)<0,-1, (-x*x*x+3*kappa*kappa*x)/(2*kappa*kappa*kappa) ));
      return(temp);
    }
    if(robu==2){
      temp=ifelse(x>10,1,ifelse(x<(-10),-1,(exp(x)-exp(-x))/(exp(x)+exp(-x))));
      return(temp);
    }
    if(robu==3){
      return(2*atan(x)/pi);
    }
    if(robu==5){ #robu=5 stands for |x|
      return(sign(x));
    }
  }
  
  #obtain hatmurobust
  hatQmurobust=function(s,n,m,t,X,hmu,beta,robu){
    hatu1=sum(Ker(t-s,hmu)*(t-s));
    hatu2=sum(Ker(t-s,hmu)*(t-s)*(t-s));
    return(  sum(Ker(t-s,hmu)*(hatu2-hatu1*(t-s))*rho(X-beta,robu)) );
  }
  searchminimumrobust=function(s,n,m,t,X,hmu,candidatebeta,robu){
    len=length(candidatebeta);
    ans=array(0,len);
    for(cand in 1:len){
      ans[cand]=hatQmurobust(s,n,m,t,X,hmu,candidatebeta[cand],robu);
    }
    return(which.min(ans));
  }
  obtainhatmurobust=function(s,n,m,t,X,hmu,robu){
    maxX=min(max(X),10);#be careful with this
    minX=max(min(X),-10);
    searchlen=5;
    candidatebeta=c(0:(searchlen-1))/(searchlen-1)*(maxX-minX)+minX;
    while((maxX-minX)>0.001){#be careful with this
      l0=searchminimumrobust(s,n,m,t,X,hmu,candidatebeta,robu);
      if(l0==1){l0=l0+1;}
      if(l0==searchlen){l0=l0-1;}
      maxX=candidatebeta[l0+1];
      minX=candidatebeta[l0-1];
      candidatebeta=c(0:(searchlen-1))/(searchlen-1)*(maxX-minX)+minX;
    }
    return((maxX+minX)/2);
  }
  
  #first half of X for training; second half of X for testing
  test=function(n,m,t,X,hmu,robu){
    testerror=0;
    for(i in ((n/2)+1):n){
      for(j in 1:m){
        testerror=testerror+rho(obtainhatmurobust(t[i,j],(n/2),m,t[1:(n/2),],X[1:(n/2),],hmu,robu)-X[i,j],robu);
      }
    }
    return(testerror);
  }
  
  #obtain hatCrobust
  EKer=function(s,h,a){
    if(a==0){
      if(s<h){
        return( 0.5+(70/81)*((s/h)-(3/4)*((s/h)**4)+(3/7)*((s/h)**7)-(1/10)*((s/h))**(10)) );
      }
      if(s>(1-h)){return(EKer(1-s,h,a));}
      return(1);
    }
    if(a==1){
      if(s<h){
        return( h*(70/81)*((1/2)-(3/5)+(3/8)-(1/11)) - h*(70/81)*((1/2)*((s/h)**2)-(3/5)*((s/h)**5)+(3/8)*((s/h)**8)-(1/11)*((s/h))**(11)) );
      }
      if(s>(1-h)){return(-EKer(1-s,h,a));}
      return(0);
    }
    if(a==2){
      if(s<h){
        return( (35/486)*h*h + h*h*(70/81)*((1/3)*((s/h)**3)-(1/2)*((s/h)**6)+(1/3)*((s/h)**9)-(1/12)*((s/h))**(12)) );
      }
      if(s>(1-h)){return(EKer(1-s,h,a));}
      return((35/243)*h*h);
    }
  }
  
  
  obtainhatCrobust=function(tau1,tau2,n,m,t,psiX,hc){
    #S00=sum(rowSums(Ker(t-tau1,hc))*rowSums(Ker(t-tau2,hc))-rowSums(Ker(t-tau1,hc)*Ker(t-tau2,hc)));
    #S10=sum(rowSums(Ker(t-tau1,hc)*(t-tau1))*rowSums(Ker(t-tau2,hc))-rowSums(Ker(t-tau1,hc)*(t-tau1)*Ker(t-tau2,hc)));
    #S01=sum(rowSums(Ker(t-tau1,hc))*rowSums(Ker(t-tau2,hc)*(t-tau2))-rowSums(Ker(t-tau1,hc)*Ker(t-tau2,hc)*(t-tau2)));
    #S20=sum(rowSums(Ker(t-tau1,hc)*(t-tau1)*(t-tau1))*rowSums(Ker(t-tau2,hc))-rowSums(Ker(t-tau1,hc)*(t-tau1)*(t-tau1)*Ker(t-tau2,hc)));
    #S11=sum(rowSums(Ker(t-tau1,hc)*(t-tau1))*rowSums(Ker(t-tau2,hc)*(t-tau2))-rowSums(Ker(t-tau1,hc)*(t-tau1)*Ker(t-tau2,hc)*(t-tau2)));
    #S02=sum(rowSums(Ker(t-tau1,hc))*rowSums(Ker(t-tau2,hc)*(t-tau2)*(t-tau2))-rowSums(Ker(t-tau1,hc)*Ker(t-tau2,hc)*(t-tau2)*(t-tau2)));
    S00=EKer(tau1,hc,0)*EKer(tau2,hc,0);
    S10=EKer(tau1,hc,1)*EKer(tau2,hc,0);
    S01=EKer(tau1,hc,0)*EKer(tau2,hc,1);
    S20=EKer(tau1,hc,2)*EKer(tau2,hc,0);
    S11=EKer(tau1,hc,1)*EKer(tau2,hc,1);
    S02=EKer(tau1,hc,0)*EKer(tau2,hc,2);
    R00=sum(rowSums(Ker(t-tau1,hc)*psiX)*rowSums(Ker(t-tau2,hc)*psiX)-rowSums(Ker(t-tau1,hc)*Ker(t-tau2,hc)*psiX*psiX));
    R10=sum(rowSums(Ker(t-tau1,hc)*psiX*(t-tau1))*rowSums(Ker(t-tau2,hc)*psiX)-rowSums(Ker(t-tau1,hc)*psiX*(t-tau1)*Ker(t-tau2,hc)*psiX));
    R01=sum(rowSums(Ker(t-tau1,hc)*psiX)*rowSums(Ker(t-tau2,hc)*psiX*(t-tau2))-rowSums(Ker(t-tau1,hc)*psiX*Ker(t-tau2,hc)*psiX*(t-tau2)));
    ini_order=(R00*(S20*S02-S11*S11)-R10*(S10*S02-S01*S11)+R01*(S10*S11-S01*S20))/(S00*(S20*S02-S11*S11)-S10*(S10*S02-S01*S11)+S01*(S10*S11-S01*S20))
    return(ini_order/(n*m*(m-1)));
  }
  
  #train C: first half of X; test: second half of X 
  testC=function(n,m,t,psiX,hc){
    testerror=0;
    for(i in (((n/2)+1):n)){
      for(j1 in 1:m){
        for(j2 in 1:m){
          testerror=testerror+(obtainhatCrobust(t[i,j1],t[i,j2],n/2,m,t[1:(n/2),],psiX[1:(n/2),],hc)-psiX[i,j1]*psiX[i,j2])**2;
        }
      }
    }
    return(testerror);
  }
  
  ##################################################
  #generate score
  xi=array(0,c(cut,n));
  #normal score
  if(normalscore==1){
    for(k in 1:cut){
      xi[k,]=rnorm(n,0,k);
    }
  }
  #cauchy score
  if(cauchyscore==1){
    for(k in 1:cut){
      xi[k,]=rcauchy(n,0,k);
    }
  }
  #t score
  if(tscore==1){
    for(k in 1:cut){
      xi[k,]=rt(n,k);
    }
  }
  #sln score
  if(slnscore==1){
    for(k in 1:cut){
      xi[k,]=rlnorm(n,0,k)*(floor(runif(n,0,2))-0.5)*2;
    }
  }
  #asym score
  if(asymscore==1){
    for(k in 1:cut){
      xi[k,]=rbeta(n,k+k,k)-2/3;
    }
  }
  #generate data
  X=array(0,c(n,m));
  t=array(0,c(n,m));
  for(i in 1:n){
    for(j in 1:m){
      t[i,j]=runif(1,0,1);
      for(k in 1:cut){
        X[i,j]=X[i,j]+xi[k,i]*phi(k,t[i,j]);
      }
    }
  }
  #add 5% modifications
  if(addmodi05==1){
    for(i in 1:n){
      for(j in 1:m){
        modi=runif(1,0,1);
        if(modi<0.05){
          X[i,j]=rnorm(1,10,0.1);
        }
      }
    }
  }
  #add 10% modifications
  if(addmodi10==1){
    for(i in 1:n){
      for(j in 1:m){
        modi=runif(1,0,1);
        if(modi<0.1){
          X[i,j]=rnorm(1,10,0.1);
        }
      }
    }
  }
  #add 15% modifications
  if(addmodi15==1){
    for(i in 1:n){
      for(j in 1:m){
        modi=runif(1,0,1);
        if(modi<0.15){
          X[i,j]=rnorm(1,10,0.1);
        }
      }
    }
  }
  #add 20% modifications
  if(addmodi20==1){
    for(i in 1:n){
      for(j in 1:m){
        modi=runif(1,0,1);
        if(modi<0.2){
          X[i,j]=rnorm(1,10,0.1);
        }
      }
    }
  }
  
  #select kappa
  lenkappa=3;
  candidatekappa=c(0.001,0.01,0.1);
  testerrorkappa=array(0,lenkappa);
  for(candikappa in 1:lenkappa){
    kappa=candidatekappa[candikappa];
    trainn=n/2;
    testn=n/2;
    trainm=m;
    testm=m;
    traint=t[1:trainn,];
    testt=t[(trainn+1):n,];
    trainX=X[1:trainn,];
    testX=X[(trainn+1):n,];
    #select hmu;
    candilen=5;
    candidatehmu=exp(seq(log(0.2),log(0.5),length.out=candilen));
    testerrorlist=array(0,candilen);
    for(candi in 1:candilen){
      testerrorlist[candi]=test(trainn,trainm,traint,trainX,candidatehmu[candi],robu);
    }
    hmu=candidatehmu[which.min(testerrorlist)];
    itresult_mu_square=0;
    for(i in 1:testn){
    for(j in 1:testm){
      temp=obtainhatmurobust(testt[i,j],trainn,trainm,traint,trainX,hmu,robu);
      testerrorkappa[candikappa]=testerrorkappa[candikappa]+(temp-testX[i,j])**2/l;
    }}
    #switch train and test
    traint=t[(trainn+1):n,];
    testt=t[1:trainn,];
    trainX=X[(trainn+1):n,];
    testX=X[1:trainn,];
    #select hmu;
    candilen=5;
    candidatehmu=exp(seq(log(0.2),log(0.5),length.out=candilen));
    testerrorlist=array(0,candilen);
    for(candi in 1:candilen){
      testerrorlist[candi]=test(trainn,trainm,traint,trainX,candidatehmu[candi],robu);
    }
    hmu=candidatehmu[which.min(testerrorlist)];
    itresult_mu_square=0;
    for(i in 1:testn){
    for(j in 1:testm){
        temp=obtainhatmurobust(testt[i,j],trainn,trainm,traint,trainX,hmu,robu);
        testerrorkappa[candikappa]=testerrorkappa[candikappa]+(temp-testX[i,j])**2/l;
    }}
  }
  kappa=candidatekappa[which.min(testerrorkappa)];
  #now kappa is selected
  
  #true value of mu
  truemu=array(0,l);
  if((asymscore==1)&&(robu!=0)){
    Fmurobust=function(menton,Xtanchork,robu,beta){
      return(sum(rho(Xtanchork-beta,robu)));
    }
    searchminimumFrobust=function(menton,Xtanchork,robu,candidatebeta){
      len=length(candidatebeta);
      ans=array(0,len);
      for(cand in 1:len){
        ans[cand]=Fmurobust(menton,Xtanchork,robu,candidatebeta[cand]);
      }
      return(which.min(ans));
    }
    obtainFmurobust=function(menton,Xtanchork,robu){
      maxX=min(max(Xtanchork),10);
      minX=max(min(Xtanchork),-10);
      searchlen=5;
      candidatebeta=c(0:(searchlen-1))/(searchlen-1)*(maxX-minX)+minX;
      while((maxX-minX)>0.01){
        l0=searchminimumFrobust(menton,Xtanchork,robu,candidatebeta);
        if(l0==1){l0=l0+1;}
        if(l0==searchlen){l0=l0-1;}
        maxX=candidatebeta[l0+1];
        minX=candidatebeta[l0-1];
        candidatebeta=c(0:(searchlen-1))/(searchlen-1)*(maxX-minX)+minX;
      }
      return((maxX+minX)/2);
    }
    mentoitermax=100;
    mentoresult=array(0,c(mentoitermax,l));
    for(it in 1:mentoitermax){
      #cat("mentoit=",it,"\n");
      menton=10**4;
      mentoxi=array(0,c(cut,menton));
      Xtanchork=array(0,c(menton,l));
      for(k in 1:cut){
        mentoxi[k,]=rbeta(menton,k+k,k)-2/3;
      }
      for(i in 1:menton){
        for(j in 1:l){
          for(k in 1:cut){
            Xtanchork[i,j]=Xtanchork[i,j]+mentoxi[k,i]*phi(k,tanchor[j]);
          }
        }
      }
      for(j in 1:l){
        mentoresult[it,j]=obtainFmurobust(menton,Xtanchork[,j],robu);
      }
    }
    for(j in 1:l){
      truemu[j]=mean(mentoresult[,j]);
    }
  }
  
  #estimate hatmu
  #select hmu;
  candilen=5;
  candidatehmu=exp(seq(log(0.2),log(0.5),length.out=candilen))
  testerrorlist=array(0,candilen);
  for(candi in 1:candilen){
    testerrorlist[candi]=test(n,m,t,X,candidatehmu[candi],robu);
  }
  hmu=candidatehmu[which.min(testerrorlist)];
  cat("hmu=",hmu,"\n");
  #robust mu
  hatmurobust=array(0,l);
  for(k in 1:l){
    hatmurobust[k]=obtainhatmurobust(tanchor[k],n,m,t,X,hmu,robu);
  }
  itresult_mu_square=0;
  itresult_mu_robust=0;
  for(k in 1:l){
    itresult_mu_square=itresult_mu_square+((hatmurobust[k]-truemu[k])**2)/l;
    itresult_mu_robust=itresult_mu_robust+rho(hatmurobust[k]-truemu[k],robu)/l;
  }
  
  #true value of C
  trueC=array(0,c(l,l));
  mentoitermax=100;
  mentoCresult=array(0,c(mentoitermax,l,l));
  for(it in 1:mentoitermax){
    #cat("C mentoit=",it,"\n");
    menton=10**4;
    mentoxi=array(0,c(cut,menton));
    Xtanchork=array(0,c(menton,l));
    if(normalscore==1){
      for(k in 1:cut){
        mentoxi[k,]=rnorm(menton,0,k);
      }
    }
    if(cauchyscore==1){
      for(k in 1:cut){
        mentoxi[k,]=rcauchy(menton,0,k);
      }
    }
    if(tscore==1){
      for(k in 1:cut){
        mentoxi[k,]=rt(menton,k);
      }
    }
    if(slnscore==1){
      for(k in 1:cut){
        mentoxi[k,]=rlnorm(n,0,k)*(floor(runif(menton,0,2))-0.5)*2;
      }
    }
    if(asymscore==1){
      for(k in 1:cut){
        mentoxi[k,]=rbeta(menton,k+k,k)-2/3;
      }
    }
    for(i in 1:menton){
      for(j in 1:l){
        for(k in 1:cut){
          Xtanchork[i,j]=Xtanchork[i,j]+mentoxi[k,i]*phi(k,tanchor[j]);
        }      
      }
    }
    for(j1 in 1:l){
      for(j2 in 1:l){
        mentoCresult[it,j1,j2]=mean(psi(Xtanchork[,j1]-truemu[j1],robu)*psi(Xtanchork[,j2]-truemu[j2],robu));
      }
    }
  }
  #calculate trueC
  for(j1 in 1:l){
    for(j2 in 1:l){
      trueC[j1,j2]=mean(mentoCresult[,j1,j2]);
    }
  }
  
  #estimate hatC
  #rescaled X
  psiX=array(0,c(n,m));
  for(i in 1:n){
    for(j in 1:m){
      psiX[i,j]=psi(X[i,j]-obtainhatmurobust(t[i,j],n,m,t,X,hmu,robu),robu);
    }
  }
  #select hc;
  candilen=5;
  candidatehc=exp(seq(log(0.2),log(0.5),length.out=candilen));
  testerrorlist=array(0,candilen);
  for(candi in 1:candilen){
    testerrorlist[candi]=testC(n,m,t,psiX,candidatehc[candi]);
  }
  hc=candidatehc[which.min(testerrorlist)];
  cat("hc=",hc,"\n");
  #robust C
  hatCrobust=array(0,c(l,l));
  for(k1 in 1:l){
    for(k2 in 1:l){
      hatCrobust[k1,k2]=obtainhatCrobust(tanchor[k1],tanchor[k2],n,m,t,psiX,hc);
    }
  }
  itresultC=sum((hatCrobust-trueC)**2)/(l*l);
  itresultCrela=sum((hatCrobust-trueC)**2)/sum((trueC)**2);
  return(list(itresult_mu_square,itresult_mu_robust,itresultC,itresultCrela));
}

#####################################################################
###  main
#settings
robu=1;
for(scorenumber in c(1,3,4,5,6)){
for(casenumber in c(1:4)){
set.seed(1);
cut=2;#dimension of score
l=10;
tanchor=c(1:l)/l;

#same size
if(casenumber==1){n=100;m=5;}
if(casenumber==2){n=100;m=10;}
if(casenumber==3){n=200;m=5;}
if(casenumber==4){n=200;m=10;}
#score
normalscore=0;
cauchyscore=0;
tscore=0;
slnscore=0;
asymscore=0;
addmodi05=0;
addmodi10=0;
addmodi15=0;
addmodi20=0;
if(scorenumber==1){normalscore=1;}
if(scorenumber==2){cauchyscore=1;}
if(scorenumber==3){tscore=1;}
if(scorenumber==4){slnscore=1;}
if(scorenumber==5){asymscore=1;addmodi10=1;}
if(scorenumber==6){asymscore=1;addmodi20=1;}

#main function
itermax=100;
cl <- makeCluster(50);
registerDoParallel(cl);
result <- foreach(it=c(1:itermax), .combine='c') %dopar% mento(it);
stopCluster(cl);

itresult_mu_square=array(0,itermax);
itresultCrela=array(0,itermax);
for(it in 1:itermax){
  itresult_mu_square[it]=result[[4*it-3]];
  itresultCrela[it]=result[[4*it]];
}
cat(robu,scorenumber,n,m,mean(itresult_mu_square),sd(itresult_mu_square)/sqrt(itermax),mean(itresultCrela),sd(itresultCrela)/sqrt(itermax),"\n")
output=c(robu,scorenumber,n,m,mean(itresult_mu_square),sd(itresult_mu_square)/sqrt(itermax),mean(itresultCrela),sd(itresultCrela)/sqrt(itermax))
write.csv(output,paste0("~/robu1_tuned_kappa_",scorenumber,"_",n,"_",m,sep=""));
}}





