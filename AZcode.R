library("ggplot2")
library("parallel") 
library("foreach")
library("doParallel")
library("pROC")
library("MASS")


load("~/Project 3 RobustFDA/hippo.RData")
###########################wash data
d=3;
iname=names(table(metainfo[,4]));
n=length(iname);
gname=names(table(metainfo[,5]));
#make Y
Y=array(0,c(n,20,d,d));
t=array(-10,c(n,20));
gin=array(0,n);
m=array(0,n);
for(l in 1:977){
    thisid=metainfo[l,4];
    thisi=0;
    for(i in 1:n){
        if(thisid==iname[i]){
            thisi=i;
            break;
        }
    }
    if(thisi==0){
        print(l);
        break;
    }
    m[thisi]=m[thisi]+1;
    Y[thisi,m[thisi],,]=dti[,,l];
    t[thisi,m[thisi]]=metainfo[l,7];
    gin[thisi]=as.numeric(metainfo[l,5]=="CN");
}
#calculate FA
X=array(0,c(n,20));
for(i in 1:n){
  for(j in 1:m[i]){
    value=eigen(Y[i,j,,],only.values=1)$values;
    X[i,j]=sqrt(1.5*sum((value-mean(value))**2)/sum(value*value));
  }
}
#select subject with age in [65,85]
#time domain
l=100;
tbegin=65;
tend=85;
tanchor=seq(tbegin,tend,length.out=l);
#gin=2 means the omitted
for(i in 1:n){
  for(j in 1:m[i]){
    if(t[i,j]<tbegin){gin[i]=2};
    if(t[i,j]>tend){gin[i]=2};
  }
}
#record the washed data
orin=n;rm(n);
oriY=Y;rm(Y);
oriX=X;rm(X);
orit=t;rm(t);
orim=m;rm(m);

###################################################################################################
#estimation
#kernel
Ker=function(s,h){
  temp=ifelse((s-h)>0,0,ifelse((s+h)<0,0,((1-(abs(s/h))**3)**3)*70/(h*81)))
  return(temp);
}
#robust function
rho=function(x,robu){
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
}
psi=function(x,robu){
  if(robu==0){
    return(2*x)
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
}

#obtain hatmurobust
hatQmurobust=function(s,n,m,t,X,hmu,beta,robu){
    hatu1=0;
    hatu2=0;
    for(i in 1:n){
        for(j in 1:m[i]){
            hatu1=hatu1+Ker(t[i,j]-s,hmu)*(t[i,j]-s);
            hatu2=hatu2+Ker(t[i,j]-s,hmu)*(t[i,j]-s)*(t[i,j]-s);
        }
    }
    temp=0;
    for(i in 1:n){
        for(j in 1:m[i]){
            temp=temp+Ker(t[i,j]-s,hmu)*(hatu2-hatu1*(t[i,j]-s))*rho(X[i,j]-beta,robu);
        }
    }
    return(temp);
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
    maxX=max(X,na.rm=1);
    minX=min(X,na.rm=1);
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
#train: first half of X; test: second half of X 
test=function(n,m,t,X,hmu,robu,ordering=initialordering){
    trainindex=ordering[1:floor(n/2)];
    testindex=ordering[(floor(n/2)+1):n];
    testerror=0;
    for(i in testindex){
        for(j in 1:m[i]){
            testerror=testerror+rho(obtainhatmurobust(t[i,j],(n/2),m[trainindex],t[trainindex,],X[trainindex,],hmu,robu)-X[i,j],robu);
        }
    }
    return(testerror);
}
#obtain hatCrobust
obtainhatCrobust=function(s1,s2,n,m,t,rawC,hc){
    S00=0;
    S01=0;
    S10=0;
    S11=0;
    S20=0;
    S02=0;
    R00=0;
    R01=0;
    R10=0;
    for(i in 1:n){
        for(j1 in 1:m[i]){
            for(j2 in 1:m[i]){
                if((abs(t[i,j1]-s1)<hc)&&(abs(t[i,j2]-s2)<hc)&&(j1!=j2)){
                    S00=S00+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc);
                    S10=S10+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc)*((t[i,j1]-s1)/hc);
                    S01=S01+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc)*((t[i,j2]-s2)/hc);
                    S11=S11+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc)*((t[i,j1]-s1)/hc)*((t[i,j2]-s2)/hc);
                    S20=S20+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc)*((t[i,j1]-s1)/hc)*((t[i,j1]-s1)/hc);
                    S02=S02+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc)*((t[i,j2]-s2)/hc)*((t[i,j2]-s2)/hc);
                    R00=R00+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc)*rawC[i,j1,j2];
                    R10=R10+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc)*((t[i,j1]-s1)/hc)*rawC[i,j1,j2];
                    R01=R01+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc)*((t[i,j2]-s2)/hc)*rawC[i,j1,j2];
                }
            }
        }
    }
    S20021111=S20*S02-S11*S11;
    S10020111=S10*S02-S01*S11;
    S10110120=S10*S11-S01*S20;
    if((S20021111*S00-S10020111*S10+S10110120*S01)==0){
        return(0);
    }
    else{
        return(((S20021111*R00-S10020111*R10+S10110120*R01)/(S20021111*S00-S10020111*S10+S10110120*S01)));
    }
}
#train C: first half of X; test: second half of X 
testC=function(n,m,t,rawC,hc,ordering=initialordering){
    trainindex=ordering[1:floor(n/2)];
    testindex=ordering[(floor(n/2)+1):n];
    testerror=0;
    for(i in testindex){
        for(j1 in 1:m[i]){
            for(j2 in 1:m[i]){
                testerror=testerror+(obtainhatCrobust(t[i,j1],t[i,j2],n/2,m,t[trainindex,],rawC[trainindex,,],hc)-rawC[i,j1,j2])**2;
            }
        }
    }
    return(testerror);
}


#set robust loss function
robu=2;
#AZCN, 0 for AZ and 1 for CN 
for(AZCN in 0:1){
n=length(which(gin==AZCN));
X=oriX[which(gin==AZCN),];
t=orit[which(gin==AZCN),];
m=orim[which(gin==AZCN)];
maxm=max(m);
set.seed(97);
initialordering=sample(c(1:n),n);
#plot all
M=sum(m);
Y=array(0,c(M,5));
k=0;
for(i in 1:n){
    for(j in 1:m[i]){
        k=k+1;
        Y[k,1]=as.numeric(t[i,j]);
        Y[k,2]=X[i,j];
        Y[k,3]=i;
    }
}
Y=as.data.frame(Y);
gg=ggplot(Y)+geom_line(aes(x=V1,y=V2,group=V3),color="grey",size=0.5)+
    geom_point(aes(x=V1,y=V2),color="grey",size=0.8)+
    scale_x_continuous(limits=c(tbegin,tend))+scale_y_continuous(limits = c(0,1))+xlab("")+ylab("")
#record, 0 for AZ and 1 for CN 
ggsave(gg,file=paste0(getwd(),"/AD_plot_",AZCN,".pdf",sep="")); 

##########################################################################
rescaleX=X;
#estimate mu
candilen=3;
candidatehmu=(tend-tbegin)*exp(seq(log(0.4),log(0.6),length.out=candilen));
testerrorlist=array(0,candilen);
for(candi in 1:candilen){
    #cat("hmu candi=",candi,"\n");
    testerrorlist[candi]=test(n,m,t,rescaleX,candidatehmu[candi],robu);
}
hmu=candidatehmu[which.min(testerrorlist)];
#cat("hmu=",hmu,"\n");
hatmurobust=array(0,l);
for(k in 1:l){
    #cat("hatmu,k=",k,"\n");
    hatmurobust[k]=obtainhatmurobust(tanchor[k],n,m,t,rescaleX,hmu,robu);
}
#record, 0 for AZ and 1 for CN
save(hmu,file=paste0(getwd(),"/hmu_",AZCN,"_robu",robu,".RData",sep=""))
save(hatmurobust,file=paste0(getwd(),"/hatmurobust_",AZCN,"_robu",robu,".RData",sep=""))

##########################################################################
#estimate C
rawC=array(0,c(n,maxm,maxm));
centerX=array(0,c(n,maxm));
psicenterX=array(0,c(n,maxm));
for(i in 1:n){
    for(j in 1:m[i]){
        centerX[i,j]=rescaleX[i,j]-obtainhatmurobust(t[i,j],n,m,t,rescaleX,hmu,robu);
        psicenterX[i,j]=psi(centerX[i,j],robu);
    }
}
for(i in 1:n){
    for(j1 in 1:m[i]){
        for(j2 in 1:m[i]){
            rawC[i,j1,j2]=psicenterX[i,j1]*psicenterX[i,j2];
        }
    }
}
#select hc;
candilen=3;
candidatehc=(tend-tbegin)*exp(seq(log(0.4),log(0.6),length.out=candilen));
testerrorlist=array(0,candilen);
for(candi in 1:candilen){
    #cat("hc candi=",candi,"\n");
    testerrorlist[candi]=testC(n,m,t,rawC,candidatehc[candi]);
}
hc=candidatehc[which.min(testerrorlist)];
#cat("hc=",hc,"\n");
#calculate hatCrobust
#parallel
parallelobtainhatCrobust=function(k1k2){
    k1=floor((k1k2-1)/l)+1;
    k2=k1k2-(k1-1)*l;
    if(k1<=k2){
        return(obtainhatCrobust(tanchor[k1],tanchor[k2],n,m,t,rawC,hc));
    }
    if(k1>k2){
        return(0);
    }
}
cl <- makeCluster(100);
registerDoParallel(cl);
temp <- foreach(k1k2=c(1:(l*l)), .combine='c') %dopar% parallelobtainhatCrobust(k1k2);
stopCluster(cl);
hatCrobust=array(0,c(l,l));
for(k1k2 in 1:(l*l)){
    k1=floor((k1k2-1)/l)+1;
    k2=k1k2-(k1-1)*l;
    hatCrobust[k1,k2]=temp[k1k2];
}
for(k1 in 1:l){
    for(k2 in 1:l){
        if(k1>k2){
            hatCrobust[k1,k2]=hatCrobust[k2,k1];
        }
    }
}
#eigen-analysis
temp=eigen(hatCrobust);
eigenvalues=temp$values/(l/(tend-tbegin));
eigenvectors=temp$vectors*sqrt(l/(tend-tbegin));# eigenvectors[,1] is the first eigenvector
for(k in 1:l){
  eigenvectors[,k]=eigenvectors[,k]*sign(sum(eigenvectors[,k]))
}

##########################################################################
#score estimation
cut=2;
hatxi_ni=array(0,c(cut,n));
for(k in 1:cut){
for(i in 1:n){
    count=0;
    inte=0;
    for(j in 1:m[i]){
      if((t[i,j]>tbegin)&&(t[i,j]<tend)){
        count=count+1;
        inte=inte+psicenterX[i,j]*eigenvectors[ceiling((t[i,j]-tbegin)*l/(tend-tbegin)),k];
      }}
    hatxi_ni[k,i]=inte/count;
}}

#record, 0 for AZ and 1 for CN 
eigenvalues=eigenvalues/sum(eigenvalues[which(eigenvalues>0)]);
save(hc,file=paste0(getwd(),"/hc_",AZCN,"_robu",robu,".RData",sep=""))
save(eigenvalues,file=paste0(getwd(),"/eigenvalues_",AZCN,"_robu",robu,".RData",sep=""))
save(eigenvectors,file=paste0(getwd(),"/eigenvectors_",AZCN,"_robu",robu,".RData",sep=""))
save(hatxi_ni,file=paste0(getwd(),"/hatxi_ni_",AZCN,"_robu",robu,".RData",sep=""))
}#end of for AZCN 0:1

# load(file=paste0(getwd(),"/eigenvalues_",0,"_robu",3,".RData",sep=""))
# eigenvalues[1:2]/sum(eigenvalues[eigenvalues>0])
  
#load 
AZCN=0;
load(file=paste0(getwd(),"/hatmurobust_",AZCN,"_robu",robu,".RData",sep=""))
load(file=paste0(getwd(),"/eigenvectors_",AZCN,"_robu",robu,".RData",sep=""))
load(file=paste0(getwd(),"/hatxi_ni_",AZCN,"_robu",robu,".RData",sep=""))
hatmurobustAZ=hatmurobust;
eigenvectorsAZ=eigenvectors;
hatxi_niAZ=hatxi_ni;
rm(hatmurobust);rm(eigenvectors);rm(hatxi_ni);
AZCN=1;
load(file=paste0(getwd(),"/hatmurobust_",AZCN,"_robu",robu,".RData",sep=""))
load(file=paste0(getwd(),"/eigenvectors_",AZCN,"_robu",robu,".RData",sep=""))
load(file=paste0(getwd(),"/hatxi_ni_",AZCN,"_robu",robu,".RData",sep=""))
hatmurobustCN=hatmurobust;
eigenvectorsCN=eigenvectors;
hatxi_niCN=hatxi_ni;
rm(hatmurobust);rm(eigenvectors);rm(hatxi_ni);

#plot hatmurobustAZ and hatmurobustCN 
Y1=hatmurobustAZ;
Y2=hatmurobustCN;
YY=as.data.frame(cbind(tanchor,Y1,Y2));
gg=ggplot(YY)+
  geom_line(aes(x=tanchor,y=Y1),color="red",size=1,show.legend=FALSE)+
  geom_line(aes(x=tanchor,y=Y2),color="black",size=1,show.legend=FALSE)+
  scale_x_continuous(limits=c(tbegin,tend))+scale_y_continuous(limits=c(4,7))+
  xlab("")+ylab("")
ggsave(gg,file=paste0(getwd(),"/AD_hatmurobust_robu",robu,".pdf",sep=""));

#plot eigenvectorsAZ and eigenvectorsCN
Y=as.data.frame(cbind(tanchor[1:l],eigenvectorsAZ[,1],eigenvectorsAZ[,2],eigenvectorsCN[,1],eigenvectorsCN[,2]))
gg=ggplot(Y)+
    geom_line(aes(x=V1,y=V2),col="red",size=1,show.legend=FALSE)+
    geom_line(aes(x=V1,y=V3),col="red",linetype="dashed",size=1,show.legend=FALSE)+
    geom_line(aes(x=V1,y=V4),col="black",size=1,show.legend=FALSE)+
    geom_line(aes(x=V1,y=V5),col="black",linetype="dashed",size=1,show.legend=FALSE)+
    scale_x_continuous(limits=c(tbegin,tend))+xlab("")+ylab("")
ggsave(gg,file=paste0(getwd(),"/AD_hateigenfunction_robu",robu,".pdf",sep=""));

#plot hatxi_niAZ and hatxi_niCN, the scores are multiplied by 100
nAZ=dim(hatxi_niAZ)[2];
nCN=dim(hatxi_niCN)[2];
#plot score with labels
group=c(array("0",nAZ),array("1",nCN));
Y <- data.frame(
  V1=c(hatxi_niAZ[1,],hatxi_niCN[1,])*100,
  V2=c(hatxi_niAZ[2,],hatxi_niCN[2,])*100,
  group=group,
  stringsAsFactors = FALSE 
)
gg=ggplot(Y,aes(x=V1,y=V2,color=group,shape=group))+geom_point(show.legend=FALSE)+
  scale_color_manual(values=c("0"="red", "1"="black"))+
  scale_shape_manual(values=c("0"=16, "1"=17))+
    xlab("")+ylab("")
ggsave(gg,file=paste0(getwd(),"/AD_hatxi_robu",robu,".pdf",sep=""));





  
