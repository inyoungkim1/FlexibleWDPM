// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
 
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
arma::vec  cpp1(double inputn, double inputp, double inputQ, double inputtype){
     int n=floor(inputn);
     int p=floor(inputp);
     int type=floor(inputtype);
     vec zv(1);
     double Q=1;
     mat Xmat(n,p),Candimat(1,1),probmat(1,1),mumat(1,1); 
     double alpha=0.05,M=2,Mu0=-0.2,Sz=3,Sn=0.1,fi=0.96,deflat=1,v1=0,W=0.01,temp=0,probsuc=0;
     v1=(1-alpha)*Sz;
     double hp=R::rnorm(0,1);
     vec h(n);
     for(int t = 0; t<n; ++t){
     if(t==0){
     h[0]=hp*fi+R::rnorm(0,sqrt(Sn));
     }
     if(t>0){
     h[t]=h[t-1]*fi+R::rnorm(0,sqrt(Sn)); 
     }
     }
     if(type>2){
     if(type==3){  
     Q=pow(3,p);
     }
     if(type==4){
     Q=inputQ;
     }
     Candimat.ones(Q,p);
     if(type==3){
     double tempz1=0,tempz2=0;
     int tempz3=0;
     vec tempv6(3);
     tempv6[0]=-3;
     tempv6[1]=0;
     tempv6[2]=3;
     for(int j=0; j<p; j++){
     tempz1=pow(3,j);  
     for(int t=0; t<Q; t++){
     tempz2=t;
     tempz3=floor(tempz2/tempz1);
     Candimat(t,j)=tempv6[tempz3%3];
     }
     }
     }
     if(type==4){
     for(int j = 0; j<p; ++j){
     for(int q = 0; q<Q; ++q){
     Candimat(q,j)=R::runif(-5,5);   
     }
     }
     }
     }
     if(type>1){
     probmat.ones(150,Q);
     mumat.ones(150,Q);
     for(int q = 0; q<Q; ++q){
     deflat=1;
     for(int l = 0; l<150; ++l){
     if (l==0){
     probmat(l,q)=R::rbeta(1,M);
     mumat(l,q)=R::rnorm(Mu0,sqrt(v1));
     deflat=1-probmat(l,q);
     }
     if (l>0){
     temp=R::rbeta(1,M);
     probmat(l,q)=deflat*temp;
     mumat(l,q)=R::rnorm(Mu0,sqrt(v1));
     deflat=deflat*(1-temp);
     }
     }
     }     
     }  
     vec y(n),distance(Q),a1(1),a2(1),index=linspace<vec>(1,Q,Q),prob(1),candiv(1),zs(n),mu(n),Z(n),
     tempv(1),stempv(1),sortindex(Q),findv(1),qmat(Q),pmat(Q),wmat(Q),qvec(Q),qvec1(Q);
     y.fill(0);
     qvec.fill(1);
     qvec1.fill(1);
     if(type>2){
     for(int t = 0; t<n; ++t){
     if(t==0){
     for (int lag = 0;  lag< p; ++lag){
     Xmat(0,lag)=R::rnorm(0,0.1);
     }
     }
     if(t>0 && t<p){
     for (int lag = 0;  lag< t; ++lag){
     Xmat(t,lag)=y[t-lag-1];
     }
     for (int lag = t;  lag< p; ++lag){
     Xmat(t,lag)=Xmat(t-1,lag-1);
     }
     }
     if(t>=p){
     for (int lag = 0;  lag< p; ++lag){
     Xmat(t,lag)=y[t-lag-1];
     }
     }
     a1=conv_to< colvec >::from(Xmat.submat(t,0,t,p-1));
     for(int j = 0; j<Q; ++j){
     a2=conv_to< colvec >::from(Candimat.submat(j,0,j,p-1)); 
     distance[j]=sqrt(sum((a1-a2)%(a1-a2)));
     distance[j]=exp(2*log(distance[j]));
     }
     tempv=distance+index/100000000;
     stempv=sort(tempv);
     for(int j = 0; j<Q; ++j){
     findv=index.elem(find(tempv==stempv[j]));
     sortindex[j]=findv[0];
     }
     tempv=1/stempv;
     qvec=qvec1;
     qmat[0]=1;
     for(int j = 0; j<Q; ++j){
     pmat[j]=tempv[j]/sum(tempv.subvec(j,Q-1));
     qvec[j]=1-pmat[j];
     if(j>0){
     qmat[j]=prod(qvec.subvec(0,j-1));
     }
     wmat[j]=(1-0.7*(1-pmat[j]))*qmat[j]*exp(j*log(0.7));
     }
     prob=wmat;
     candiv=sortindex;  
     for (int q = 0;  q< Q; ++q) {
     if (q==0) {
     deflat=1.0;
     }
     if (q>0) {
     deflat=deflat-prob[q-1];
     }
     probsuc=prob[q]/deflat;
     if (R::runif(0,1)<probsuc) {
     Z[t]=candiv[q];
     break;}
     }     
     prob=conv_to< colvec >::from(probmat.submat(0,Z[t]-1,149,Z[t]-1));
     candiv=conv_to< colvec >::from(mumat.submat(0,Z[t]-1,149,Z[t]-1));    
     for (int j = 0;  j< 150; ++j) {
     if (j==0) {
     deflat=1.0;
     }
     if (j>0) {
     deflat=deflat-prob[j-1];
     }
     probsuc=prob[j]/deflat;
     if (R::runif(0,1)<probsuc) {
     mu[t]=candiv[j];
     break;}
     }
     zs[t]=R::rnorm(mu[t],sqrt(alpha*Sz));
     y[t]=h[t]+zs[t];  
     y[t]=exp(y[t]/2);
     if(R::runif(0,1)<W){
     y[t]=0;
     }
     if(R::runif(0,1)<0.5){
     y[t]=-y[t];
     } 
     }
     }
     if(type==2){
     for(int t = 0; t<n; ++t){ 
     prob=conv_to< colvec >::from(probmat.submat(0,0,149,0));
     candiv=conv_to< colvec >::from(mumat.submat(0,0,149,0));    
     for (int j = 0;  j< 150; ++j) {
     if (j==0) {
     deflat=1.0;
     }
     if (j>0) {
     deflat=deflat-prob[j-1];
     }
     probsuc=prob[j]/deflat;
     if (R::runif(0,1)<probsuc) {
     mu[t]=candiv[j];
     break;}
     }
     zs[t]=R::rnorm(mu[t],sqrt(alpha*Sz));
     y[t]=h[t]+zs[t];  
     y[t]=exp(y[t]/2);
     if(R::runif(0,1)<W){
     y[t]=0;
     }
     if(R::runif(0,1)<0.5){
     y[t]=-y[t];
     }
     }
     }
     if(type==1){
     vec mixmean=NumericVector::create(-10.12999, -3.97281, -8.56686, 
     2.77786, 0.61942, 1.79518, -1.08819);
     mixmean=mixmean-1.2704;
     vec mixprob=NumericVector::create(0.0073, 0.10556, 0.00002, 
     0.04395, 0.34001, 0.24566, 0.25750);
     vec mixsig=NumericVector::create(5.79596, 2.61369,5.17950, 
     0.16735 ,0.64009 ,0.34023, 1.26261);
     prob=mixprob;
     candiv=mixmean;  
     vec varv(n);
     for(int t = 0; t<n; ++t){   
     for (int j = 0;  j< 7; ++j) {
     if (j==0) {
     deflat=1.0;
     }
     if (j>0) {
     deflat=deflat-prob[j-1];
     }
     probsuc=prob[j]/deflat;
     if (R::runif(0,1)<probsuc) {
     mu[t]=candiv[j];
     varv[t]=mixsig[j];
     break;}
     }
     zs[t]=R::rnorm(mu[t],sqrt(varv[t]));
     y[t]=h[t]+zs[t];  
     y[t]=exp(y[t]/2);
     if(R::runif(0,1)<0.5){
     y[t]=-y[t];
     } 
     }
     }
     zv=y;
     return zv;
     }

// [[Rcpp::export]]
arma::vec  cpp2(double inputn, double inputp, double inputQ, double inputtype, double method){
     int n=floor(inputn);
     int p=floor(inputp);
     int type=floor(inputtype);
     double Q=1;
     mat Xmat(n,p),Candimat(1,1),probmat(1,1),mumat(1,1); 
     double alpha=0.05,M=2,Mu0=-0.2,Sz=3,Sn=0.1,fi=0.96,deflat=1,v1=0,W=0.01,temp=0,probsuc=0;
     v1=(1-alpha)*Sz;
     double hp=R::rnorm(0,1);
     vec h(n);
     for(int t = 0; t<n; ++t){
     if(t==0){
     h[0]=hp*fi+R::rnorm(0,sqrt(Sn));
     }
     if(t>0){
     h[t]=h[t-1]*fi+R::rnorm(0,sqrt(Sn)); 
     }
     }
     if(type==1){  
     Q=pow(3,p);
     }
     if(type==2){
     Q=inputQ;
     }
     
     Candimat.ones(Q,p);
     if(type==1){
     double tempz1=0,tempz2=0;
     int tempz3=0;
     vec tempv6(3);
     tempv6[0]=-3;
     tempv6[1]=0;
     tempv6[2]=3;
     for(int j=0; j<p; j++){
     tempz1=pow(3,j);  
     for(int t=0; t<Q; t++){
     tempz2=t;
     tempz3=floor(tempz2/tempz1);
     Candimat(t,j)=tempv6[tempz3%3];
     }
     }
     }
     if(type==2){
     for(int j = 0; j<p; ++j){
     for(int q = 0; q<Q; ++q){
     Candimat(q,j)=R::runif(-5,5);   
     }
     }
     }
     probmat.ones(150,Q);
     mumat.ones(150,Q);
     for(int q = 0; q<Q; ++q){
     deflat=1;
     for(int l = 0; l<150; ++l){
     if (l==0){
     probmat(l,q)=R::rbeta(1,M);
     mumat(l,q)=R::rnorm(Mu0,sqrt(v1));
     deflat=1-probmat(l,q);
     }
     if (l>0){
     temp=R::rbeta(1,M);
     probmat(l,q)=deflat*temp;
     mumat(l,q)=R::rnorm(Mu0,sqrt(v1));
     deflat=deflat*(1-temp);
     }
     }
     }
     vec y(n),distance(Q),a1(1),a2(1),index=linspace<vec>(1,Q,Q),prob(1),candiv(1),zs(n),mu(n),Z(n),
     tempv(1),stempv(1),sortindex(Q),findv(1),qmat(Q),pmat(Q),wmat(Q),qvec(Q),qvec1(Q);
     y.fill(0);
     qvec.fill(1);
     qvec1.fill(1);
     vec psiv(Q),gamv(Q);
     psiv.fill(3);
     gamv.fill(1);
     for(int j = 0; j<Q; ++j){
     if(method<2){
     gamv[j]=R::rgamma(0.1,0.01);
     }
     if(method==3){
     psiv[j]=R::rgamma(2,1);
     }
     if(method==4){
     psiv[j]=R::rexp(2.42);
     }
     if(method==5){
     psiv[j]=R::rcauchy(0,0.003);
     }
     }
     double tempdist;
     for(int t = 0; t<n; ++t){
     if(t==0){
     for (int lag = 0;  lag< p; ++lag){
     Xmat(0,lag)=R::rnorm(0,0.1);
     }
     }
     if(t>0 && t<p){
     for (int lag = 0;  lag< t; ++lag){
     Xmat(t,lag)=y[t-lag-1];
     }
     for (int lag = t;  lag< p; ++lag){
     Xmat(t,lag)=Xmat(t-1,lag-1);
     }
     }
     if(t>=p){
     for (int lag = 0;  lag< p; ++lag){
     Xmat(t,lag)=y[t-lag-1];
     }
     }
     a1=conv_to< colvec >::from(Xmat.submat(t,0,t,p-1));
     for(int j = 0; j<Q; ++j){ 
     a2=conv_to< colvec >::from(Candimat.submat(j,0,j,p-1)); 
     tempdist=sum((a1-a2)%(a1-a2));
     wmat[j]=log(gamv[j])-psiv[j]*tempdist;
     }
     wmat=wmat-max(wmat);
     wmat=exp(wmat);
     prob=wmat/sum(wmat);
     for (int q = 0;  q< Q; ++q) {
     if (q==0) {
     deflat=1.0;
     }
     if (q>0) {
     deflat=deflat-prob[q-1];
     }
     probsuc=prob[q]/deflat;
     if (R::runif(0,1)<probsuc) {
     Z[t]=q+1;
     break;}
     }  
     prob=conv_to< colvec >::from(probmat.submat(0,Z[t]-1,149,Z[t]-1));
     candiv=conv_to< colvec >::from(mumat.submat(0,Z[t]-1,149,Z[t]-1));    
     for (int j = 0;  j< 150; ++j) {
     if (j==0) {
     deflat=1.0;
     }
     if (j>0) {
     deflat=deflat-prob[j-1];
     }
     probsuc=prob[j]/deflat;
     if (R::runif(0,1)<probsuc) {
     mu[t]=candiv[j];
     break;}
     }
     zs[t]=R::rnorm(mu[t],sqrt(alpha*Sz));
     y[t]=h[t]+zs[t];  
     y[t]=exp(y[t]/2);
     if(R::runif(0,1)<W){
     y[t]=0;
     }
     if(R::runif(0,1)<0.5){
     y[t]=-y[t];
     } 
     }
     return y;
     }

// [[Rcpp::export]]
Rcpp::List cpp3(arma::vec y,double maxiter,double Nb){ 
vec datay=y;  
int n=y.n_elem;
double logc=-10;
for(int t = 0; t<n; ++t){
y[t]=log(y[t]*y[t]+exp(logc));
} 
vec S(n);
vec seven=linspace<vec>(1,7,7);
for(int i = 0; i<n; ++i){
S[i]=seven[i%7];
}
double fi=0.86,sigmaeta=0.2,mu0=0.3;
vec h(n),m(n),C(n);
vec mixprob(7),mixmean(7),mixvar(7),prob(7);
mixprob[0]=0.0073;
mixprob[1]=0.10556;
mixprob[2]=0.00002;
mixprob[3]=0.04395;
mixprob[4]=0.34001;
mixprob[5]=0.24566;
mixprob[6]=0.2575;
mixmean[0]=-10.12999;
mixmean[1]=-3.97281;
mixmean[2]=-8.56686;
mixmean[3]=2.77786;
mixmean[4]=0.61942;
mixmean[5]=1.79518;
mixmean[6]=-1.08819;
mixmean=mixmean-1.2704;
mixvar[0]=5.79596;
mixvar[1]=2.61369;
mixvar[2]=5.17950;
mixvar[3]=0.16735;
mixvar[4]=0.64009;
mixvar[5]=0.34023;
mixvar[6]=1.26261;
prob.fill(0);
double a=1,R=1,f=1,Q=1,A=1,e=1,V=1,C0=sigmaeta/(1-fi*fi),mean1=1,var1=1,state=0;
double deflat=1,probsuc=1,temp=0,lb=0,ub=0,newfi=0;
double logratio=0;
vec tempv3(1),tempv4(1);
vec hstar(1);
int para=0;
vec vart(n),mt(n);
double temp2=0,temp3=0;
double Nm=maxiter/10.0;
vec fiT(Nm),SnT(Nm),Mu0T(Nm);
Nb=Nb-1;
double dt=0;
double temps=0;
int draw=0;
vec hT(n);
hT.fill(0);
double reptime=0;
double nmax=maxiter+2000;
double fimean=0.0,snmeaninv=0.0,nscore=0.0;
double tempnb=0.0;

for(int iter = 0; iter < nmax; ++iter){
dt=mu0*(1-fi);
state=S[0];
V=mixvar[state-1];
R=fi*fi*C0+sigmaeta;
f=mixmean[state-1];
Q=R+V;
A=R/Q;
e=y[0]-f;
m[0]=A*e;
C[0]=R-A*A*Q;
for(int i = 1; i <n; ++i){
state=S[i];
V=mixvar[state-1];
a=fi*m[i-1]+dt;
R=fi*fi*C[i-1]+sigmaeta;
f=a+mixmean[state-1];
Q=R+V;
A=R/Q;
e=y[i]-f;
m[i]=a+A*e;
C[i]=R-A*A*Q;
}
h[n-1]=R::rnorm(m[n-1],sqrt(C[n-1]));
for(int i = 1; i < n; ++i){
var1=1/(fi*fi/sigmaeta+1/C[n-i-1]);
mean1=var1*(fi*h[n-i]/sigmaeta+m[n-i-1]/C[n-i-1]-dt*fi/sigmaeta);
h[n-i-1]=R::rnorm(mean1,sqrt(var1));
}
for(int i = 0; i<n; ++i){
for (int j = 0;  j< 7; ++j){
temp=y[i]-h[i];
prob[j]=mixprob[j]*R::dnorm(temp,mixmean[j],sqrt(mixvar[j]),0);
}
prob=prob/sum(prob);
for (int j = 0;  j< 7; ++j) {
if (j==0) {
deflat=1.0;
}
if (j!=0) {
deflat=deflat-prob[j-1];
}
probsuc=prob[j]/deflat;
if (R::runif(0,1)<probsuc) {
S[i]=j+1;
break;}
}
}
para=0;
if(R::runif(0,1)<0.5){
para=1;
}
if(para==1){
hstar=h-mu0;
for (int i = 0; i < n; ++i){
vart[i]=mixvar[S[i]-1];
mt[i]=(y[i]-hstar[i]-mixmean[S[i]-1])/vart[i];
}
var1=1/sum(1/vart);
mean1=sum(mt)*var1;
mu0=R::rnorm(mean1,sqrt(var1));
logratio=0;
tempv3=hstar.subvec(1,n-1);
tempv4=hstar.subvec(0,n-2);
mean1=sum(tempv3%tempv4)/sum(tempv4%tempv4);
var1=sigmaeta/sum(tempv4%tempv4);
if(var1<0.0001){
var1=0.0001;
}
lb=R::pnorm(-0.99,mean1,sqrt(var1),1,0);
ub=R::pnorm(0.99,mean1,sqrt(var1),1,0);
newfi=R::qnorm(R::runif(lb,ub),mean1,sqrt(var1),1,0);
logratio=log(R::dnorm(fi,mean1,sqrt(var1),0)/(ub-lb))-log(R::dnorm(newfi,mean1,sqrt(var1),0)/(ub-lb));
lb=R::pnorm(-0.99,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.99,0.95,sqrt(0.5),1,0);
logratio=logratio+log(R::dnorm(newfi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio-(hstar[0]-mu0)*(hstar[0]-mu0)*(1-newfi*newfi)/2/sigmaeta+0.5*log(1-newfi*newfi);
logratio=logratio-log(R::dnorm(fi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio+(hstar[0]-mu0)*(hstar[0]-mu0)*(1-fi*fi)/2/sigmaeta-0.5*log(1-fi*fi);
if(log(R::runif(0,1))<logratio){
fi=newfi;
}

if(iter>maxiter-1){
fi=fimean;
}

temp2=n;
temp2=temp2/2+2.5;
temp3=0.025+hstar[0]*hstar[0]*(1-fi*fi)/2;
for(int i = 1; i<n; ++i){
temp3=temp3+0.5*(hstar[i]-hstar[i-1]*fi)*(hstar[i]-hstar[i-1]*fi);
}
sigmaeta=R::rgamma(temp2,1/temp3);
sigmaeta=1/sigmaeta;

if(iter>maxiter-1){
nscore=nscore+snmeaninv*snmeaninv*R::dgamma(snmeaninv,temp2,1/temp3,0);
}

h=hstar+mu0;
}
if(para==0){
var1=sigmaeta/((n-1)*(1-fi)*(1-fi)+(1-fi*fi));
mean1=var1*((1-fi*fi)/sigmaeta*h[0]+(1-fi)/sigmaeta*sum(h.subvec(1,n-1)-fi*h.subvec(0,n-2)));
mu0=R::rnorm(mean1,sqrt(var1));
logratio=0;
tempv3=h.subvec(1,n-1);
tempv3=tempv3-mu0;
tempv4=h.subvec(0,n-2);
tempv4=tempv4-mu0;
mean1=sum(tempv3%tempv4)/sum(tempv4%tempv4);
var1=sigmaeta/sum(tempv4%tempv4);
if(var1<0.0001){
var1=0.0001;
}
lb=R::pnorm(-0.99,mean1,sqrt(var1),1,0);
ub=R::pnorm(0.99,mean1,sqrt(var1),1,0);
newfi=R::qnorm(R::runif(lb,ub),mean1,sqrt(var1),1,0);
logratio=log(R::dnorm(fi,mean1,sqrt(var1),0)/(ub-lb))-log(R::dnorm(newfi,mean1,sqrt(var1),0)/(ub-lb));
lb=R::pnorm(-0.99,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.99,0.95,sqrt(0.5),1,0);
logratio=logratio+log(R::dnorm(newfi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio-(hstar[0]-mu0)*(hstar[0]-mu0)*(1-newfi*newfi)/2/sigmaeta+0.5*log(1-newfi*newfi);
logratio=logratio-log(R::dnorm(fi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio+(hstar[0]-mu0)*(hstar[0]-mu0)*(1-fi*fi)/2/sigmaeta-0.5*log(1-fi*fi);
if(log(R::runif(0,1))<logratio){
fi=newfi;
}

if(iter>maxiter-1){
fi=fimean;
}

temp2=n;
temp2=temp2/2+2.5;
temp3=0.025+(h[0]-mu0)*(h[0]-mu0)*(1-fi*fi)/2;
for(int i = 1; i<n; ++i){
temp3=temp3+0.5*((h[i]-mu0)-(h[i-1]-mu0)*fi)*((h[i]-mu0)-(h[i-1]-mu0)*fi);
}
sigmaeta=R::rgamma(temp2,1/temp3);
sigmaeta=1/sigmaeta;
if(iter>maxiter-1){
nscore=nscore+snmeaninv*snmeaninv*R::dgamma(snmeaninv,temp2,1/temp3,0);
}

}
if(iter%10==0 && iter<maxiter){
fiT[draw]=fi;
SnT[draw]=sigmaeta;
Mu0T[draw]=mu0;
draw=draw+1;
if(iter>Nb){
hT=hT+h;
reptime=reptime+1;
}
}
if(iter<maxiter && iter>maxiter-2){
tempnb=floor(Nb/10.0);
snmeaninv=mean(SnT.subvec(tempnb,Nm-1));
fimean=mean(fiT.subvec(tempnb,Nm-1));
snmeaninv=1/snmeaninv;
} 
}
Nb=floor(Nb/10.0);
hT=hT/reptime;

double Sn=mean(SnT.subvec(Nb,Nm-1));
fi=mean(fiT.subvec(Nb,Nm-1));

vec fiTsub(1);
fiTsub=fiT.subvec(Nb,Nm-1);
double tempsd=stddev(fiTsub);
double GN=Nm-Nb,mscore=0.0;
double mscore1=0.0,mscore2=0.0;
tempsd=tempsd*exp(-0.2*log(GN));
for(int j = 0; j<GN; ++j){ 
mscore=mscore+R::dnorm(fiTsub[j],fi,tempsd,0);
}
mscore=mscore/GN;
mscore1=log(mscore);
nscore=nscore/2000.0;
mscore2=log(nscore);
nscore=mscore1+mscore2;


double mu1=mean(Mu0T.subvec(Nb,Nm-1));
double score=0;
vec LPS(n), MLV(n);
LPS.fill(0);
MLV.fill(0.0);
double mlscore=0.0;
for(int i = 0; i<10000; ++i){
for(int t = 0; t<n; ++t){
if(t==0){
h[0]=R::rnorm(mu1,sqrt(Sn/(1-fi*fi)));
}
if(t>0){
h[t]=mu1+fi*(h[t-1]-mu1)+R::rnorm(0,sqrt(Sn));
}
for (int j = 0;  j< 7; ++j) {
if (j==0) {
deflat=1.0;
}
if (j!=0) {
deflat=deflat-mixprob[j-1];
}
probsuc=mixprob[j]/deflat;
if (R::runif(0,1)<probsuc) {
temps=j;
break;}
}
LPS[t]=LPS[t]+R::dnorm(y[t]-h[t],mixmean[temps],sqrt(mixvar[temps]),0);
if(datay[t]!=0){
MLV[t]=MLV[t]+R::dnorm(y[t]-h[t],mixmean[temps],sqrt(mixvar[temps]),0); 
}
}
}



LPS=LPS/10000;
MLV=MLV/10000;
LPS=LPS.subvec(10,n-1);
for(int t = 0; t<n; ++t){
if(MLV[t]!=0){
mlscore=mlscore+log(MLV[t]);  
}  
}

mlscore=mlscore+nscore;
lb=R::pnorm(-0.995,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.995,0.95,sqrt(0.5),1,0);
mlscore=mlscore+log(R::dnorm(fi,0.95,sqrt(0.5),0)/(ub-lb));
mlscore=mlscore+log(snmeaninv*snmeaninv*R::dgamma(snmeaninv,2.5,1/0.025,0));


score=mean(log(LPS));
return List::create(Named("fi")=fi,Named("Sn")=Sn,Named("Mu")=mu1, 
Named("fiT")=fiT,Named("SnT")=SnT,Named("MuT")=Mu0T,Named("LPS")=score,Named("MLL")=mlscore,Named("hT")=hT);
}

// [[Rcpp::export]]
Rcpp::List cpp4(arma::vec y,int p,int Q,double maxiter,double Nb){ 
vec datay=y;  
double alpha=0.05,Sz0=0.05,psiub=4;
double logc=-20;
double minp=exp(-50);
double rob1=0.00001;
if(p<3){
rob1=0.001;
}
int n1=y.n_elem;
int n=n1-p;
mat Xmat(n,p);
for(int t = 0; t<n; ++t){
for (int lag = 0;  lag< p; ++lag){
Xmat(t,lag)=y[t+p-lag-1];
}
}
y=y.subvec(p,n1-1);
for(int t = 0; t<n; ++t){
y[t]=log(y[t]*y[t]+exp(logc));
} 
double tempz1=0,tempz2=0,Qmax=pow(3,p);
int tempz3=0;
vec tempv6(3);
tempv6[0]=-3;
tempv6[1]=0;
tempv6[2]=3;
mat candimat(Qmax,p);
for(int j=0; j<p; j++){
tempz1=pow(3,j);  
for(int t=0; t<Qmax; t++){
tempz2=t;
tempz3=floor(tempz2/tempz1);
candimat(t,j)=tempv6[tempz3%3];
}
}
vec tempv4(1),tempv5(1); 
double psi=0.8;
mat distance(n,Qmax),distance1(n,Qmax);
vec index=linspace<vec>(1,Qmax,Qmax),check(Qmax),findv1(1);
int loc=0;
for(int j = 0; j<Qmax; ++j){
for(int i = 0; i<n; ++i){
tempv4=conv_to< colvec >::from(Xmat.submat(i,0,i,p-1));
tempv5=conv_to< colvec >::from(candimat.submat(j,0,j,p-1)); 
distance1(i,j)=sqrt(sum((tempv4-tempv5)%(tempv4-tempv5)))+rob1;
distance(i,j)=log(distance1(i,j));
}
check[j]=mean(distance1.col(j));
}
if(Q>Qmax){
Q=Qmax;
}
if(Q<Qmax){
tempv4=sort(check);
tempv4=tempv4.subvec(0,Q-1);
for(int q = 0; q<Q; ++q){
findv1=index.elem(find(check==tempv4[q]));
loc=findv1[0]-1;
for(int i = 0; i<n; ++i){
distance(i,q)=distance(i,loc);
}
}
}
distance=distance.submat(0,0,n-1,Q-1);
double temp=0;
mat sortdistance(n,Q),sortindex(n,Q);
vec Z(n),countz(Q),S(n),stempv(1);
countz.fill(0);
index=linspace<vec>(1,Q,Q);
for(int i = 0; i<n; ++i){
tempv4=conv_to< colvec >::from(distance.submat(i,0,i,Q-1));
tempv4=tempv4+index/100000000;
stempv=sort(tempv4);
for(int j = 0; j<Q; ++j){
sortdistance(i,j)=stempv[j];
temp=sortdistance(i,j);
findv1=index.elem(find(tempv4==temp));
sortindex(i,j)=findv1[0];
}
Z[i]=sortindex(i,0);
countz[Z[i]-1]=countz[Z[i]-1]+1;
S[i]=R::rpois(6);
}
double tempj=1.0;
vec tempv9(1);
double rv=0.7;
mat qmat(n,Q),pmat(n,Q),wmat(n,Q),wmatn(n,Q),odwmat(n,Q),odwmatn(n,Q),marwmat(n,Q),marwmatn(n,Q);
vec qvec(Q),qvec1(Q);
qvec.fill(1);
qvec1.fill(1);
for(int i = 0; i<n; ++i){
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-psi*tempv5;
qvec=qvec1;
qmat(i,0)=1;
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
tempv9=tempv5.subvec(j,Q-1);
tempv9=tempv9-tempv9[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
pmat(i,j)=tempv9[0];
qvec[j]=1-pmat(i,j);
if(j>0){
qmat(i,j)=prod(qvec.subvec(0,j-1));
}
tempj=j;
wmat(i,j)=(1-rv*(1-pmat(i,j)))*qmat(i,j)*exp(j*log(rv));
odwmat(i,loc-1)=wmat(i,j);
marwmat(i,loc-1)=qmat(i,j)*(1/(tempj+1.0)-(1-pmat(i,j))/(tempj+2.0));
if(tempv9[0]==1){
for(int pw = j+1; pw<Q; ++pw){  
loc=sortindex(i,pw);
wmat(i,pw)=0;
odwmat(i,loc-1)=0;
marwmat(i,loc-1)=0;
}
break;
}
}
}
double fi=0.86,sigmaeta=0.08;
vec h(n);
for(int t = 0; t<n; ++t){
if(t==0){
h[0]=R::rnorm(0,1)*fi+R::rnorm(0,sqrt(sigmaeta));
}
if(t>0){
h[t]=h[t-1]*fi+R::rnorm(0,sqrt(sigmaeta)); 
}
}
S=10000*Z+S;  
int od=1,len=0;
vec tempa(1),Sod(n);
uvec findv2(1);
for (int i = 0; i < n; ++i) {
if(i==0){
if(S[i]==0){
Sod[i]=0;
}
if(S[i]>0){
Sod[i]=1;
od=od+1;
}
}
if(i>0){
if(S[i]==0){
Sod[i]=0;
}
if(S[i]>0){
tempa=S.subvec(0,i-1);
findv2=find(tempa==S[i]);
len=findv2.n_elem;
if(len==0){
Sod[i]=od;
od=od+1;
}
if(len>0){
tempv4=Sod.subvec(0,i-1);
findv1=tempv4.elem(find(tempa==S[i]));
Sod[i]=findv1[0];
}
}
}
}
S=Sod;
Sod=unique(S);
int k=Sod.n_elem;
vec mu(n),mustar(k+1);
for(int i = 1; i<k+1; ++i){
mustar[i]=R::rnorm(0.4,0.1);
}
mustar[0]=logc;
for (int i = 0; i<n; ++i){
mu[i]=mustar[S[i]];
}
double sigmasq=2.5,W=0.01,M=4,mu0=0;
h=ones<vec>(n);
temp=0;
countz=ones<vec>(1);
findv1=ones<vec>(1);
loc=0;
len=0;
k=0;
od=1;
tempv4=ones<vec>(1);
tempv5=ones<vec>(1);
tempa=ones<vec>(1);
Sod=ones<vec>(n);
rv=1;
index=linspace<vec>(1,Q,Q);
vec tempv1(n),tempv2(n); 
double temp1=0,temp2=0,temp3=0,temp4=0,temp5=0,temp6=0,V=Sz0,C0=sigmaeta/(1-fi*fi),
mean1=0,var1=0;
vec countzp(1),countzus(1),locv(1),sp0(1),sp1(1),sp(1),pos(1),ord(1),maxS(1),Zuniq(1);
vec q0v(Q);
q0v.fill(0);
int zt=0,st=0,locs=0,kz=0;
uvec findv(1);
vec index1=linspace<vec>(0,n-1,n);
mat basecountcluster(n,Q),basecluster(n,Q),baseclusterS0(n,Q),
baseclusterS(n,Q),baseconfigS(n,Q),basemu(n,Q);
int position=0,n2=0,lmu=0,addnew=0,nz=0;
double probsum=0,deflat=1,probsuc=0,newmu=1;
vec newS(n),tempS(1),location(1),count(1),tempmustar(1),prob(1),
trueS(1),lamz(1),countz1(1),counts(1);
vec S0v(n+1);
S0v.fill(0);
vec tempmean(1),tempvar(1);
vec zprob(1),sz(1),scount(1),zcluster(1);
mat Q0k(Q,n),zprobmat(1,1);
Q0k.fill(0);
double sm=1e-320;
double newpsi=1,logratio=0;
vec tempv(1),gamv(n);
double a0=1,d0=0.2,templogb=0,tempinvb=1;
vec mix(1),coef(1),mixgam(1);
vec hstar(1),mustar1(1);
vec tempv3(1);
double newfi=0,lb=0,ub=0;
int para=0;
vec x=linspace<vec>(-30,20,1000);
int nx=1000;
mat predmat(nx,Q);
predmat.fill(0);
double reptime=0;
vec tempv7(Q),tempv8(Q);
int draw=0;
double Nm=maxiter/10.0;
vec fiT(Nm),SnT(Nm),SzT(Nm),MT(Nm),WT(Nm),Mu0T(Nm),psiT(Nm);
Nb=Nb-1;
double tempsum=0,tempbr=0;
vec hT(n);
hT.fill(0);
mat mumat(Nm,n);
mat zmat(n,Q);
zmat.fill(0.0);
double reptime1=0.0;
double nmax=maxiter+2000;
double fimean=0.0,snmeaninv=0.0,nscore=0.0;
double tempnb=0.0;

for(int iter = 0; iter <nmax; ++iter){
temp=0;
countz=q0v;
Zuniq=countz;
for(int i = 0; i<n; ++i){
zt=Z[i];
countz[zt-1]=countz[zt-1]+1;
}
for(int j = 0; j<Q; ++j){
if(countz[j]>0){
temp=temp+1;
Zuniq[j]=temp;
}
}
kz=temp;
countzp=q0v.subvec(0,kz-1);
countzus=countzp;
pos=countzp;
ord=countzp;
maxS=countzp;
basecountcluster.fill(0);
basecluster.fill(-1);
baseclusterS.fill(-1);
baseclusterS0.fill(-1);
baseconfigS.fill(0);
basemu.fill(logc);
for(int i = 0; i<n; ++i){
loc=Zuniq[Z[i]-1]-1;
st=S[i];
baseclusterS0(ord[loc],loc)=st;
baseconfigS(ord[loc],loc)=i+1;
if(st==0){
baseclusterS(ord[loc],loc)=0;
basecountcluster(0,loc)=basecountcluster(0,loc)+1;
basecluster(0,loc)=0;
}
if(st>0){
countzp[loc]=countzp[loc]+1;
if(ord[loc]==0){
baseclusterS(0,loc)=1;
maxS[loc]=maxS[loc]+1;
countzus[loc]=countzus[loc]+1;
pos[loc]=pos[loc]+1; 
basecluster(pos[loc],loc)=st;
basemu(pos[loc],loc)=mustar[st];
basecountcluster(pos[loc],loc)=basecountcluster(pos[loc],loc)+1;
}
if(ord[loc]>0){
sp0=baseclusterS0.col(loc);
sp0=sp0.subvec(0,ord[loc]-1);
findv=find(sp0==st);
len=findv.n_elem;
if(len==0){
countzus[loc]=countzus[loc]+1;
pos[loc]=pos[loc]+1;
basecluster(pos[loc],loc)=st;
basemu(pos[loc],loc)=mustar[st];
basecountcluster(pos[loc],loc)=basecountcluster(pos[loc],loc)+1;
maxS[loc]=maxS[loc]+1;
baseclusterS(ord[loc],loc)=maxS[loc];
}
if(len>0){
sp1=baseclusterS.col(loc);
sp1=sp1.subvec(0,ord[loc]-1);
findv1=sp1.elem(find(sp0==st));
baseclusterS(ord[loc],loc)=findv1[0];
sp=basecluster.col(loc);
locv=index1.elem(find(sp==st));
locs=locv[0];
basecountcluster(locs,loc)=basecountcluster(locs,loc)+1;
}
}
}
ord[loc]=ord[loc]+1;
}
countz1=ones<vec>(kz);
od=-1;
for(int i = 0; i < Q; ++i) {
if(countz[i]>0){
od=od+1;
countz1[od]=countz[i];
}
}
temp2=a0+k-kz;
mix=ones<vec>(kz+1);
coef=mix-1;
mix[0]=lgamma(temp2);
temp3=M+1;
temp4=d0;
temp4=temp4-log(R::rbeta(temp3,countz1[0]));
coef[0]=countz1[0];
coef[1]=1;
if(kz>1){
for(int i = 1; i<kz; ++i){ 
temp4=temp4-log(R::rbeta(temp3,countz1[i]));
mix[i]=lgamma(temp2+i);
tempv=coef.subvec(0,i);
coef.subvec(1,i+1)=tempv;
coef[0]=0;
coef.subvec(0,i)=coef.subvec(0,i)+tempv*countz1[i];
coef.subvec(0,i+1)=coef.subvec(0,i+1)/max(coef.subvec(0,i+1));
}
}
mix[kz]=lgamma(temp2+kz);
templogb=log(temp4);
mix[0]=mix[0]+log(coef[0])-temp2*templogb;
temp5=mix[0];
mixgam=mix; 
tempinvb=1/temp4;   
mixgam[0]=R::rgamma(temp2,tempinvb);
for(int i = 1; i<kz+1; ++i){ 
temp1=temp2+i;
mix[i]=mix[i]+log(coef[i])-temp1*templogb;
if(mix[i]>temp5){
temp5=mix[i];
}
mixgam[i]=R::rgamma(temp1,tempinvb);
}
mix=mix-temp5;
mix=exp(mix);
mix=mix/sum(mix);
for (int q = 0;  q< kz+1; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-mix[q-1];
}
probsuc=mix[q]/deflat;
if (R::runif(0,1)<probsuc) {
M=mixgam[q];
break;
}
}
temp=alpha*sigmasq;
if(S[0]>0){
V=temp;
}
temp2=fi*fi*C0+sigmaeta;
temp3=mu[0];
temp4=temp2+V;
temp5=temp2/temp4;
temp6=y[0]-temp3;
tempv1[0]=temp5*temp6;
tempv2[0]=temp2-temp5*temp5*temp4;
for(int i = 1; i < n; ++i){
V=Sz0;
if(S[i]>0) {
V=temp;
}
temp1=fi*tempv1[i-1];
temp2=fi*fi*tempv2[i-1]+sigmaeta;
temp3=temp1+mu[i];
temp4=temp2+V;
temp5=temp2/temp4;
temp6=y[i]-temp3;
tempv1[i]=temp1+temp5*temp6;
tempv2[i]=temp2-temp5*temp5*temp4;
}
h[n-1]=R::rnorm(tempv1[n-1],sqrt(tempv2[n-1]));
for(int i = 1; i < n; ++i){
var1=1/(fi*fi/sigmaeta+1/tempv2[n-i-1]);
mean1=var1*(fi*h[n-i]/sigmaeta+tempv1[n-i-1]/tempv2[n-i-1]);
h[n-i-1]=R::rnorm(mean1,sqrt(var1));
}
for(int l = 0; l<kz; ++l){
nz=countz1[l];
tempS=conv_to< colvec >::from(baseclusterS.submat(0,l,nz-1,l));
location=conv_to< colvec >::from(baseconfigS.submat(0,l,nz-1,l));
count=conv_to< colvec >::from(basecountcluster.submat(0,l,countzus[l],l));
tempmustar=conv_to< colvec >::from(basemu.submat(0,l,countzus[l],l));  
lmu=countzus[l]+1;
n2=countzp[l]; 
lamz=ones<vec>(nz);
for(int i = 0; i<nz; ++i){
addnew=0;
st=tempS[i];
position=location[i]-1;
count[st]=count[st]-1;
if(st>0){
n2=n2-1;
}
if (count[st]==0 && st>0){
count.shed_row(st);
tempmustar=tempmustar.elem(find(tempmustar!=tempmustar[st]));
tempS.elem(find(tempS>st))=tempS.elem(find(tempS>st))-1;
lmu=lmu-1;
}
prob=ones<vec>(lmu+1);
prob[0]=W*exp(-0.5*(y[position]-h[position]-logc)*(y[position]-h[position]-logc)/Sz0)/sqrt(Sz0);
probsum=prob[0];
if(lmu>1){
for (int j = 1;  j< lmu; ++j) {
prob[j]=(1-W)*count[j]/sqrt(alpha*sigmasq)/(n2-1+M)*
exp(-0.5*(y[position]-h[position]-tempmustar[j])*(y[position]-h[position]-tempmustar[j])/alpha/sigmasq);
probsum=probsum+prob[j];
}
}
prob[lmu]=(1-W)*M/sqrt(sigmasq)/(n2-1+M)*
exp(-(y[position]-h[position]-mu0)*(y[position]-h[position]-mu0)/sigmasq/2);  
probsum=probsum+prob[lmu];
prob=prob/probsum;
for (int j = 0;  j< lmu+1; ++j) {
if (j==0) {
deflat=1.0;
}
if (j!=0) {
deflat=deflat-prob[j-1];
}
probsuc=prob[j]/deflat;
if (R::runif(0,1)<probsuc) {
tempS[i]=j;
break;
}
}
if(tempS[i]>0){
n2=n2+1;
if(tempS[i]>lmu-1){
addnew=1;
}
}
if(addnew==1){
var1=1/(1 / (alpha * sigmasq) + 1 / ((1 - alpha) * sigmasq));
mean1=var1*((y[position]-h[position])/alpha/sigmasq+mu0/(1-alpha)/sigmasq);
newmu=R::rnorm(mean1,sqrt(var1));
tempv4=tempmustar.subvec(0,lmu-1);
tempmustar=ones<vec>(lmu+1);
tempmustar.subvec(0,lmu-1)=tempv4;
tempmustar[lmu]=newmu;
tempv4=count;
count=ones<vec>(lmu+1);
count.subvec(0,lmu-1)=tempv4;
count[lmu]=1;
lamz[i]=newmu;
lmu=lmu+1;
}
if(addnew==0){
st=tempS[i];
count[st]=count[st]+1;
lamz[i]=tempmustar[st];
}
}
od=0;
trueS=ones<vec>(nz);
trueS.fill(1);
position=location[0]-1;
if(tempS[0]==0){
newS[position]=0; 
}  
if(tempS[0]>0){
newS[position]=1+10000*l;
od=od+1;
trueS[0]=od;
}
if(nz>1){
for(int i = 1; i<nz; ++i){
position=location[i]-1;  
if(tempS[i]==0){
newS[position]=0;
}
if(tempS[i]>0){
tempv5=lamz.subvec(0,i);
findv1=trueS.elem(find(tempv5==lamz[i]));
len=findv1.n_elem;
if(len==1){
od=od+1;
trueS[i]=od;
}
if(len>1){
trueS[i]=findv1[0];
}
newS[position]=trueS[i]+10000*l;
}
}       
}
}
tempS=newS;
od=1;
len=0;
for (int i = 0; i < n; ++i) {
if(i==0){
if(tempS[i]==0){
Sod[i]=0;
}
if(tempS[i]>0){
Sod[i]=1;
od=od+1;
}
}
if(i>0){
if(tempS[i]==0){
Sod[i]=0;
}
if(tempS[i]>0){
tempa=tempS.subvec(0,i-1);
findv=find(tempa==tempS[i]);
len=findv.n_elem;
if(len==0){
Sod[i]=od;
od=od+1;
}
if(len>0){
tempv5=Sod.subvec(0,i-1);
findv1=tempv5.elem(find(tempa==tempS[i]));
Sod[i]=findv1[0];
}
}
}
}
S=Sod;
counts=S0v.subvec(0,od-1);  
for (int i = 0; i < n; ++i) {
counts[S[i]]=counts[S[i]]+1;
}
n2=n-counts[0];
k=od-1;
if(iter<maxiter){
if(iter>Nb && iter%10==9){
hT=hT+h;
tempv7.fill(0);
tempv8.fill(0);
od=-1;
for (int q = 0;  q< Q; ++q){
if(countz[q]>0){
od=od+1;
tempv7[q]=countzus[od];
tempv8[q]=countzp[od];
}
}
for(int i = 0; i<nx; ++i){
for (int q = 0;  q< Q; ++q){
temp6=W*R::dnorm(x[i],logc,sqrt(Sz0),0);
if(tempv7[q]>0){ 
for (int l = 1;  l< tempv7[q]+1; ++l){
temp6=temp6+(1-W)*basecountcluster(l,q)/(tempv8[q]+M)*
R::dnorm(x[i],basemu(l,q),sqrt(alpha*sigmasq),0);
}
}
temp6=temp6+(1-W)*M/(tempv7[q]+M)*R::dnorm(x[i],mu0,sqrt(sigmasq),0);
predmat(i,q)=predmat(i,q)+temp6;      
}          
}
reptime=reptime+1;
}
}
if(Q>1){
if(k>0.5){ 
scount=counts.subvec(1,k); 
sz=S0v.subvec(0,k-1);
zcluster=sz;
zprobmat=Q0k.submat(0,0,Q-1,k-1);
tempv4=q0v;
for(int i = 0; i<n; ++i){
if(S[i]>0){
sz[S[i]-1]=Z[i];
zprob=conv_to< colvec >::from(odwmat.submat(i,0,i,Q-1));
zprob=zprob+sm;
zprobmat.col(S[i]-1)=zprobmat.col(S[i]-1)+log(zprob);
tempv4[Z[i]-1]=tempv4[Z[i]-1]+1;
}
}
for(int j = 0; j<k; ++j){
tempv4[sz[j]-1]=tempv4[sz[j]-1]-scount[j];
for(int q = 0; q<Q; ++q){
zprobmat(q,j)=zprobmat(q,j)+lgamma(M+tempv4[q])-lgamma(M+tempv4[q]+scount[j]);
}
zprob=conv_to< colvec >::from(zprobmat.submat(0,j,Q-1,j));
zprob=zprob-max(zprob);
zprob=exp(zprob);
zprob=zprob/sum(zprob);
for (int q = 0;  q< Q; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-zprob[q-1];
}
probsuc=zprob[q]/deflat;
if (R::runif(0,1)<probsuc) {
zcluster[j]=q+1;
break;}
}
tempv4[zcluster[j]-1]=tempv4[zcluster[j]-1]+scount[j];
}
} 
for(int i = 0; i<n; ++i){
if(S[i]==0){
zprob=conv_to< colvec >::from(odwmat.submat(i,0,i,Q-1));
for (int q = 0;  q< Q; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-zprob[q-1];
}
probsuc=zprob[q]/deflat;
if (R::runif(0,1)<probsuc) {
Z[i]=q+1;
break;
}
}
}
if(S[i]>0){
Z[i]=zcluster[S[i]-1];
}
}
countz=q0v;
for(int i = 0; i<n; ++i){
countz[Z[i]-1]=countz[Z[i]-1]+1;
tempv=conv_to< colvec >::from(sortindex.submat(i,0,i,(Q-1)));
locv=index.elem(find(tempv==Z[i]));
loc=locv[0];
temp1=1-pmat(i,loc-1);
if(loc<Q){
rv=R::qbeta(R::runif(0,1)*R::pbeta(temp1,loc,2,1,0),loc,2,1,0)/temp1;
}
if(loc==Q){
rv=exp(log(R::runif(0,1))/Q);
}
gamv[i]=rv;
if(rv==0){
gamv[i]=0.00001;
}
}
psiub=log(4.0/psi);
newpsi=R::qnorm(R::runif(0,R::pnorm(psiub,-0.001,0.05,1,0)),-0.001,0.05,1,0);
newpsi=psi*exp(newpsi);
logratio=0;
if(n<601){
logratio=-R::dnorm(log(psi),0.7,0.05,1)+R::dnorm(log(newpsi),0.7,0.05,1);
}
if(n>600){
logratio=-R::dnorm(log(psi),0.08,0.05,1)+R::dnorm(log(newpsi),0.08,0.05,1);
}
for(int i = 0; i<n; ++i){
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-newpsi*tempv5;
tempv9=tempv5-tempv5[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
qvec=qvec1;
qmat(i,0)=1;
rv=gamv[i];
tempsum=0;
tempbr=0;
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
pmat(i,j)=tempv9[j]/(1-tempsum);
qvec[j]=1-pmat(i,j);
if(j>0){
qmat(i,j)=prod(qvec.subvec(0,j-1));
}
tempj=j;
wmatn(i,j)=(1-rv*(1-pmat(i,j)))*qmat(i,j)*exp(j*log(rv));
odwmatn(i,loc-1)=wmatn(i,j);
marwmatn(i,loc-1)=qmat(i,j)*(1/(tempj+1.0)-(1-pmat(i,j))/(tempj+2.0));
tempsum=tempsum+tempv9[j];
if(tempsum==1 || qmat(i,j)<minp){
tempsum=j+1;
tempbr=1;
break;
}
}
if(tempbr>0){
for(int pw = tempsum; pw<Q; ++pw){  
loc=sortindex(i,pw);
wmatn(i,pw)=0;
odwmatn(i,loc-1)=0;
marwmatn(i,loc-1)=0;
}
}
logratio=logratio+log(marwmatn(i,Z[i]-1))-log(marwmat(i,Z[i]-1));
}
if(log(R::runif(0,1))<logratio){
psi=newpsi;
wmat=wmatn;
odwmat=odwmatn;
marwmat=marwmatn;
} 
}

if(iter<maxiter){
if(iter%10==2 && iter>Nb){
for(int i = 0; i < n; ++i) {
zmat(i,Z[i]-1)=zmat(i,Z[i]-1)+1;
}
reptime1=reptime1+1;
}
}

mustar=S0v.subvec(0,k);
mustar[0]=logc;
if(k>0.5){
tempmean=S0v.subvec(0,k-1);
tempvar=tempmean;
for(int i = 0; i < n; ++i) {
st=S[i];
if (st>0){
tempmean[st-1]=(y[i]-h[i])*(1-alpha)/(alpha+counts[st]*(1-alpha))+tempmean[st-1];
}
}
for(int j = 0; j < k; ++j) {
tempmean[j]=tempmean[j]+mu0*alpha/(alpha+counts[j+1]*(1-alpha));
tempvar[j]=sigmasq/(counts[j+1]/alpha+1/(1-alpha)); 
mustar[j+1]=R::rnorm(tempmean[j],sqrt(tempvar[j]));
}
}
for(int i = 0; i < n; ++i) {
mu[i]=mustar[S[i]];
if(iter%10==0 && iter<maxiter){
mumat(draw,i)=mu[i];
}
}

para=0;
if(R::runif(0,1)<0.5){
para=1;
}
if(para==1){
hstar=h+mu0;
mustar1=mustar.subvec(1,k)-mu0;
var1=sigmaeta/((n-1)*(1-fi)*(1-fi)+(1-fi*fi));
mean1=var1*((1-fi*fi)/sigmaeta*hstar[0]+(1-fi)/sigmaeta*sum(hstar.subvec(1,n-1)-fi*hstar.subvec(0,n-2)));
temp6=mu0;
mu0=R::rnorm(mean1,sqrt(var1));
temp2=n2+k;
temp2=temp2/2+1;
temp3=3;
for(int i = 0; i<n; ++i){
if(S[i]>0){
temp3=temp3+(y[i]-h[i]-mu[i])*(y[i]-h[i]-mu[i])/2/alpha;
}
}
temp3=temp3+sum(mustar1%mustar1)/2/(1-alpha);
sigmasq=R::rgamma(temp2,1/temp3);
sigmasq=1/sigmasq;
logratio=0;
tempv3=hstar.subvec(1,n-1);
tempv3=tempv3-mu0;
tempv4=hstar.subvec(0,n-2);
tempv4=tempv4-mu0;
mean1=sum(tempv3%tempv4)/sum(tempv4%tempv4);
var1=sigmaeta/sum(tempv4%tempv4);
if(var1<0.0001){
var1=0.0001;
}
lb=R::pnorm(-0.995,mean1,sqrt(var1),1,0);
ub=R::pnorm(0.995,mean1,sqrt(var1),1,0);
newfi=R::qnorm(R::runif(lb,ub),mean1,sqrt(var1),1,0);
logratio=log(R::dnorm(fi,mean1,sqrt(var1),0)/(ub-lb))-log(R::dnorm(newfi,mean1,sqrt(var1),0)/(ub-lb));
lb=R::pnorm(-0.995,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.995,0.95,sqrt(0.5),1,0);
logratio=logratio+log(R::dnorm(newfi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio-(hstar[0]-mu0)*(hstar[0]-mu0)*(1-newfi*newfi)/2/sigmaeta+0.5*log(1-newfi*newfi);
logratio=logratio-log(R::dnorm(fi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio+(hstar[0]-mu0)*(hstar[0]-mu0)*(1-fi*fi)/2/sigmaeta-0.5*log(1-fi*fi);
if(log(R::runif(0,1))<logratio){
fi=newfi;
}
if(iter>maxiter-1){
fi=fimean;
}


temp2=n;
temp2=temp2/2+2.5;
temp3=0.025;
for(int i = 1; i<n; ++i){
temp3=temp3+0.5*((hstar[i]-mu0)-(hstar[i-1]-mu0)*fi)*((hstar[i]-mu0)-(hstar[i-1]-mu0)*fi);
}
sigmaeta=R::rgamma(temp2,1/temp3);
sigmaeta=1/sigmaeta;

if(iter>maxiter-1){
nscore=nscore+snmeaninv*snmeaninv*R::dgamma(snmeaninv,temp2,1/temp3,0);
}

h=hstar-mu0;
mustar.subvec(1,k)=mustar1+mu0;
for(int i = 0; i<n; ++i){
if(S[i]>0){
mu[i]=mu[i]-temp6+mu0;
}
}
}
if(para==0){
mean1=(sum(mustar)-mustar[0])/k;
var1=(1-alpha)*sigmasq/k;
mu0=R::rnorm(mean1,sqrt(var1));
temp2=n2+k;
temp2=temp2/2+1;
temp3=3;
for(int i = 0; i<n; ++i){
if(S[i]>0){
temp3=temp3+(y[i]-h[i]-mu[i])*(y[i]-h[i]-mu[i])/2/alpha;
}
}
temp3=temp3+sum((mustar.subvec(1,k)-mu0)%(mustar.subvec(1,k)-mu0))/2/(1-alpha);
sigmasq=R::rgamma(temp2,1/temp3);
sigmasq=1/sigmasq;
logratio=0;
tempv3=h.subvec(1,n-1);
tempv4=h.subvec(0,n-2);
mean1=sum(tempv3%tempv4)/sum(tempv4%tempv4);
var1=sigmaeta/sum(tempv4%tempv4);
if(var1<0.0001){
var1=0.0001;
}
lb=R::pnorm(-0.995,mean1,sqrt(var1),1,0);
ub=R::pnorm(0.995,mean1,sqrt(var1),1,0);
newfi=R::qnorm(R::runif(lb,ub),mean1,sqrt(var1),1,0);
logratio=log(R::dnorm(fi,mean1,sqrt(var1),0)/(ub-lb))-log(R::dnorm(newfi,mean1,sqrt(var1),0)/(ub-lb));
lb=R::pnorm(-0.995,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.995,0.95,sqrt(0.5),1,0);
logratio=logratio+log(R::dnorm(newfi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio-h[0]*h[0]*(1-newfi*newfi)/2/sigmaeta+0.5*log(1-newfi*newfi);
logratio=logratio-log(R::dnorm(fi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio+h[0]*h[0]*(1-fi*fi)/2/sigmaeta-0.5*log(1-fi*fi);
if(log(R::runif(0,1))<logratio){
fi=newfi;
}

if(iter>maxiter-1){
fi=fimean;
}


temp2=n;
temp2=temp2/2+2.5;
temp3=0.025;
for(int i = 1; i<n; ++i){
temp3=temp3+0.5*(h[i]-h[i-1]*fi)*(h[i]-h[i-1]*fi);
}
sigmaeta=R::rgamma(temp2,1/temp3);
sigmaeta=1/sigmaeta;

if(iter>maxiter-1){
nscore=nscore+snmeaninv*snmeaninv*R::dgamma(snmeaninv,temp2,1/temp3,0);
}

}
temp2=n-n2;
temp2=temp2+0.1;
temp3=n2;
temp3=temp3+0.9;
W=R::rbeta(temp2,temp3);
if(iter%10==0 && iter<maxiter){
fiT[draw]=fi;
SnT[draw]=sigmaeta;
SzT[draw]=sigmasq;
MT[draw]=M;
WT[draw]=W;
Mu0T[draw]=mu0;
psiT[draw]=psi;
draw=draw+1;
}
if(iter<maxiter && iter>maxiter-2){
tempnb=floor(Nb/10.0);
snmeaninv=mean(SnT.subvec(tempnb,Nm-1));
fimean=mean(fiT.subvec(tempnb,Nm-1));
snmeaninv=1/snmeaninv;
}  
}
for(int iter = 0; iter < Nm; ++iter){  
newpsi=R::qnorm(R::runif(0,R::pnorm(log(5.0/psi),-0.5,1,1,0)),-0.5,1,1,0);
newpsi=psi*exp(newpsi);
logratio=5*log(newpsi)-5*log(psi)+2*(psi-newpsi);
for(int i = 0; i<1; ++i){
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-newpsi*tempv5;
tempv9=tempv5-tempv5[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
qvec=qvec1;
qmat(i,0)=1;
tempsum=0;
tempbr=0;
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
pmat(i,j)=tempv9[j]/(1-tempsum);
qvec[j]=1-pmat(i,j);
if(j>0){
qmat(i,j)=prod(qvec.subvec(0,j-1));
}
tempj=j;
marwmatn(i,loc-1)=qmat(i,j)*(1/(tempj+1.0)-(1-pmat(i,j))/(tempj+2.0));
tempsum=tempsum+tempv9[j];
if(tempsum==1 || qmat(i,j)<minp){
tempsum=j+1;
tempbr=1;
break;
}
}
if(tempbr>0){
for(int pw = tempsum; pw<Q; ++pw){  
loc=sortindex(i,pw);
marwmatn(i,loc-1)=0;
}
}
logratio=logratio+log(marwmatn(i,Z[i]-1))-log(marwmat(i,Z[i]-1));
}
if(log(R::runif(0,1))<logratio){
psi=newpsi;
marwmat=marwmatn;
}
psiT[iter]=psi;
}
predmat=predmat/reptime;
zmat=zmat/reptime1;

hT=hT/reptime;
hT=exp(hT/2.0);
for(int q = 0; q<Q; ++q){
tempv=predmat.col(q);
temp=sum(tempv.subvec(1,999))+sum(tempv.subvec(0,998));
temp=temp/2.0;
temp=temp*50.0/999;
tempv=tempv/temp;
predmat.col(q)=tempv;
}
Nb=floor(Nb/10.0);
psi=mean(psiT.subvec(Nb,Nm-1));
for(int i = 0; i<n; ++i){
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-psi*tempv5;
tempv9=tempv5-tempv5[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
qvec=qvec1;
qmat(i,0)=1;
tempsum=0;
tempbr=0;
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
pmat(i,j)=tempv9[j]/(1-tempsum);
qvec[j]=1-pmat(i,j);
if(j>0){
qmat(i,j)=prod(qvec.subvec(0,j-1));
}
temp=j;
odwmat(i,loc-1)=qmat(i,j)*(1/(temp+1.0)-(1-pmat(i,j))/(temp+2.0));
tempsum=tempsum+tempv9[j];
if(tempsum==1 || qmat(i,j)<minp){
tempsum=j+1;
tempbr=1;
break;
}
}
if(tempbr>0){
for(int pw = tempsum; pw<Q; ++pw){  
loc=sortindex(i,pw);
odwmat(i,loc-1)=0;
}
}
}
temp1=mean(SnT.subvec(Nb,Nm-1));
temp2=mean(fiT.subvec(Nb,Nm-1));
temp3=0;

vec fiTsub(1);
fiTsub=fiT.subvec(Nb,Nm-1);
double tempsd=stddev(fiTsub);
double GN=Nm-Nb,mscore=0.0;
double mscore1=0.0,mscore2=0.0;

tempsd=tempsd*exp(-0.2*log(GN));
for(int j = 0; j<GN; ++j){ 
mscore=mscore+R::dnorm(fiTsub[j],temp2,tempsd,0);
}
mscore=mscore/GN;
mscore1=log(mscore);
nscore=nscore/2000.0;
mscore2=log(nscore);
nscore=mscore1+mscore2;
vec pred0(1),pred1(1),predv(1),odwv(1),odwv1(1),LPS(n),MLV(n);
LPS.fill(0);
MLV.fill(0);
double err=0,mlscore=0.0;
double step=50.0/999;
for(int i = 0; i<10000; ++i){
if(err>1){
break;
}
for(int t = 0; t<n; ++t){
if(err>1){
break;
}
if(t==0){
h[0]=R::rnorm(0,sqrt(temp1/(1-temp2*temp2)));
}
if(t>0){
h[t]=temp2*h[t-1]+R::rnorm(0,sqrt(temp1));
}
temp4=y[t]-h[t];
if(temp4<x[0]){
err=1000;
}
if(temp4>x[nx-1]){
err=1000;
}
if(err==0){
loc=floor((temp4-x[0])/step);
pred0=conv_to< colvec >::from(predmat.submat(loc,0,loc,Q-1));
pred1=conv_to< colvec >::from(predmat.submat(loc+1,0,loc+1,Q-1));
predv=(temp4-x[loc])/step*(pred1-pred0)+pred0;
odwv=conv_to< colvec >::from(odwmat.submat(t,0,t,Q-1));
odwv1=conv_to< colvec >::from(zmat.submat(t,0,t,Q-1));
LPS[t]=LPS[t]+sum(predv%odwv);
MLV[t]=MLV[t]+sum(predv%odwv1);
}
}
}
if(err<1){
LPS=LPS/10000;
MLV=MLV/10000;
LPS=LPS.subvec(10-p,n-1);
temp3=mean(log(LPS));
for(int t = 10-p; t<n; ++t){
if(datay[t]!=0){
mlscore=mlscore+log(MLV[t]);
}
}
mlscore=mlscore+nscore;
}
if(err>1){
temp3=datum::nan;
mlscore=datum::nan;
}
sigmasq=mean(SzT.subvec(Nb,Nm-1));
W=mean(WT.subvec(Nb,Nm-1));
mu0=mean(Mu0T.subvec(Nb,Nm-1));
M=mean(MT.subvec(Nb,Nm-1));

lb=R::pnorm(-0.995,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.995,0.95,sqrt(0.5),1,0);
mlscore=mlscore+log(R::dnorm(temp2,0.95,sqrt(0.5),0)/(ub-lb));
mlscore=mlscore+log(snmeaninv*snmeaninv*R::dgamma(snmeaninv,2.5,1/0.025,0));

return List::create(
Named("fi")=temp2,Named("Sn")=temp1,Named("mu0")=mu0,Named("Sz")=sigmasq,Named("M")=M,Named("W")=W,Named("psi")=psi,
Named("fiT")=fiT, Named("SnT")=SnT,Named("Mu0T")=Mu0T,Named("SzT")=SzT,Named("MT")=MT,Named("WT")=WT,Named("psiT")=psiT,
Named("LPS")=temp3,Named("MLL")=mlscore,Named("hT")=hT,Named("mumat")=mumat);
} 

// [[Rcpp::export]]
Rcpp::List cpp5(arma::vec y,int p,int Q,double maxiter,double Nb){ 
vec datay=y;    
double alpha=0.05,Sz0=0.05,psiub=4;
double logc=-20;
double minp=exp(-50);
int n1=y.n_elem;
int n=n1-p;
mat Xmat(n,p);
for(int t = 0; t<n; ++t){
for (int lag = 0;  lag< p; ++lag){
Xmat(t,lag)=y[t+p-lag-1];
}
}
y=y.subvec(p,n1-1);
for(int t = 0; t<n; ++t){
y[t]=log(y[t]*y[t]+exp(logc));
} 
double Qmax=Q;
mat candimat(Qmax,p);
for(int j=0; j<p; j++){  
for(int t=0; t<Qmax; t++){
candimat(t,j)=R::runif(-5,5);
}
}
vec tempv4(1),tempv5(1); 
double psi=0.8;
mat distance(n,Qmax),distance1(n,Qmax);
vec index=linspace<vec>(1,Qmax,Qmax),check(Qmax),findv1(1);
int loc=0;
for(int j = 0; j<Qmax; ++j){
for(int i = 0; i<n; ++i){
tempv4=conv_to< colvec >::from(Xmat.submat(i,0,i,p-1));
tempv5=conv_to< colvec >::from(candimat.submat(j,0,j,p-1)); 
distance1(i,j)=sqrt(sum((tempv4-tempv5)%(tempv4-tempv5)));
distance(i,j)=log(distance1(i,j));
}
check[j]=mean(distance1.col(j));
}
if(Q>Qmax){
Q=Qmax;
}
if(Q<Qmax){
tempv4=sort(check);
tempv4=tempv4.subvec(0,Q-1);
for(int q = 0; q<Q; ++q){
findv1=index.elem(find(check==tempv4[q]));
loc=findv1[0]-1;
for(int i = 0; i<n; ++i){
distance(i,q)=distance(i,loc);
}
}
}
distance=distance.submat(0,0,n-1,Q-1);
double temp=0;
mat sortdistance(n,Q),sortindex(n,Q);
vec Z(n),countz(Q),S(n),stempv(1);
countz.fill(0);
index=linspace<vec>(1,Q,Q);
for(int i = 0; i<n; ++i){
tempv4=conv_to< colvec >::from(distance.submat(i,0,i,Q-1));
tempv4=tempv4+index/100000000;
stempv=sort(tempv4);
for(int j = 0; j<Q; ++j){
sortdistance(i,j)=stempv[j];
temp=sortdistance(i,j);
findv1=index.elem(find(tempv4==temp));
sortindex(i,j)=findv1[0];
}
Z[i]=sortindex(i,0);
countz[Z[i]-1]=countz[Z[i]-1]+1;
S[i]=R::rpois(6);
}
double tempj=1.0;
vec tempv9(1);
double rv=0.7;
mat qmat(n,Q),pmat(n,Q),wmat(n,Q),wmatn(n,Q),odwmat(n,Q),odwmatn(n,Q),marwmat(n,Q),marwmatn(n,Q);
vec qvec(Q),qvec1(Q);
qvec.fill(1);
qvec1.fill(1);
for(int i = 0; i<n; ++i){
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-psi*tempv5;
qvec=qvec1;
qmat(i,0)=1;
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
tempv9=tempv5.subvec(j,Q-1);
tempv9=tempv9-tempv9[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
pmat(i,j)=tempv9[0];
qvec[j]=1-pmat(i,j);
if(j>0){
qmat(i,j)=prod(qvec.subvec(0,j-1));
}
tempj=j;
wmat(i,j)=(1-rv*(1-pmat(i,j)))*qmat(i,j)*exp(j*log(rv));
odwmat(i,loc-1)=wmat(i,j);
marwmat(i,loc-1)=qmat(i,j)*(1/(tempj+1.0)-(1-pmat(i,j))/(tempj+2.0));
if(tempv9[0]==1){
for(int pw = j+1; pw<Q; ++pw){  
loc=sortindex(i,pw);
wmat(i,pw)=0;
odwmat(i,loc-1)=0;
marwmat(i,loc-1)=0;
}
break;
}
}
}
double fi=0.86,sigmaeta=0.08;
vec h(n);
for(int t = 0; t<n; ++t){
if(t==0){
h[0]=R::rnorm(0,1)*fi+R::rnorm(0,sqrt(sigmaeta));
}
if(t>0){
h[t]=h[t-1]*fi+R::rnorm(0,sqrt(sigmaeta)); 
}
}
S=10000*Z+S;  
int od=1,len=0;
vec tempa(1),Sod(n);
uvec findv2(1);
for (int i = 0; i < n; ++i) {
if(i==0){
if(S[i]==0){
Sod[i]=0;
}
if(S[i]>0){
Sod[i]=1;
od=od+1;
}
}
if(i>0){
if(S[i]==0){
Sod[i]=0;
}
if(S[i]>0){
tempa=S.subvec(0,i-1);
findv2=find(tempa==S[i]);
len=findv2.n_elem;
if(len==0){
Sod[i]=od;
od=od+1;
}
if(len>0){
tempv4=Sod.subvec(0,i-1);
findv1=tempv4.elem(find(tempa==S[i]));
Sod[i]=findv1[0];
}
}
}
}
S=Sod;
Sod=unique(S);
int k=Sod.n_elem;
vec mu(n),mustar(k+1);
for(int i = 1; i<k+1; ++i){
mustar[i]=R::rnorm(0.4,0.1);
}
mustar[0]=logc;
for (int i = 0; i<n; ++i){
mu[i]=mustar[S[i]];
}
double sigmasq=2.5,W=0.01,M=4,mu0=0;
h=ones<vec>(n);
temp=0;
countz=ones<vec>(1);
findv1=ones<vec>(1);
loc=0;
len=0;
k=0;
od=1;
tempv4=ones<vec>(1);
tempv5=ones<vec>(1);
tempa=ones<vec>(1);
Sod=ones<vec>(n);
rv=1;
index=linspace<vec>(1,Q,Q);
vec tempv1(n),tempv2(n); 
double temp1=0,temp2=0,temp3=0,temp4=0,temp5=0,temp6=0,V=Sz0,C0=sigmaeta/(1-fi*fi),
mean1=0,var1=0;
vec countzp(1),countzus(1),locv(1),sp0(1),sp1(1),sp(1),pos(1),ord(1),maxS(1),Zuniq(1);
vec q0v(Q);
q0v.fill(0);
int zt=0,st=0,locs=0,kz=0;
uvec findv(1);
vec index1=linspace<vec>(0,n-1,n);
mat basecountcluster(n,Q),basecluster(n,Q),baseclusterS0(n,Q),
baseclusterS(n,Q),baseconfigS(n,Q),basemu(n,Q);
int position=0,n2=0,lmu=0,addnew=0,nz=0;
double probsum=0,deflat=1,probsuc=0,newmu=1;
vec newS(n),tempS(1),location(1),count(1),tempmustar(1),prob(1),
trueS(1),lamz(1),countz1(1),counts(1);
vec S0v(n+1);
S0v.fill(0);
vec tempmean(1),tempvar(1);
vec zprob(1),sz(1),scount(1),zcluster(1);
mat Q0k(Q,n),zprobmat(1,1);
Q0k.fill(0);
double sm=1e-320;
double newpsi=1,logratio=0;
vec tempv(1),gamv(n);
double a0=1,d0=0.2,templogb=0,tempinvb=1;
vec mix(1),coef(1),mixgam(1);
vec hstar(1),mustar1(1);
vec tempv3(1);
double newfi=0,lb=0,ub=0;
int para=0;
vec x=linspace<vec>(-30,20,1000);
int nx=1000;
mat predmat(nx,Q);
predmat.fill(0);
double reptime=0;
vec tempv7(Q),tempv8(Q);
int draw=0;
double Nm=maxiter/10.0;
vec fiT(Nm),SnT(Nm),SzT(Nm),MT(Nm),WT(Nm),Mu0T(Nm),psiT(Nm);
Nb=Nb-1;
double tempsum=0,tempbr=0;
vec hT(n);
hT.fill(0);
double stp=0.2;
mat mumat(Nm,n);
mat zmat(n,Q);
zmat.fill(0.0);
double reptime1=0.0;
double nmax=maxiter+2000;
double fimean=0.0,snmeaninv=0.0,nscore=0.0;
double tempnb=0.0;

for(int iter = 0; iter < nmax; ++iter){
temp=0;
countz=q0v;
Zuniq=countz;
for(int i = 0; i<n; ++i){
zt=Z[i];
countz[zt-1]=countz[zt-1]+1;
}
for(int j = 0; j<Q; ++j){
if(countz[j]>0){
temp=temp+1;
Zuniq[j]=temp;
}
}
kz=temp;
countzp=q0v.subvec(0,kz-1);
countzus=countzp;
pos=countzp;
ord=countzp;
maxS=countzp;
basecountcluster.fill(0);
basecluster.fill(-1);
baseclusterS.fill(-1);
baseclusterS0.fill(-1);
baseconfigS.fill(0);
basemu.fill(logc);
for(int i = 0; i<n; ++i){
loc=Zuniq[Z[i]-1]-1;
st=S[i];
baseclusterS0(ord[loc],loc)=st;
baseconfigS(ord[loc],loc)=i+1;
if(st==0){
baseclusterS(ord[loc],loc)=0;
basecountcluster(0,loc)=basecountcluster(0,loc)+1;
basecluster(0,loc)=0;
}
if(st>0){
countzp[loc]=countzp[loc]+1;
if(ord[loc]==0){
baseclusterS(0,loc)=1;
maxS[loc]=maxS[loc]+1;
countzus[loc]=countzus[loc]+1;
pos[loc]=pos[loc]+1; 
basecluster(pos[loc],loc)=st;
basemu(pos[loc],loc)=mustar[st];
basecountcluster(pos[loc],loc)=basecountcluster(pos[loc],loc)+1;
}
if(ord[loc]>0){
sp0=baseclusterS0.col(loc);
sp0=sp0.subvec(0,ord[loc]-1);
findv=find(sp0==st);
len=findv.n_elem;
if(len==0){
countzus[loc]=countzus[loc]+1;
pos[loc]=pos[loc]+1;
basecluster(pos[loc],loc)=st;
basemu(pos[loc],loc)=mustar[st];
basecountcluster(pos[loc],loc)=basecountcluster(pos[loc],loc)+1;
maxS[loc]=maxS[loc]+1;
baseclusterS(ord[loc],loc)=maxS[loc];
}
if(len>0){
sp1=baseclusterS.col(loc);
sp1=sp1.subvec(0,ord[loc]-1);
findv1=sp1.elem(find(sp0==st));
baseclusterS(ord[loc],loc)=findv1[0];
sp=basecluster.col(loc);
locv=index1.elem(find(sp==st));
locs=locv[0];
basecountcluster(locs,loc)=basecountcluster(locs,loc)+1;
}
}
}
ord[loc]=ord[loc]+1;
}
countz1=ones<vec>(kz);
od=-1;
for(int i = 0; i < Q; ++i) {
if(countz[i]>0){
od=od+1;
countz1[od]=countz[i];
}
}
temp2=a0+k-kz;
mix=ones<vec>(kz+1);
coef=mix-1;
mix[0]=lgamma(temp2);
temp3=M+1;
temp4=d0;
temp4=temp4-log(R::rbeta(temp3,countz1[0]));
coef[0]=countz1[0];
coef[1]=1;
if(kz>1){
for(int i = 1; i<kz; ++i){ 
temp4=temp4-log(R::rbeta(temp3,countz1[i]));
mix[i]=lgamma(temp2+i);
tempv=coef.subvec(0,i);
coef.subvec(1,i+1)=tempv;
coef[0]=0;
coef.subvec(0,i)=coef.subvec(0,i)+tempv*countz1[i];
coef.subvec(0,i+1)=coef.subvec(0,i+1)/max(coef.subvec(0,i+1));
}
}
mix[kz]=lgamma(temp2+kz);
templogb=log(temp4);
mix[0]=mix[0]+log(coef[0])-temp2*templogb;
temp5=mix[0];
mixgam=mix; 
tempinvb=1/temp4;   
mixgam[0]=R::rgamma(temp2,tempinvb);
for(int i = 1; i<kz+1; ++i){ 
temp1=temp2+i;
mix[i]=mix[i]+log(coef[i])-temp1*templogb;
if(mix[i]>temp5){
temp5=mix[i];
}
mixgam[i]=R::rgamma(temp1,tempinvb);
}
mix=mix-temp5;
mix=exp(mix);
mix=mix/sum(mix);
for (int q = 0;  q< kz+1; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-mix[q-1];
}
probsuc=mix[q]/deflat;
if (R::runif(0,1)<probsuc) {
M=mixgam[q];
break;
}
}
temp=alpha*sigmasq;
if(S[0]>0){
V=temp;
}
temp2=fi*fi*C0+sigmaeta;
temp3=mu[0];
temp4=temp2+V;
temp5=temp2/temp4;
temp6=y[0]-temp3;
tempv1[0]=temp5*temp6;
tempv2[0]=temp2-temp5*temp5*temp4;
for(int i = 1; i < n; ++i){
V=Sz0;
if(S[i]>0) {
V=temp;
}
temp1=fi*tempv1[i-1];
temp2=fi*fi*tempv2[i-1]+sigmaeta;
temp3=temp1+mu[i];
temp4=temp2+V;
temp5=temp2/temp4;
temp6=y[i]-temp3;
tempv1[i]=temp1+temp5*temp6;
tempv2[i]=temp2-temp5*temp5*temp4;
}
h[n-1]=R::rnorm(tempv1[n-1],sqrt(tempv2[n-1]));
for(int i = 1; i < n; ++i){
var1=1/(fi*fi/sigmaeta+1/tempv2[n-i-1]);
mean1=var1*(fi*h[n-i]/sigmaeta+tempv1[n-i-1]/tempv2[n-i-1]);
h[n-i-1]=R::rnorm(mean1,sqrt(var1));
}
for(int l = 0; l<kz; ++l){
nz=countz1[l];
tempS=conv_to< colvec >::from(baseclusterS.submat(0,l,nz-1,l));
location=conv_to< colvec >::from(baseconfigS.submat(0,l,nz-1,l));
count=conv_to< colvec >::from(basecountcluster.submat(0,l,countzus[l],l));
tempmustar=conv_to< colvec >::from(basemu.submat(0,l,countzus[l],l));  
lmu=countzus[l]+1;
n2=countzp[l]; 
lamz=ones<vec>(nz);
for(int i = 0; i<nz; ++i){
addnew=0;
st=tempS[i];
position=location[i]-1;
count[st]=count[st]-1;
if(st>0){
n2=n2-1;
}
if (count[st]==0 && st>0){
count.shed_row(st);
tempmustar=tempmustar.elem(find(tempmustar!=tempmustar[st]));
tempS.elem(find(tempS>st))=tempS.elem(find(tempS>st))-1;
lmu=lmu-1;
}
prob=ones<vec>(lmu+1);
prob[0]=W*exp(-0.5*(y[position]-h[position]-logc)*(y[position]-h[position]-logc)/Sz0)/sqrt(Sz0);
probsum=prob[0];
if(lmu>1){
for (int j = 1;  j< lmu; ++j) {
prob[j]=(1-W)*count[j]/sqrt(alpha*sigmasq)/(n2-1+M)*
exp(-0.5*(y[position]-h[position]-tempmustar[j])*(y[position]-h[position]-tempmustar[j])/alpha/sigmasq);
probsum=probsum+prob[j];
}
}
prob[lmu]=(1-W)*M/sqrt(sigmasq)/(n2-1+M)*
exp(-(y[position]-h[position]-mu0)*(y[position]-h[position]-mu0)/sigmasq/2);  
probsum=probsum+prob[lmu];
prob=prob/probsum;
for (int j = 0;  j< lmu+1; ++j) {
if (j==0) {
deflat=1.0;
}
if (j!=0) {
deflat=deflat-prob[j-1];
}
probsuc=prob[j]/deflat;
if (R::runif(0,1)<probsuc) {
tempS[i]=j;
break;
}
}
if(tempS[i]>0){
n2=n2+1;
if(tempS[i]>lmu-1){
addnew=1;
}
}
if(addnew==1){
var1=1/(1 / (alpha * sigmasq) + 1 / ((1 - alpha) * sigmasq));
mean1=var1*((y[position]-h[position])/alpha/sigmasq+mu0/(1-alpha)/sigmasq);
newmu=R::rnorm(mean1,sqrt(var1));
tempv4=tempmustar.subvec(0,lmu-1);
tempmustar=ones<vec>(lmu+1);
tempmustar.subvec(0,lmu-1)=tempv4;
tempmustar[lmu]=newmu;
tempv4=count;
count=ones<vec>(lmu+1);
count.subvec(0,lmu-1)=tempv4;
count[lmu]=1;
lamz[i]=newmu;
lmu=lmu+1;
}
if(addnew==0){
st=tempS[i];
count[st]=count[st]+1;
lamz[i]=tempmustar[st];
}
}
od=0;
trueS=ones<vec>(nz);
trueS.fill(1);
position=location[0]-1;
if(tempS[0]==0){
newS[position]=0; 
}  
if(tempS[0]>0){
newS[position]=1+10000*l;
od=od+1;
trueS[0]=od;
}
if(nz>1){
for(int i = 1; i<nz; ++i){
position=location[i]-1;  
if(tempS[i]==0){
newS[position]=0;
}
if(tempS[i]>0){
tempv5=lamz.subvec(0,i);
findv1=trueS.elem(find(tempv5==lamz[i]));
len=findv1.n_elem;
if(len==1){
od=od+1;
trueS[i]=od;
}
if(len>1){
trueS[i]=findv1[0];
}
newS[position]=trueS[i]+10000*l;
}
}       
}
}
tempS=newS;
od=1;
len=0;
for (int i = 0; i < n; ++i) {
if(i==0){
if(tempS[i]==0){
Sod[i]=0;
}
if(tempS[i]>0){
Sod[i]=1;
od=od+1;
}
}
if(i>0){
if(tempS[i]==0){
Sod[i]=0;
}
if(tempS[i]>0){
tempa=tempS.subvec(0,i-1);
findv=find(tempa==tempS[i]);
len=findv.n_elem;
if(len==0){
Sod[i]=od;
od=od+1;
}
if(len>0){
tempv5=Sod.subvec(0,i-1);
findv1=tempv5.elem(find(tempa==tempS[i]));
Sod[i]=findv1[0];
}
}
}
}
S=Sod;
counts=S0v.subvec(0,od-1);  
for (int i = 0; i < n; ++i) {
counts[S[i]]=counts[S[i]]+1;
}
n2=n-counts[0];
k=od-1;
if(iter<maxiter){
if(iter>Nb && iter%10==9){
hT=hT+h;
tempv7.fill(0);
tempv8.fill(0);
od=-1;
for (int q = 0;  q< Q; ++q){
if(countz[q]>0){
od=od+1;
tempv7[q]=countzus[od];
tempv8[q]=countzp[od];
}
}
for(int i = 0; i<nx; ++i){
for (int q = 0;  q< Q; ++q){
temp6=W*R::dnorm(x[i],logc,sqrt(Sz0),0);
if(tempv7[q]>0){ 
for (int l = 1;  l< tempv7[q]+1; ++l){
temp6=temp6+(1-W)*basecountcluster(l,q)/(tempv8[q]+M)*
R::dnorm(x[i],basemu(l,q),sqrt(alpha*sigmasq),0);
}
}
temp6=temp6+(1-W)*M/(tempv7[q]+M)*R::dnorm(x[i],mu0,sqrt(sigmasq),0);
predmat(i,q)=predmat(i,q)+temp6;      
}          
}
reptime=reptime+1;
}
}
if(Q>1){
if(k>0.5){ 
scount=counts.subvec(1,k); 
sz=S0v.subvec(0,k-1);
zcluster=sz;
zprobmat=Q0k.submat(0,0,Q-1,k-1);
tempv4=q0v;
for(int i = 0; i<n; ++i){
if(S[i]>0){
sz[S[i]-1]=Z[i];
zprob=conv_to< colvec >::from(odwmat.submat(i,0,i,Q-1));
zprob=zprob+sm;
zprobmat.col(S[i]-1)=zprobmat.col(S[i]-1)+log(zprob);
tempv4[Z[i]-1]=tempv4[Z[i]-1]+1;
}
}
for(int j = 0; j<k; ++j){
tempv4[sz[j]-1]=tempv4[sz[j]-1]-scount[j];
for(int q = 0; q<Q; ++q){
zprobmat(q,j)=zprobmat(q,j)+lgamma(M+tempv4[q])-lgamma(M+tempv4[q]+scount[j]);
}
zprob=conv_to< colvec >::from(zprobmat.submat(0,j,Q-1,j));
zprob=zprob-max(zprob);
zprob=exp(zprob);
zprob=zprob/sum(zprob);
for (int q = 0;  q< Q; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-zprob[q-1];
}
probsuc=zprob[q]/deflat;
if (R::runif(0,1)<probsuc) {
zcluster[j]=q+1;
break;}
}
tempv4[zcluster[j]-1]=tempv4[zcluster[j]-1]+scount[j];
}
} 
for(int i = 0; i<n; ++i){
if(S[i]==0){
zprob=conv_to< colvec >::from(odwmat.submat(i,0,i,Q-1));
for (int q = 0;  q< Q; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-zprob[q-1];
}
probsuc=zprob[q]/deflat;
if (R::runif(0,1)<probsuc) {
Z[i]=q+1;
break;
}
}
}
if(S[i]>0){
Z[i]=zcluster[S[i]-1];
}
}
countz=q0v;
for(int i = 0; i<n; ++i){
countz[Z[i]-1]=countz[Z[i]-1]+1;
tempv=conv_to< colvec >::from(sortindex.submat(i,0,i,(Q-1)));
locv=index.elem(find(tempv==Z[i]));
loc=locv[0];
temp1=1-pmat(i,loc-1);
if(loc<Q){
rv=R::qbeta(R::runif(0,1)*R::pbeta(temp1,loc,2,1,0),loc,2,1,0)/temp1;
}
if(loc==Q){
rv=exp(log(R::runif(0,1))/Q);
}
gamv[i]=rv;
if(rv==0){
gamv[i]=0.00001;
}
}
newpsi=R::qnorm(R::runif(R::pnorm(0,psi,stp,1,0),R::pnorm(psiub,psi,stp,1,0)),psi,stp,1,0);
logratio=0;
logratio=logratio+3*log(newpsi)-3*log(psi)+2*(psi-newpsi);;
lb=R::pnorm(0,psi,stp,1,0);
ub=R::pnorm(psiub,psi,stp,1,0);
logratio=logratio+log(ub-lb);
lb=R::pnorm(0,newpsi,stp,1,0);
ub=R::pnorm(psiub,newpsi,stp,1,0);
logratio=logratio-log(ub-lb);
for(int i = 0; i<n; ++i){
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-newpsi*tempv5;
tempv9=tempv5-tempv5[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
qvec=qvec1;
qmat(i,0)=1;
rv=gamv[i];
tempsum=0;
tempbr=0;
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
pmat(i,j)=tempv9[j]/(1-tempsum);
qvec[j]=1-pmat(i,j);
if(j>0){
qmat(i,j)=prod(qvec.subvec(0,j-1));
}
tempj=j;
wmatn(i,j)=(1-rv*(1-pmat(i,j)))*qmat(i,j)*exp(j*log(rv));
odwmatn(i,loc-1)=wmatn(i,j);
marwmatn(i,loc-1)=qmat(i,j)*(1/(tempj+1.0)-(1-pmat(i,j))/(tempj+2.0));
tempsum=tempsum+tempv9[j];
if(tempsum==1 || qmat(i,j)<minp){
tempsum=j+1;
tempbr=1;
break;
}
}
if(tempbr>0){
for(int pw = tempsum; pw<Q; ++pw){  
loc=sortindex(i,pw);
wmatn(i,pw)=0;
odwmatn(i,loc-1)=0;
marwmatn(i,loc-1)=0;
}
}
logratio=logratio+log(marwmatn(i,Z[i]-1))-log(marwmat(i,Z[i]-1));
}
if(log(R::runif(0,1))<logratio){
psi=newpsi;
wmat=wmatn;
odwmat=odwmatn;
marwmat=marwmatn;
} 
}

if(iter<maxiter){
if(iter%10==2 && iter>Nb){
for(int i = 0; i < n; ++i) {
zmat(i,Z[i]-1)=zmat(i,Z[i]-1)+1;
}
reptime1=reptime1+1;
}
}


mustar=S0v.subvec(0,k);
mustar[0]=logc;
if(k>0.5){
tempmean=S0v.subvec(0,k-1);
tempvar=tempmean;
for(int i = 0; i < n; ++i) {
st=S[i];
if (st>0){
tempmean[st-1]=(y[i]-h[i])*(1-alpha)/(alpha+counts[st]*(1-alpha))+tempmean[st-1];
}
}
for(int j = 0; j < k; ++j) {
tempmean[j]=tempmean[j]+mu0*alpha/(alpha+counts[j+1]*(1-alpha));
tempvar[j]=sigmasq/(counts[j+1]/alpha+1/(1-alpha)); 
mustar[j+1]=R::rnorm(tempmean[j],sqrt(tempvar[j]));
}
}
for(int i = 0; i < n; ++i) {
mu[i]=mustar[S[i]];
if(iter%10==0 && iter<maxiter){
mumat(draw,i)=mu[i];
}
}
para=0;
if(R::runif(0,1)<0.5){
para=1;
}
if(para==1){
hstar=h+mu0;
mustar1=mustar.subvec(1,k)-mu0;
var1=sigmaeta/((n-1)*(1-fi)*(1-fi)+(1-fi*fi));
mean1=var1*((1-fi*fi)/sigmaeta*hstar[0]+(1-fi)/sigmaeta*sum(hstar.subvec(1,n-1)-fi*hstar.subvec(0,n-2)));
temp6=mu0;
mu0=R::rnorm(mean1,sqrt(var1));
temp2=n2+k;
temp2=temp2/2+1;
temp3=3;
for(int i = 0; i<n; ++i){
if(S[i]>0){
temp3=temp3+(y[i]-h[i]-mu[i])*(y[i]-h[i]-mu[i])/2/alpha;
}
}
temp3=temp3+sum(mustar1%mustar1)/2/(1-alpha);
sigmasq=R::rgamma(temp2,1/temp3);
sigmasq=1/sigmasq;
logratio=0;
tempv3=hstar.subvec(1,n-1);
tempv3=tempv3-mu0;
tempv4=hstar.subvec(0,n-2);
tempv4=tempv4-mu0;
mean1=sum(tempv3%tempv4)/sum(tempv4%tempv4);
var1=sigmaeta/sum(tempv4%tempv4);
if(var1<0.0001){
var1=0.0001;
}
lb=R::pnorm(-0.995,mean1,sqrt(var1),1,0);
ub=R::pnorm(0.995,mean1,sqrt(var1),1,0);
newfi=R::qnorm(R::runif(lb,ub),mean1,sqrt(var1),1,0);
logratio=log(R::dnorm(fi,mean1,sqrt(var1),0)/(ub-lb))-log(R::dnorm(newfi,mean1,sqrt(var1),0)/(ub-lb));
lb=R::pnorm(-0.995,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.995,0.95,sqrt(0.5),1,0);
logratio=logratio+log(R::dnorm(newfi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio-(hstar[0]-mu0)*(hstar[0]-mu0)*(1-newfi*newfi)/2/sigmaeta+0.5*log(1-newfi*newfi);
logratio=logratio-log(R::dnorm(fi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio+(hstar[0]-mu0)*(hstar[0]-mu0)*(1-fi*fi)/2/sigmaeta-0.5*log(1-fi*fi);
if(log(R::runif(0,1))<logratio){
fi=newfi;
}

if(iter>maxiter-1){
fi=fimean;
}

temp2=n;
temp2=temp2/2+2.5;
temp3=0.025;
for(int i = 1; i<n; ++i){
temp3=temp3+0.5*((hstar[i]-mu0)-(hstar[i-1]-mu0)*fi)*((hstar[i]-mu0)-(hstar[i-1]-mu0)*fi);
}
sigmaeta=R::rgamma(temp2,1/temp3);
sigmaeta=1/sigmaeta;


if(iter>maxiter-1){
nscore=nscore+snmeaninv*snmeaninv*R::dgamma(snmeaninv,temp2,1/temp3,0);
}


h=hstar-mu0;
mustar.subvec(1,k)=mustar1+mu0;
for(int i = 0; i<n; ++i){
if(S[i]>0){
mu[i]=mu[i]-temp6+mu0;
}
}
}
if(para==0){
mean1=(sum(mustar)-mustar[0])/k;
var1=(1-alpha)*sigmasq/k;
mu0=R::rnorm(mean1,sqrt(var1));
temp2=n2+k;
temp2=temp2/2+1;
temp3=3;
for(int i = 0; i<n; ++i){
if(S[i]>0){
temp3=temp3+(y[i]-h[i]-mu[i])*(y[i]-h[i]-mu[i])/2/alpha;
}
}
temp3=temp3+sum((mustar.subvec(1,k)-mu0)%(mustar.subvec(1,k)-mu0))/2/(1-alpha);
sigmasq=R::rgamma(temp2,1/temp3);
sigmasq=1/sigmasq;
logratio=0;
tempv3=h.subvec(1,n-1);
tempv4=h.subvec(0,n-2);
mean1=sum(tempv3%tempv4)/sum(tempv4%tempv4);
var1=sigmaeta/sum(tempv4%tempv4);
if(var1<0.0001){
var1=0.0001;
}
lb=R::pnorm(-0.995,mean1,sqrt(var1),1,0);
ub=R::pnorm(0.995,mean1,sqrt(var1),1,0);
newfi=R::qnorm(R::runif(lb,ub),mean1,sqrt(var1),1,0);
logratio=log(R::dnorm(fi,mean1,sqrt(var1),0)/(ub-lb))-log(R::dnorm(newfi,mean1,sqrt(var1),0)/(ub-lb));
lb=R::pnorm(-0.995,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.995,0.95,sqrt(0.5),1,0);
logratio=logratio+log(R::dnorm(newfi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio-h[0]*h[0]*(1-newfi*newfi)/2/sigmaeta+0.5*log(1-newfi*newfi);
logratio=logratio-log(R::dnorm(fi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio+h[0]*h[0]*(1-fi*fi)/2/sigmaeta-0.5*log(1-fi*fi);
if(log(R::runif(0,1))<logratio){
fi=newfi;
}

if(iter>maxiter-1){
fi=fimean;
}


temp2=n;
temp2=temp2/2+2.5;
temp3=0.025;
for(int i = 1; i<n; ++i){
temp3=temp3+0.5*(h[i]-h[i-1]*fi)*(h[i]-h[i-1]*fi);
}
sigmaeta=R::rgamma(temp2,1/temp3);
sigmaeta=1/sigmaeta;
}

if(iter>maxiter-1){
nscore=nscore+snmeaninv*snmeaninv*R::dgamma(snmeaninv,temp2,1/temp3,0);
}

temp2=n-n2;
temp2=temp2+0.1;
temp3=n2;
temp3=temp3+0.9;
W=R::rbeta(temp2,temp3);
if(iter%10==0 && iter<maxiter){
fiT[draw]=fi;
SnT[draw]=sigmaeta;
SzT[draw]=sigmasq;
MT[draw]=M;
WT[draw]=W;
Mu0T[draw]=mu0;
psiT[draw]=psi;
draw=draw+1;
}
if(iter<maxiter && iter>maxiter-2){
tempnb=floor(Nb/10.0);
snmeaninv=mean(SnT.subvec(tempnb,Nm-1));
fimean=mean(fiT.subvec(tempnb,Nm-1));
snmeaninv=1/snmeaninv;
}  
}
for(int iter = 0; iter < Nm; ++iter){  
newpsi=R::qnorm(R::runif(0,R::pnorm(log(5.0/psi),-0.5,1,1,0)),-0.5,1,1,0);
newpsi=psi*exp(newpsi);
logratio=5*log(newpsi)-5*log(psi)+2*(psi-newpsi);
for(int i = 0; i<1; ++i){
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-newpsi*tempv5;
tempv9=tempv5-tempv5[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
qvec=qvec1;
qmat(i,0)=1;
tempsum=0;
tempbr=0;
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
pmat(i,j)=tempv9[j]/(1-tempsum);
qvec[j]=1-pmat(i,j);
if(j>0){
qmat(i,j)=prod(qvec.subvec(0,j-1));
}
tempj=j;
marwmatn(i,loc-1)=qmat(i,j)*(1/(tempj+1.0)-(1-pmat(i,j))/(tempj+2.0));
tempsum=tempsum+tempv9[j];
if(tempsum==1 || qmat(i,j)<minp){
tempsum=j+1;
tempbr=1;
break;
}
}
if(tempbr>0){
for(int pw = tempsum; pw<Q; ++pw){  
loc=sortindex(i,pw);
marwmatn(i,loc-1)=0;
}
}
logratio=logratio+log(marwmatn(i,Z[i]-1))-log(marwmat(i,Z[i]-1));
}
if(log(R::runif(0,1))<logratio){
psi=newpsi;
marwmat=marwmatn;
}
psiT[iter]=psi;
}
predmat=predmat/reptime;
zmat=zmat/reptime1;


hT=hT/reptime;
hT=exp(hT/2.0);
for(int q = 0; q<Q; ++q){
tempv=predmat.col(q);
temp=sum(tempv.subvec(1,999))+sum(tempv.subvec(0,998));
temp=temp/2.0;
temp=temp*50.0/999;
tempv=tempv/temp;
predmat.col(q)=tempv;
}


Nb=floor(Nb/10.0);
psi=mean(psiT.subvec(Nb,Nm-1));
for(int i = 0; i<n; ++i){
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-psi*tempv5;
tempv9=tempv5-tempv5[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
qvec=qvec1;
qmat(i,0)=1;
tempsum=0;
tempbr=0;
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
pmat(i,j)=tempv9[j]/(1-tempsum);
qvec[j]=1-pmat(i,j);
if(j>0){
qmat(i,j)=prod(qvec.subvec(0,j-1));
}
temp=j;
odwmat(i,loc-1)=qmat(i,j)*(1/(temp+1.0)-(1-pmat(i,j))/(temp+2.0));
tempsum=tempsum+tempv9[j];
if(tempsum==1 || qmat(i,j)<minp){
tempsum=j+1;
tempbr=1;
break;
}
}
if(tempbr>0){
for(int pw = tempsum; pw<Q; ++pw){  
loc=sortindex(i,pw);
odwmat(i,loc-1)=0;
}
}
}

temp1=mean(SnT.subvec(Nb,Nm-1));
temp2=mean(fiT.subvec(Nb,Nm-1));
temp3=0;

vec fiTsub(1);
fiTsub=fiT.subvec(Nb,Nm-1);
double tempsd=stddev(fiTsub);
double GN=Nm-Nb,mscore=0.0;
double mscore1=0.0,mscore2=0.0;

tempsd=tempsd*exp(-0.2*log(GN));
for(int j = 0; j<GN; ++j){ 
mscore=mscore+R::dnorm(fiTsub[j],temp2,tempsd,0);
}
mscore=mscore/GN;
mscore1=log(mscore);
nscore=nscore/2000.0;
mscore2=log(nscore);
nscore=mscore1+mscore2;


vec pred0(1),pred1(1),predv(1),odwv(1),odwv1(1),LPS(n),MLV(n);
LPS.fill(0);
MLV.fill(0);
double err=0,mlscore=0.0;
double step=50.0/999;
for(int i = 0; i<10000; ++i){
if(err>1){
break;
}
for(int t = 0; t<n; ++t){
if(err>1){
break;
}
if(t==0){
h[0]=R::rnorm(0,sqrt(temp1/(1-temp2*temp2)));
}
if(t>0){
h[t]=temp2*h[t-1]+R::rnorm(0,sqrt(temp1));
}
temp4=y[t]-h[t];
if(temp4<x[0]){
err=1000;
}
if(temp4>x[nx-1]){
err=1000;
}
if(err==0){
loc=floor((temp4-x[0])/step);
pred0=conv_to< colvec >::from(predmat.submat(loc,0,loc,Q-1));
pred1=conv_to< colvec >::from(predmat.submat(loc+1,0,loc+1,Q-1));
predv=(temp4-x[loc])/step*(pred1-pred0)+pred0;
odwv=conv_to< colvec >::from(odwmat.submat(t,0,t,Q-1));
odwv1=conv_to< colvec >::from(zmat.submat(t,0,t,Q-1));
LPS[t]=LPS[t]+sum(predv%odwv);
MLV[t]=MLV[t]+sum(predv%odwv1);
}
}
}
if(err<1){
LPS=LPS/10000;
MLV=MLV/10000;
LPS=LPS.subvec(10-p,n-1);
temp3=mean(log(LPS));
for(int t = 10-p; t<n; ++t){
if(datay[t]!=0){
mlscore=mlscore+log(MLV[t]);
}
}
mlscore=mlscore+nscore;
}
if(err>1){
temp3=datum::nan;
mlscore=datum::nan;
}
sigmasq=mean(SzT.subvec(Nb,Nm-1));
W=mean(WT.subvec(Nb,Nm-1));
mu0=mean(Mu0T.subvec(Nb,Nm-1));
M=mean(MT.subvec(Nb,Nm-1));

lb=R::pnorm(-0.995,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.995,0.95,sqrt(0.5),1,0);
mlscore=mlscore+log(R::dnorm(temp2,0.95,sqrt(0.5),0)/(ub-lb));
mlscore=mlscore+log(snmeaninv*snmeaninv*R::dgamma(snmeaninv,2.5,1/0.025,0));
return List::create(
Named("fi")=temp2,Named("Sn")=temp1,Named("mu0")=mu0,Named("Sz")=sigmasq,Named("M")=M,Named("W")=W,Named("psi")=psi,
Named("fiT")=fiT, Named("SnT")=SnT,Named("Mu0T")=Mu0T,Named("SzT")=SzT,Named("MT")=MT,Named("WT")=WT,Named("psiT")=psiT,
Named("LPS")=temp3,Named("MLL")=mlscore,Named("hT")=hT,Named("mumat")=mumat);
} 

// [[Rcpp::export]]
Rcpp::List cpp6(arma::vec y,int p,int Q, int wdp,double maxiter,double Nb){
vec datay=y;  
double alpha=0.05,Sz0=0.05;
double logc=-20;
double caup=0.003*0.003;
int n1=y.n_elem;
int n=n1-p;
mat Xmat(n,p);
for(int t = 0; t<n; ++t){
for (int lag = 0;  lag< p; ++lag){
Xmat(t,lag)=y[t+p-lag-1];
}
}
y=y.subvec(p,n1-1);
for(int t = 0; t<n; ++t){
y[t]=log(y[t]*y[t]+exp(logc));
} 
double tempz1=0,tempz2=0,Qmax=pow(3,p);
int tempz3=0;
vec tempv6(3);
tempv6[0]=-3;
tempv6[1]=0;
tempv6[2]=3;
mat candimat(Qmax,p);
for(int j=0; j<p; j++){
tempz1=pow(3,j);  
for(int t=0; t<Qmax; t++){
tempz2=t;
tempz3=floor(tempz2/tempz1);
candimat(t,j)=tempv6[tempz3%3];
}
}
vec tempv4(1),tempv5(1); 
double psi=3;
mat distance(n,Qmax),distance1(n,Qmax);
vec index=linspace<vec>(1,Qmax,Qmax),check(Qmax),findv1(1);
int loc=0;
for(int j = 0; j<Qmax; ++j){
for(int i = 0; i<n; ++i){
tempv4=conv_to< colvec >::from(Xmat.submat(i,0,i,p-1));
tempv5=conv_to< colvec >::from(candimat.submat(j,0,j,p-1)); 
distance1(i,j)=sqrt(sum((tempv4-tempv5)%(tempv4-tempv5)));
}
check[j]=mean(distance1.col(j));
}
if(Q>Qmax){
Q=Qmax;
}
if(Q<Qmax){
tempv4=sort(check);
tempv4=tempv4.subvec(0,Q-1);
for(int q = 0; q<Q; ++q){
findv1=index.elem(find(check==tempv4[q]));
loc=findv1[0]-1;
for(int i = 0; i<n; ++i){
distance(i,q)=distance1(i,loc);
}
}
}
if(Q==Qmax){
distance=distance1;
}
distance=distance.submat(0,0,n-1,Q-1);
double temp=0;
mat sortdistance(n,Q),sortindex(n,Q),sortdistance1(n,Q);
vec Z(n),countz(Q),S(n),stempv(1);
countz.fill(0);
index=linspace<vec>(1,Q,Q);
for(int i = 0; i<n; ++i){
tempv4=conv_to< colvec >::from(distance.submat(i,0,i,Q-1));
tempv4=tempv4+index/100000000;
stempv=sort(tempv4);
for(int j = 0; j<Q; ++j){
sortdistance(i,j)=stempv[j];
sortdistance1(i,j)=-log(sortdistance(i,j))/psi;
temp=sortdistance(i,j);
findv1=index.elem(find(tempv4==temp));
sortindex(i,j)=findv1[0];
}
Z[i]=sortindex(i,0);
countz[Z[i]-1]=countz[Z[i]-1]+1;
S[i]=R::rpois(6);
}

mat wmat(n,Q),wmatn(n,Q);
double fi=0.86,sigmaeta=0.08;
vec h(n);
for(int t = 0; t<n; ++t){
if(t==0){
h[0]=R::rnorm(0,1)*fi+R::rnorm(0,sqrt(sigmaeta));
}
if(t>0){
h[t]=h[t-1]*fi+R::rnorm(0,sqrt(sigmaeta)); 
}
}
S=10000*Z+S;  
int od=1,len=0;
vec tempa(1),Sod(n);
uvec findv2(1);
for (int i = 0; i < n; ++i) {
if(i==0){
if(S[i]==0){
Sod[i]=0;
}
if(S[i]>0){
Sod[i]=1;
od=od+1;
}
}
if(i>0){
if(S[i]==0){
Sod[i]=0;
}
if(S[i]>0){
tempa=S.subvec(0,i-1);
findv2=find(tempa==S[i]);
len=findv2.n_elem;
if(len==0){
Sod[i]=od;
od=od+1;
}
if(len>0){
tempv4=Sod.subvec(0,i-1);
findv1=tempv4.elem(find(tempa==S[i]));
Sod[i]=findv1[0];
}
}
}
}
S=Sod;
mat odwmat(n,Q);
for(int i = 0; i<n; ++i){
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
odwmat(i,loc-1)=wmat(i,j);
}
}
Sod=unique(S);
int k=Sod.n_elem;
vec mu(n),mustar(k+1);
for(int i = 1; i<k+1; ++i){
mustar[i]=R::rnorm(0.4,0.1);
}
mustar[0]=logc;
for (int i = 0; i<n; ++i){
mu[i]=mustar[S[i]];
}
double sigmasq=2.5,W=0.01,M=4,mu0=0;
sortdistance=sortdistance1;
h=ones<vec>(n);
temp=0;
countz=ones<vec>(1);
findv1=ones<vec>(1);
loc=0;
len=0;
k=0;
od=1;
tempv4=ones<vec>(1);
tempv5=ones<vec>(1);
tempa=ones<vec>(1);
Sod=ones<vec>(n);
index=linspace<vec>(1,Q,Q);
vec tempv1(n),tempv2(n); 
double temp1=0,temp2=0,temp3=0,temp4=0,temp5=0,temp6=0,V=Sz0,C0=sigmaeta/(1-fi*fi),
mean1=0,var1=0;
vec countzp(1),countzus(1),locv(1),sp0(1),sp1(1),sp(1),pos(1),ord(1),maxS(1),Zuniq(1);
vec q0v(Q);
q0v.fill(0);
int zt=0,st=0,locs=0,kz=0;
uvec findv(1);
vec index1=linspace<vec>(0,n-1,n);
mat basecountcluster(n,Q),basecluster(n,Q),baseclusterS0(n,Q),
baseclusterS(n,Q),baseconfigS(n,Q),basemu(n,Q);
int position=0,n2=0,lmu=0,addnew=0,nz=0;
double probsum=0,deflat=1,probsuc=0,newmu=1;
vec newS(n),tempS(1),location(1),count(1),tempmustar(1),prob(1),
trueS(1),lamz(1),countz1(1),counts(1);
vec S0v(n+1);
S0v.fill(0);
vec tempmean(1),tempvar(1);
vec zprob(1),sz(1),scount(1),zcluster(1);
mat Q0k(Q,n),zprobmat(1,1);
Q0k.fill(0);
double sm=1e-320;
double newpsi=3.1,logratio=0;
vec tempv(1);
double a0=1,d0=0.2,templogb=0,tempinvb=1;
vec mix(1),coef(1),mixgam(1);
vec hstar(1),mustar1(1);
vec tempv3(1);
double newfi=0,lb=0,ub=0;
int para=0;
vec x=linspace<vec>(-30,20,1000);
int nx=1000;
mat predmat(nx,Q);
predmat.fill(0);
double reptime=0;
vec tempv7(Q),tempv8(Q);
int draw=0;
double logL1=0,logL2=0;
mat KS(n,Q),SS(n,Q);
vec sumv(n);
vec psiv(Q);
psiv.fill(1.2);
for(int t = 0; t<n; ++t){
for(int q = 0; q<Q; ++q){
KS(t,q)=exp(-psiv[q]*distance(t,q));
}
sumv[t]=sum(KS.row(t));
}
vec sumpv(1);
int change=0;
vec psivt(Q);
psivt.fill(0);
double Nm=maxiter/10.0;
vec fiT(Nm),SnT(Nm),SzT(Nm),MT(Nm),WT(Nm),Mu0T(Nm),psiT(Nm);
Nb=Nb-1;
wmat.fill(1.0/n);
vec probv1(Q),probv2(Q);
double repsi=0;
mat mumat(Nm,n);
vec hT(n);
hT.fill(0);
mat zmat(n,Q);
zmat.fill(0.0);
double reptime1=0.0;
double nmax=maxiter+2000;
double fimean=0.0,snmeaninv=0.0,nscore=0.0;
double tempnb=0.0;



for(int iter = 0; iter < nmax; ++iter){ 
temp=0;
countz=q0v;
Zuniq=countz;
for(int i = 0; i<n; ++i){
zt=Z[i];
countz[zt-1]=countz[zt-1]+1;
}
for(int j = 0; j<Q; ++j){
if(countz[j]>0){
temp=temp+1;
Zuniq[j]=temp;
}
}
kz=temp;
countzp=q0v.subvec(0,kz-1);
countzus=countzp;
pos=countzp;
ord=countzp;
maxS=countzp;
basecountcluster.fill(0);
basecluster.fill(-1);
baseclusterS.fill(-1);
baseclusterS0.fill(-1);
baseconfigS.fill(0);
basemu.fill(logc);
for(int i = 0; i<n; ++i){
loc=Zuniq[Z[i]-1]-1;
st=S[i];
baseclusterS0(ord[loc],loc)=st;
baseconfigS(ord[loc],loc)=i+1;
if(st==0){
baseclusterS(ord[loc],loc)=0;
basecountcluster(0,loc)=basecountcluster(0,loc)+1;
basecluster(0,loc)=0;
}
if(st>0){
countzp[loc]=countzp[loc]+1;
if(ord[loc]==0){
baseclusterS(0,loc)=1;
maxS[loc]=maxS[loc]+1;
countzus[loc]=countzus[loc]+1;
pos[loc]=pos[loc]+1; 
basecluster(pos[loc],loc)=st;
basemu(pos[loc],loc)=mustar[st];
basecountcluster(pos[loc],loc)=basecountcluster(pos[loc],loc)+1;
}
if(ord[loc]>0){
sp0=baseclusterS0.col(loc);
sp0=sp0.subvec(0,ord[loc]-1);
findv=find(sp0==st);
len=findv.n_elem;
if(len==0){
countzus[loc]=countzus[loc]+1;
pos[loc]=pos[loc]+1;
basecluster(pos[loc],loc)=st;
basemu(pos[loc],loc)=mustar[st];
basecountcluster(pos[loc],loc)=basecountcluster(pos[loc],loc)+1;
maxS[loc]=maxS[loc]+1;
baseclusterS(ord[loc],loc)=maxS[loc];
}
if(len>0){
sp1=baseclusterS.col(loc);
sp1=sp1.subvec(0,ord[loc]-1);
findv1=sp1.elem(find(sp0==st));
baseclusterS(ord[loc],loc)=findv1[0];
sp=basecluster.col(loc);
locv=index1.elem(find(sp==st));
locs=locv[0];
basecountcluster(locs,loc)=basecountcluster(locs,loc)+1;
}
}
}
ord[loc]=ord[loc]+1;
}
countz1=ones<vec>(kz);
od=-1;
for(int i = 0; i < Q; ++i) {
if(countz[i]>0){
od=od+1;
countz1[od]=countz[i];
}
}
temp2=a0+k-kz;
mix=ones<vec>(kz+1);
coef=mix-1;
mix[0]=lgamma(temp2);
temp3=M+1;
temp4=d0;
temp4=temp4-log(R::rbeta(temp3,countz1[0]));
coef[0]=countz1[0];
coef[1]=1;
if(kz>1){
for(int i = 1; i<kz; ++i){ 
temp4=temp4-log(R::rbeta(temp3,countz1[i]));
mix[i]=lgamma(temp2+i);
tempv=coef.subvec(0,i);
coef.subvec(1,i+1)=tempv;
coef[0]=0;
coef.subvec(0,i)=coef.subvec(0,i)+tempv*countz1[i];
coef.subvec(0,i+1)=coef.subvec(0,i+1)/max(coef.subvec(0,i+1));
}
}
mix[kz]=lgamma(temp2+kz);
templogb=log(temp4);
mix[0]=mix[0]+log(coef[0])-temp2*templogb;
temp5=mix[0];
mixgam=mix; 
tempinvb=1/temp4;   
mixgam[0]=R::rgamma(temp2,tempinvb);
for(int i = 1; i<kz+1; ++i){ 
temp1=temp2+i;
mix[i]=mix[i]+log(coef[i])-temp1*templogb;
if(mix[i]>temp5){
temp5=mix[i];
}
mixgam[i]=R::rgamma(temp1,tempinvb);
}
mix=mix-temp5;
mix=exp(mix);
mix=mix/sum(mix);
for (int q = 0;  q< kz+1; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-mix[q-1];
}
probsuc=mix[q]/deflat;
if (R::runif(0,1)<probsuc) {
M=mixgam[q];
break;
}
}
temp=alpha*sigmasq;
if(S[0]>0){
V=temp;
}
temp2=fi*fi*C0+sigmaeta;
temp3=mu[0];
temp4=temp2+V;
temp5=temp2/temp4;
temp6=y[0]-temp3;
tempv1[0]=temp5*temp6;
tempv2[0]=temp2-temp5*temp5*temp4;
for(int i = 1; i < n; ++i){
V=Sz0;
if(S[i]>0) {
V=temp;
}
temp1=fi*tempv1[i-1];
temp2=fi*fi*tempv2[i-1]+sigmaeta;
temp3=temp1+mu[i];
temp4=temp2+V;
temp5=temp2/temp4;
temp6=y[i]-temp3;
tempv1[i]=temp1+temp5*temp6;
tempv2[i]=temp2-temp5*temp5*temp4;
}
h[n-1]=R::rnorm(tempv1[n-1],sqrt(tempv2[n-1]));
for(int i = 1; i < n; ++i){
var1=1/(fi*fi/sigmaeta+1/tempv2[n-i-1]);
mean1=var1*(fi*h[n-i]/sigmaeta+tempv1[n-i-1]/tempv2[n-i-1]);
h[n-i-1]=R::rnorm(mean1,sqrt(var1));
}
for(int l = 0; l<kz; ++l){
nz=countz1[l];
tempS=conv_to< colvec >::from(baseclusterS.submat(0,l,nz-1,l));
location=conv_to< colvec >::from(baseconfigS.submat(0,l,nz-1,l));
count=conv_to< colvec >::from(basecountcluster.submat(0,l,countzus[l],l));
tempmustar=conv_to< colvec >::from(basemu.submat(0,l,countzus[l],l));  
lmu=countzus[l]+1;
n2=countzp[l]; 
lamz=ones<vec>(nz);
for(int i = 0; i<nz; ++i){
addnew=0;
st=tempS[i];
position=location[i]-1;
count[st]=count[st]-1;
if(st>0){
n2=n2-1;
}
if (count[st]==0 && st>0){
count.shed_row(st);
tempmustar=tempmustar.elem(find(tempmustar!=tempmustar[st]));
tempS.elem(find(tempS>st))=tempS.elem(find(tempS>st))-1;
lmu=lmu-1;
}
prob=ones<vec>(lmu+1);
prob[0]=W*exp(-0.5*(y[position]-h[position]-logc)*(y[position]-h[position]-logc)/Sz0)/sqrt(Sz0);
probsum=prob[0];
if(lmu>1){
for (int j = 1;  j< lmu; ++j) {
prob[j]=(1-W)*count[j]/sqrt(alpha*sigmasq)/(n2-1+M)*
exp(-0.5*(y[position]-h[position]-tempmustar[j])*(y[position]-h[position]-tempmustar[j])/alpha/sigmasq);
probsum=probsum+prob[j];
}
}
prob[lmu]=(1-W)*M/sqrt(sigmasq)/(n2-1+M)*
exp(-(y[position]-h[position]-mu0)*(y[position]-h[position]-mu0)/sigmasq/2);  
probsum=probsum+prob[lmu];
prob=prob/probsum;
for (int j = 0;  j< lmu+1; ++j) {
if (j==0) {
deflat=1.0;
}
if (j!=0) {
deflat=deflat-prob[j-1];
}
probsuc=prob[j]/deflat;
if (R::runif(0,1)<probsuc) {
tempS[i]=j;
break;
}
}
if(tempS[i]>0){
n2=n2+1;
if(tempS[i]>lmu-1){
addnew=1;
}
}
if(addnew==1){
var1=1/(1 / (alpha * sigmasq) + 1 / ((1 - alpha) * sigmasq));
mean1=var1*((y[position]-h[position])/alpha/sigmasq+mu0/(1-alpha)/sigmasq);
newmu=R::rnorm(mean1,sqrt(var1));
tempv4=tempmustar.subvec(0,lmu-1);
tempmustar=ones<vec>(lmu+1);
tempmustar.subvec(0,lmu-1)=tempv4;
tempmustar[lmu]=newmu;
tempv4=count;
count=ones<vec>(lmu+1);
count.subvec(0,lmu-1)=tempv4;
count[lmu]=1;
lamz[i]=newmu;
lmu=lmu+1;
}
if(addnew==0){
st=tempS[i];
count[st]=count[st]+1;
lamz[i]=tempmustar[st];
}
}
od=0;
trueS=ones<vec>(nz);
trueS.fill(1);
position=location[0]-1;
if(tempS[0]==0){
newS[position]=0; 
}  
if(tempS[0]>0){
newS[position]=1+10000*l;
od=od+1;
trueS[0]=od;
}
if(nz>1){
for(int i = 1; i<nz; ++i){
position=location[i]-1;  
if(tempS[i]==0){
newS[position]=0;
}
if(tempS[i]>0){
tempv5=lamz.subvec(0,i);
findv1=trueS.elem(find(tempv5==lamz[i]));
len=findv1.n_elem;
if(len==1){
od=od+1;
trueS[i]=od;
}
if(len>1){
trueS[i]=findv1[0];
}
newS[position]=trueS[i]+10000*l;
}
}       
}
}
tempS=newS;
od=1;
len=0;
for (int i = 0; i < n; ++i) {
if(i==0){
if(tempS[i]==0){
Sod[i]=0;
}
if(tempS[i]>0){
Sod[i]=1;
od=od+1;
}
}
if(i>0){
if(tempS[i]==0){
Sod[i]=0;
}
if(tempS[i]>0){
tempa=tempS.subvec(0,i-1);
findv=find(tempa==tempS[i]);
len=findv.n_elem;
if(len==0){
Sod[i]=od;
od=od+1;
}
if(len>0){
tempv5=Sod.subvec(0,i-1);
findv1=tempv5.elem(find(tempa==tempS[i]));
Sod[i]=findv1[0];
}
}
}
}
S=Sod;
counts=S0v.subvec(0,od-1);  
for (int i = 0; i < n; ++i) {
counts[S[i]]=counts[S[i]]+1;
}
n2=n-counts[0];
k=od-1;
if(iter<maxiter){
if(iter>Nb && iter%10==9){
hT=hT+h;
tempv7.fill(0);
tempv8.fill(0);
od=-1;
for (int q = 0;  q< Q; ++q){
if(countz[q]>0){
od=od+1;
tempv7[q]=countzus[od];
tempv8[q]=countzp[od];
}
}
for(int i = 0; i<nx; ++i){
for (int q = 0;  q< Q; ++q){
temp6=W*R::dnorm(x[i],logc,sqrt(Sz0),0);
if(tempv7[q]>0){ 
for (int l = 1;  l< tempv7[q]+1; ++l){
temp6=temp6+(1-W)*basecountcluster(l,q)/(tempv8[q]+M)*
R::dnorm(x[i],basemu(l,q),sqrt(alpha*sigmasq),0);
}
}
temp6=temp6+(1-W)*M/(tempv7[q]+M)*R::dnorm(x[i],mu0,sqrt(sigmasq),0);
predmat(i,q)=predmat(i,q)+temp6;      
}          
}
reptime=reptime+1;
}
}
if(Q>1){
if(k>0.5){ 
scount=counts.subvec(1,k); 
sz=S0v.subvec(0,k-1);
zcluster=sz;
zprobmat=Q0k.submat(0,0,Q-1,k-1);
tempv4=q0v;
for(int i = 0; i<n; ++i){
if(S[i]>0){
sz[S[i]-1]=Z[i];
zprob=conv_to< colvec >::from(wmat.submat(i,0,i,Q-1));
zprob=zprob+sm;
zprobmat.col(S[i]-1)=zprobmat.col(S[i]-1)+log(zprob);
tempv4[Z[i]-1]=tempv4[Z[i]-1]+1;
}
}
for(int j = 0; j<k; ++j){
tempv4[sz[j]-1]=tempv4[sz[j]-1]-scount[j];
for(int q = 0; q<Q; ++q){
zprobmat(q,j)=zprobmat(q,j)+lgamma(M+tempv4[q])-lgamma(M+tempv4[q]+scount[j]);
}
zprob=conv_to< colvec >::from(zprobmat.submat(0,j,Q-1,j));
zprob=zprob-max(zprob);
zprob=exp(zprob);
zprob=zprob/sum(zprob);
for (int q = 0;  q< Q; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-zprob[q-1];
}
probsuc=zprob[q]/deflat;
if (R::runif(0,1)<probsuc) {
zcluster[j]=q+1;
break;}
}
tempv4[zcluster[j]-1]=tempv4[zcluster[j]-1]+scount[j];
}
} 
for(int i = 0; i<n; ++i){
if(S[i]==0){
zprob=conv_to< colvec >::from(wmat.submat(i,0,i,Q-1));
for (int q = 0;  q< Q; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-zprob[q-1];
}
probsuc=zprob[q]/deflat;
if (R::runif(0,1)<probsuc) {
Z[i]=q+1;
break;
}
}
}
if(S[i]>0){
Z[i]=zcluster[S[i]-1];
}
}



if(wdp<4){
logL1=0;
logL2=0;
logratio=0;
newpsi=psi*exp(R::rnorm(0,0.1));
for(int t = 0; t<n; ++t){
for(int q = 0; q<Q; ++q){
probv1[q]=-psi*distance(t,q);
probv2[q]=-newpsi*distance(t,q);
}
probv1=probv1-max(probv1);
probv1=exp(probv1);
probv1=probv1/sum(probv1);
wmat.row(t)=conv_to< rowvec >::from(probv1);
probv2=probv2-max(probv2);
probv2=exp(probv2);
probv2=probv2/sum(probv2);
wmatn.row(t)=conv_to< rowvec >::from(probv2);
logL1=logL1+log(wmat(t,Z[t]-1));
logL2=logL2+log(wmatn(t,Z[t]-1));
}
logL1=logL1+R::dnorm(log(psi),1.03,0.02,1);
logL2=logL2+R::dnorm(log(newpsi),1.03,0.02,1); 
logratio=exp(logL2-logL1);
if(logratio>R::runif(0,1)){
psi=newpsi;
wmat=wmatn;
}
}
if(wdp>3){
sumpv=sumv;
logL1=0;
logL2=0;
logratio=0;
change=0;
for(int j = 0; j<Q; ++j){
if(wdp<6){
newpsi=R::qnorm(R::runif(R::pnorm(0,psiv[j],0.01,1,0),R::pnorm(4,psiv[j],0.01,1,0)),psiv[j],0.01,1,0);
}
if(wdp==6){
newpsi=R::qnorm(R::runif(R::pnorm(0,psiv[j],0.01,1,0),R::pnorm(2,psiv[j],0.01,1,0)),psiv[j],0.01,1,0);
}
logL1=0;
logL2=0;
logratio=0;
for(int t = 0; t<n; ++t){  
SS(t,j)=exp(-newpsi*distance(t,j));
sumv[t]=sumpv[t]+SS(t,j)-KS(t,j);
logL1=logL1+log(KS(t,Z[t]-1)/sumpv[t]);
if(j==Z[t]-1){
logL2=logL2+log(SS(t,Z[t]-1)/sumv[t]);
}
if(j!=Z[t]-1){
logL2=logL2+log(KS(t,Z[t]-1)/sumv[t]);
}
}
if(wdp<6){
logL1=logL1-log(R::pnorm(4,psiv[j],0.01,1,0)-R::pnorm(0,psiv[j],0.01,1,0));
logL2=logL2-log(R::pnorm(4,newpsi,0.01,1,0)-R::pnorm(0,newpsi,0.01,1,0));
}
if(wdp==6){
logL1=logL1-log(R::pnorm(2,psiv[j],0.01,1,0)-R::pnorm(0,psiv[j],0.01,1,0));
logL2=logL2-log(R::pnorm(2,newpsi,0.01,1,0)-R::pnorm(0,newpsi,0.01,1,0));
}
if(wdp==4){
logL1=logL1-psiv[j]+log(psiv[j]);
logL2=logL2-newpsi+log(newpsi);
}
if(wdp==5){
logL1=logL1-0.413*psiv[j];
logL2=logL2-0.413*newpsi;
}
if(wdp==6){
logL1=logL1-log(caup+psiv[j]*psiv[j]);
logL2=logL2-log(caup+newpsi*newpsi);
}
logratio=exp(logL2-logL1);
if(logratio>R::runif(0,1)){
psiv[j]=newpsi;
KS.col(j)=SS.col(j);
sumpv=sumv;
change=1;
}
if (change==0){
SS.col(j)=KS.col(j);
sumv=sumpv;
}
}
for(int t = 0; t<n; ++t){
for(int q = 0; q<Q; ++q){
wmat(t,q)=exp(-psiv[q]*distance(t,q));
}
wmat.row(t)=wmat.row(t)/sum(wmat.row(t));
}
}
countz=q0v;
for(int i = 0; i<n; ++i){
countz[Z[i]-1]=countz[Z[i]-1]+1;
}
}

if(iter<maxiter){
if(iter%10==2 && iter>Nb){
for(int i = 0; i < n; ++i) {
zmat(i,Z[i]-1)=zmat(i,Z[i]-1)+1;
}
reptime1=reptime1+1;
}
}

mustar=S0v.subvec(0,k);
mustar[0]=logc;
if(k>0.5){
tempmean=S0v.subvec(0,k-1);
tempvar=tempmean;
for(int i = 0; i < n; ++i) {
st=S[i];
if (st>0){
tempmean[st-1]=(y[i]-h[i])*(1-alpha)/(alpha+counts[st]*(1-alpha))+tempmean[st-1];
}
}
for(int j = 0; j < k; ++j) {
tempmean[j]=tempmean[j]+mu0*alpha/(alpha+counts[j+1]*(1-alpha));
tempvar[j]=sigmasq/(counts[j+1]/alpha+1/(1-alpha)); 
mustar[j+1]=R::rnorm(tempmean[j],sqrt(tempvar[j]));
}
}
for(int i = 0; i < n; ++i) {
mu[i]=mustar[S[i]];
if(iter%10==0 && iter<maxiter){
mumat(draw,i)=mu[i];
}
}
para=0;
if(R::runif(0,1)<0.5){
para=1;
}
if(para==1){
hstar=h+mu0;
mustar1=mustar.subvec(1,k)-mu0;
var1=sigmaeta/((n-1)*(1-fi)*(1-fi)+(1-fi*fi));
mean1=var1*((1-fi*fi)/sigmaeta*hstar[0]+(1-fi)/sigmaeta*sum(hstar.subvec(1,n-1)-fi*hstar.subvec(0,n-2)));
temp6=mu0;
mu0=R::rnorm(mean1,sqrt(var1));
temp2=n2+k;
temp2=temp2/2+1;
temp3=3;
for(int i = 0; i<n; ++i){
if(S[i]>0){
temp3=temp3+(y[i]-h[i]-mu[i])*(y[i]-h[i]-mu[i])/2/alpha;
}
}
temp3=temp3+sum(mustar1%mustar1)/2/(1-alpha);
sigmasq=R::rgamma(temp2,1/temp3);
sigmasq=1/sigmasq;
logratio=0;
tempv3=hstar.subvec(1,n-1);
tempv3=tempv3-mu0;
tempv4=hstar.subvec(0,n-2);
tempv4=tempv4-mu0;
mean1=sum(tempv3%tempv4)/sum(tempv4%tempv4);
var1=sigmaeta/sum(tempv4%tempv4);
if(var1<0.0001){
var1=0.0001;
}
lb=R::pnorm(-0.995,mean1,sqrt(var1),1,0);
ub=R::pnorm(0.995,mean1,sqrt(var1),1,0);
newfi=R::qnorm(R::runif(lb,ub),mean1,sqrt(var1),1,0);
logratio=log(R::dnorm(fi,mean1,sqrt(var1),0)/(ub-lb))-log(R::dnorm(newfi,mean1,sqrt(var1),0)/(ub-lb));
lb=R::pnorm(-0.995,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.995,0.95,sqrt(0.5),1,0);
logratio=logratio+log(R::dnorm(newfi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio-(hstar[0]-mu0)*(hstar[0]-mu0)*(1-newfi*newfi)/2/sigmaeta+0.5*log(1-newfi*newfi);
logratio=logratio-log(R::dnorm(fi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio+(hstar[0]-mu0)*(hstar[0]-mu0)*(1-fi*fi)/2/sigmaeta-0.5*log(1-fi*fi);
if(log(R::runif(0,1))<logratio){
fi=newfi;
}
if(iter>maxiter-1){
fi=fimean;
}

temp2=n;
temp2=temp2/2+2.5;
temp3=0.025;
for(int i = 1; i<n; ++i){
temp3=temp3+0.5*((hstar[i]-mu0)-(hstar[i-1]-mu0)*fi)*((hstar[i]-mu0)-(hstar[i-1]-mu0)*fi);
}
sigmaeta=R::rgamma(temp2,1/temp3);
sigmaeta=1/sigmaeta;

if(iter>maxiter-1){
nscore=nscore+snmeaninv*snmeaninv*R::dgamma(snmeaninv,temp2,1/temp3,0);
}


h=hstar-mu0;
mustar.subvec(1,k)=mustar1+mu0;
for(int i = 0; i<n; ++i){
if(S[i]>0){
mu[i]=mu[i]-temp6+mu0;
}
}
}
if(para==0){
mean1=(sum(mustar)-mustar[0])/k;
var1=(1-alpha)*sigmasq/k;
mu0=R::rnorm(mean1,sqrt(var1));
temp2=n2+k;
temp2=temp2/2+1;
temp3=3;
for(int i = 0; i<n; ++i){
if(S[i]>0){
temp3=temp3+(y[i]-h[i]-mu[i])*(y[i]-h[i]-mu[i])/2/alpha;
}
}
temp3=temp3+sum((mustar.subvec(1,k)-mu0)%(mustar.subvec(1,k)-mu0))/2/(1-alpha);
sigmasq=R::rgamma(temp2,1/temp3);
sigmasq=1/sigmasq;
logratio=0;
tempv3=h.subvec(1,n-1);
tempv4=h.subvec(0,n-2);
mean1=sum(tempv3%tempv4)/sum(tempv4%tempv4);
var1=sigmaeta/sum(tempv4%tempv4);
if(var1<0.0001){
var1=0.0001;
}
lb=R::pnorm(-0.995,mean1,sqrt(var1),1,0);
ub=R::pnorm(0.995,mean1,sqrt(var1),1,0);
newfi=R::qnorm(R::runif(lb,ub),mean1,sqrt(var1),1,0);
logratio=log(R::dnorm(fi,mean1,sqrt(var1),0)/(ub-lb))-log(R::dnorm(newfi,mean1,sqrt(var1),0)/(ub-lb));
lb=R::pnorm(-0.995,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.995,0.95,sqrt(0.5),1,0);
logratio=logratio+log(R::dnorm(newfi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio-(hstar[0]-mu0)*(hstar[0]-mu0)*(1-newfi*newfi)/2/sigmaeta+0.5*log(1-newfi*newfi);
logratio=logratio-log(R::dnorm(fi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio+(hstar[0]-mu0)*(hstar[0]-mu0)*(1-fi*fi)/2/sigmaeta-0.5*log(1-fi*fi);
if(log(R::runif(0,1))<logratio){
fi=newfi;
}


if(iter>maxiter-1){
fi=fimean;
}

temp2=n;
temp2=temp2/2+2.5;
temp3=0.025;
for(int i = 1; i<n; ++i){
temp3=temp3+0.5*(h[i]-h[i-1]*fi)*(h[i]-h[i-1]*fi);
}
sigmaeta=R::rgamma(temp2,1/temp3);
sigmaeta=1/sigmaeta;
}

if(iter>maxiter-1){
nscore=nscore+snmeaninv*snmeaninv*R::dgamma(snmeaninv,temp2,1/temp3,0);
}

temp2=n-n2;
temp2=temp2+0.1;
temp3=n2;
temp3=temp3+0.9;
W=R::rbeta(temp2,temp3);
if(iter%10==0 && iter<maxiter){
fiT[draw]=fi;
SnT[draw]=sigmaeta;
SzT[draw]=sigmasq;
MT[draw]=M;
WT[draw]=W;
Mu0T[draw]=mu0;
psiT[draw]=psi;
if(wdp>3 && iter>Nb){
psivt=psivt+psiv;
repsi=repsi+1;
}
draw=draw+1;
}
if(iter<maxiter && iter>maxiter-2){
tempnb=floor(Nb/10.0);
snmeaninv=mean(SnT.subvec(tempnb,Nm-1));
fimean=mean(fiT.subvec(tempnb,Nm-1));
snmeaninv=1/snmeaninv;
} 
}

predmat=predmat/reptime;
zmat=zmat/reptime1;
hT=hT/reptime;
hT=exp(hT/2.0);
for(int q = 0; q<Q; ++q){
tempv=predmat.col(q);
temp=sum(tempv.subvec(1,999))+sum(tempv.subvec(0,998));
temp=temp/2.0;
temp=temp*50.0/999;
tempv=tempv/temp;
predmat.col(q)=tempv;
}
Nb=floor(Nb/10.0);
if(wdp<4){
psi=mean(psiT.subvec(Nb,Nm-1));
psivt.fill(psi);
}
if(wdp>3){
psivt=psivt/repsi;
}
for(int t = 0; t<n; ++t){
for(int q = 0; q<Q; ++q){
odwmat(t,q)=exp(-psiv[q]*distance(t,q));
}
odwmat.row(t)=odwmat.row(t)/sum(odwmat.row(t));
}


temp1=mean(SnT.subvec(Nb,Nm-1));
temp2=mean(fiT.subvec(Nb,Nm-1));
temp3=0;

vec fiTsub(1);
fiTsub=fiT.subvec(Nb,Nm-1);
double tempsd=stddev(fiTsub);
double GN=Nm-Nb,mscore=0.0;
double mscore1=0.0,mscore2=0.0;

tempsd=tempsd*exp(-0.2*log(GN));
for(int j = 0; j<GN; ++j){ 
mscore=mscore+R::dnorm(fiTsub[j],temp2,tempsd,0);
}
mscore=mscore/GN;
mscore1=log(mscore);
nscore=nscore/2000.0;
mscore2=log(nscore);
nscore=mscore1+mscore2;

vec pred0(1),pred1(1),predv(1),odwv(1),odwv1(1),LPS(n),MLV(n);
LPS.fill(0);
MLV.fill(0);
double err=0,mlscore=0.0;
double step=50.0/999;
for(int i = 0; i<10000; ++i){
if(err>1){
break;
}
for(int t = 0; t<n; ++t){
if(err>1){
break;
}
if(t==0){
h[0]=R::rnorm(0,sqrt(temp1/(1-temp2*temp2)));
}
if(t>0){
h[t]=temp2*h[t-1]+R::rnorm(0,sqrt(temp1));
}
temp4=y[t]-h[t];
if(temp4<x[0]){
err=1000;
}
if(temp4>x[nx-1]){
err=1000;
}
if(err==0){
loc=floor((temp4-x[0])/step);
pred0=conv_to< colvec >::from(predmat.submat(loc,0,loc,Q-1));
pred1=conv_to< colvec >::from(predmat.submat(loc+1,0,loc+1,Q-1));
predv=(temp4-x[loc])/step*(pred1-pred0)+pred0;
odwv=conv_to< colvec >::from(odwmat.submat(t,0,t,Q-1));
odwv1=conv_to< colvec >::from(zmat.submat(t,0,t,Q-1));
LPS[t]=LPS[t]+sum(predv%odwv);
MLV[t]=MLV[t]+sum(predv%odwv1);
}
}
}
if(err<1){
LPS=LPS/10000;
MLV=MLV/10000;
LPS=LPS.subvec(10-p,n-1);
temp3=mean(log(LPS));
for(int t = 10-p; t<n; ++t){
if(datay[t]!=0){
mlscore=mlscore+log(MLV[t]);
}
}
mlscore=mlscore+nscore;
}
if(err>1){
temp3=datum::nan;
mlscore=datum::nan;
}
sigmasq=mean(SzT.subvec(Nb,Nm-1));
W=mean(WT.subvec(Nb,Nm-1));
mu0=mean(Mu0T.subvec(Nb,Nm-1));
M=mean(MT.subvec(Nb,Nm-1));

lb=R::pnorm(-0.995,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.995,0.95,sqrt(0.5),1,0);
mlscore=mlscore+log(R::dnorm(temp2,0.95,sqrt(0.5),0)/(ub-lb));
mlscore=mlscore+log(snmeaninv*snmeaninv*R::dgamma(snmeaninv,2.5,1/0.025,0));
return List::create(
Named("fi")=temp2,Named("Sn")=temp1,Named("mu0")=mu0,Named("Sz")=sigmasq,Named("M")=M,Named("W")=W,Named("psi")=psi,
Named("fiT")=fiT, Named("SnT")=SnT,Named("Mu0T")=Mu0T,Named("SzT")=SzT,Named("MT")=MT,Named("WT")=WT,Named("psiT")=psiT,
Named("LPS")=temp3,Named("MLL")=mlscore,Named("hT")=hT,Named("mumat")=mumat);
} 

// [[Rcpp::export]]
Rcpp::List cpp7(arma::vec Macf, arma::vec Szacf, arma::vec Snacf, arma::vec Wacf, 
arma::vec psiacf, arma::vec fiacf, arma::vec Mu0acf, double len, double tsize){
double temp=2/sqrt(tsize);
double a1=0,a2=0,a3=0,a4=0,a5=0,a6=0,a7=0,a8=0;
double b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,b7=0;
for(int i = 0; i<1001; ++i){
if(b1<0.5){
if(Snacf[i]<temp){
b1=1;
}
if(b1<0.5){
a1=a1+Snacf[i];
}
}
if(b2<0.5){
if(fiacf[i]<temp){
b2=1;
}
if(b2<0.5){
a2=a2+fiacf[i];
}
}
if(b3<0.5){
if(Mu0acf[i]<temp){
b3=1;
}
if(b3<0.5){
a3=a3+Mu0acf[i];
}
}
}
if(len>9.5){
for(int i = 0; i<1001; ++i){
if(b4<0.5){
if(Szacf[i]<temp){
b4=1;
}
if(b4<0.5){
a4=a4+Szacf[i];
}
}
if(b5<0.5){
if(Wacf[i]<temp){
b5=1;
}
if(b5<0.5){
a5=a5+Wacf[i];
}
}
if(b6<0.5){
if(Macf[i]<temp){
b6=1;
}
if(b6<0.5){
a6=a6+Macf[i];
}
}
if(b7<0.5){
if(psiacf[i]<temp){
b7=1;
}
if(b7<0.5){
a7=a7+psiacf[i];
}
}
}
}
a1=2*a1-1;
a2=2*a2-1;
a3=2*a3-1;
a4=2*a4-1;
a5=2*a5-1;
a6=2*a6-1;
a7=2*a7-1;
if(len>9.5){
a8=datum::nan;
}
if(len<9.5){
a6=datum::nan;
a4=datum::nan;
a7=datum::nan;
a5=datum::nan;
a8=a3;
a3=datum::nan;
}
return List::create(Named("M")=a6,Named("Sz")=a4,Named("fi")=a2,Named("Sn")=a1,
Named("W")=a5,Named("Mu0")=a3,Named("psi")=a7,Named("Mu")=a8);
}

// [[Rcpp::export]]
Rcpp::List fitgarch(arma::vec y,int p,int Q,double method,double maxiter){ 
if(method<1.5){
p=3;  
Q=1;
}
int n1=y.n_elem;
int n=n1-p;
mat Xmat(n,p);
for(int t = 0; t<n; ++t){
for (int lag = 0;  lag< p; ++lag){
Xmat(t,lag)=y[t+p-lag-1];
}
}
double y0=y[p-1];
y=y.subvec(p,n1-1);
double tempz1=0,tempz2=0,Qmax=pow(3,p);
int tempz3=0;
vec tempv6(3);
tempv6[0]=-1.0;
tempv6[1]=0.0;
tempv6[2]=1.0;
mat candimat(Qmax,p),candimat1(Qmax,p);
for(int j=0; j<p; j++){
tempz1=pow(3,j);  
for(int t=0; t<Qmax; t++){
tempz2=t;
tempz3=floor(tempz2/tempz1);
candimat(t,j)=tempv6[tempz3%3];
}
}
vec tempv4(1),tempv5(1); 
mat distance(n,Qmax),distance1(n,Qmax),distance2(n,Qmax);
vec index=linspace<vec>(1,Qmax,Qmax),check(Qmax),findv1(1);
int loc=0;
for(int j = 0; j<Qmax; ++j){
for(int i = 0; i<n; ++i){
tempv4=conv_to< colvec >::from(Xmat.submat(i,0,i,p-1));
tempv5=conv_to< colvec >::from(candimat.submat(j,0,j,p-1)); 
distance1(i,j)=sqrt(sum((tempv4-tempv5)%(tempv4-tempv5)))+0.00001;
distance(i,j)=log(distance1(i,j));
}
check[j]=mean(distance1.col(j));
}
distance2=distance;
candimat1=candimat;
if(Q>Qmax){
Q=Qmax;
}

if(Q<Qmax){
tempv4=sort(check);
tempv4=tempv4.subvec(0,Q-1);
for(int q = 0; q<Q; ++q){
findv1=index.elem(find(check==tempv4[q]));
loc=findv1[0]-1;
candimat.row(q)=candimat1.row(loc);
for(int i = 0; i<n; ++i){
distance(i,q)=distance2(i,loc);
}
}
}
distance=distance.submat(0,0,n-1,Q-1);
candimat=candimat.submat(0,0,Q-1,p-1);
double temp=0;
mat sortdistance(n,Q),sortindex(n,Q);
vec Z(n),countz(Q),S(n),stempv(1);
countz.fill(0);
index=linspace<vec>(1,Q,Q);
for(int i = 0; i<n; ++i){
tempv4=conv_to< colvec >::from(distance.submat(i,0,i,Q-1));
tempv4=tempv4+index/100000000;
stempv=sort(tempv4);
for(int j = 0; j<Q; ++j){
sortdistance(i,j)=stempv[j];
temp=sortdistance(i,j);
findv1=index.elem(find(tempv4==temp));
sortindex(i,j)=findv1[0];
}
Z[i]=sortindex(i,0);
countz[Z[i]-1]=countz[Z[i]-1]+1;
S[i]=R::rpois(6);
}
double psi=2,tempj=1.0,newpsi=1,caup=0.003*0.003,am0=10,dm0=3;
vec tempv9(1);
mat qmat(n,Q),pmat(n,Q),wmat(n,Q),wmatn(n,Q),odwmat(n,Q),odwmatn(n,Q),marwmat(n,Q),marwmatn(n,Q);
odwmatn.fill(0);
vec qvec(Q),qvec1(Q);
qvec.fill(1);
qvec1.fill(1);
vec tempq(Q+1);
double tempsum=0,tempbr=0,logL1=0,logL2=0;
vec sumv(n),psiv(Q),probv1(Q),probv2(Q);
psiv.fill(0.2);
mat KS(n,Q),SS(n,Q);
for(int t = 0; t<n; ++t){
for(int q = 0; q<Q; ++q){
KS(t,q)=exp(-psiv[q]*distance(t,q));
}
sumv[t]=sum(KS.row(t));
}
vec sumpv(1);
int change=0;
if(method<2.5 && method>1.5){
for(int iter = 0; iter<1; ++iter){
for(int i = 0; i<n; ++i){  
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-psi*tempv5;
tempv9=tempv5-tempv5[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
tempsum=0;
tempbr=0;
for(int j = 0; j<Q; ++j){
tempj=j;
loc=sortindex(i,j);
pmat(i,j)=tempv9[j]/(1-tempsum);
if(tempj==0){
tempq[0]=1;
}
tempq[j+1]=tempq[j]*(1-pmat(i,j))*(tempj+1)/(tempj+2);
odwmat(i,loc-1)=tempq[j]-tempq[j+1];
tempsum=tempsum+tempv9[j];
if(tempsum==1){
tempsum=j+1;
tempbr=1;
break;
}
}
if(tempbr>0){
for(int pw = tempsum; pw<Q; ++pw){  
loc=sortindex(i,pw);
odwmat(i,loc-1)=0;
}
}
}
odwmatn=odwmatn+odwmat;
}
}
if(method>2.5){
for(int t = 0; t<n; ++t){
for(int q = 0; q<Q; ++q){
probv1[q]=exp(-psiv[q]*distance(t,q));
}
probv1=probv1-max(probv1);
probv1=exp(probv1);
probv1=probv1/sum(probv1);
wmat.row(t)=conv_to< rowvec >::from(probv1);
}
odwmat=wmat;
wmatn=wmat; 
}
vec h(n),hnew(n);
double a0=0.01,a1=0.2,beta=0.5,h0=1,tempy=0.0;
double wg=0,rh1=0,rh2=0,tau=0,wgp=0,rh1p=0,rh2p=0,taup=0;
double stpw=0.015,stp1=0.15,stp2=0.15,stpt=0.3;
a0=wgp*wgp;
a1=exp(rh1p)/(1+exp(rh1p))/(1+exp(rh2p));
beta=exp(rh1p+rh2p)/(1+exp(rh1p))/(1+exp(rh2p));
h0=exp(taup);
for(int t = 0; t<n; ++t){
if(t==0){
h[t]=a0+h0*beta+a1*y0*y0;
}
if(t>0){
h[t]=a0+h[t-1]*beta+a1*y[t-1]*y[t-1];
}
}
S=10000*Z+S;  
int od=1,len=0;
vec tempa(1),Sod(n);
uvec findv2(1);
for (int i = 0; i < n; ++i) {
if(i==0){
if(S[i]==0){
Sod[i]=0;
}
if(S[i]>0){
Sod[i]=1;
od=od+1;
}
}
if(i>0){
if(S[i]==0){
Sod[i]=0;
}
if(S[i]>0){
tempa=S.subvec(0,i-1);
findv2=find(tempa==S[i]);
len=findv2.n_elem;
if(len==0){
Sod[i]=od;
od=od+1;
}
if(len>0){
tempv4=Sod.subvec(0,i-1);
findv1=tempv4.elem(find(tempa==S[i]));
Sod[i]=findv1[0];
}
}
}
}
S=Sod;
Sod=unique(S);
int k=Sod.n_elem;
vec mustar(k+1),lam(n),lamp(n);
lam.fill(1);
for(int i = 1; i<k+1; ++i){
mustar[i]=R::rnorm(0,0.1);
}
mustar[0]=-20;
for (int i = 0; i<n; ++i){
lamp[i]=mustar[S[i]];
}
double M=4.0;
temp=0;
countz=ones<vec>(1);
findv1=ones<vec>(1);
loc=0;
len=0;
k=0;
od=1;
tempv4=ones<vec>(1);
tempv5=ones<vec>(1);
tempa=ones<vec>(1);
Sod=ones<vec>(n);
index=linspace<vec>(1,Q,Q);
vec tempv1(n),tempv2(n); 
double temp1=0,temp2=0,temp3=0,temp4=0,temp5=0;
vec countzp(1),countzus(1),locv(1),sp0(1),sp1(1),sp(1),pos(1),ord(1),maxS(1),Zuniq(1);
vec q0v(Q);
q0v.fill(0);
int zt=0,st=0,locs=0,kz=0;
uvec findv(1);
vec index1=linspace<vec>(0,n-1,n);
mat basecountcluster(n,Q),basecluster(n,Q),baseclusterS0(n,Q),
baseclusterS(n,Q),baseconfigS(n,Q),basemu(n,Q);
int position=0,n2=0,lmu=0,addnew=0,nz=0;
double probsum=0,deflat=1,probsuc=0,newmu=1;
vec newS(n),tempS(1),location(1),count(1),tempmustar(1),prob(1),
trueS(1),lamz(1),countz1(1),counts(1);
vec S0v(n+1);
S0v.fill(0);
vec tempmean(1),tempvar(1);
vec zprob(1),sz(1),scount(1),zcluster(1);
mat Q0k(Q,n),zprobmat(1,1);
Q0k.fill(0);
double sm=1e-320;
double logratio=0;
vec tempv(1),gamv(n);
double templogb=0,tempinvb=1;
vec mix(1),coef(1),mixgam(1);
vec hstar(1),mustar1(1);
vec tempv3(1);
vec tempv7(Q),tempv8(Q);
if(method>0.5){
if(maxiter<15000){
cout << "#of iterations are set to be 15000. Use maxiter between 15000 and 40000." << endl;
maxiter=15000;
}
if(maxiter>40000){
cout << "#of iterations are set to be 40000. Use maxiter between 15000 and 40000." << endl;
maxiter=40000;
}
}
double Nm=50000+maxiter+1000;
vec A0T(50000),A1T(50000),mu0T(50000),BT(50000),HT(50000),MT(Nm),psiT(Nm),SzT(Nm),aaa(Nm);
double reptime4=0;
mat zmat(n,Q);
zmat.fill(0.0);
if(method>3.5){
psiT=q0v;
}
double zord=0.0;
vec x=linspace<vec>(-40,40,4000);
double temp6=0;
mat predmat(4000,Q);
predmat.fill(0);
double mu0=-0.1,score2=0,score3=0;
vec scorepv(3000),scorev(5000),scorev1(5000);
mat lammat(n,3000),mumat(n,5000);
double reptime1=0,reptime2=0,reptime3=0,mu0mean=0;
vec lamv(1),tauv(3000),rh1v(3000),rh2v(3000),wgv(3000),mu0v(3000);
double MLL=0,MLR=0;
vec h1t(n);
h1t.fill(0);
vec posv1(n),posv2(n),ratv1(n),ratv2(n),posv1n(n),posv2n(n),ratv1n(n),ratv2n(n),predv(n);
predv.fill(0);
vec sdy(n),sdyn(n);
double sigmasq=1,alpha=0.05,mean1=0,var1=1;
vec nill(1);
nill.fill(0);
mat nill1(1,1);
nill1.fill(0);
if(method<0.5){
Nm=50000;
}
double vfix=7.0;
double temp4n=0,temp5n=0,tempy1=0;
vec predv1(n);
predv1.fill(0);
double wgp1=0,rh1p1=0,rh2p1=0,taup1=0;
for(int iter = 0; iter < Nm; ++iter){
if(iter>49999){  
temp=0;
countz=q0v;
Zuniq=countz;
for(int i = 0; i<n; ++i){
zt=Z[i];
countz[zt-1]=countz[zt-1]+1;
}
for(int j = 0; j<Q; ++j){
if(countz[j]>0){
temp=temp+1;
Zuniq[j]=temp;
}
}
kz=temp;
countzp=q0v.subvec(0,kz-1);
countzus=countzp;
pos=countzp;
ord=countzp;
maxS=countzp;
basecountcluster.fill(0);
basecluster.fill(-1);
baseclusterS.fill(-1);
baseclusterS0.fill(-1);
baseconfigS.fill(0);
basemu.fill(-20);
for(int i = 0; i<n; ++i){
loc=Zuniq[Z[i]-1]-1;
st=S[i];
baseclusterS0(ord[loc],loc)=st;
baseconfigS(ord[loc],loc)=i+1;
if(st==0){
baseclusterS(ord[loc],loc)=0;
basecountcluster(0,loc)=basecountcluster(0,loc)+1;
basecluster(0,loc)=0;
}
if(st>0){
countzp[loc]=countzp[loc]+1;
if(ord[loc]==0){
baseclusterS(0,loc)=1;
maxS[loc]=maxS[loc]+1;
countzus[loc]=countzus[loc]+1;
pos[loc]=pos[loc]+1; 
basecluster(pos[loc],loc)=st;
basemu(pos[loc],loc)=mustar[st];
basecountcluster(pos[loc],loc)=basecountcluster(pos[loc],loc)+1;
}
if(ord[loc]>0){
sp0=baseclusterS0.col(loc);
sp0=sp0.subvec(0,ord[loc]-1);
findv=find(sp0==st);
len=findv.n_elem;
if(len==0){
countzus[loc]=countzus[loc]+1;
pos[loc]=pos[loc]+1;
basecluster(pos[loc],loc)=st;
basemu(pos[loc],loc)=mustar[st];
basecountcluster(pos[loc],loc)=basecountcluster(pos[loc],loc)+1;
maxS[loc]=maxS[loc]+1;
baseclusterS(ord[loc],loc)=maxS[loc];
}
if(len>0){
sp1=baseclusterS.col(loc);
sp1=sp1.subvec(0,ord[loc]-1);
findv1=sp1.elem(find(sp0==st));
baseclusterS(ord[loc],loc)=findv1[0];
sp=basecluster.col(loc);
locv=index1.elem(find(sp==st));
locs=locv[0];
basecountcluster(locs,loc)=basecountcluster(locs,loc)+1;
}
}
}
ord[loc]=ord[loc]+1;
}
countz1=ones<vec>(kz);
od=-1;
for(int i = 0; i < Q; ++i) {
if(countz[i]>0){
od=od+1;
countz1[od]=countz[i];
}
}
}
if(iter>49999 && iter<60000){
if(iter%10==0){
tempv7.fill(0);
tempv8.fill(0);
od=-1;
for (int q = 0;  q< Q; ++q){
if(countz[q]>0){
od=od+1;
tempv7[q]=countzus[od];
tempv8[q]=countzp[od];
}
}
for(int i = 0; i<4000; ++i){
for (int q = 0;  q< Q; ++q){
temp6=0;
if(tempv7[q]>0){ 
for (int l = 1;  l< tempv7[q]+1; ++l){
temp6=temp6+basecountcluster(l,q)/(tempv8[q]+M)*R::dnorm(x[i],basemu(l,q),sqrt(alpha*sigmasq),0);
}
}
temp6=temp6+M/(tempv7[q]+M)*R::dnorm(x[i],mu0,sqrt(sigmasq),0);
predmat(i,q)=predmat(i,q)+temp6;      
}          
}
}
}
if(iter>59999 && iter<60001){
predmat=predmat/1000.0;
for(int q = 0; q<Q; ++q){
tempv=predmat.col(q);
temp=sum(tempv.subvec(1,3999))+sum(tempv.subvec(0,3998));
temp=temp/2.0;
temp=temp*80.0/3999.0;
tempv=tempv/temp;
predmat.col(q)=tempv;
}
}
if(iter>59999 && iter<65000){
for(int t = 0; t<n; ++t){ 
temp4=posv1[t];
temp4n=posv1n[t];
temp5=posv2[t];
temp5n=posv2n[t];
if(temp4<(-0.1)){
predv[t]=predv[t]+predmat(0,Z[t]-1);
predv[t]=predv[t]+predmat(3999,Z[t]-1);
if(sdy[t]<(-39.999)){
predv1[t]=predv1[t]+predmat(0,Z[t]-1);
}
if(sdy[t]>39.999){
predv1[t]=predv1[t]+predmat(3999,Z[t]-1);
}
}
if(temp4>-0.1){
predv[t]=predv[t]+ratv1[t]*predmat(temp4,Z[t]-1)+ratv2[t]*predmat(temp5,Z[t]-1);
predv1[t]=predv1[t]+ratv1[t]*predmat(temp4,Z[t]-1)+ratv2[t]*predmat(temp5,Z[t]-1);
predv[t]=predv[t]+ratv1n[t]*predmat(temp4n,Z[t]-1)+ratv2n[t]*predmat(temp5n,Z[t]-1);
}
}
}
if(iter>49999){
temp2=am0+k-kz;
mix=ones<vec>(kz+1);
coef=mix-1;
mix[0]=lgamma(temp2);
temp3=M+1;
temp4=dm0;
temp4=temp4-log(R::rbeta(temp3,countz1[0]));
coef[0]=countz1[0];
coef[1]=1;
if(kz>1){
for(int i = 1; i<kz; ++i){ 
temp4=temp4-log(R::rbeta(temp3,countz1[i]));
mix[i]=lgamma(temp2+i);
tempv=coef.subvec(0,i);
coef.subvec(1,i+1)=tempv;
coef[0]=0;
coef.subvec(0,i)=coef.subvec(0,i)+tempv*countz1[i];
coef.subvec(0,i+1)=coef.subvec(0,i+1)/max(coef.subvec(0,i+1));
}
}
mix[kz]=lgamma(temp2+kz);
templogb=log(temp4);
mix[0]=mix[0]+log(coef[0])-temp2*templogb;
temp5=mix[0];
mixgam=mix; 
tempinvb=1/temp4;   
mixgam[0]=R::rgamma(temp2,tempinvb);
for(int i = 1; i<kz+1; ++i){ 
temp1=temp2+i;
mix[i]=mix[i]+log(coef[i])-temp1*templogb;
if(mix[i]>temp5){
temp5=mix[i];
}
mixgam[i]=R::rgamma(temp1,tempinvb);
}
mix=mix-temp5;
mix=exp(mix);
mix=mix/sum(mix);
for (int q = 0;  q< kz+1; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-mix[q-1];
}
probsuc=mix[q]/deflat;
if (R::runif(0,1)<probsuc) {
M=mixgam[q];
break;
}
}
}
if(iter>49999){
for(int l = 0; l<kz; ++l){
nz=countz1[l];
tempS=conv_to< colvec >::from(baseclusterS.submat(0,l,nz-1,l));
location=conv_to< colvec >::from(baseconfigS.submat(0,l,nz-1,l));
count=conv_to< colvec >::from(basecountcluster.submat(0,l,countzus[l],l));
tempmustar=conv_to< colvec >::from(basemu.submat(0,l,countzus[l],l));  
lmu=countzus[l]+1;
n2=countzp[l]; 
lamz=ones<vec>(nz);
for(int i = 0; i<nz; ++i){
addnew=0;
st=tempS[i];
position=location[i]-1;
tempy=sdy[position];

count[st]=count[st]-1;
if(st>0){
n2=n2-1;
}
if (count[st]==0 && st>0){
count.shed_row(st);
tempmustar=tempmustar.elem(find(tempmustar!=tempmustar[st]));
tempS.elem(find(tempS>st))=tempS.elem(find(tempS>st))-1;
lmu=lmu-1;
}
prob=ones<vec>(lmu+1);
prob[0]=0;
probsum=prob[0];
if(lmu>1){
for (int j = 1;  j< lmu; ++j) {
prob[j]=count[j]/(n2-1+M)*R::dnorm(tempy,tempmustar[j],sqrt(alpha*sigmasq),0);
probsum=probsum+prob[j];  
}
}
prob[lmu]=M/(n2-1+M)*R::dnorm(tempy,mu0,sqrt(sigmasq),0); 
probsum=probsum+prob[lmu];
prob=prob/probsum;
for (int j = 0;  j< lmu+1; ++j) {
if (j==0) {
deflat=1.0;
}
if (j!=0) {
deflat=deflat-prob[j-1];
}
probsuc=prob[j]/deflat;
if (R::runif(0,1)<probsuc) {
tempS[i]=j;
break;
}
}
if(tempS[i]>0){
n2=n2+1;
if(tempS[i]>lmu-1){
addnew=1;
}
}
if(addnew==1){
var1=1/(1 / (alpha * sigmasq) + 1 / ((1 - alpha) * sigmasq));
mean1=var1*((tempy)/alpha/sigmasq+mu0/(1-alpha)/sigmasq);
newmu=R::rnorm(mean1,sqrt(var1));
tempv4=tempmustar.subvec(0,lmu-1);
tempmustar=ones<vec>(lmu+1);
tempmustar.subvec(0,lmu-1)=tempv4;
tempmustar[lmu]=newmu;
tempv4=count;
count=ones<vec>(lmu+1);
count.subvec(0,lmu-1)=tempv4;
count[lmu]=1;
lamz[i]=newmu;
lmu=lmu+1;
}
if(addnew==0){
st=tempS[i];
count[st]=count[st]+1;
lamz[i]=tempmustar[st];
}
}
od=0;
trueS=ones<vec>(nz);
trueS.fill(1);
position=location[0]-1;
if(tempS[0]==0){
newS[position]=0; 
}  
if(tempS[0]>0){
newS[position]=1+10000*l;
od=od+1;
trueS[0]=od;
}
if(nz>1){
for(int i = 1; i<nz; ++i){
position=location[i]-1;  
if(tempS[i]==0){
newS[position]=0;
}
if(tempS[i]>0){
tempv5=lamz.subvec(0,i);
findv1=trueS.elem(find(tempv5==lamz[i]));
len=findv1.n_elem;
if(len==1){
od=od+1;
trueS[i]=od;
}
if(len>1){
trueS[i]=findv1[0];
}
newS[position]=trueS[i]+10000*l;
}
}       
}
}
tempS=newS;
od=1;
len=0;
for (int i = 0; i < n; ++i) {
if(i==0){
if(tempS[i]==0){
Sod[i]=0;
}
if(tempS[i]>0){
Sod[i]=1;
od=od+1;
}
}
if(i>0){
if(tempS[i]==0){
Sod[i]=0;
}
if(tempS[i]>0){
tempa=tempS.subvec(0,i-1);
findv=find(tempa==tempS[i]);
len=findv.n_elem;
if(len==0){
Sod[i]=od;
od=od+1;
}
if(len>0){
tempv5=Sod.subvec(0,i-1);
findv1=tempv5.elem(find(tempa==tempS[i]));
Sod[i]=findv1[0];
}
}
}
}
S=Sod;
counts=S0v.subvec(0,od-1);  
for (int i = 0; i < n; ++i) {
counts[S[i]]=counts[S[i]]+1;
}
n2=n-counts[0];
k=od-1;
}
if(iter>49999 && Q>1){
if(k>0.5){ 
scount=counts.subvec(1,k); 
sz=S0v.subvec(0,k-1);
zcluster=sz;
zprobmat=Q0k.submat(0,0,Q-1,k-1);
tempv4=q0v;
for(int i = 0; i<n; ++i){
if(S[i]>0){
sz[S[i]-1]=Z[i];
zprob=conv_to< colvec >::from(odwmat.submat(i,0,i,Q-1));
zprob=zprob+sm;
zprobmat.col(S[i]-1)=zprobmat.col(S[i]-1)+log(zprob);
tempv4[Z[i]-1]=tempv4[Z[i]-1]+1;
}
}
for(int j = 0; j<k; ++j){
tempv4[sz[j]-1]=tempv4[sz[j]-1]-scount[j];
for(int q = 0; q<Q; ++q){
zprobmat(q,j)=zprobmat(q,j)+lgamma(M+tempv4[q])-lgamma(M+tempv4[q]+scount[j]);
}
zprob=conv_to< colvec >::from(zprobmat.submat(0,j,Q-1,j));
zprob=zprob-max(zprob);
zprob=exp(zprob);
zprob=zprob/sum(zprob);
for (int q = 0;  q< Q; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-zprob[q-1];
}
probsuc=zprob[q]/deflat;
if (R::runif(0,1)<probsuc) {
zcluster[j]=q+1;
break;}
}
tempv4[zcluster[j]-1]=tempv4[zcluster[j]-1]+scount[j];
}
} 
for(int i = 0; i<n; ++i){
if(S[i]==0){
zprob=conv_to< colvec >::from(odwmat.submat(i,0,i,Q-1));
for (int q = 0;  q< Q; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-zprob[q-1];
}
probsuc=zprob[q]/deflat;
if (R::runif(0,1)<probsuc) {
Z[i]=q+1;
break;
}
}
}
if(S[i]>0){
Z[i]=zcluster[S[i]-1];
}
if(iter>54999){
zmat(i,Z[i]-1)=zmat(i,Z[i]-1)+1;
}
}
countz=q0v;
for(int i = 0; i<n; ++i){
countz[Z[i]-1]=countz[Z[i]-1]+1;
}
}
if(iter>49999 && method>3.5){
sumpv=sumv;
logL1=0;
logL2=0;
logratio=0;
change=0;
for(int j = 0; j<Q; ++j){
newpsi=R::qnorm(R::runif(R::pnorm(0,psi,1,1,0),R::pnorm(4,psi,1,1,0)),psi,1,1,0);
logL1=0;
logL2=0;
logratio=0;
for(int t = 0; t<n; ++t){  
SS(t,j)=exp(-newpsi*distance(t,j));
sumv[t]=sumpv[t]+SS(t,j)-KS(t,j);
logL1=logL1+log(KS(t,Z[t]-1)/sumpv[t]);
if(j==Z[t]-1){
logL2=logL2+log(SS(t,Z[t]-1)/sumv[t]);
}
if(j!=Z[t]-1){
logL2=logL2+log(KS(t,Z[t]-1)/sumv[t]);
}
}
logL1=logL1-log(R::pnorm(4,psiv[j],1,1,0)-R::pnorm(0,psiv[j],1,1,0));
logL2=logL2-log(R::pnorm(4,newpsi,1,1,0)-R::pnorm(0,newpsi,1,1,0));
if(method<4.5){
logL1=logL1-psiv[j]+log(psiv[j]);
logL2=logL2-newpsi+log(newpsi);
}
if(method>4.5 && method<5.5){
logL1=logL1-0.413*psiv[j];
logL2=logL2-0.413*newpsi;
}
if(method>5.5 && method<6.5){
logL1=logL1-log(caup+psiv[j]*psiv[j]);
logL2=logL2-log(caup+newpsi*newpsi);
}
logratio=exp(logL2-logL1);
if(logratio>R::runif(0,1)){
psiv[j]=newpsi;
KS.col(j)=SS.col(j);
sumpv=sumv;
change=1;
}
if (change==0){
SS.col(j)=KS.col(j);
sumv=sumpv;
}
}
if(iter>54999){
psiT=psiT+psiv;
zord=zord+1;
}
for(int t = 0; t<n; ++t){
for(int q = 0; q<Q; ++q){
probv1[q]=exp(-psiv[q]*distance(t,q));
}
probv1=probv1-max(probv1);
probv1=exp(probv1);
probv1=probv1/sum(probv1);
odwmat.row(t)=conv_to< rowvec >::from(probv1);
}
}
if(iter<50000){
for(int i = 0; i < n; ++i) {
temp1=(y[i]/sqrt(h[i])-mu0)*(y[i]/sqrt(h[i])-mu0)/2.0+vfix/2.0;
temp1=1.0/temp1;
lam[i]= R::rgamma((vfix+1)/2.0,temp1);
}
}
if(iter>49999){
mustar=S0v.subvec(0,k);
mustar[0]=-20;
if(k>0.5){
tempmean=S0v.subvec(0,k-1);
tempvar=tempmean;
for(int i = 0; i < n; ++i) {
st=S[i];
if (st>0){
tempmean[st-1]=(sdy[i])*(1-alpha)/(alpha+counts[st]*(1-alpha))+tempmean[st-1];
}
}
for(int j = 0; j < k; ++j) {
tempmean[j]=tempmean[j]+mu0*alpha/(alpha+counts[j+1]*(1-alpha));
tempvar[j]=sigmasq/(counts[j+1]/alpha+1/(1-alpha)); 
mustar[j+1]=R::rnorm(tempmean[j],sqrt(tempvar[j]));
}
}
for(int i = 0; i < n; ++i) {
lam[i]=mustar[S[i]];
}
if(iter<55000){
mumat.col(reptime3)=lam;
reptime3=reptime3+1;
}  
}
if(iter<50000){
wg=R::rnorm(wgp,stpw);
rh1=R::rnorm(rh1p,stp1);
rh2=R::rnorm(rh2p,stp2);
tau=R::rnorm(taup,stpt);
a0=wg*wg;
a1=exp(rh1)/(1+exp(rh1))/(1+exp(rh2));
beta=exp(rh1+rh2)/(1+exp(rh1))/(1+exp(rh2));
h0=exp(tau);
logratio=0;
score2=0;
score3=0;
for(int t = 0; t<n; ++t){
if(t==0){
hnew[t]=a0+h0*beta+a1*y0*y0;
}
if(t>0){
hnew[t]=a0+hnew[t-1]*beta+a1*y[t-1]*y[t-1];
}
temp3=(y[t]/sqrt(hnew[t])-mu0)*(y[t]/sqrt(hnew[t])-mu0)/2*lam[t]+0.5*log(hnew[t]);
score2=score2+temp3;
temp3=(y[t]/sqrt(h[t])-mu0)*(y[t]/sqrt(h[t])-mu0)/2*lam[t]+0.5*log(h[t]);
score3=score3+temp3;
}
score2=score2+wg*wg/8+(rh1-2)*(rh1-2)/8+(rh2-2)*(rh2-2)/8+(tau+2)*(tau+2)/8;
score3=score3+wgp*wgp/8+(rh1p-2)*(rh1p-2)/8+(rh2p-2)*(rh2p-2)/8+(taup+2)*(taup+2)/8;
logratio=score3-score2;
if(iter<45000){
if(log(R::runif(0,1))<logratio){
wgp=wg;
rh1p=rh1;
rh2p=rh2;
taup=tau;
h=hnew;
score3=score2;
}
a0=wgp*wgp;
a1=exp(rh1p)/(1+exp(rh1p))/(1+exp(rh2p));
beta=exp(rh1p+rh2p)/(1+exp(rh1p))/(1+exp(rh2p));
h0=exp(taup);
A0T[iter]=a0;
A1T[iter]=a1;
BT[iter]=beta;
HT[iter]=h0;
mu0T[iter]=mu0;
if(iter>29999 && iter%5==0){
h1t=h1t+h; 
lammat.col(reptime1)=lam;
scorepv[reptime1]=score3;
wgv[reptime1]=wgp;
rh1v[reptime1]=rh1p;
rh2v[reptime1]=rh2p;
tauv[reptime1]=taup;
mu0v[reptime1]=mu0;
reptime1=reptime1+1;
}
if(iter>44998){
wgp=mean(wgv);
rh1p=mean(rh1v);
rh2p=mean(rh2v);
taup=mean(tauv);
mu0mean=mean(mu0v);
wgp1=wgp;
rh1p1=rh1p;
rh2p1=rh2p1;
taup1=taup;
a0=wgp*wgp;
a1=exp(rh1p)/(1+exp(rh1p))/(1+exp(rh2p));
beta=exp(rh1p+rh2p)/(1+exp(rh1p))/(1+exp(rh2p));
h0=exp(taup);
h1t=h1t/3000.0;
for(int t = 0; t<n; ++t){ 
sdy[t]=y[t]/sqrt(h1t[t])-mu0mean;
sdyn[t]=-sdy[t];
posv1[t]=-1.0;
posv2[t]=-1.0;
ratv1[t]=-1.0;
ratv2[t]=-1.0;
posv1n[t]=-1.0;
posv2n[t]=-1.0;
ratv1n[t]=-1.0;
ratv2n[t]=-1.0;
if(sdy[t]>-40 && sdy[t]<40){
temp4=sdy[t];  
posv1[t]=floor((temp4-x[0])*3999/80.0)+0.0;
posv2[t]=posv1[t]+1.0;
temp5=posv1[t];
temp5=(temp4-x[temp5])*3999/80.0;
ratv2[t]=temp5;
ratv1[t]=1-temp5;
}
if(sdyn[t]>-40 && sdyn[t]<40){
temp4n=sdyn[t];  
posv1n[t]=floor((temp4n-x[0])*3999/80.0)+0.0;
posv2n[t]=posv1n[t]+1.0;
temp5n=posv1n[t];
temp5n=(temp4n-x[temp5n])*3999/80.0;
ratv2n[t]=temp5n;
ratv1n[t]=1-temp5n;
}
}
for(int t = 0; t<n; ++t){
if(t==0){
hnew[t]=a0+h0*beta+a1*y0*y0;
}
if(t>0){
hnew[t]=a0+hnew[t-1]*beta+a1*y[t-1]*y[t-1];
}
}
h=hnew;
for(int j = 0; j<3000; ++j){
lamv=lammat.col(j);
logratio=0;
for(int t = 0; t<n; ++t){
logratio=logratio-(y[t]/sqrt(hnew[t])-mu0v[j])*(y[t]/sqrt(hnew[t])-mu0v[j])/2*lam[t];
logratio=logratio-0.5*log(hnew[t]);
}
logratio=logratio-wgp*wgp/8-(rh1p-2)*(rh1p-2)/8-(rh2p-2)*(rh2p-2)/8-(taup+2)*(taup+2)/8;
logratio=logratio+R::dnorm(wgp,wgv[j],stpw,1)+R::dnorm(rh1p,rh1v[j],stp1,1)+
R::dnorm(rh2p,rh2v[j],stp2,1)+R::dnorm(taup,tauv[j],stpt,1);
scorepv[j]=scorepv[j]+logratio;
if(scorepv[j]>1){
scorepv[j]=1;
}
}
}
}
if(iter>44999){
scorev[reptime2]=logratio;
if(scorev[reptime2]>1){
scorev[reptime2]=1;
}
scorev1[reptime2]=R::dnorm(mu0mean,temp2,temp1,1);
reptime2=reptime2+1;
}
temp1=1.0;
temp2=0.1;
for(int t = 0; t<n; ++t){
temp1=temp1+lam[t];
temp2=temp2+y[t]/sqrt(h[t])*lam[t];
}
temp1=1.0/temp1;
temp2=temp2*temp1;
temp1=sqrt(temp1);
mu0=R::qnorm(R::runif(0,R::pnorm(0,temp2,temp1,1,0)),temp2,temp1,1,0);
if(iter>49998){
mu0=0;
lam=lamp;
}
}
if(iter>49999){
mean1=(sum(mustar)-mustar[0])/k;
var1=(1-alpha)*sigmasq/k;
mu0=R::rnorm(mean1,sqrt(var1));
aaa[reptime4]=mu0;
reptime4=reptime4+1;
temp2=n2+k;
temp2=temp2/2.0+1.0;
temp3=3;
for(int i = 0; i<n; ++i){
if(S[i]>0){
temp3=temp3+(sdy[i]-lam[i])*(sdy[i]-lam[i])/2/alpha;
}
}
temp3=temp3+sum((mustar.subvec(1,k)-mu0)%(mustar.subvec(1,k)-mu0))/2/(1-alpha);
sigmasq=R::rgamma(temp2,1/temp3);
sigmasq=1/sigmasq;    
}
if(iter>49999){
MT[iter]=M;
SzT[iter]=sigmasq; 
}
}
MLL=log(mean(exp(scorepv)))-log(mean(exp(scorev)))+log(mean(exp(scorev1)));
MLR=MLL;
predv=predv/5000.0;
predv1=predv1/5000.0;
for(int t = 0; t<n; ++t){
temp1=0;  
for(int q = 0; q<Q; ++q){
if(zmat(t,q)>temp1){
temp1=zmat(t,q);  
Z[t]=q+1;
}
}
}
temp1=0;
if(method>3.5){
psiT=psiT/zord;  
psiv=psiT;
}
double minp=1e-50;
if(method<3.5 && method>1.5){
psiT=psiT.subvec(0,9999);
if(method<2.5){
for(int iter = 0; iter < 10000; ++iter){  
newpsi=R::qnorm(R::runif(0,R::pnorm(log(5.0/psi),-0.5,1,1,0)),-0.5,1,1,0);
newpsi=psi*exp(newpsi);
logratio=5*log(newpsi)-5*log(psi)+2*(psi-newpsi);
for(int i = 0; i<1; ++i){
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-newpsi*tempv5;
tempv9=tempv5-tempv5[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
qvec=qvec1;
qmat(i,0)=1;
tempsum=0;
tempbr=0;
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
pmat(i,j)=tempv9[j]/(1-tempsum);
qvec[j]=1-pmat(i,j);
if(j>0){
qmat(i,j)=prod(qvec.subvec(0,j-1));
}
tempj=j;
marwmatn(i,loc-1)=qmat(i,j)*(1/(tempj+1.0)-(1-pmat(i,j))/(tempj+2.0));
tempsum=tempsum+tempv9[j];
if(tempsum==1 || qmat(i,j)<minp){
tempsum=j+1;
tempbr=1;
break;
}
}
if(tempbr>0){
for(int pw = tempsum; pw<Q; ++pw){  
loc=sortindex(i,pw);
marwmatn(i,loc-1)=0;
}
}
logratio=logratio+log(marwmatn(i,Z[i]-1))-log(marwmat(i,Z[i]-1));
}
if(log(R::runif(0,1))<logratio){
psi=newpsi;
marwmat=marwmatn;
}
psiT[iter]=psi;
}
}
if(method>2.5){
for(int iter = 0; iter < 10000; ++iter){  
logL1=0;
logL2=0;
logratio=0;
newpsi=psi*exp(R::rnorm(0,0.02));
for(int t = 0; t<n; ++t){
for(int q = 0; q<Q; ++q){
probv1[q]=-psi*distance(t,q);
probv2[q]=-newpsi*distance(t,q);
}
probv1=probv1-max(probv1);
probv1=exp(probv1);
probv1=probv1/sum(probv1);
wmat.row(t)=conv_to< rowvec >::from(probv1);
probv2=probv2-max(probv2);
probv2=exp(probv2);
probv2=probv2/sum(probv2);
wmatn.row(t)=conv_to< rowvec >::from(probv2);
logL1=logL1+log(wmat(t,Z[t]-1));
logL2=logL2+log(wmatn(t,Z[t]-1));
}
logL1=logL1+R::dnorm(log(psi),1.03,0.02,1);
logL2=logL2+R::dnorm(log(newpsi),1.03,0.02,1); 
logratio=exp(logL2-logL1);
if(logratio>R::runif(0,1)){
psi=newpsi;
wmat=wmatn;
}
psiT[iter]=psi;
}
}
}
reptime2=0;
int rpn=floor(maxiter/10.0);
for(int j = 0; j < 50000; ++j){ 
if(j%10==0){
A0T[reptime2]=A0T[j];
A1T[reptime2]=A1T[j];
BT[reptime2]=BT[j];
HT[reptime2]=HT[j];
mu0T[reptime2]=mu0T[j];
reptime2=reptime2+1;
}
}
A0T=A0T.subvec(499,498+rpn);  
A1T=A1T.subvec(499,498+rpn); 
BT=BT.subvec(499,498+rpn); 
HT=HT.subvec(499,498+rpn); 
mu0T=mu0T.subvec(499,498+rpn); 
reptime1=0;
for(int iter = 0; iter<maxiter; ++iter){
if(iter%10==0){
aaa[reptime1]=aaa[iter];  
reptime1=reptime1+1;
}    
}  
aaa=aaa.subvec(0,rpn-1);
if(method<0.5){
MT=nill;
psiT=nill;
SzT=nill;
predmat=nill1;
mumat=nill1;
}
if(method>0.5 && method<1.5){
psiT=nill;
}
if(method>0.5 && method<6.5){
if(method>1.5 && method<3.5){
psiT=psiT.subvec(5000,5000+rpn-1);
}
reptime2=0;
for(int iter = 50000; iter < Nm; ++iter){ 
if(iter%10==0){
MT[reptime2]=MT[iter];
SzT[reptime2]=SzT[iter];
reptime2=reptime2+1;
}
}
MT=MT.subvec(99,99+rpn-1);
SzT=SzT.subvec(99,99+rpn-1);
}
if(method>0.5){
mumat=mumat.submat(0,0,n-1,rpn-1);  
}  
a0=mean(A0T);
a1=mean(A1T);
beta=mean(BT);
mu0=mean(mu0T);
h0=mean(HT);
vec h2t(n);
for(int t = 0; t<n; ++t){
if(t==0){
h2t[t]=a0+h0*beta+a1*y0*y0;
}
if(t>0){
h2t[t]=a0+h2t[t-1]*beta+a1*y[t-1]*y[t-1];
}
if(t>9-p){
if(method>0.5){
MLR=MLR+log(predv1[t])-log(sqrt(h2t[t]));  
}  
if(method<0.5){
MLR=MLR+R::dt(sdy[t],vfix,1)-log(sqrt(h2t[t]));
}  
if(y[t]!=0){
tempy1=y[t];
if(y[t]<0){
tempy1=-tempy1;
}
if(method>0.5){
MLL=MLL+log(predv[t])-log(sqrt(h2t[t]))-log(2*tempy1/(tempy1*tempy1+exp(-20)));
}
if(method<0.5){
MLL=MLL+log(2)+R::dt(sdy[t],vfix,1)-log(sqrt(h2t[t]))-log(2*tempy1/(tempy1*tempy1+exp(-20)));
}
}
}
}

MLL=MLL+R::dnorm(wgp1,0,2,1)+R::dnorm(rh1p1,2,2,1)+R::dnorm(rh2p1,2,2,1)+R::dnorm(taup1,-2,2,1);

if(method>1.5 && method<2.5){
psi=mean(psiT);
}
if(method>2.5 && method<3.5){
psiv.fill(psi);
}

mat cdfmat(4000,Q);
cdfmat.fill(0);
if(method>0.5){
for(int q = 0; q < Q; ++q){ 
for(int j = 1; j < 4000; ++j){ 
cdfmat(j,q)=cdfmat(j-1,q)+(predmat(j,q)+predmat(j-1,q))*40/3999.0;
}
cdfmat(3999,q)=1.0;
}
}
vec ypred(p+10);
ypred.fill(0);
vec predictr(10000),distancev(Q),sortdistancev(Q),sortindexv(Q),pmatv(Q);
double temph=h1t[n-1],tempz=0;
double ub=3999,lb=0,mid=0,findx=0,err=0.0;
int loclb=0,locub=0,locmid=0;
vec dense(1);
for(int iter = 0; iter < 10000; ++iter){
for(int j = 0; j < p; ++j){
ypred[j]=y[n-p+j];  
} 
temph=h1t[n-1];
for(int j = 0; j < 10; ++j){  
temph=a0+a1*(ypred[p+j-1])*(ypred[p+j-1])+beta*temph;
tempv5=ypred.subvec(j,j+p-1);
for(int q = 0; q<Q; ++q){
tempv4=conv_to< colvec >::from(candimat.submat(q,0,q,p-1));
distancev[q]=sqrt(sum((tempv4-tempv5)%(tempv4-tempv5)))+0.00001;
distancev[q]=log(distancev[q]);
}
tempz=1;
if(method>1.5){
if(method>2.5){
for(int q = 0; q<Q; ++q){
probv1[q]=exp(-psiv[q]*distancev[q]);
}
probv1=probv1-max(probv1);
probv1=exp(probv1);
probv1=probv1/sum(probv1);
}
if(method<2.5){
index=linspace<vec>(1,Q,Q);
tempv4=distancev;
tempv4=tempv4+index/100000000;
stempv=sort(tempv4);
for(int q = 0; q<Q; ++q){
sortdistancev[q]=stempv[q];
temp=sortdistancev[q];
findv1=index.elem(find(tempv4==temp));
sortindexv[q]=findv1[0];
}
tempv5=sortdistancev;
tempv5=-psi*tempv5;
tempv9=tempv5-tempv5[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
tempsum=0;
tempbr=0;
for(int q = 0; q<Q; ++q){
tempj=q;
loc=sortindexv[q];
pmatv[q]=tempv9[q]/(1-tempsum);
if(tempj==0){
tempq[0]=1;
}
tempq[q+1]=tempq[q]*(1-pmatv[q])*(tempj+1)/(tempj+2);
probv1[loc-1]=tempq[q]-tempq[q+1];
tempsum=tempsum+tempv9[q];
if(tempsum==1){
tempsum=q+1;
tempbr=1;
break;
}
}
if(tempbr>0){
for(int pw = tempsum; pw<Q; ++pw){  
loc=sortindexv[pw];
probv1[loc-1]=0;
}
}
}
for (int q = 0;  q< Q; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-zprob[q-1];
}
probsuc=zprob[q]/deflat;
if (R::runif(0,1)<probsuc) {
tempz=q+1;
break;
}
}
}
if(method>0.5){
dense=cdfmat.col(tempz-1);
temp1=R::runif(0,1);  
ub=3999;
lb=0;
for (int tr = 0; tr < 12; ++tr){
findx=0;  
mid=(ub+lb)*0.5;
locmid=floor(mid);
if(dense[locmid]>temp1){
ub=mid;
}
if(dense[locmid]<temp1){
lb=mid;
}
if(dense[locmid]<temp1){
temp2=dense[locmid+1];
if(temp2>temp1){
findx=1;
loclb=locmid;
locub=locmid+1;
}
}
if(dense[locmid]>temp1){
temp2=dense[locmid-1];
if(temp2<temp1){
findx=1;
loclb=locmid-1;
locub=locmid;
}
}
if(findx>0.5){
break;      
}
}
err=x[loclb]+(temp1-dense[loclb])/(dense[locub]-dense[loclb])*80/3999.0;
}
if(method<0.5){
err=R::rt(vfix);
}
ypred[p+j]=(err+mu0)*sqrt(temph);
}
predictr[iter]=exp(sum(ypred.subvec(p,p+9))/100.0);
}
return List::create(
Named("A0T")=A0T,Named("A1T")=A1T,Named("BT")=BT,Named("HT")=HT,Named("SzT")=SzT,Named("MT")=MT,Named("psiT")=psiT,
Named("Mu0T")=aaa,
Named("MLL")=MLL,Named("MLR")=MLR, Named("h")=h1t,Named("lammat")=mumat,Named("predictr")=predictr);
}

// [[Rcpp::export]]
arma::vec  simdatag(double n, double p, double Q, double method){ 
vec output(1);
for(int iter = 0; iter<10; ++iter){
if(method<1.5){
Q=1;
p=1;
}
double Qmax=pow(3,p);
if(Q>Qmax){
Q=Qmax;
}
mat Candimat(Qmax,p);
double tempz1=0,tempz2=0;
int tempz3=0;
vec tempv6(3);
tempv6[0]=-3;
tempv6[1]=0;
tempv6[2]=3;
for(int j=0; j<p; j++){
tempz1=pow(3,j);  
for(int t=0; t<Qmax; t++){
tempz2=t;
tempz3=floor(tempz2/tempz1);
Candimat(t,j)=tempv6[tempz3%3];
}
}
Candimat=Candimat.submat(0,0,Q-1,p-1);
double a0=0.03,al1=0.1,beta=0.82,h0=0.03,v=7.0,M=6,mu0=0,alpha=0.05,Sz=1.5;
double deflat=1,temp=0,probsuc=0,v1=(1-alpha)*Sz,tempdist=0,check1=0,check2=0;
vec h(n); 
mat Xmat(n,p),probmat(150,Q),mumat(150,Q);   
for(int q = 0; q<Q; ++q){
deflat=1;
for(int l = 0; l<150; ++l){
if (l==0){
probmat(l,q)=R::rbeta(1,M);
mumat(l,q)=R::rnorm(0,sqrt(v1));
deflat=1-probmat(l,q);
}
if (l>0){
temp=R::rbeta(1,M);
probmat(l,q)=deflat*temp;
mumat(l,q)=R::rnorm(0,sqrt(v1));
deflat=deflat*(1-temp);
}
}
}
vec y(n),distance(Q),a1(1),a2(1),index=linspace<vec>(1,Q,Q),prob(1),candiv(1),zs(n),mu(n),Z(n),
tempv(1),stempv(1),sortindex(Q),findv(1),qmat(Q),pmat(Q),wmat(Q),qvec(Q),qvec1(Q);
y.fill(0);
qvec.fill(1);
qvec1.fill(1);
vec psiv(Q),gamv(Q);
psiv.fill(2);
gamv.fill(1);
Z.fill(1);
for(int j = 0; j<Q; ++j){
if(method>3.5 && method<4.5){
psiv[j]=R::rgamma(2,1);
}
if(method>4.5 && method<5.5){
psiv[j]=R::rexp(2.42);
}
if(method>5.5 && method<6.5){
psiv[j]=R::rcauchy(0,0.003);
}
}
for(int t = 0; t<n; ++t){  
if(t==0){
for (int lag = 0;  lag< p; ++lag){
Xmat(0,lag)=R::rnorm(0,1);
}
}
if(t>0 && t<p){
for (int lag = 0;  lag< t; ++lag){
Xmat(t,lag)=y[t-lag-1];
}
for (int lag = t;  lag< p; ++lag){
Xmat(t,lag)=Xmat(t-1,lag-1);
}
}
if(t>=p){
for (int lag = 0;  lag< p; ++lag){
Xmat(t,lag)=y[t-lag-1];
}
}
if(method>2.5){
a1=conv_to< colvec >::from(Xmat.submat(t,0,t,p-1));
for(int j = 0; j<Q; ++j){ 
a2=conv_to< colvec >::from(Candimat.submat(j,0,j,p-1)); 
tempdist=sum((a1-a2)%(a1-a2));
wmat[j]=log(gamv[j])-psiv[j]*tempdist;
}
wmat=wmat-max(wmat);
wmat=exp(wmat);
prob=wmat/sum(wmat);
for (int q = 0;  q< Q; ++q) {
if (q==0) {
deflat=1.0;
}
if (q>0) {
deflat=deflat-prob[q-1];
}
probsuc=prob[q]/deflat;
if (R::runif(0,1)<probsuc) {
Z[t]=q+1;
break;}
}
}
if(method<2.5 && method>1.5){
a1=conv_to< colvec >::from(Xmat.submat(t,0,t,p-1));
for(int j = 0; j<Q; ++j){
a2=conv_to< colvec >::from(Candimat.submat(j,0,j,p-1)); 
distance[j]=sqrt(sum((a1-a2)%(a1-a2)));
distance[j]=exp(2*log(distance[j]));
}
tempv=distance+index/100000000;
stempv=sort(tempv);
for(int j = 0; j<Q; ++j){
findv=index.elem(find(tempv==stempv[j]));
sortindex[j]=findv[0];
}
tempv=1/stempv;
qvec=qvec1;
qmat[0]=1;
for(int j = 0; j<Q; ++j){
pmat[j]=tempv[j]/sum(tempv.subvec(j,Q-1));
qvec[j]=1-pmat[j];
if(j>0){
qmat[j]=prod(qvec.subvec(0,j-1));
}
wmat[j]=(1-0.7*(1-pmat[j]))*qmat[j]*exp(j*log(0.7));
}
prob=wmat;
candiv=sortindex;  
for (int q = 0;  q< Q; ++q) {
if (q==0) {
deflat=1.0;
}
if (q>0) {
deflat=deflat-prob[q-1];
}
probsuc=prob[q]/deflat;
if (R::runif(0,1)<probsuc) {
Z[t]=candiv[q];
break;}
}  
}
if(method>0.5){
prob=conv_to< colvec >::from(probmat.submat(0,Z[t]-1,149,Z[t]-1));
candiv=conv_to< colvec >::from(mumat.submat(0,Z[t]-1,149,Z[t]-1)); 
for (int j = 0;  j< 150; ++j) {
if (j==0) {
deflat=1.0;
}
if (j>0) {
deflat=deflat-prob[j-1];
}
probsuc=prob[j]/deflat;
if (R::runif(0,1)<probsuc) {
mu[t]=candiv[j];
break;}
}
}
if(t==0){
h[0]=a0+beta*h0;
}
if(t>0){
h[t]=a0+al1*y[t-1]*y[t-1]+beta*h[t-1];
}     
if(method>0.5){  
y[t]=sqrt(h[t])*(mu0+R::rnorm(mu[t],sqrt(alpha*Sz)));
}
if(method<0.5){
y[t]=sqrt(h[t])*(mu0+R::rt(v)); 
}
}
check1=max(y);
check2=min(y);
if(check1<180 && check2>-90){
output=y;
break;
}
}
return output;
}

// [[Rcpp::export]]
Rcpp::List cpp8(arma::vec A0acf, arma::vec A1acf, arma::vec Bacf, arma::vec Hacf, 
arma::vec Szacf, arma::vec Macf, arma::vec Muacf, arma::vec psiacf,double tsize){
double temp=2/sqrt(tsize);
double a1=0,a2=0,a3=0,a4=0,a5=0,a6=0,a7=0,a8=0;
double b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,b7=0,b8=0;
for(int i = 0; i<1001; ++i){
  
if(b1<0.5){
if(A0acf[i]<temp){
b1=1;
}
if(b1<0.5){
a1=a1+A0acf[i];
}
}

if(b2<0.5){
if(A1acf[i]<temp){
b2=1;
}
if(b2<0.5){
a2=a2+A1acf[i];
}
}

if(b3<0.5){
if(Bacf[i]<temp){
b3=1;
}
if(b3<0.5){
a3=a3+Bacf[i];
}
}

if(b4<0.5){
if(Hacf[i]<temp){
b4=1;
}
if(b4<0.5){
a4=a4+Hacf[i];
}
}

if(b5<0.5){
if(Szacf[i]<temp){
b5=1;
}
if(b5<0.5){
a5=a5+Szacf[i];
}
}

if(b6<0.5){
if(Macf[i]<temp){
b6=1;
}
if(b6<0.5){
a6=a6+Macf[i];
}
}


if(b7<0.5){
if(Muacf[i]<temp){
b7=1;
}
if(b7<0.5){
a7=a7+Muacf[i];
}
}

if(b8<0.5){
if(psiacf[i]<temp){
b8=1;
}
if(b8<0.5){
a8=a8+psiacf[i];
}
}
}

a1=2*a1-1;
a2=2*a2-1;
a3=2*a3-1;
a4=2*a4-1;
a5=2*a5-1;
a6=2*a6-1;
a7=2*a7-1;
a8=2*a8-1;

return List::create(Named("a0")=a1,Named("a1")=a2,Named("beta")=a3,Named("h0")=a4,
Named("Sz")=a5,Named("M")=a6,Named("Mu0")=a7,Named("psi")=a8);
}


// [[Rcpp::export]]
Rcpp::List BTR0(arma::vec predv, arma::mat check, double n2, double n, double holdstart){
vec hold(200),Mnv(200);
double Mn=100.0,C0=0.01,holdp=-1.0,profit=0,loc1=0,loc2=0;
double btr=0.0;
hold.fill(0);
for(int j = 0; j<n2; ++j){
if(j==0){
holdp=holdstart;
}
if(j>0){
holdp=hold[j-1];
}
if(holdp>0.5){
loc1=10*j+n;
loc2=loc1-1;
profit=Mn*predv[j]-Mn*check(loc1,0)/check(loc2,1)+C0;
if(profit<0){
hold[j]=0;
}
if(profit>0){
hold[j]=1;
}
}
if(holdp<0.5){
loc1=10*j+n;
loc2=loc1-1;
profit=Mn*predv[j]*check(loc2,1)/check(loc1,0)-Mn-C0;
if(profit>0){
hold[j]=1;
}
if(profit<0){
hold[j]=0;
}
}
loc1=10*j+n+10;
loc2=loc1-10;
if(hold[j]>0.5 && holdp>0.5){
Mn=Mn*check(loc1,1)/check(loc2,1);
}
if(hold[j]>0.5 && holdp<0.5){
Mn=Mn*check(loc1,1)/check(loc2+1,0)-C0;
}
if(hold[j]<0.5 && holdp>0.5){
Mn=Mn*check(loc2+1,0)/check(loc2,1)-C0;
}
if(hold[j]<0.5 && holdp<0.5){
Mn=Mn;
}
Mnv[j]=Mn;
}
btr=Mn-100.0;
return List::create(Named("btr")=btr,Named("hold")=hold,Named("Mnv")=Mnv);
}
// [[Rcpp::export]]
arma::vec mixsim(double n1,double n2){
double n=n1+n2;
vec output(1),price(n+1);
price.fill(100);
for(int iter = 0; iter<20; ++iter){
double a0=0.03,al1=0.1,beta=0.82,h0=0.03;
double deflat=1,temps=0,probsuc=0,check1=0,check2=0;
vec y(n),h(n),err(n),prob(5);
prob[0]=0.24;
prob[1]=0.37;
prob[2]=0.2;
prob[3]=0.095;
prob[4]=0.095;
y.fill(0);
for(int t = 0; t<n; ++t){  
for (int j = 0;  j< 5; ++j) {
if (j==0) {
deflat=1.0;
}
if (j!=0) {
deflat=deflat-prob[j-1];
}
probsuc=prob[j]/deflat;
if (R::runif(0,1)<probsuc) {
temps=j;
break;}
}
if(temps==0){
err[t]=R::rnorm(-0.08,0.4);
}
if(temps==1){
err[t]=R::rnorm(-0.5,1);
}
if(temps==2){
err[t]=R::rnorm(1,0.5);
}
if(temps==3){
err[t]=R::rgamma(2,1);
}
if(temps==4){
err[t]=R::rgamma(2,1);
err[t]=-err[t];
}

if(t==0){
h[0]=a0+beta*h0;
}
if(t>0){
h[t]=a0+al1*y[t-1]*y[t-1]+beta*h[t-1];
}     

y[t]=sqrt(h[t])*err[t];

}
check1=max(y);
check2=min(y);
double finalcheck=10;
if(check1<180 && check2>-90){
for(int t = 1; t<n+1; ++t){  
price[t]=price[t-1]*exp(y[t-1]/100);
}
finalcheck=price[n]/price[n2];
}
if(finalcheck<1.5 && finalcheck>0.7){
output=y;
break;
}
}
return output;
}

// [[Rcpp::export]]
Rcpp::List BTR(arma::vec predv, arma::mat check, double n2, double n, double holdstart){
vec hold(n2),Mnv(n2);
double Mn=100.0,C0=0.01,holdp=-1.0,profit=0,loc1=0,loc2=0;
double btr=0.0;
hold.fill(0);
Mnv.fill(0);
for(int j = 0; j<n2; ++j){
if(j==0){
holdp=holdstart;
}
if(j>0){
holdp=hold[j-1];
}
if(holdp>0.5){
loc1=10*j+1;
loc2=loc1-1;
profit=Mn*predv[j]-Mn*check(loc1,0)/check(loc2,1)+C0;
if(profit<0){
hold[j]=0;
}
if(profit>0){
hold[j]=1;
}
}
if(holdp<0.5){
loc1=10*j+1;
loc2=loc1-1;
profit=Mn*predv[j]*check(loc2,1)/check(loc1,0)-Mn-C0;
if(profit>0){
hold[j]=1;
}
if(profit<0){
hold[j]=0;
}
}
loc1=10*j+10;
loc2=loc1-10;
if(hold[j]>0.5 && holdp>0.5){
Mn=Mn*check(loc1,1)/check(loc2,1);
}
if(hold[j]>0.5 && holdp<0.5){
Mn=Mn*check(loc1,1)/check(loc2+1,0)-C0;
}
if(hold[j]<0.5 && holdp>0.5){
Mn=Mn*check(loc2+1,0)/check(loc2,1)-C0;
}
if(hold[j]<0.5 && holdp<0.5){
Mn=Mn;
}
Mnv[j]=Mn;
}
btr=Mn-100.0;
return List::create(Named("btr")=btr,Named("hold")=hold,Named("Mnv")=Mnv);
}