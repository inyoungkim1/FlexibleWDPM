fitdata=function(y,p=1,Q=3,model='pm',Nm=50000,Nb=10000,trace=F){
ptm=proc.time()
validname=c("svpm","svdpm","svwdpmp","svwdpmpr","svewdpm","svwdpmg","svwdpme","svwdpmh",
            "gpm","gdpm","gwdpmp","gewdpm","gwdpmg","gwdpme","gwdpmh")
id=0
for(i in 1:length(validname)){
if(validname[i]==model){
id=i 
}
}
if(id<0.5){
stop("invalid model name")}

fitresult=list(model=model)

if(id>1.5 && id<2.5){
Q=1
p=1
}
if(id>0.5 && id<1.5){
result=cpp3(y,Nm,Nb)
}
if(id>1.5 && id<3.5){
result=cpp4(y,p,Q,Nm,Nb)
}
if(id>3.5 && id<4.5){
result=cpp5(y,p,Q,Nm,Nb)
}
if(id>4.5 && id<8.5){
id1=id-2  
result=cpp6(y,p,Q,id1,Nm,Nb)
}
if(id>8.5 && id<15.5){
id1=id-9  
result=fitgarch(y,p,Q,id1,Nm)
}

if(id>0.5 && id<8.5){
len=length(result)
a1=floor(Nb/10)+1
a2=floor(Nm/10)
tsize=a2-a1+1
if(len>9.5){
Macf=as.vector(acf(result$MT[a1:a2],lag.max=1000,plot=F)$acf)
Szacf=as.vector(acf(result$SzT[a1:a2],lag.max=1000,plot=F)$acf)
Snacf=as.vector(acf(result$SnT[a1:a2],lag.max=1000,plot=F)$acf)
Wacf=as.vector(acf(result$WT[a1:a2],lag.max=1000,plot=F)$acf)
fiacf=as.vector(acf(result$fiT[a1:a2],lag.max=1000,plot=F)$acf)
Mu0acf=as.vector(acf(result$Mu0T[a1:a2],lag.max=1000,plot=F)$acf)
if(Q>1){
psiacf=as.vector(acf(result$psiT[a1:a2],lag.max=1000,plot=F)$acf)
}
if(Q==1){
psiacf=c(1)
}
}
if(len<9.5){
Macf=c(1)
Szacf=c(1)
Snacf=as.vector(acf(result$SnT[a1:a2],lag.max=1000,plot=F)$acf)
Wacf=c(1)
psiacf=c(1)
fiacf=as.vector(acf(result$fiT[a1:a2],lag.max=1000,plot=F)$acf)
Mu0acf=as.vector(acf(result$MuT[a1:a2],lag.max=1000,plot=F)$acf)
}
effi_list=cpp7(Macf, Szacf, Snacf, Wacf, psiacf, fiacf, Mu0acf, len,tsize)
fitresult=c(fitresult,list(LPS=result$LPS, MLL=result$MLL,
Sn=result$Sn,effi_Sn=effi_list$Sn, fi=result$fi, effi_fi=effi_list$fi))
if(model=="svpm"){
fitresult=c(fitresult,Mu=result$Mu, effi_mu=effi_list$Mu)
}
if(model!="svpm"){
fitresult=c(fitresult,mu0=result$mu0,effi_Mu0=effi_list$Mu0,M=result$M,effi_M=effi_list$M,
Sz=result$Sz,effi_Sz=effi_list$Sz,W=result$W,effi_W=effi_list$W)
if(model=="svwdpmp" || model=="svwdpmpr" || model=="svewdpm"){
fitresult=c(fitresult,psi=result$psi,effi_psi=effi_list$psi)
}
}
if(trace==T){
SnT=list(SnT=result$SnT)
fiT=list(fiT=result$fiT)
if(model=="svpm"){
MuT=list(MuT=result$MuT)
fitresult=c(fitresult,SnT,fiT,MuT)
}
if(model!="svpm"){
Mu0T=list(Mu0T=result$Mu0T)
SzT=list(SzT=result$SzT)
MT=list(MT=result$MT)
WT=list(WT=result$WT)
fitresult=c(fitresult,SnT,fiT,Mu0T,SzT,MT,WT)
}
if(model=="svwdpmp" || model=="svwdpmpr" || model=="svewdpm"){
psiT=list(psiT=result$psiT)
fitresult=c(fitresult,psiT)
}
hT=list(hT=result$hT)
fitresult=c(fitresult,hT)
}
if(model!="pm" && trace==T){
mumat=list(mumat=result$mumat)
fitresult= c(fitresult,mumat)
}
}
if(id>8.5){
fitresult=c(fitresult,list(MLL=result$MLL,MLR=result$MLR,
a0=mean(result$A0T),a1=mean(result$A1T),beta=mean(result$BT),h0=mean(result$HT)))
if(id>9.5){
fitresult=c(fitresult,list(Sz=mean(result$SzT),Mu0=mean(result$Mu0T),M=mean(result$MT)))
}
if(id>10.5 && id<12.5){
fitresult=c(fitresult,list(psi=mean(result$psiT)))  
}  

A0acf=as.vector(acf(result$A0T,lag.max=1000,plot=F)$acf)
A1acf=as.vector(acf(result$A1T,lag.max=1000,plot=F)$acf)
Bacf=as.vector(acf(result$BT,lag.max=1000,plot=F)$acf)
Hacf=as.vector(acf(result$HT,lag.max=1000,plot=F)$acf)


Szacf=rep(-10,1001)
Macf=rep(-10,1001)
Muacf=rep(-10,1001)
psiacf=rep(-10,1001)

if(id>9.5){
Szacf=as.vector(acf(result$SzT,lag.max=1000,plot=F)$acf)
Macf=as.vector(acf(result$MT,lag.max=1000,plot=F)$acf)
Muacf=as.vector(acf(result$Mu0T,lag.max=1000,plot=F)$acf)
}
if(id>10.5 && id<12.5){
psiacf=as.vector(acf(result$psiT,lag.max=1000,plot=F)$acf)  
}  
tsize=Nm/10

re1=cpp8(A0acf, A1acf, Bacf, Hacf, Szacf, Macf, Muacf, psiacf, tsize)
fitresult=c(fitresult,list(effi_a0=re1$a0,effi_a1=re1$a1,effi_beta=re1$beta,effi_h0=re1$h0))
if(id>9.5){
fitresult=c(fitresult,list(effi_Sz=re1$Sz,effi_M=re1$M,effi_Mu0=re1$Mu0))
}
if(id>10.5 && id<12.5){
fitresult=c(fitresult,list(effi_psi=re1$psi))  
}
if(trace==T){
fitresult=c(fitresult,list(A0T=result$A0T,A1T=result$A1T,BT=result$BT,HT=result$HT))
if(id>9.5){
fitresult=c(fitresult,list(SzT=result$SzT,MT=result$MT,Mu0T=result$Mu0T))
}
if(id>10.5 && id<12.5){
fitresult=c(fitresult,list(psiT=result$psiT))
}  
fitresult=c(fitresult,list(h=result$h))  
if(id>9.5){
fitresult=c(fitresult,list(mumat=result$lammat))
}
}    
}

ptm=proc.time()-ptm
fitresult=c(fitresult,cpu_time=as.vector(ptm)[3])
return (fitresult)
}