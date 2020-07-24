simdata=function(n=500,p=1,Q=1,model="pm",random=F){
if(model!="pm" && model!="dpm" && model!="wdpmp" && model!="wdpme" && model!="wdpmh" && model!="wdpmd"){
stop("invalid model name")
}
if(n<100){
stop("Sample size should be at least 100.")
}
if(n>10000){
stop("To avoid slow simulation, sample size should not be larger than 10000.")
}
spind=1
if(model=="pm" || model=="dpm"){
spind=0
if(model=="pm"){
y=cpp1(n,1,1,1)
}
if(model=="dpm"){
y=cpp1(n,1,1,2)
}
}
if(spind>0.5){
if(p>6){
stop("To avoid slow simulation, p should not be larger than 6.")
}
if(Q>800){
stop("To avoid slow simulation, Q should not be larger than 800.")
}
if(p<1){
p=1
}
if(Q<1){
Q=3
}
if(model=="wdpmp" && random==F){
y=cpp1(n,p,Q,3)
}
if(model=="wdpmp" && random==T){
y=cpp1(n,p,Q,4)
}
if(model=="wdpmd" && random==F){
y=cpp2(n,p,Q,1,1)
}
if(model=="wdpmd" && random==T){
y=cpp2(n,p,Q,2,1)
}
if(model=="wdpme" && random==F){
y=cpp2(n,p,Q,1,2)
}
if(model=="wdpme" && random==T){
y=cpp2(n,p,Q,2,2)
}
if(model=="wdpmh" && random==F){
y=cpp2(n,p,Q,1,5)
}
if(model=="wdpmh" && random==T){
y=cpp2(n,p,Q,2,5)
}
}
return (y)
}
  
  



