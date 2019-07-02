
REF=readRDS('../TUMOR_REF.RDS')
.norm_exp<-function(x){
    y=x/sum(x)
    y=y*1000000
    return(y)
    }
REF=apply(REF,2,.norm_exp)

colnames(REF)=c('A','E','M','N','OL','OP','Q','NK','TA','TO','TP1','TP2','T')



###############
NUM=100
ALLR=matrix(0,ncol=NUM,nrow=13)
EXP=matrix(0,ncol=NUM,nrow=nrow(REF))

set.seed(1)
i=1
while(i<=NUM){
    this_r=runif(13)
    this_r=this_r/sum(this_r)
    this_exp=REF %*% this_r
    ALLR[,i]=this_r
    EXP[,i]=this_exp
    if(i%%10==1){print(i)}
i=i+1
}

rownames(ALLR)=colnames(REF)
rownames(EXP)=rownames(REF)
###############


