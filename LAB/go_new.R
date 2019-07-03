source('SCD_orig.R')

REF=read.table('Reference_expression.txt',sep='\t',header=T,row.names=1)
#REF[,4] = REF[,4] + REF[,5]
#REF=REF[,c(1,2,3, 4,6,7)]

REF=log(REF+1,10)
REF=apply(REF, 2, .norm_exp)
hist(REF[,1])

colnames(REF)=c('ASTRO','NEURON','OPC','NFOL','MOL','MICRO','ENDO')



###############
NUM=100
ALLR=matrix(0,ncol=NUM,nrow=ncol(REF))
EXP=matrix(0,ncol=NUM,nrow=nrow(REF))

set.seed(1)
i=1
while(i<=NUM){
    this_r=rnorm(ncol(REF))
    this_r=pnorm(this_r)
  
    #this_r=runif(this_r)
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

EXP=apply(EXP, 2, .norm_exp)


#####################
getRanMap<-function(x){
  oldx=x
  x[order(x)]=sort(rnorm(length(x)))
  y=x-min(x)
  y[which(oldx==0)]=0
  return(y)
  }


set.seed(1)
REXP=apply(EXP,2,getRanMap)

REXP=REXP
  set.seed(123)
  addNOI=function(x){
     M=mean(x)
     y=x+M/3 * rnorm(length(x))#(runif(length(x))*2-1)
     return(y)
  }

  REXP=t(apply(REXP,1,addNOI))
  REXP[which( REXP<0)]=0
  rownames( REXP)=rownames(REXP)
  colnames(REXP)=colnames(REXP)

REXP=apply(REXP, 2, .norm_exp)

colnames(REXP)=paste0('NOI_',c(1:ncol(REXP)))
colnames(EXP)=paste0('RAN_',c(1:ncol(EXP)))


COR=cor(REXP, EXP, method='spearman')


plot(REXP[,1] , EXP[,1])
################
#rownames(REXP)=toupper(rownames(REXP))
write.table(REXP, file='REXP_mix.txt',sep='\t',row.names=T,col.names=T,quote=F)





VAR=apply(REF,1,var)
VG=which(VAR>=-(sort(-VAR)[2000]))
VREF=REF[VG,]
write.table(VREF, file='VREF_sig.txt',sep='\t',row.names=T,col.names=T,quote=F)




source('SCD_orig.R')
OUT=SCD(REXP, VREF)
#plot(OUT$l, col=OUT$col, pch=16)

RATIO=OUT$out
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))

SCORE=0
i=1
while(i<=ncol(CORMAT)){
   SCORE=SCORE+CORMAT[i,i]
   i=i+1
}
print(SCORE)

LOSS=0
i=1
while(i<=ncol(CORMAT)){
   LOSS=LOSS+sum(CORMAT[which(c(1:ncol(CORMAT))!=i),i])
   i=i+1
}
print(LOSS)

print(SCORE-LOSS)









CB=read.table('CIBERSORT.Output_Job18.txt',header=T,row.names=1,sep='\t')
RATIO=t(CB[,c(1:(ncol(CB)-3))])
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))

SCORE=0
i=1
while(i<=ncol(CORMAT)){
   SCORE=SCORE+CORMAT[i,i]
   i=i+1
}
print(SCORE)

LOSS=0
i=1
while(i<=ncol(CORMAT)){
   LOSS=LOSS+sum(CORMAT[which(c(1:ncol(CORMAT))!=i),i])
   i=i+1
}
print(LOSS)

print(SCORE-LOSS)
















sc_mat=read.table('Zeisel_exp_sc_mat.txt',sep='\t',header=T,row.names=1)
sc_mat=log(sc_mat+1,10)
sc_mat=apply(sc_mat, 2, .norm_exp)



TAG=read.table('Zeisel_exp_sc_mat_cluster_merged.txt',header=T,sep='\t')
table(TAG[,2])
TAG[,2]=as.character(TAG[,2])
TAG[which(TAG[,2]=='astrocytes_ependymal'),2]='ASTRO'
TAG[which(TAG[,2]=='endothelial-mural'),2]='ENDO'
TAG[which(TAG[,2]=='microglia'),2]='MICRO'
TAG[which(TAG[,2]=='neurons'),2]='NEURON'
TAG[which(TAG[,2]=='oligodendrocytes'),2]='OLIGO'



SC.REF=.generate_ref(sc_mat, TAG, min_cell=1)
SC.REF=apply(SC.REF, 2, .norm_exp)

SC.REF=SC.REF[,c(1,4,5,3,2)]

VAR=apply(SC.REF, 1, var)
VG=which(rank(-VAR) <= 2000  )
V.SC.REF=SC.REF[VG,]

#V.SC.REF=V.SC.REF[,c(1,4,5,3,2)]


TMP=sc_mat[VG,]
colnames(TMP)=TAG[,2]

write.table(V.SC.REF, file='V.SC.REF_sig.txt',sep='\t',row.names=T,col.names=T,quote=F)
write.table(TMP, file='V.SC.REF_scmat_sig.txt',sep='\t',row.names=T,col.names=T,quote=F)






source('SCD_orig.R')
OUT=SCD(REXP, V.SC.REF, COMBAT=FALSE)
#OUT=SCD(REXP, V.SC.REF, COMBAT=FALSE)
#plot(OUT$l, col=OUT$col, pch=16)
RATIO=OUT$out
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]
sum(CORMAT)-(CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7])



source('SCD_orig.R')
OUT=SCD(REXP, V.SC.REF, COMBAT=TRUE)
#OUT=SCD(REXP, V.SC.REF, COMBAT=FALSE)
#plot(OUT$l, col=OUT$col, pch=16)
RATIO=OUT$out
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]
sum(CORMAT)-(CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7])


CB=read.table('CIBERSORT.Output_Job19.txt',header=T,row.names=1,sep='\t')
RATIO=t(CB[,c(1:(ncol(CB)-3))])
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]
sum(CORMAT)-(CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7])




CB=read.table('CIBERSORTx_Job7_Adjusted.txt',header=T,row.names=1,sep='\t')
RATIO=t(CB[,c(1:(ncol(CB)-3))])
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]
sum(CORMAT)-(CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7])



















































.simple_combine <- function(exp_sc_mat1, exp_sc_mat2){    
    exp_sc_mat=exp_sc_mat1
    exp_ref_mat=exp_sc_mat2
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    OUT=list()
    OUT$exp_sc_mat1=exp_sc_mat
    OUT$exp_sc_mat2=exp_ref_mat
    OUT$combine=cbind(exp_sc_mat,exp_ref_mat)
    return(OUT)
    }

COM=.simple_combine(REXP,VREF)$combine
#COM=as.data.frame(COM)

this_com= COM[,c(1,c((ncol(REXP)+1): ncol(COM) ))]
colnames(this_com)[1]='NOI'
this_com=as.data.frame(this_com)
fit=lm(NOI ~ ., data=this_com) 
this_coef=fit$coefficients
this_ratio=.norm_one(this_coef[c(2:length(this_coef))])








EXP=REXP
REF=VREF
LR=1
N=20
method='spearman'


REF=REF
    EXP=EXP
    method=method
    LR=LR
    N=N
    L=c()
    
    REF=apply(REF,2,.norm_exp)
    EXP=apply(EXP,2,.norm_exp)    
      
    ###########
    COR0=.scdcor(EXP, REF,method=method)
    NCOR0=apply(COR0, 2, .norm_one)
    
    
    MAT.LIST=list()
    
    i=1
    print(i)
    COR1=COR0
    NCOR1=NCOR0
    colnames(NCOR1)=colnames(EXP)
    rownames(NCOR1)=colnames(REF)
    MAT.LIST=c(MAT.LIST, list(NCOR1))
    
    this_exp = REF %*% NCOR1
    cor.mat=.scdcor(this_exp ,EXP, method=method)
    L=c(L, mean(cor.mat))

    
    EXP1=REF %*% NCOR1
    COR2=.scdcor(EXP1, REF, method=method)
    NCOR2=apply(COR2, 2, .norm_one)


    predict.NCOR0=NCOR0
    this_ref=1
    while(this_ref<=nrow(NCOR0)){
        this_fit=lm(NCOR0[this_ref,]~NCOR2[this_ref,])
        predict.NCOR0[this_ref,]=predict(this_fit)
        this_ref=this_ref+1
        }
    predict.NCOR0=apply(predict.NCOR0, 2, .norm_one)
    NCOR1 = NCOR1  + (NCOR0-predict.NCOR0) * LR









