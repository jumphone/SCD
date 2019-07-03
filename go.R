source('SCD.R')

REF=read.table('Reference_expression.txt',sep='\t',header=T,row.names=1)
REF[,4] = REF[,4] + REF[,5]
REF=REF[,c(1,2,4,6,7)]

REF=log(REF+1,10)
REF=apply(REF, 2, .norm_exp)
hist(REF[,1])

colnames(REF)=c('ASTRO','NEURON','OLIGO','MICRO','ENDO')




###############
NUM=100
ALLR=matrix(0,ncol=NUM,nrow=ncol(REF))
EXP=matrix(0,ncol=NUM,nrow=nrow(REF))

set.seed(1)
i=1
while(i<=NUM){
    this_r=runif(ncol(REF))
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
     y=x+M/3*(runif(length(x))*2-1)
     return(y)
  }

  REXP=t(apply(REXP,1,addNOI))
  REXP[which( REXP<0)]=0
  rownames( REXP)=rownames(REXP)
  colnames(REXP)=colnames(REXP)

REXP=apply(REXP,2,.norm_exp) 

################
#rownames(REXP)=toupper(rownames(REXP))
write.table(REXP, file='REXP_mix.txt',sep='\t',row.names=T,col.names=T,quote=F)
#####################################


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

VAR=apply(SC.REF, 1, var)
VG=which(rank(-VAR) <= 2000  )
V.SC.REF=SC.REF[VG,]

V.SC.REF=V.SC.REF[,c(1,4,5,3,2)]


TMP=sc_mat[VG,]
colnames(TMP)=TAG[,2]

write.table(V.SC.REF, file='V.SC.REF_sig.txt',sep='\t',row.names=T,col.names=T,quote=F)
write.table(TMP, file='V.SC.REF_scmat_sig.txt',sep='\t',row.names=T,col.names=T,quote=F)
#####################################


OUT=SCD(REXP, V.SC.REF, N=20, method='spearman')
plot(OUT$l, col=OUT$col, pch=16)


CORMAT=cor(t(OUT$out), t(ALLR))


pdf('RESULT_SCD.pdf',width=7,height=7)
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))
dev.off()



CB=read.table('CIBERSORT.Output_Job12.txt',header=T,row.names=1,sep='\t')

CORMAT=cor(CB[,c(1:(ncol(CB)-3))], t(ALLR))
pdf('RESULT_CB.pdf',width=7,height=7)
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))
dev.off()




CBX=read.table('CIBERSORTx_Job5_Adjusted.txt',header=T,row.names=1,sep='\t')



CORMAT=cor(CBX[,c(1:(ncol(CBX)-3))], t(ALLR))
pdf('RESULT_CBX.pdf',width=7,height=7)
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))
dev.off()












































.get_cor  <- function(exp_sc_mat, exp_ref_mat, method='kendall',CPU=4, print_step=10, gene_check=FALSE){
    #method = "pearson", "kendall", "spearman"
    #exp_sc_mat: single-cell gene expression matrix; row is gene; col is sample; should have row name and col name
    #exp_ref_mat: reference gene expression matrix; row is gene; col is sample; should have row name and col name
    
    library(parallel)
    #Step 1. get overlapped genes
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    
    #Step 2. calculate prob
    SINGLE <- function(i){
        library('pcaPP')
        exp_sc = as.array(exp_sc_mat[,i])
        log_p_sc_given_ref_list=rep(0,length(colname_ref))
        j=1
        while(j<=length(colname_ref)){
            exp_ref = as.array(exp_ref_mat[,j])
            log_p_sc_given_ref=cor(exp_sc,exp_ref, method=method)
            log_p_sc_given_ref_list[j]=log_p_sc_given_ref
            j=j+1}
        ################################
        #if(i%%print_step==1){print(i)}
        return(log_p_sc_given_ref_list)
        }
    #######################################
    cl= makeCluster(CPU,outfile='')
    RUN = parLapply(cl=cl,1:length(exp_sc_mat[1,]), SINGLE)
    stopCluster(cl)
    #RUN = mclapply(1:length(colname_sc), SINGLE, mc.cores=CPU)
    LOG_P_SC_GIVEN_REF = c()
    for(log_p_sc_given_ref_list in RUN){
        LOG_P_SC_GIVEN_REF=cbind(LOG_P_SC_GIVEN_REF, log_p_sc_given_ref_list)}
    #######################################
    rownames(LOG_P_SC_GIVEN_REF)=colname_ref
    colnames(LOG_P_SC_GIVEN_REF)=colname_sc
    return(LOG_P_SC_GIVEN_REF)
    }



.norm_one <- function(x, cutoff=0){
    used_index=which(x>cutoff)
    non_used_index=which(x<=cutoff)
    y=x
    #y[used_index]=x[used_index]/sum(x[used_index])
    y[used_index]=scale(x[used_index]) 
    y[used_index]=pnorm(y[used_index])
    y[used_index]=y[used_index]/sum(y[used_index])
    y[non_used_index]=0
    return(y)
    }

.norm_one <- function(x){
    y=scale(x) 
    y=pnorm(y)
    y=y/sum(y)

    return(y)
    }





.scdcor  <- function(exp_sc_mat, exp_ref_mat, method='spearman'){
    #method = "pearson", "kendall", "spearman"
    method=method

    #Step 1. get overlapped genes
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    this_cor=cor(exp_ref_mat,exp_sc_mat,method=method)    
    return(this_cor)
    }
##






COR0=.scdcor(EXP, REF,method='spearman')
NCOR0=apply(COR0, 2, .norm_one)
cor.test(NCOR0[,1],ALLR[,1])


COR1=COR0
NCOR1=NCOR0
tmp=0

i=1
while(i<=50){
    print(i)
    EXP1=REF %*% NCOR1
    COR2=.scdcor(EXP1, REF,method='spearman')
    COR1= COR0 + COR1 - COR2
    NCOR1=apply(COR1, 2, .norm_one)
    this_exp = REF %*% NCOR1
    tmp=this_cor
    i=i+1
}



plot(NCOR1[,1],ALLR[,1])

CORMAT=cor(t(NCOR1),t(ALLR))

library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))






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
     y=x+M/3*(runif(length(x))*2-1)
     return(y)
  }

  REXP=t(apply(REXP,1,addNOI))
  REXP[which( REXP<0)]=0
  rownames( REXP)=rownames(REXP)
  colnames(REXP)=colnames(REXP)

REXP=apply(REXP,2,.norm_exp) 



EXP=REXP



COR0=.scdcor(EXP, REF,method='spearman')
NCOR0=apply(COR0, 2, .norm_one)
cor.test(NCOR0[,1],ALLR[,1])


COR1=COR0
NCOR1=NCOR0

i=1
while(i<=50){
    print(i)
    EXP1=REF %*% NCOR1
    COR2=.scdcor(EXP1, REF,method='spearman')
    COR1= COR0 + COR1 - COR2
    NCOR1=apply(COR1, 2, .norm_one)
    this_exp = REF %*% NCOR1
    i=i+1
}

OUT=NCOR1
colnames(OUT)=colnames(EXP)
rownames(OUT)=colnames(REF)











COR0=.scdcor(EXP, REF,method='spearman')
NCOR0=apply(COR0, 2, .norm_one)
cor.test(NCOR0[,1],ALLR[,1])


COR1=COR0
NCOR1=NCOR0

i=1
while(i<=50){
    print(i)
    EXP1=REF %*% NCOR1
    COR2=.scdcor(EXP1, REF,method='spearman')
    COR1= COR0 + COR1 - COR2
    NCOR1=apply(COR1, 2, .norm_one)
    this_exp = REF %*% NCOR1
    i=i+1
}

OUT=NCOR1
colnames(OUT)=colnames(EXP)
rownames(OUT)=colnames(REF)

















this_out=SCD(EXP,REF,N=50)

CORMAT=cor(t(this_out),t(ALLR))

library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))

