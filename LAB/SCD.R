.norm_exp<-function(x){
    y=x/sum(x)
    y=y*1000000
    return(y)
    }

.scdcor  <- function(exp_sc_mat, exp_ref_mat, method='spearman'){
    #method = "pearson", "kendall", "spearman"
    method=method
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

.norm_one <- function(x){
    y=scale(x) 
    y=pnorm(y)
    y=y/sum(y)

    return(y)
    }


.generate_ref <- function(exp_sc_mat, TAG, min_cell=1, refnames=FALSE){
    NewRef=c()
    TAG[,2]=as.character(TAG[,2])
    if(refnames==FALSE){
        refnames=names(table(TAG[,2]))}
        else{refnames=refnames}
    outnames=c()
    for(one in refnames){
        this_col=which(TAG[,2]==one)
        if(length(this_col)>= min_cell){
            outnames=c(outnames,one)
            if(length(this_col) >1){
                #this_new_ref=apply(exp_sc_mat[,this_col],1,sum)
                this_new_ref=apply(exp_sc_mat[,this_col],1,mean)
                }
                else{this_new_ref = exp_sc_mat[,this_col]}
            NewRef=cbind(NewRef,this_new_ref)
            }
        }
    rownames(NewRef)=rownames(exp_sc_mat)
    colnames(NewRef)=outnames
    if(length(NewRef[1,])==1){
        NewRef=cbind(NewRef[,1], NewRef[,1])
        rownames(NewRef)=rownames(exp_sc_mat)
        colnames(NewRef)=c(outnames,outnames)
        }
    NewRef=apply(NewRef, 2, .norm_exp)
    return(NewRef)
    }



SCD <- function(EXP, REF, LR=0.1, N=10, method='spearman'){

    ######
    REF=REF
    LR=LR
    RREF=apply(REF,2,rank)
    SREF=t(apply(RREF,1,scale))
    PSREF=pnorm(SREF)
    XREF=REF
    
    i=1
    while(i<=ncol(REF)){
        this_r_max= PSREF[,i] - apply(PSREF[,which(c(1:ncol(REF)) !=i)], 1, max)
        #this_r=RREF[,i] - apply(RREF[,which(c(1:ncol(REF)) !=i)], 1, max)
        #this_r[which(this_r>=0.5)]= 2**this_r[which(this_r>=0.5)]
        #this_r[which(this_r>=0)]=  this_r[which(this_r>=0)]
        this_r=rep(0, nrow(REF))
        this_r[which(this_r_max>0)] = this_r_max[which(this_r_max>0)]
        #this_r[which(this_r_min<0)]= this_r_min[which(this_r_min<0)]
        this_r=this_r#/sum(this_r) * nrow(REF)
        XREF[, i] = REF[,i] * this_r
    	i=i+1
    }
    
    EXP=EXP
    method=method
    N=N
    L=c()
    
    REF=apply(REF,2,.norm_exp)
    XREF=apply(XREF,2,.norm_exp)
    colnames(XREF)=colnames(REF)
    rownames(XREF)=rownames(REF)
    EXP=apply(EXP,2,.norm_exp)    
      
    ###########
    COR0=.scdcor(EXP, XREF,method=method)
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
    this_exp=apply(this_exp,2,.norm_exp)
    cor.mat=.scdcor(this_exp ,EXP, method=method)
    L=c(L, mean(cor.mat))
    
    
    ###############
    i=2
    while(i<=N){
        print(i)
        EXP1=XREF %*% NCOR1
        COR2=.scdcor(EXP1, XREF, method=method)
        NCOR2=apply(COR2, 2, .norm_one)
        #COR1= COR1 + (COR0 - COR2)*LR
        #NCOR1=apply(COR1, 2, .norm_one)
        #NCOR1= NCOR1 * NCOR0 / NCOR2
        tmp_exp = REF %*% NCOR2
        tmp_exp=apply(tmp_exp,2,.norm_exp)
        cor.mat=.scdcor(tmp_exp ,EXP, method=method)

        tmp_mean=mean(cor.mat)
        if(tmp_mean > L[1]){
           NCOR1= NCOR1 + (NCOR2 - NCOR0) * LR
           }else{
           NCOR1= NCOR1 + (NCOR0 - NCOR2) * LR
           }

        #NCOR1= NCOR1 + (NCOR0-NCOR2) * LR
        NCOR1=apply(NCOR1, 2, .norm_one)

        colnames(NCOR1)=colnames(EXP)
        rownames(NCOR1)=colnames(REF)
        MAT.LIST=c(MAT.LIST, list(NCOR1))
        
        this_exp = REF %*% NCOR1
        this_exp=apply(this_exp,2,.norm_exp)
        cor.mat=.scdcor(this_exp ,EXP, method=method)
        L=c(L, mean(cor.mat))
        i=i+1
    }
    ###############
    
    max.index=which(L==max(L))[1]
    
    
    OUT=MAT.LIST[[max.index]]
    
    EXP.OUT=REF %*% OUT
    
    RESULT=list()
    RESULT$out=OUT
    RESULT$exp=EXP.OUT
    RESULT$l=L
    RESULT$col=rep('black',length(L))
    RESULT$col[max.index]='red'
    RESULT$max.index=max.index
    RESULT$mat.list=MAT.LIST   
    return(RESULT)
    }



