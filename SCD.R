
.norm_exp<-function(x){
    y=x/sum(x)
    y=y*1000000
    return(y)
    }



.norm_one <- function(x){
    y=scale(x) 
    y=pnorm(y)
    y=y/sum(y)

    return(y)
    }

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



.combat <- function(EXP, BATCH){
    library(sva)
    library(limma)
    pheno = data.frame(batch=as.matrix(BATCH))
    edata = EXP
    batch = pheno$batch
    modcombat = model.matrix(~1, data=pheno)
    combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
    return(combat_edata)
    }


SCD <- function(EXP, REF, COMBAT=TRUE){
    ##############################
    print('Start!')
    print(Sys.time())
    ##############################
    ######
    
    REF=REF
    EXP=EXP
    COMBAT=COMBAT

    NREF=apply(REF,2,.norm_exp)
    NEXP=apply(EXP,2,.norm_exp)
    colnames(NREF)=colnames(REF)
    rownames(NREF)=rownames(REF)
    colnames(NEXP)=colnames(EXP)
    rownames(NEXP)=rownames(EXP)


    ##########
    if(COMBAT==TRUE){

        ALL=.simple_combine(NREF, NEXP)$combine
        #print(dim(ALL))
        BATCH=c(rep('REF',ncol(NREF)),rep('EXP',ncol(NEXP)))
        #print(BATCH)
        ALL.combat=.combat(ALL, BATCH)
    
        NREF=ALL.combat[,c(1:ncol(REF))]
        NEXP=ALL.combat[,c((ncol(REF)+1):ncol(ALL))]
        }
    ############
      

    COM=.simple_combine(NEXP,NREF)$combine
    #COM=as.data.frame(COM)
    OUT=c()
    i=1
    while(i<=ncol(EXP)){
        this_com= COM[,c(i,c((ncol(EXP)+1): ncol(COM) ))]
        colnames(this_com)[1]='NOI'
        this_com=as.data.frame(this_com)
        fit=lm(NOI ~ ., data=this_com) 
        this_coef=fit$coefficients
        this_ratio=.norm_one(this_coef[c(2:length(this_coef))])
        OUT=cbind(OUT,this_ratio)
        i=i+1
    }
    rownames(OUT)=colnames(REF)
    colnames(OUT)=colnames(EXP)


    ##################################
    RESULT=list()
    RESULT$exp=EXP
    RESULT$ref=REF
    RESULT$combat=COMBAT
    
    RESULT$out=OUT
    
    if(COMBAT==TRUE){
        RESULT$combat.exp=ALL
        RESULT$combat.batch=BATCH
        RESULT$combat.adjexp=ALL.combat
        }
    ##############################
    print('Finished!')
    print(Sys.time())
    ##############################
    return(RESULT)
    }



