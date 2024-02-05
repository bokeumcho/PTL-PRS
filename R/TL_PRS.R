##cor=bim_sum_stats[which(LDblocks2[[1]]==i),]; num=which(i==unique(LDblocks2[[1]]));nsnp=nrow(bim_sum_stats);temp.file=paste0(tempfile,"_step1")
block_calculation2<-function(cor,num,nsnp,temp.file, lr_list, iter){
  temp_file=paste0(temp.file,"_block_",num)
  # write.table(cor$V2,file=temp_file,col.names=F,row.names=F,quote=F)
  # cmd = paste0("plink --bfile ",train_file," --extract ",temp_file,   " --recodeA  --out ", temp_file,"_Geno.txt")
  # system(cmd)
  
  Gtemp = try(readSomeSnp(cor$V2),silent=T) #ref
  cat('read Gtemp block ',num,'\n')
  # Gtemp=try(as.data.frame(fread(paste0(temp_file,"_Geno.txt.raw"),header=T)),silent=T)
  # if (file.exists(temp_file)) {file.remove(temp_file)}
  # if (file.exists(paste0(temp_file,"_Geno.txt.nosex"))) {file.remove(paste0(temp_file,"_Geno.txt.nosex"))}
  # if (file.exists(paste0(temp_file,"_Geno.txt.log"))) {file.remove(paste0(temp_file,"_Geno.txt.log"))}
  # if (file.exists(paste0(temp_file,"_Geno.txt.raw"))) {file.remove(paste0(temp_file,"_Geno.txt.raw"))}

  if (class(Gtemp)=="try-error"){
    return(NULL)
    #GG=diag(nrow(cor));colnames(GG)=paste0(cor$V2,"_",cor$V5) 
    #geno_info=as.data.frame(t(sapply(colnames(GG),   split_SNPandA1   ) )) ;colnames(geno_info)=c("SNP","A1")
    #geno_info$mean=NA; geno_info$maf=NA; geno_info$sd=NA
  }else{
    Gtemp1 = t(Gtemp)

    #handling NA (mean imputation)
    Gtemp1 <- apply(Gtemp1, 2, function(x){
      #x[is.na(x)] <- mean(x, na.rm=TRUE)
      x[which(x==9)] <- mean(x[which(x!=9)])
      return(x)
    })

    Gtemp2 = Gtemp1[,apply(Gtemp1, 2, sd, na.rm = TRUE)!=0] #4727 -> 4710 snps, n*p
    #Gtemp2 <- fread('output/HDL/G_ref')

    #cat('dimension of Gtemp2 and cor',dim(Gtemp2), dim(cor))

    GG = cor(as.matrix(Gtemp2))
    geno_info = as.data.frame(cor[,c('V2','V5')]);colnames(geno_info)=c("SNP","A1")
    geno_info = geno_info[apply(Gtemp1, 2, sd, na.rm = TRUE)!=0,]    
    geno_info$mean=colMeans(as.matrix(Gtemp2),na.rm=T); geno_info$maf=geno_info$mean/2; geno_info$sd=sqrt(2*geno_info$maf*(1-geno_info$maf))

    # GG=cor(as.matrix(Gtemp[,7:ncol(Gtemp)]))
    # geno_info=as.data.frame(t(sapply(colnames(Gtemp)[7:ncol(Gtemp)],   split_SNPandA1   ) )) ;colnames(geno_info)=c("SNP","A1")
    # geno_info$mean=colMeans(as.matrix(Gtemp[,7:ncol(Gtemp)]),na.rm=T); geno_info$maf=geno_info$mean/2; geno_info$sd=sqrt(2*geno_info$maf*(1-geno_info$maf))
  }

  list1=which(geno_info$sd==0)
  if (length(list1)>0){
    geno_info=geno_info[-list1,]
    GG=GG[-list1,-list1]
  }
  if (nrow(geno_info)==0){
    return(NULL)
  } else {
    geno_info$order=1:nrow(geno_info)
    geno_info2=merge(cor[,c("V1","V2","V5","Beta2","cor")],geno_info, by.x="V2",by.y="SNP",sort=F)
    # No need of allele-flipping again
    # flag_nomatch=which(geno_info2$A1 != geno_info2$V5)
    # if (length(flag_nomatch)>0){
    #   geno_info2$Beta2[flag_nomatch]=-geno_info2$Beta2[flag_nomatch]
    #   geno_info2$cor[flag_nomatch]=-geno_info2$cor[flag_nomatch]
    # }
    GG2=as.matrix(GG[geno_info2$order,geno_info2$order])
    gy=geno_info2$cor
    betatemp=geno_info2$Beta2*geno_info2$sd
    u0=gy-GG2%*%betatemp
    beta.all=cbind(u0, betatemp)
    for (factor1 in lr_list){
      k=1
      betatemp=beta.all[,2]
      u0=beta.all[,1]
      while (k<=iter){
        ##betanew=c()
        learningrate=1/nsnp*factor1
        if (learningrate>1){learningrate=1}
        ##print(learningrate)
        for (j in 1:length(betatemp)){
          beta_old=betatemp[j]
          betatemp[j]=(learningrate*u0[j]+beta_old)/ 1
          u0=u0-GG2[,j]*(betatemp[j]-beta_old)
        }
        beta.all=cbind(beta.all,betatemp)
        k=k+1
      } 
    }
    geno_info2=cbind(geno_info2,beta.all)
    return(geno_info2)
  } 
}##function end




##PRStr_calculation2(sum_stats_target, train_file, sum_stats, LDblocks, cluster=cluster,temp.file=paste0(tempfile,"_step1"))
##temp.file=paste0(tempfile,"_step1")
PRStr_calculation2<-function(sum_stats_target, ref_file, sum_stats, LDblocks, cluster=NULL,temp.file, lr_list, iter){
  possible.LDblocks <- c("EUR.hg19", "AFR.hg19", "ASN.hg19", 
                         "EUR.hg38", "AFR.hg38", "ASN.hg38") 
  if(!is.null(LDblocks)) {
    if(is.character(LDblocks) && length(LDblocks) == 1) {
      if(LDblocks %in% possible.LDblocks) {
        LDblocks <- data.table::fread(system.file(paste0("data/Berisa.",  LDblocks, ".bed"),  package="lassosum"), header=T)
      } else {
        stop(paste("I cannot recognize this LDblock. Specify one of", 
                   paste(possible.LDblocks, collapse=", ")))
      }
    }
    if(is.factor(LDblocks)) LDblocks <- as.integer(LDblocks)
    if(is.vector(LDblocks)) stopifnot(length(LDblocks) == length(cor)) else 
      if(is.data.frame(LDblocks) || is.data.table(LDblocks)) {
        LDblocks <- as.data.frame(LDblocks)
        stopifnot(ncol(LDblocks) == 3)
        stopifnot(all(LDblocks[,3] >= LDblocks[,2]))
        LDblocks[,1] <- as.character(sub("^chr", "", LDblocks[,1], ignore.case = T))
      }
  } else {
    stop(paste0("LDblocks must be specified. Specify one of ", 
                paste(possible.LDblocks, collapse=", "), 
                ". Alternatively, give an integer vector defining the blocks, ", 
                "or a .bed file with three columns read as a data.frame."))
  }
  
  ref.bim <- fread(paste0(ref_file, ".bim"))
  ref.bim$V1 <- as.character(sub("^chr", "", ref.bim$V1, ignore.case = T))
  ref.bim$order=1:nrow(ref.bim)
  bim_sum_stats=merge(ref.bim, sum_stats_target,by.x="V2",by.y="SNP",order=F)
  bim_sum_stats=bim_sum_stats[order(bim_sum_stats$order),]
  
# remove duplicates
  bim_sum_stats = bim_sum_stats[!duplicated(bim_sum_stats$V2),]

  bim_sum_stats$Beta2=NA
  flag1=which(bim_sum_stats$V5==bim_sum_stats$A1)
  if (length(flag1)>0){  bim_sum_stats$Beta2[flag1]=bim_sum_stats$Beta[flag1]}
  flag2=which(bim_sum_stats$V6==bim_sum_stats$A1)
  if (length(flag2)>0){  bim_sum_stats$Beta2[flag2]=-bim_sum_stats$Beta[flag2];  bim_sum_stats$cor[flag2]=-bim_sum_stats$cor[flag2];}
  
  # fix the allele of ref 
  # ref.bim.flip = ref.bim
  # flag3=which(ref.bim.flip$V2 %in% bim_sum_stats[flag2,'V2'][[1]])
  
  # ref.bim.flip$A1 = NA
  # temp = ref.bim.flip[flag3,'V6']
  # ref.bim.flip[flag3,'V6'] = ref.bim.flip[flag3,'V5']
  # ref.bim.flip[flag3,'V5'] = temp

  # write.table(ref.bim.flip, file='/media/leelabsg-storage1/bokeum/TLPRS/complication/data/1000G/EURID_1000G_plink_no3alleles_flip.bim',col.names=F,row.names=F,quote=F)
  
  ## do not flip allele
  # bim_sum_stats$Beta2=bim_sum_stats$Beta

  bim_sum_stats=bim_sum_stats[which(! is.na(bim_sum_stats$Beta2)),c("V2","V1","V4","V5","V6","order","Beta2","cor")] 
  
  ref.extract <- rep(FALSE, nrow(ref.bim))
  ref.extract[bim_sum_stats$order] <- TRUE
  
  if(!is.null(LDblocks)) {
      LDblocks2 <- splitgenome2(CHR = ref.bim$V1[ ref.extract], 
                              POS = ref.bim$V4[ ref.extract],
                              ref.CHR = LDblocks[,1], 
                              ref.breaks = LDblocks[,3])
      # Assumes base 1 for the 3rd column of LDblocks (like normal bed files)
  } 
  
  if(is.null(cluster)) {
  	results.list <- lapply(unique(LDblocks2[[1]]), function(i) {
    		block_calculation2(cor=bim_sum_stats[which(LDblocks2[[1]]==i),], num=which(i==unique(LDblocks2[[1]])),nsnp=nrow(bim_sum_stats),temp.file, lr_list, iter)
  	})
  } else {
  	results.list <-  parallel::parLapplyLB(cluster,unique(LDblocks2[[1]]), function(i) {
    		block_calculation2(cor=bim_sum_stats[which(LDblocks2[[1]]==i),], num=which(i==unique(LDblocks2[[1]])),nsnp=nrow(bim_sum_stats),temp.file, lr_list, iter)
  	})
  }
  
  results.list<-do.call("rbind", results.list)

  return(results.list)
}







######################The function used summary statistics for training############### 
##ped_file=ped.file;Covar_name="";Y_name=kword;Ytype="C"; train_file=train.bfile;test_file=test.bfile;sum_stats_file=beta.file;LDblocks="EUR.hg19"
##ped.file,"",kword, Ytype="C",train.bfile,test.bfile,beta.file,target_sumstats_file,LDblocks="EUR.hg19",tempfile
TL_PRS<-function(ped_val_file,Covar_name,Y_name, Ytype="C",ref_file,sum_stats_file,target_sumstats_file, LDblocks="EUR.hg19",outfile,cluster=NULL, target_only=FALSE, lr_list, iter){
	tempfile=outfile
	# out1=PRStr_main_check(ped_val_file,Covar_name,Y_name, Ytype,ref_file,sum_stats_file,LDblocks)
	#if (out1!=0){stop(out1)}

  sourceCpp("/media/leelabsg-storage1/bokeum/TLPRS/source/BedFileReader.cpp")
  setClass("BedFileReader", representation( pointer = "externalptr" ) )
  BedFileReader_method <- function(name) {paste( "BedFileReader", name, sep = "__" ) }
  setMethod( "$", "BedFileReader", function(x, name ) {function(...) .Call(BedFileReader_method(name),x@pointer, ... )})
  setMethod("initialize","BedFileReader", function(.Object, ...) {
    .Object@pointer <-.Call(BedFileReader_method("new"), ... )
    .Object})

  BedFileReader <- new( "BedFileReader", paste0(ref_file,".fam"), paste0(ref_file,".bim"), paste0(ref_file,".bed"))

  result = try(BedFileReader$snp_index_func(), silent=TRUE)
  
  # preprocess summary stats
	sum_stats=data.frame(fread(sum_stats_file))
	if (ncol(sum_stats)==4){ 
		if (sum(colnames(sum_stats) %in% c("V1","V2","V3","V4"))==4){
			colnames(sum_stats)=c("SNP","CHR","A1","Beta")
		}
	} 
	sum_stats=sum_stats[,c("SNP","CHR","A1","Beta")] 
	sum_stats_file=paste0(tempfile,"_original_sum_stats.txt")
	write.table(sum_stats, file=sum_stats_file,col.names=T,row.names=F,quote=F)
	
  # target sum stats
  sum_stats_target=fread(target_sumstats_file)[1:1000,]
  
  if (target_only) snps_top100 = sum_stats_target[order(sum_stats_target$p),][1:100,];

  sum_stats_target=merge(sum_stats,sum_stats_target,by="SNP",sort=F)

  if (target_only){
    colnames(sum_stats_target)[4] = 'A1'
    sum_stats_target = merge(sum_stats_target, snps_top100, by=c("SNP",'A1','beta','N','p'), all=TRUE, sort=FALSE)

    sum_stats_target$Beta[is.na(sum_stats_target$Beta)] = 0
    sum_stats_target$A1.x[is.na(sum_stats_target$A1.x)] = sum_stats_target[is.na(sum_stats_target$A1.x), 'A1'] 
  }

  if (sum(sum_stats_target$p<=1E-320)>0){ sum_stats_target$p[sum_stats_target$p<=1E-320]=1E-320}

  if (!is.numeric(sum_stats_target$beta)) sum_stats_target$beta = as.numeric(sum_stats_target$beta)
  if (!is.numeric(sum_stats_target$p)) sum_stats_target$p = as.numeric(sum_stats_target$p)

  sum_stats_target$cor=lassosum::p2cor(p = sum_stats_target$p, n = median(sum_stats_target$N,na.rm=T), sign=sum_stats_target$beta)
  
  # allele-flipping
  flag=which(sum_stats_target$A1.x !=sum_stats_target$A1.y)
  if (length(flag)>0){sum_stats_target$cor[flag]=-sum_stats_target$cor[flag]}
  sum_stats_target=sum_stats_target[,c("SNP","A1.x","Beta","cor")];colnames(sum_stats_target)[2]="A1";
  gc()

  # block-wise gradient descent
	beta_list=as.data.frame(PRStr_calculation2(sum_stats_target, ref_file, sum_stats, LDblocks, cluster=cluster,temp.file=paste0(tempfile,"_step1"), lr_list, iter))
	beta_list=as.data.frame(beta_list[,-c(6,10)]) #del A1&order; idx changed due to V1
	colnames(beta_list)[1:3]=c("SNP","CHR","A1")
	write.table(beta_list,file=paste0(tempfile,"_beta.candidates.txt"),row.names=F,quote=F,col.names=T)
  
  #tempfile=outfile
  #beta_list = as.data.frame(fread(paste0(tempfile,"_beta.candidates.txt")))
	
  out1=PRStr_tuning(beta_list, ped_val_file, Covar_name, Y_name, Ytype, lr_list, iter)

  if (file.exists(paste0(tempfile,"_original_sum_stats.txt"))) {file.remove(paste0(tempfile,"_original_sum_stats.txt"))}
  if (file.exists(paste0(tempfile,"_step0.train.PRS.nosex"))) {file.remove(paste0(tempfile,"_step0.train.PRS.nosex"))}
  if (file.exists(paste0(tempfile,"_step0.train.PRS.log"))) {file.remove(paste0(tempfile,"_step0.train.PRS.log"))}
  if (file.exists(paste0(tempfile,"_step0.train.PRS.profile"))) {file.remove(paste0(tempfile,"_step0.train.PRS.profile"))}  
	write.table(out1$best.beta,file=paste0(tempfile,"_best.beta.txt"),row.names=F,quote=F,col.names=T)
	write.table(out1$best.params,file=paste0(tempfile,"_best.params.txt"),row.names=F,quote=F,col.names=T)

	return(out1)
}

PTL_PRS<-function(Covar_name,Y_name, Ytype="C",sum_stats_file,target_sumstats_file, subprop, ref_file_ps, LDblocks="EUR.hg19",outfile,cluster=NULL, target_sumstats_train_file=NULL, target_sumstats_val_file=NULL, ps){
	tempfile=outfile
	# out1=PRStr_main_check(ped_val_file,Covar_name,Y_name, Ytype,ref_file,sum_stats_file,LDblocks)
	#if (out1!=0){stop(out1)}

	sum_stats=data.frame(fread(sum_stats_file))
	if (ncol(sum_stats)==4){ 
		if (sum(colnames(sum_stats) %in% c("V1","V2","V3","V4"))==4){
			colnames(sum_stats)=c("SNP","CHR","A1","Beta")
		}
	} 
	sum_stats=sum_stats[,c("SNP","CHR","A1","Beta")] 
	sum_stats_file=paste0(tempfile,"_original_sum_stats.txt")
	write.table(sum_stats, file=sum_stats_file,col.names=T,row.names=F,quote=F)

if (ps){
  ## pseudo summ generation
  target_sumstats <- pseudo_sum(target_sumstats_file, subprop, ref_file_ps, tempfile)
  target_sumstats_train <- target_sumstats[[1]]
  target_sumstats_val <- target_sumstats[[2]]
  rm(target_sumstats)
}

else{
  target_sumstats_train <- fread(target_sumstats_train_file)
  target_sumstats_val <- fread(target_sumstats_val_file)
}

for(group in c('train','val')){
    sum_stats_target=get(paste0('target_sumstats_',group))
    
    sum_stats_target=merge(sum_stats,sum_stats_target,by="SNP", sort=F)
	
      # allele flipping (only cor, not Beta, needs flipping)
      flag=which(sum_stats_target$A1.x !=sum_stats_target$A1.y)
      if (length(flag)>0){sum_stats_target$cor[flag]=-sum_stats_target$cor[flag]}
      # CLUMPING
      # sum(sum_stats_target$cor>=0.2)
      
    if(group=='train'){
      #Beta from source summ stat !!
      sum_stats_target=sum_stats_target[,c("SNP","A1.x","Beta","cor","N")]
    }
    else{
      sum_stats_target=sum_stats_target[,c("SNP","A1.x","cor","N")]
    }
    colnames(sum_stats_target)[2]="A1"; 
    assign(paste0('sum_stats_target_',group), sum_stats_target)
    gc()
  }

  # 1699 blocks in total
	beta_list=as.data.frame(PRStr_calculation2(sum_stats_target_train, ref_file, sum_stats, LDblocks, cluster=cluster,temp.file=paste0(tempfile,"_step1"), lr_list, iter))
	beta_list=as.data.frame(beta_list[,-c(6,10)]) #del A1&order; idx changed due to V1
	colnames(beta_list)[1:3]=c("SNP","CHR","A1")
	write.table(beta_list,file=paste0(tempfile,"_beta.candidates.txt"),row.names=F,quote=F,col.names=T)
  #Beta.all2 = as.data.frame(fread(paste0(tempfile,"_noflip_beta.candidates.txt")))
	
  out1=PRStr_tuning_pv(beta_list, sum_stats_target_train, sum_stats_target_val, sum_stats_target_test, maf_file, Covar_name, Y_name, Ytype, lr_list, iter)

  	if (file.exists(paste0(tempfile,"_original_sum_stats.txt"))) {file.remove(paste0(tempfile,"_original_sum_stats.txt"))}
  	if (file.exists(paste0(tempfile,"_step0.train.PRS.nosex"))) {file.remove(paste0(tempfile,"_step0.train.PRS.nosex"))}
  	if (file.exists(paste0(tempfile,"_step0.train.PRS.log"))) {file.remove(paste0(tempfile,"_step0.train.PRS.log"))}
  	if (file.exists(paste0(tempfile,"_step0.train.PRS.profile"))) {file.remove(paste0(tempfile,"_step0.train.PRS.profile"))}  
	write.table(out1$best.beta,file=paste0(tempfile,"_best.beta.txt"),row.names=F,quote=F,col.names=T)
	write.table(out1$best.params,file=paste0(tempfile,"_best.params.txt"),row.names=F,quote=F,col.names=T)

	return(out1)
}

