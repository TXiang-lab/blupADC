
#calcualte corrected phenotype
cal_corrected_phe<-function(target_trait_name=NULL,
							phe_col_names=NULL,
							genetic_effect_number=NULL,
							genetic_effect_name="Id",
							dmu_result_path=NULL,
							residual_average=TRUE,
							dmu_module="dmuai",
							phe_path=NULL,  #没有列名 
							phe_name=NULL,
							analysis_model=NULL, #在自定义模型的文件中，输出多个随机效应的值
							relationship_name=NULL,
							pedigree_path=NULL, #计算debv需要提供系谱
							pedigree_name=NULL, #计算debv需要提供系谱
							provided_prior_file_path=NULL,
							provided_prior_file_name=NULL,							
							cal_debv=FALSE,
							cal_reliability=FALSE,
							debv_id=NULL, #计算有基因型个体的debv
							output_ebv_path=NULL,
							output_ebv_name=NULL,
							return_result=FALSE
							){
library(data.table)
O0O0OOO0O0O0OOOOO0OO=length(target_trait_name)
phe=fread(paste0(phe_path,"/",phe_name),header=F,data.table=F)
colnames(phe)=phe_col_names

if(is.null(genetic_effect_number)){

if(analysis_model%in%c("User_define","GBLUP_A","GBLUP_AD")){genetic_effect_number=3}
if(analysis_model=="SSBLUP_A"|analysis_model=="PBLUP_A"){genetic_effect_number=4}
}
#针对 User_define 模型
if(analysis_model%in%c("User_define","GBLUP_A","GBLUP_AD")){

SOL=fread(paste0(dmu_result_path,"/SOL"),heade=F,data.table=F)
SOL=SOL[,c(1,4,2,5,8,9)];
colnames(SOL)=c("OOO0O0OOOOOOOOO0O0OO","Random_number","Trait","Id","EBV","SE")
tmp_SOL=SOL
OOOOOOO0O0OOO0OOOOOO=unique(SOL[(SOL$OOO0O0OOOOOOOOO0O0OO%in%genetic_effect_number)&(SOL$Random_number%in%(1:length(relationship_name))),c("Id")])
O0O0OOO0O0OOOOO0O0OO=matrix(NA,nrow=length(OOOOOOO0O0OOO0OOOOOO),ncol=1+(2+length(relationship_name))*O0O0OOO0O0O0OOOOO0OO)
O0O0OOO0O0OOOOO0O0OO=as.data.frame(O0O0OOO0O0OOOOO0O0OO,stringsAsFactors=F)
O0O0OOO0O0OOOOO0O0OO$V1=OOOOOOO0O0OOO0OOOOOO
colnames(O0O0OOO0O0OOOOO0O0OO)[1]="Id"
RESIDUAL=fread(paste0(dmu_result_path,"/RESIDUAL"),heade=F,data.table=F)

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(target_trait_name)){
residual=RESIDUAL[,c(1,1+OOO0OOO0O0OOO0O0OOO0OOOO+2*O0O0OOO0O0O0OOOOO0OO)];colnames(residual)=c("Pos","RESIDUAL")
residual[residual=="-9999"]=NA
# 重复记录的性状，对应多个残差，将残差取平均值
if(sum(duplicated(residual[,1]))>=1){
cat("Due to one animal have multiple records for a trait, the residuals will be calcualted on average \n")
attach(residual)
residual=aggregate(residual, by=list(Pos), FUN=mean)
residual=residual[,-1]
detach(residual)
}
residual[,1]=phe[as.numeric(residual[,"Pos"]),genetic_effect_name]
OOOOO0OOOOOOOOOOO0OO=match(O0O0OOO0O0OOOOO0O0OO[,1],residual[,1])

for(OOO0OOO0O0OOO0O0OOO0OOOOO in 1:length(relationship_name)){
additive=SOL[(SOL$OOO0O0OOOOOOOOO0O0OO%in%genetic_effect_number)&(SOL$Trait%in%OOO0OOO0O0OOO0O0OOO0OOOO)&(SOL$Random_number%in%OOO0OOO0O0OOO0O0OOO0OOOOO),c("Id","EBV")]
OOOOOOOOOOO0O0OOO0OO=match(O0O0OOO0O0OOOOO0O0OO[,1],additive[,1])
O0O0OOO0O0OOOOO0O0OO[,(OOO0OOO0O0OOO0O0OOO0OOOO-1)*(length(relationship_name)+2)+OOO0OOO0O0OOO0O0OOO0OOOOO+1]=additive[OOOOOOOOOOO0O0OOO0OO,"EBV"]
O0O0OOO0O0OOOOO0O0OO[,(OOO0OOO0O0OOO0O0OOO0OOOO-1)*(length(relationship_name)+2)+OOO0OOO0O0OOO0O0OOO0OOOOO+2]=residual[OOOOO0OOOOOOOOOOO0OO,"RESIDUAL"]
O0O0OOO0O0OOOOO0O0OO[,(OOO0OOO0O0OOO0O0OOO0OOOO-1)*(length(relationship_name)+2)+OOO0OOO0O0OOO0O0OOO0OOOOO+3]=O0O0OOO0O0OOOOO0O0OO[,(OOO0OOO0O0OOO0O0OOO0OOOO-1)*(length(relationship_name)+2)+OOO0OOO0O0OOO0O0OOO0OOOOO+2]+O0O0OOO0O0OOOOO0O0OO[,(OOO0OOO0O0OOO0O0OOO0OOOO-1)*(length(relationship_name)+2)+OOO0OOO0O0OOO0O0OOO0OOOOO+1]
colnames(O0O0OOO0O0OOOOO0O0OO)[(OOO0OOO0O0OOO0O0OOO0OOOO-1)*(length(relationship_name)+2)+OOO0OOO0O0OOO0O0OOO0OOOOO+1]=paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],"_R",OOO0OOO0O0OOO0O0OOO0OOOOO,"_Value")
colnames(O0O0OOO0O0OOOOO0O0OO)[(OOO0OOO0O0OOO0O0OOO0OOOO-1)*(length(relationship_name)+2)+OOO0OOO0O0OOO0O0OOO0OOOOO+2]=paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],"_Residual")
colnames(O0O0OOO0O0OOOOO0O0OO)[(OOO0OOO0O0OOO0O0OOO0OOOO-1)*(length(relationship_name)+2)+OOO0OOO0O0OOO0O0OOO0OOOOO+3]=paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],"_R",OOO0OOO0O0OOO0O0OOO0OOOOO,"_Value_Plus_Residual")

}
#O0O0OOO0O0OOOOO0O0OO[,OOO0OOO0O0OOO0O0OOO0OOOO*(length(relationship_name)+2)+1]= residual[OOOOO0OOOOOOOOOOO0OO,"RESIDUAL"]
#O0O0OOO0O0OOOOO0O0OO[,OOO0OOO0O0OOO0O0OOO0OOOO*(length(relationship_name)+2)+2]= O0O0OOO0O0OOOOO0O0OO[,OOO0OOO0O0OOO0O0OOO0OOOO*(length(relationship_name)+2)+0]+O0O0OOO0O0OOOOO0O0OO[,OOO0OOO0O0OOO0O0OOO0OOOO*(length(relationship_name)+2)+1]
#colnames(O0O0OOO0O0OOOOO0O0OO)[OOO0OOO0O0OOO0O0OOO0OOOO*(length(relationship_name)+2)+1]=paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],"_Residual")
#colnames(O0O0OOO0O0OOOOO0O0OO)[(OOO0OOO0O0OOO0O0OOO0OOOO-1)*(length(relationship_name)+2)+OOO0OOO0O0OOO0O0OOO0OOOOO+2]=paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],"_R",OOO0OOO0O0OOO0O0OOO0OOOOO,"_Value_Plus_Residual")
}


# 常规模型-计算EBV
}else{
SOL=fread(paste0(dmu_result_path,"/SOL"),heade=F,data.table=F)
SOL=SOL[,c(1,2,5,8,9)];
colnames(SOL)=c("OOO0O0OOOOOOOOO0O0OO","Trait","Id","EBV","SE")
tmp_SOL=SOL
OOOOOOO0O0OOO0OOOOOO=unique(SOL[(SOL$OOO0O0OOOOOOOOO0O0OO%in%genetic_effect_number),c("Id")])
O0O0OOO0O0OOOOO0O0OO=matrix(NA,nrow=length(OOOOOOO0O0OOO0OOOOOO),ncol=1+4*O0O0OOO0O0O0OOOOO0OO)
O0O0OOO0O0OOOOO0O0OO=as.data.frame(O0O0OOO0O0OOOOO0O0OO,stringsAsFactors=F)
O0O0OOO0O0OOOOO0O0OO$V1=OOOOOOO0O0OOO0OOOOOO
colnames(O0O0OOO0O0OOOOO0O0OO)[1]="Id"



RESIDUAL=fread(paste0(dmu_result_path,"/RESIDUAL"),heade=F,data.table=F)



for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(target_trait_name)){

additive=SOL[(SOL$OOO0O0OOOOOOOOO0O0OO%in%genetic_effect_number)&(SOL$Trait%in%OOO0OOO0O0OOO0O0OOO0OOOO),c("Id","EBV")]
residual=RESIDUAL[,c(1,1+OOO0OOO0O0OOO0O0OOO0OOOO+2*O0O0OOO0O0O0OOOOO0OO)];colnames(residual)=c("Pos","RESIDUAL")
residual[residual=="-9999"]=NA
# 重复记录的性状，对应多个残差，将残差取平均值
if(sum(duplicated(residual[,1]))>=1){
cat("Due to one animal have multiple records for a trait, the residuals will be calcualted on average \n")
attach(residual)
residual=aggregate(residual, by=list(Pos), FUN=mean)
residual=residual[,-1]
detach(residual)
}

residual[,1]=phe[as.numeric(residual[,"Pos"]),genetic_effect_name]
OOOOOOOOOOO0O0OOO0OO=match(O0O0OOO0O0OOOOO0O0OO[,1],additive[,1])
OOOOO0OOOOOOOOOOO0OO=match(O0O0OOO0O0OOOOO0O0OO[,1],residual[,1])

O0O0OOO0O0OOOOO0O0OO[,1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*4+1]=   additive[OOOOOOOOOOO0O0OOO0OO,"EBV"]
O0O0OOO0O0OOOOO0O0OO[,1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*4+2]= residual[OOOOO0OOOOOOOOOOO0OO,"RESIDUAL"]
O0O0OOO0O0OOOOO0O0OO[,1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*4+3]=O0O0OOO0O0OOOOO0O0OO[,1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*3+OOO0OOO0O0OOO0O0OOO0OOOO]+O0O0OOO0O0OOOOO0O0OO[,1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*3+OOO0OOO0O0OOO0O0OOO0OOOO+1]

colnames(O0O0OOO0O0OOOOO0O0OO)[1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*4+1]=paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],"_EBV")
colnames(O0O0OOO0O0OOOOO0O0OO)[1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*4+2]=paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],"_Res")
colnames(O0O0OOO0O0OOOOO0O0OO)[1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*4+3]=paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],"_EBV_Plus_Residual")
colnames(O0O0OOO0O0OOOOO0O0OO)[1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*4+4]=paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],"_dEBV")
}


#计算debv 
if(cal_debv==TRUE){
if(is.null(pedigree_path)&is.null(pedigree_name)){
stop("Please provide OOO0OOO0O0OOO0O0OOO0OOOOOO path and name for calculating dEBV!")
}
cat("Output dEBV......\n")
debv_result=OOO0O0OOO0OOOOOOOOO0(target_trait_name=target_trait_name,
				   genetic_effect_number=genetic_effect_number,
				   genetic_effect_name=genetic_effect_name,
				   dmu_result_path=dmu_result_path,
				   pedigree_path=pedigree_path,
				   pedigree_name=pedigree_name,							
				   provided_prior_file_path=provided_prior_file_path,
				   provided_prior_file_name=provided_prior_file_name,
				   debv_id=debv_id #计算有基因型个体的debv
				   )

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(target_trait_name)){
O0O0OOO0O0OOOOO0O0OO[,1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*4+OOO0OOO0O0OOO0O0OOO0OOOO+3]=as.numeric(debv_result[match(O0O0OOO0O0OOOOO0O0OO[,1],debv_result[,1]),OOO0OOO0O0OOO0O0OOO0OOOO+1])
}
original_O0O0OOO0O0OOOOO0O0OO=O0O0OOO0O0OOOOO0O0OO #带有NA的
}
}

if(cal_reliability==TRUE){

Sum_reliability=NULL

if(dmu_module=="dmu5"){stop("dmu5 module coundn't be used to estimate reliability,please use dmuai module!")}
if(!file.exists("PAROUT")){stop("PAROUT file is needed to estimate reliability,please use dmuai module!")}
cat("Output trait reliability result! \n")
O0O0O0OOOOO0OOOOO0O0<-data.frame(read.table("PAROUT",header=F))
se=tmp_SOL[tmp_SOL[,1]%in%genetic_effect_number,]
se=se[,c("Id","Trait","EBV","SE")]
se=se[!se[,4]%in%0,]  #去除为0的SE

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(target_trait_name)){
trait_se=se[se[,2]%in%OOO0OOO0O0OOO0O0OOO0OOOO,]
trait_var_Additive=O0O0O0OOOOO0OOOOO0O0[O0O0O0OOOOO0OOOOO0O0[,1]%in%1&O0O0O0OOOOO0OOOOO0O0[,2]%in%OOO0OOO0O0OOO0O0OOO0OOOO&O0O0O0OOOOO0OOOOO0O0[,3]%in%OOO0OOO0O0OOO0O0OOO0OOOO,4]
reliability=round(mean(1-(trait_se[,"SE"]^2)/trait_var_Additive),4)
Sum_reliability=rbind(Sum_reliability,c(target_trait_name[[OOO0OOO0O0OOO0O0OOO0OOOO]],reliability))
}

colnames(Sum_reliability)=c("Trait","Reliability")
write.table(Sum_reliability,paste("Trait_reliability.txt",sep=""),row.names=F,col.names=T,quote=F)
}


if(!is.null(output_ebv_path)&!is.null(output_ebv_name)){
setwd(output_ebv_path)
write.table(O0O0OOO0O0OOOOO0O0OO,paste0("colnames_",output_ebv_name,".txt"),quote=F,row.names=F,sep=" ")

O0O0OOO0O0OOOOO0O0OO[is.na(O0O0OOO0O0OOOOO0O0OO)]="-9999"
write.table(O0O0OOO0O0OOOOO0O0OO,paste0(output_ebv_name,".txt"),quote=F,row.names=F,col.names=F,sep=" ")
}

if(return_result==TRUE){return(original_O0O0OOO0O0OOOOO0O0OO)}
}


#calcualte corrected phenotype
OOO0O0OOO0OOOOOOOOO0_old<-function(target_trait_name=NULL,
				   genetic_effect_number=NULL,
				   genetic_effect_name="Id",
				   dmu_result_path=NULL,
				   pedigree_path=NULL,
				   pedigree_name=NULL,							
				   provided_prior_file_path=NULL,
				   provided_prior_file_name=NULL,
				   debv_id=NULL #计算有基因型个体的debv
				   ){

O0OOO0O0OOO0OOOOO0O0=calculate_reliability(return_individual_result=TRUE,
					    target_trait_name=target_trait_name,
					    genetic_effect_number=genetic_effect_number,
						provided_prior_file_path=provided_prior_file_path,
						provided_prior_file_name=provided_prior_file_name,
						dmu_result_path=dmu_result_path,
						debv_id=NULL)

OOOOOOO0OOOOO0O0OOO0=matrix(nrow=length(debv_id),ncol=length(target_trait_name)+1)
OOOOOOO0OOOOO0O0OOO0=as.data.frame(OOOOOOO0OOOOO0O0OOO0,stringsAsFactors=F)
OOOOOOO0OOOOO0O0OOO0[,1]=debv_id

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(target_trait_name)){

ebv_r2=O0OOO0O0OOO0OOOOO0O0$ebv_r2[[OOO0OOO0O0OOO0O0OOO0OOOO]]
h2=O0OOO0O0OOO0OOOOO0O0$h2[[OOO0OOO0O0OOO0O0OOO0OOOO]]
OOO0OOO0O0OOO0O0OOO0OOOOOO=fread(paste0(pedigree_path,"/",pedigree_name),header=F,data.table=F)

OOO0OOO0O0OOO0O0OOO0OOOOOO=OOO0OOO0O0OOO0O0OOO0OOOOOO[,1:3]
colnames(OOO0OOO0O0OOO0O0OOO0OOOOOO)=c("ID","Sire","Dam")
p.varSNP=0.5

debv_id=data.frame(debv_id,stringsAsFactors=F);

colnames(debv_id)[1]="Id"

OOOOO0O0OOOOO0OOO0OO <- merge(x=ebv_r2,y=OOO0OOO0O0OOO0O0OOO0OOOOOO,by.x=1,by.y=1)
  OOOOO0OOO0OOOOO0O0O0 <- merge(x=OOOOO0O0OOOOO0OOO0OO,y=debv_id,by.x=1,by.y=1)
  colnames(x=OOOOO0OOO0OOOOO0O0O0) <- c("ID","ID_EBV","ID_Rel","SireID","DamID") 
  Rel.sire <- merge(x=OOOOO0O0OOOOO0OOO0OO,y=OOOOO0OOO0OOOOO0O0O0,by.x=1,by.y=4)
  Rel.sire <- Rel.sire[,-4:-5]  
  colnames(Rel.sire)[1:3] <- c('SireID','Sire_EBV','Sire_R2')
  Rel.dam <- merge(x=OOOOO0O0OOOOO0OOO0OO,y=Rel.sire,by.x=1,by.y=7)
  Rel.dam <- Rel.dam[,-4:-5]  
  colnames(Rel.dam)[1:3] <- c('DamID','Dam_EBV','Dam_R2')
  data <- Rel.dam[,c('ID','ID_EBV','ID_Rel','SireID','Sire_EBV','Sire_R2','DamID','Dam_EBV','Dam_R2')]
  dEBV <- data[,-c(1,4,7)]
  dEBV$h2 <- h2
  dEBV$p.varSNP <- 1-p.varSNP
  dEBV=as.matrix(dEBV)
  dEBV=apply(dEBV,2,as.numeric)
  Deregress <- t(apply(X=dEBV,MARGIN=1,FUN=OOOOOOO0O0OOOOOOOOOO))
  Deregress <- cbind.data.frame(Id=data[,1],Deregress)

  OOOOOOO0OOOOO0O0OOO0[na.omit(match(OOOOOOO0OOOOO0O0OOO0[,1],Deregress[,1])),OOO0OOO0O0OOO0O0OOO0OOOO+1]=Deregress[na.omit(match(OOOOOOO0OOOOO0O0OOO0[,1],Deregress[,1])),"dEBV"]
}

colnames(OOOOOOO0OOOOO0O0OOO0)=c("Id",paste0(target_trait_name,"_dEBV"))
return(OOOOOOO0OOOOO0O0OOO0)

}

#calcualte corrected phenotype
OOO0O0OOO0OOOOOOOOO0<-function(target_trait_name=NULL,
				   genetic_effect_number=NULL,
				   genetic_effect_name="Id",
				   dmu_result_path=NULL,
				   pedigree_path=NULL,
				   pedigree_name=NULL,							
				   provided_prior_file_path=NULL,
				   provided_prior_file_name=NULL,
				   p_var=0.5,
				   debv_id=NULL #计算有基因型个体的debv
				   ){

O0OOO0O0OOO0OOOOO0O0=calculate_reliability(return_individual_result=TRUE,
					    target_trait_name=target_trait_name,
					    genetic_effect_number=genetic_effect_number,
						provided_prior_file_path=provided_prior_file_path,
						provided_prior_file_name=provided_prior_file_name,
						dmu_result_path=dmu_result_path,
						debv_id=NULL)

OOO0OOO0O0OOO0O0OOO0OOOOOO=fread(paste0(pedigree_path,"/",pedigree_name),header=F,data.table=F)
OOO0OOO0O0OOO0O0OOO0OOOOOO=OOO0OOO0O0OOO0O0OOO0OOOOOO[,1:3]
colnames(OOO0OOO0O0OOO0O0OOO0OOOOOO)=c("ID","Sire","Dam")
p_var=p_var


ebv_r2=data.frame(O0OOO0O0OOO0OOOOO0O0$ebv_r2[[1]],stringsAsFactors=F)

h2=O0OOO0O0OOO0OOOOO0O0$h2[[1]]

debv=debv_method1(ebv=ebv_r2,ped=OOO0OOO0O0OOO0O0OOO0OOOOOO,h2=h2,p_var=p_var,debv_id=debv_id)[,c("Id","dEBV")]

if(length(target_trait_name)>=2){
for(OOO0OOO0O0OOO0O0OOO0OOOO in 2:length(target_trait_name)){

ebv_r2=data.frame(O0OOO0O0OOO0OOOOO0O0$ebv_r2[[OOO0OOO0O0OOO0O0OOO0OOOO]],stringsAsFactors=F)
h2=O0OOO0O0OOO0OOOOO0O0$h2[[OOO0OOO0O0OOO0O0OOO0OOOO]]

debv_trait=debv_method1(ebv=ebv_r2,ped=OOO0OOO0O0OOO0O0OOO0OOOOOO,h2=h2,p_var=p_var,debv_id=debv_id)[,c("Id","dEBV")]

debv=cbind(debv,debv_trait[match(debv[,1],debv_trait[,1]),2])

}
}


return(debv)

}

debv_method1<-function(ebv,ped,h2=1,p_var=0.5,debv_id=NULL){
 OOOOOOO0O0OOOOOOOOOO <- function (ebv_mat){
    lambda = (1 - ebv_mat[7])/ebv_mat[7]
    PA = (ebv_mat[3] + ebv_mat[5])/2
    rPA = (ebv_mat[4] + ebv_mat[6])/4
    alpha = 1/(0.5 - rPA)
    delta = (0.5 - rPA)/(1 - ebv_mat[2])
    ZpZPA = lambda*(0.5*alpha - 4) + 0.5*lambda*sqrt(alpha^2 +16/delta)
    ZpZi = delta*ZpZPA + 2*lambda*(2*delta - 1)
    LHS = rbind(cbind(ZpZPA + 4*lambda, -2*lambda),cbind(-2*lambda, ZpZi + 2*lambda))
    #L1 = solve(LHS)
    RHS = LHS %*% c(PA, ebv_mat[1])
    drgi = RHS[2]/ZpZi
    rdrg = 1 - lambda/(ZpZi + lambda)
    we = (1 - ebv_mat[7])/((ebv_mat[8] + (1 - rdrg)/rdrg)*ebv_mat[7])
    ret = c(ebv_mat[1],ebv_mat[2],round(drgi,3),round(rdrg,3),round(we,3))
    names(ret) <- c("EBV", "r2EBV", "dEBV","r2dEBV","weight")
    return(ret)
  }
ped=ped[,1:3]
colnames(ebv)=c("Id","EBV","r2")
colnames(ped)=c("Id","Sire","Dam")

ped$Id_EBV=ebv[match(ped$Id,ebv$Id),"EBV"]
ped$Id_r2=ebv[match(ped$Id,ebv$Id),"r2"]
ped$Sire_EBV=ebv[match(ped$Sire,ebv$Id),"EBV"]
ped$Sire_r2=ebv[match(ped$Sire,ebv$Id),"r2"]
ped$Dam_EBV=ebv[match(ped$Dam,ebv$Id),"EBV"]
ped$Dam_r2=ebv[match(ped$Dam,ebv$Id),"r2"]


if(is.null(debv_id)){
debv_id=ebv[,1]
}

ped=ped[ped[,1]%in%debv_id,]

#过滤父母无ebv的个体
ped=ped[!is.na(ped$Id_EBV)&!is.na(ped$Dam_EBV)&!is.na(ped$Sire_EBV),]
ped=ped[!is.na(ped$Id_r2)&!is.na(ped$Sire_r2)&!is.na(ped$Dam_r2),]

ped$h2=h2
ped$p_var=p_var

Deregress <- t(apply(X=as.matrix(ped[,-c(1:3)]),MARGIN=1,FUN=OOOOOOO0O0OOOOOOOOOO))
ped=cbind(ped[,1:3],Deregress)  
  


return(ped)
}

calculate_reliability<-function(
								target_trait_name,
								genetic_effect_number=4,
								dmu_result_path=NULL,
								provided_prior_file_path=NULL,
								provided_prior_file_name=NULL,
								return_individual_result=FALSE,
								debv_id=NULL #  在计算原有所有个体reliability的基础上，进一步计算  有基因型个体的reliability 和 没有基因型个体的reliability
								){

select_id_r2=NULL								
unselect_id_r2=NULL
r2=NULL								
trait_ebv_r2=NULL
h2=NULL

setwd(dmu_result_path)

SOL=fread("SOL",heade=F,data.table=F)
SOL=SOL[,c(1,2,5,8,9)];
colnames(SOL)=c("OOO0O0OOOOOOOOO0O0OO","Trait","Id","EBV","SE")
SOL=SOL[SOL$OOO0O0OOOOOOOOO0O0OO%in%genetic_effect_number,]

if(is.null(provided_prior_file_path)&is.null(provided_prior_file_name)){

O0O0O0OOOOO0OOOOO0O0<-data.frame(read.table("PAROUT",header=F),stringsAsFactors=F)
}else{

O0O0O0OOOOO0OOOOO0O0<-data.frame(read.table(paste0(provided_prior_file_path,"/",provided_prior_file_name),header=F),stringsAsFactors=F)}

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(target_trait_name)){

trait_se=SOL[SOL$Trait%in%OOO0OOO0O0OOO0O0OOO0OOOO,]
trait_var_Additive=as.numeric(O0O0O0OOOOO0OOOOO0O0[(O0O0O0OOOOO0OOOOO0O0[,1]%in%1)&(O0O0O0OOOOO0OOOOO0O0[,2]%in%OOO0OOO0O0OOO0O0OOO0OOOO)&(O0O0O0OOOOO0OOOOO0O0[,3]%in%OOO0OOO0O0OOO0O0OOO0OOOO),4])

reliability=mean(1-(trait_se[,"SE"]^2)/trait_var_Additive)
r2=c(r2,list(reliability))
write.table(reliability,paste0("The_reliability_of_",target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],".txt"),row.names=F,col.names=F,quote=F)

ebv_r2=cbind(SOL[,"Id"],SOL[,"EBV"],(1-(trait_se[,"SE"]^2)/trait_var_Additive))
h2=c(h2,list(trait_var_Additive/sum(O0O0O0OOOOO0OOOOO0O0[,4])))
colnames(ebv_r2)=c("Id","EBV","Reliability")
trait_ebv_r2=c(trait_ebv_r2,list(ebv_r2))

if(!is.null(debv_id)){
selected_id_se=trait_se[trait_se[,1]%in%debv_id,]
unselected_id_se=trait_se[!trait_se[,1]%in%debv_id,]
selected_id_reliability=sum(1-(selected_id_se[,"SE"]^2)/trait_var_Additive)/nrow(selected_id_se)
unselected_id_reliability=sum(1-(unselected_id_se[,"SE"]^2)/trait_var_Additive)/nrow(unselected_id_se)


select_id_r2=c(select_id_r2,list(selected_id_reliability))
unselect_id_r2=c(unselect_id_r2,list(unselected_id_reliability))
}}

if(return_individual_result){return(list(ebv_r2=trait_ebv_r2,h2=h2,r2=r2,
                                                 select_id_r2=select_id_r2,unselect_id_r2=unselect_id_r2,
										trait=target_trait_name))}
}

		#
geno_check<-function(
		input_data_type=NULL, #"Plink" , "Hapmap" , "VCF" , "BLUPF90","Numeric"  "Haplotype"
		input_data_path=NULL,
		input_data_name=NULL,  #hapmap 为全名，vcf 和 plink 默认不包括后缀，为 name.vcf name.ped name.map 
		input_data_hmp=NULL,
		input_data_plink_ped=NULL,
		input_data_plink_map=NULL,
		input_data_blupf90=NULL,
		input_data_numeric=NULL,
		input_data_vcf=NULL,
		input_data_haplotype_hap=NULL,
		input_data_haplotype_map=NULL,
		input_data_haplotype_sample=NULL,
		input_data_numeric_map=NULL,
		input_data_blupf90_map=NULL,
		miss_base=NULL,
		keep_inds_set=NULL, #字符串向量
		keep_snps_set=NULL, #字符串向量
		keep_chroms_set=NULL, #数值向量
		chr_set=NULL,
		output_data_name="genotype_conversion",
		output_data_type=NULL,  # "Plink" , "Hapmap" , "VCF"  "Blupf90"
		output_data_path=NULL,
		return_result=TRUE,	
		cpu_cores=1, #并行计算所有的cpu数目
          selected_snps=1000, # 默认随机从芯片中抽取2K出来计算overlap比例
		overlap_name=TRUE,
		duplication_check=FALSE,
		duplication_threshold=0.95,#overlap比例超过阈值就认为是重复个体
		breed_check=FALSE,
		breed_record=NULL
			){
library(data.table)
input=O0OOO0O0OOO0O0O0OOOO(input_data_type=input_data_type,input_data_hmp=input_data_hmp,input_data_plink_ped=input_data_plink_ped,
                                          input_data_plink_map=input_data_plink_map,input_data_blupf90=input_data_blupf90,
								    input_data_numeric=input_data_numeric,input_data_vcf=input_data_vcf,
									input_data_haplotype_hap=input_data_haplotype_hap,input_data_haplotype_map=input_data_haplotype_map,
									input_data_haplotype_sample=input_data_haplotype_sample,
								     miss_base=miss_base)

input_data_type=input$input_data_type
miss_base=input$miss_base

sum_data=geno_format(   
          input_data_hmp=input_data_hmp,      
          input_data_name=input_data_name,  
	     input_data_path=input_data_path,		  
          input_data_type=input_data_type,  
	     input_data_plink_map=input_data_plink_map, 
	     input_data_plink_ped=input_data_plink_ped, 
	     input_data_blupf90=input_data_blupf90,
		 input_data_numeric_map=input_data_numeric_map,
		input_data_blupf90_map=input_data_blupf90_map,
		cpu_cores=cpu_cores,
		miss_base=miss_base,
		output_data_type="Numeric",
		output_data_name=NULL,
		output_data_path=NULL,		
		return_result=TRUE)		
		
data_numeric=sum_data$numeric
IND_geno=rownames(data_numeric)

OOOOOOOOOOOOO0O0O0O0=NULL
if(duplication_check==TRUE){
#reduce snps number in detecting overlap for saving time
if(ncol(data_numeric)>selected_snps){
set.seed(19980204)
selected_line=sample(1:ncol(data_numeric),selected_snps,replace=FALSE) # 随机选择1000行
data_numeric_subset=data_numeric[,selected_line]
}else{
data_numeric_subset=data_numeric
}	

message("Detecting overlap genotype......")	
OOOOOOOOOOOOO0O0O0O0=numeric_overlap_cpp(numeric=data_numeric_subset,dis_progress=TRUE,overlap_threshold=duplication_threshold,cpu_cores=cpu_cores)
OOOOOOOOOOOOO0O0O0O0=OOOOOOOOOOOOO0O0O0O0[OOOOOOOOOOOOO0O0O0O0[,3]>duplication_threshold,]

#根据OOOOOOOOOOOOO0O0O0O0的初步结果，利用所有SNP数据,再重新进行overlap计算
if(nrow(OOOOOOOOOOOOO0O0O0O0)>500){stop("Plesase select appropriate numbers of selected snps in detecting overlap ")}
if(nrow(OOOOOOOOOOOOO0O0O0O0)>=1){
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:nrow(OOOOOOOOOOOOO0O0O0O0)){
OOOOOOOOOOOOO0O0O0O0[OOO0OOO0O0OOO0O0OOO0OOOO,3]=sum(data_numeric[OOOOOOOOOOOOO0O0O0O0[OOO0OOO0O0OOO0O0OOO0OOOO,1],]==data_numeric[OOOOOOOOOOOOO0O0O0O0[OOO0OOO0O0OOO0O0OOO0OOOO,2],])/ncol(data_numeric)
}}

OOOOOOOOOOOOO0O0O0O0=OOOOOOOOOOOOO0O0O0O0[as.numeric(OOOOOOOOOOOOO0O0O0O0[,3])>=duplication_threshold,]

if(nrow(OOOOOOOOOOOOO0O0O0O0)>=1){
message(paste0("Found ",nrow(OOOOOOOOOOOOO0O0O0O0)," duplicate genotypes in your chip data"))
IND1=IND_geno[OOOOOOOOOOOOO0O0O0O0[,1]]
IND2=IND_geno[OOOOOOOOOOOOO0O0O0O0[,2]]
OOOOOOOOOOOOO0O0O0O0[,1]=IND1
OOOOOOOOOOOOO0O0O0O0[,2]=IND2
if(!is.null(output_data_path)){
setwd(output_data_path)
write.table(OOOOOOOOOOOOO0O0O0O0,"Duplicated_Genotype_Percentage.txt",quote=F,row.names=F,,col.names=F,sep="\t")}
overlap_IND=unique(c(OOOOOOOOOOOOO0O0O0O0[,1],OOOOOOOOOOOOO0O0O0O0[,2]))

keep_inds_set=IND_geno[!IND_geno%in%overlap_IND]

sum_data=geno_format(   
          input_data_hmp=input_data_hmp,      
          input_data_name=input_data_name,  
	     input_data_path=input_data_path,		  
          input_data_type=input_data_type,  
	     input_data_plink_map=input_data_plink_map, 
	     input_data_plink_ped=input_data_plink_ped, 
	     input_data_blupf90=input_data_blupf90,
		input_data_numeric_map=input_data_numeric_map,
		input_data_blupf90_map=input_data_blupf90_map,
		cpu_cores=cpu_cores,
		miss_base=miss_base,
		output_data_type=output_data_type,
		output_data_name=output_data_name,
		output_data_path=output_data_path,			
		return_result=TRUE)		
cat(paste0("Finally, there are total ",length(keep_inds_set)," individuals in your data! \n"))
		
}else{
message(paste0("Found no duplicate genotypes in your chip data"))
}

}

pca_outlier=NULL
if(breed_check==TRUE){
if(is.null(breed_record)){stop("Please provide the information of breed")}

pca_outlier=detect_pca_plot(input_data_numeric=data_numeric,
                                  Breed=breed_record)
							
}


if(return_result==TRUE){
return(list(sum_data=sum_data,duplicate_genotype=OOOOOOOOOOOOO0O0O0O0,breed=pca_outlier))
}
}




#读取完所有数据再进行合并
genotype_data_combine<-function(
               input_data_type=NULL,  # 可以使 vcf ,plink ,hapmap 中的任何一个，也可以是三种类型的混合
               input_data_path=NULL,     # 可以是多个路径，
               return_result=FALSE,               #在r中返回结果
               combine_type="intersect",           #多个基因型文件合并的类型，intersect, union 				
			miss_base=NULL,
			duplication_threshold=0.95,
			selected_snps=1000,
			cpu_cores=1, #并行计算所有的cpu数目
			overlap_id=TRUE,  #检测基因型中是否有重复个体名称
			duplication_check=TRUE, #检测基因型中是否有重复基因型			 
			output_data_name="combine_result",		
               output_data_type=NULL,   #输出的文件类型
               output_data_path=NULL			
                ){
library(data.table)
library(purrr)
#程序自动判断输入类型
input=O0OOO0O0OOO0O0O0OOOO(input_data_type=input_data_type,input_data_hmp=input_data_hmp,input_data_plink_ped=input_data_plink_ped,
                                          input_data_plink_map=input_data_plink_map,input_data_blupf90=input_data_blupf90,
								     miss_base=miss_base)

input_data_type=input$input_data_type
miss_base=input$miss_base

			
library(data.table);library(purrr)
data_ped_list=NULL;data_map_list=NULL;data_hmp_list=NULL
#获得Plink 和  hapmap所有的snps信息

data_map=NULL
#统计Plink信息，找到所有SNP的交集或并集
for(path in input_data_path){
setwd(path)
input_plink_data_name=try(list.files(getwd(),pattern = '.ped$'),silent=T)
#####读取plink数据
if(length(input_plink_data_name)!=0){
input_data_name=unique(unlist(strsplit(input_plink_data_name,split = ".ped")))
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(input_data_name)){
data_map=fread(paste0(path,"/",input_data_name[OOO0OOO0O0OOO0O0OOO0OOOO],".map"),data.table=F,header=F)
data_map_list=c(data_map_list,list(data_map))	
}}}

if(length(data_map_list)>=1){
if(combine_type=="intersect"){
data_map= reduce(data_map_list,inner_join,by="V2")
data_map=data_map[,1:4]
}else{
data_map= reduce(data_map_list,full_join,by="V2")
data_map=data_map[,1:4]
}
}

plink_iter=1;hmp_iter=1
for(path in input_data_path){
setwd(path)
input_plink_data_name=try(list.files(getwd(),pattern = '.ped$'),silent=T)
input_hapmap_data_name=try(list.files(getwd(),pattern = '.hmp.txt$'),silent=T)
#####读取plink数据
if(length(input_plink_data_name)!=0){
input_data_name=unique(unlist(strsplit(input_plink_data_name,split = ".ped")))
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(input_data_name)){
cat("Reading plink data......\n")
data_ped_subset=fread(paste0(path,"/",input_data_name[OOO0OOO0O0OOO0O0OOO0OOOO],".ped"),data.table=F,header=F,nThread=cpu_cores)
data_map_subset=fread(paste0(path,"/",input_data_name[OOO0OOO0O0OOO0O0OOO0OOOO],".map"),data.table=F,header=F)
cat(paste0("Plink data file path",plink_iter,": ",path," ;\nPlink data file name",plink_iter,": ",input_data_name[OOO0OOO0O0OOO0O0OOO0OOOO],"\n"))
plink_iter=plink_iter+1
if(combine_type=="intersect"){
pos=plink_index(match(data_map[,2],data_map_subset[,2]))[,1]
data_ped_list=c(data_ped_list,list(data_ped_subset[,pos]))
rm(data_ped_subset);gc();
}else if(combine_type=="union"){
pos=plink_index(match(data_map_subset[,2],data_map[,2]))[,1]
data_ped=data.frame(matrix(0,nrow=nrow(data_ped_subset),ncol=nrow(data_map)*2+6))
data_ped[,pos]=data_ped_subset
data_ped_list=c(data_ped_list,list(data_ped))
rm(data_ped_subset,data_ped);gc();
}
}}

#####读取hapmap数据
if(length(input_hapmap_data_name)!=0){
input_data_name=unique(unlist(strsplit(input_hapmap_data_name,split = ".hmp.txt")))
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(input_data_name)){
cat("Reading hapmap data......\n")
data_hmp=fread(paste0(path,"/",input_data_name[OOO0OOO0O0OOO0O0OOO0OOOO],".hmp.txt"),data.table=F,header=T)
cat(paste0("Hapmap data file path",hmp_iter,": ",path," ;\nHapmap data file name",hmp_iter,": ",input_data_name[OOO0OOO0O0OOO0O0OOO0OOOO],"\n"))
colnames(data_hmp)[1]="rs#"
hmp_iter=hmp_iter+1
if(!is.null(data_hmp_list)){data_hmp=data_hmp[,-c(2:11)]}
data_hmp_list=c(data_hmp_list,list(data_hmp))
}}
}
data_hmp=NULL;data_ped=NULL

#合并 plink数据 
if(length(data_map_list)>=1){
data_ped=do.call(rbind,data_ped_list)
rm(data_ped_list);gc();
data_map=data_map
cat("Combine plink data......\n")
}

#合并 hapmap数据 
if(length(data_hmp_list)>=1){
cat("Combine hapmap data......")
if(combine_type=="intersect"){
#data_hmp=data_hmp_list%>%reduce(inner_join,by="rs#")
data_hmp=reduce(data_hmp_list,inner_join,by="rs#")
}else if(combine_type=="union"){
#data_hmp=data_hmp_list%>%reduce(full_join,by="rs#")
data_hmp=reduce(data_hmp_list,full_join,by="rs#")
}}

#输出结果


#如果overlap为TRUE,那么在检测 overlap之前不输出任何结果

if(duplication_check==TRUE){
before_overlap_out_put_path=NULL;
before_overlap_out_put_name=NULL;
before_overlap_return_result=TRUE;
}else{
before_overlap_out_put_path=output_data_path
before_overlap_out_put_name=output_data_name
before_overlap_return_result=FALSE
}

if(!is.null(output_data_path)){
if(!file.exists(output_data_path)){
dir.create(output_data_path,recursive=T)
cat("Output data path is not exist,software will constructs this path  automatically \n")
}
}

#只有plink数据，没有hapmap数据
if(is.null(data_hmp)&!is.null(data_ped)){
						 
if(duplication_check==TRUE){




sum_data=genotype_data_check(
		input_data_type="Plink", #"Plink" , "Hapmap" , "VCF" , "Blupf90"
		input_data_plink_ped=data_ped,
		input_data_plink_map=data_map,
		output_data_path=output_data_path,
		output_data_name=output_data_name,
		return_result=return_result,
		output_data_type=output_data_type,
		cpu_cores=cpu_cores, #并行计算所有的cpu数目
		miss_base=miss_base,
          selected_snps=selected_snps, # 默认随机从芯片中抽取2K出来计算overlap比例
		overlap_name=TRUE,
		duplication_check=TRUE,
		duplication_threshold=duplication_threshold#overlap比例超过阈值就认为是重复个体
		)
		}else{

sum_data=geno_format(   
                           input_data_plink_ped=data_ped,
	                       input_data_plink_map=data_map,
	                     input_data_type="Plink",
					 cpu_cores=cpu_cores,
					 miss_base=miss_base,
	                     output_data_type=output_data_type,
				      output_data_name=output_data_name,
					 output_data_path=output_data_path,
	                     return_result=return_result)
		
		}
		
		
		}
		

#只有hapmap数据，没有plink数据
if(!is.null(data_hmp)&is.null(data_ped)){						 
sum_data=geno_format(   
                           input_data_hmp=data_hmp,
	                     input_data_type="Hapmap",
					 cpu_cores=cpu_cores,
	                     output_data_type=output_data_type,
				      output_data_name=before_overlap_out_put_name,
					 output_data_path=before_overlap_out_put_path,
	                     return_result=before_overlap_return_result)
						 
if(duplication_check==TRUE){
sum_data=genotype_data_check(
		input_data_type="Hapmap", #"Plink" , "Hapmap" , "VCF" , "Blupf90"
		input_data_hmp=sum_data$hmp,
		output_data_path=output_data_path,
		output_data_name=output_data_name,
		return_result=return_result,
		output_data_type=before_overlap_output_data_type,		
		cpu_cores=cpu_cores, #并行计算所有的cpu数目
          selected_snps=selected_snps, # 默认随机从芯片中抽取2K出来计算overlap比例
		overlap_name=TRUE,
		duplication_check=TRUE,
		duplication_threshold=duplication_threshold#overlap比例超过阈值就认为是重复个体
		)}}		

#既有hapmap数据，又有plink数据
if(!is.null(data_hmp)&!is.null(data_ped)){	

data_hmp_conversion=geno_format(   
                           input_data_plink_ped=data_ped,
	                       input_data_plink_map=data_map,
	                     input_data_type="Plink",
					  cpu_cores=cpu_cores,
	                     output_data_type="Hapmap",	 
	                     return_result=TRUE)
data_hmp_conversion=data_hmp_conversion$hmp
#data_hmp=list(data_hmp,data_hmp_conversion)%>%reduce(inner_join,by="rs#")
data_hmp=reduce(list(data_hmp,data_hmp_conversion),inner_join,by="rs#")
				 
sum_data=geno_format(   
                           input_data_hmp=data_hmp,
	                     input_data_type="Hapmap",
					 cpu_cores=cpu_cores,
	                     output_data_type=output_data_type,
				      output_data_name=before_overlap_out_put_name,
					 output_data_path=before_overlap_out_put_path,
	                     return_result=before_overlap_return_result)
						 
if(duplication_check==TRUE){
sum_data=genotype_data_check(
		input_data_type="Hapmap", #"Plink" , "Hapmap" , "VCF" , "Blupf90"
		input_data_hmp=sum_data$hmp,
		output_data_path=output_data_path,
		output_data_name=output_data_name,
		return_result=return_result,
		output_data_type=before_overlap_output_data_type,		
		cpu_cores=cpu_cores, #并行计算所有的cpu数目
          selected_snps=selected_snps, # 默认随机从芯片中抽取2K出来计算overlap比例
		overlap_name=TRUE,
		duplication_check=TRUE,
		duplication_threshold=duplication_threshold#overlap比例超过阈值就认为是重复个体
		)}}

if(return_result==TRUE){return(sum_data)}
		
}


geno_qc_impute<-function(
      input_data_hmp=NULL,       # for a given hapmap not read from  R
      input_data_name=NULL,      # get data from the work directory which don't need to read by R, but data must put in the input_data_path
      input_data_type=NULL,  # there are three type of data, include  hapmap 、plink 、 vcf 、Blupf90 the default is hapmap
	 input_data_plink_map=NULL,  # for a given map not read from  R
	 input_data_plink_ped=NULL,  # for a given ped not read from  R
	 input_data_blupf90=NULL,
	 input_data_path=NULL,     	# the directory of genotype data 
	 input_data_vcf=NULL,
	 input_data_numeric=NULL,
	 input_data_haplotype_hap=NULL,
	 input_data_haplotype_map=NULL,
	 input_data_haplotype_sample=NULL,
	input_data_numeric_map=NULL,
	input_data_blupf90_map=NULL,	 
	 bigmemory_cal=FALSE,
	 bigmemory_data_type="integer",
	 bigmemory_data_path=getwd(),
	 bigmemory_data_name="blupADC",
	 phased_symbol="|",
	 phased_genotype=FALSE,
	 haplotype_window_nSNP=NULL,
	 haplotype_window_kb=NULL,
	 haplotype_window_block=NULL,	 
	 output_data_path=NULL,
	 output_data_name=NULL,
	 output_data_type="Plink",
	 data_analysis_method="QC_Imputation",   # there are three data_analysis_method , Plink 、 Beagle 、 QC_Imputation	 
	 cpu_cores=1, 
      miss_base="0",          # the default missing value of genotype 
	 qc_snp_rate=0.1,
	 qc_ind_rate=0.1,
	 qc_maf=0.05,
	 qc_hwe=0.0000001,
	 extra_parameter=NULL,
	 plink_software_path=ifelse(as.character(Sys.info()["sysname"])=="Linux",system.file("extdata/bin", package = "blupADC"),system.file("extdata/bin_windows", package = "blupADC")),
	 plink_software_name="plink", 
	 chr_set=NULL, # 识别多个染色体-plink, 最多可识别多个95个 	 
	 keep_inds_set=NULL,   #字符串向量
	 keep_snps_set=NULL,   #字符串向量
	 keep_chroms_set=NULL, #数值向量
	 beagle_software_path=ifelse(as.character(Sys.info()["sysname"])=="Linux",system.file("extdata/bin", package = "blupADC"),system.file("extdata/bin_windows", package = "blupADC")),
	 beagle_software_name="beagle.5.2.jar",
	 beagle_ped_path=NULL,  
	 beagle_ped_name=NULL,
	 beagle_ref_data_path=NULL, #Using Reference For Inputaion
	 beagle_ref_data_name=NULL,
	 Java_Space=NULL,	 
	 sigle_duo_trio=NULL  # becuase if you want to use OOO0OOO0O0OOO0O0OOO0OOOOOO, trio and duos would make your computation more slowly , the standard_format is  sigle_duo_trio=c(" singlescale=20  duoscale=20  trioscale=20 ")	
){  # determin whether output the hmp.file 
	#author:  Quanshun Mei 
	#date: 2019.09.09
	#purpose: conduct QC and Imputation 
	#data: the standard genotype file, include  .vcf 、.ped .map 、 .hmp ，
	#missing: the symbols for the missing genotype
	#SNPcall, INDcall, maf, hwe: the conditions that will be used during using Plink
	#path默认的格式为  ：  eg： beaglepath_4.1="/home/qsmei/beagle_4.1"
	# the version of beagle include beagle.4.0.jar 、beagle.4.1.jar 、beagle.5.1.jar，   only beagle.4.0.jar can add OOO0OOO0O0OOO0O0OOO0OOOOOO into imputation 

if(is.null(output_data_name)){stop("Please specify your output data name!!!")}
if(is.null(output_data_path)){stop("Please specify your output data path!!!")}
if(!is.null(chr_set)){chr_set=paste0(" --chr-set ",chr_set," ")}

library(data.table)

input=O0OOO0O0OOO0O0O0OOOO(input_data_type=input_data_type,input_data_hmp=input_data_hmp,input_data_plink_ped=input_data_plink_ped,
                                          input_data_plink_map=input_data_plink_map,input_data_blupf90=input_data_blupf90,
								    input_data_numeric=input_data_numeric,input_data_vcf=input_data_vcf,
									input_data_haplotype_hap=input_data_haplotype_hap,input_data_haplotype_map=input_data_haplotype_map,
									input_data_haplotype_sample=input_data_haplotype_sample,
								     miss_base=miss_base)

input_data_type=input$input_data_type
miss_base=input$miss_base

#当提供输入数据为vcf格式，且不需要质控的时候，直接开始填充
if(!(input_data_type=="VCF"&data_analysis_method=="Imputation")){


#如果输入的基因型文件为本地Plink格式的基因型文件，直接拷贝基因型文件即可
if(input_data_type=="Plink"&!is.null(input_data_path)){
cat("Start Read The Plink Format Genotype Data \n")
system(paste0("cp -r ",input_data_path,"/",input_data_name,".ped   ",output_data_path,"/temp_plink.ped"))
system(paste0("cp -r ",input_data_path,"/",input_data_name,".map   ",output_data_path,"/temp_plink.map"))

input_data_ped=fread(paste0(output_data_path,"/temp_plink.ped"),select=1,data.table=F)
input_data_map=fread(paste0(output_data_path,"/temp_plink.map"),data.table=F)
cat("Complete Read The Plink Format Genotype Data \n")

}else if(input_data_type=="Plink"&!is.null(input_data_plink_ped)){

fwrite(data.table(input_data_plink_ped),paste0(output_data_path,"/temp_plink.ped"),quote=F,row.names=F,col.names=F,sep=" ")
fwrite(data.table(input_data_plink_map),paste0(output_data_path,"/temp_plink.map"),quote=F,row.names=F,col.names=F,sep=" ")

}else{

#将输入数据转换为 Plink格式
sum_data=geno_format(   
          input_data_hmp=input_data_hmp,      
          input_data_name=input_data_name,  
	     input_data_path=input_data_path,		  
          input_data_type=input_data_type,  
	     input_data_plink_map=input_data_plink_map, 
	     input_data_plink_ped=input_data_plink_ped, 
	     input_data_blupf90=input_data_blupf90,
		input_data_numeric=input_data_numeric,
		input_data_vcf=input_data_vcf,
	input_data_numeric_map=input_data_numeric_map,
	input_data_blupf90_map=input_data_blupf90_map,		
		miss_base=miss_base,
		cpu_cores=cpu_cores,
		output_data_type="Plink",
		output_data_name="temp_plink",
		output_data_path=output_data_path,	
		return_result=TRUE)	

input_data_ped=sum_data$ped
input_data_map=sum_data$map 
rm(sum_data);gc();

}

IND_n=nrow(input_data_ped)
SNP_n=nrow(input_data_map)
rm(input_data_ped,input_data_map);gc();


setwd(output_data_path)
output_data_path=getwd()

OOO0OOO0O0OOOOO0OOO0=as.character(Sys.info()["sysname"])

if(OOO0OOO0O0OOOOO0OOO0=="Linux"){

#针对软件路径含有空格的情况，如果路径存在空格，将软件copy 到输出路径下
if(length(grep(" ",plink_software_path))>0){
#cat("Attention: Phe path include space!")
file.copy(paste0(plink_software_path,"/",plink_software_name),getwd())
plink_software_path=getwd()
}

}else{

if(length(grep(" ",plink_software_path))>0){
#cat("Attention: Phe path include space!")
file.copy(paste0(plink_software_path,"/",plink_software_name,".exe"),getwd())
plink_software_path=getwd()
}

}


if(length(grep(" ",beagle_software_path))>0){
#cat("Attention: Phe path include space!")
file.copy(paste0(beagle_software_path,"/",beagle_software_name),getwd())
beagle_software_path=getwd()
}




#Part1:调用Plink,对Plink格式的数据进行质控
if(data_analysis_method%in%c("QC","QC_Imputation")){

output_data_name=paste0(output_data_name,"_QC")
if(!is.null(plink_software_path)){plink_software_path=paste0(plink_software_path,"/")}

temp=system(paste0(plink_software_path,plink_software_name),ignore.stdout=T)
if(temp==127){
stop("Software couldn't find plink in your computer!")
}

system(paste0(plink_software_path,plink_software_name," --noweb --file ",output_data_path,"/","temp_plink",
	        " --missing-genotype ",miss_base," --geno ",qc_snp_rate," --mind ",
			qc_ind_rate," --maf ",qc_maf," --hwe ",qc_hwe," --make-just-fam  ",chr_set," ",extra_parameter," ",
			" --allow-extra-chr --write-snplist --out ",	
			output_data_path,"/",output_data_name," --recode vcf-iid"),ignore.stdout=T)

if(file.exists(paste0(output_data_name,".fam"))){
IND_QC=fread(paste0(output_data_name,".fam"),header=F,data.table=F)
IND_QC_n=nrow(IND_QC)
}else{stop("All individuals have been filtered, please check your data carefully!")}

if(file.exists(paste0(output_data_name,".snplist"))){
SNP_QC=fread(paste0(output_data_name,".snplist"),header=F,data.table=F)
SNP_QC_n=nrow(SNP_QC)
}else{stop("All snps have been filtered, please check your data carefully!")}

cat("Quality control is done \n")
cat(paste0("Before QC ,there are ", SNP_n," SNPs; After QC ,there are ",SNP_QC_n," SNPs remained \n"))
cat(paste0("Before QC ,there are ", IND_n," indviduals; After QC ,there are ",IND_QC_n," indviduals remained \n"))



}else{
output_data_name=output_data_name
if(!is.null(plink_software_path)){plink_software_path=paste0(plink_software_path,"/")}
system(paste0(plink_software_path,plink_software_name," --noweb --file ",output_data_path,"/","temp_plink",
	        " --missing-genotype ",miss_base,chr_set," ",extra_parameter," ",
			" --allow-extra-chr  --out ",
			output_data_path,"/",output_data_name," --recode vcf-iid"),ignore.stdout=T)
}

#删除temp_plink data
system(paste0("rm -rf temp_plink.map")) 
system(paste0("rm -rf temp_plink.ped")) 


#删除copy过来的软件
if(file.exists(paste0(plink_software_name))){system(paste0("rm -rf ",plink_software_name)) }
if(file.exists(paste0(plink_software_name,".exe"))){system(paste0("rm -rf ",paste0(plink_software_name,".exe")))}



if(file.exists(paste0(output_data_name,".nosex"))){system(paste0("rm -rf ",output_data_name,".nosex")) }
#if(file.exists(paste0(output_data_name,".log"))){system(paste0("rm -rf ",output_data_name,".log")) }
file.rename(paste0(output_data_name,".log"),paste0(output_data_name,"_plink.log"))
if(file.exists(paste0(output_data_name,".snplist"))){system(paste0("rm -rf ",output_data_name,".snplist")) }

}else{

if(input_data_type=="VCF"&data_analysis_method=="Imputation"){
file.copy(paste0(input_data_path,"/",input_data_name),paste0(output_data_path,"/",output_data_name,".vcf"))
#system(paste0("cp -r ",input_data_path,"/",input_data_name,"   ",output_data_path,"/",output_data_name,".vcf"))
}

setwd(output_data_path)
output_data_path=getwd()

if(length(grep(" ",beagle_software_path))>0){
#cat("Attention: Phe path include space!")
file.copy(paste0(beagle_software_path,"/",beagle_software_name),getwd())
beagle_software_path=getwd()
}

}

if(file.exists(paste0(output_data_name,".fam"))){unlink(paste0(output_data_name,".fam"))}

#Part2:调用 Beagle 对数据进行填充
if(data_analysis_method%in%c("QC_Imputation","Imputation")){ 
beagle_software=paste0(beagle_software_path,"/",beagle_software_name)  
OOOOOOO0OOO0O0O0O0O0=paste0(output_data_path,"/",output_data_name,".vcf")
OOOOO0O0OOO0OOOOOOOO=paste0(output_data_path,"/",output_data_name,"_Imputation")
output_data_name=paste0(output_data_name,"_Imputation")



if(!is.null(beagle_ped_path)){
cat("Please make sure the version of Beagle in this analysi is 4.0 for which included pedigree_ module \n")
beagle_pedigree=paste0(" ped=",beagle_ped_path,"/",beagle_ped_name)}else{beagle_pedigree=NULL} 
				 
if(!is.null(beagle_ref_data_path)){
beagle_reference=paste0(" ref=",beagle_ref_data_path,"/",beagle_ref_data_name)}else{beagle_reference=NULL}	
			 
cat("Using Beagle to perform  imputation ......\n")

#调用 beagle

temp1=system("java",ignore.stdout = T,ignore.stderr = T)
if(temp1==127){
stop("Java doesn't install on you computer, please intall Java before!")
}

temp=system(paste0("java ",Java_Space," -jar ", beagle_software),ignore.stdout=T)
if(temp==1){
stop("Software couldn't find beagle in your computer!")
}

system(paste0("java ",Java_Space," -jar ", beagle_software,beagle_reference,beagle_pedigree,sigle_duo_trio,
			" gt=",OOOOOOO0OOO0O0O0O0O0," nthreads=",cpu_cores," out=",OOOOO0O0OOO0OOOOOOOO))

file.rename(paste0(OOOOO0O0OOO0OOOOOOOO,".log"),paste0(OOOOO0O0OOO0OOOOOOOO,"_beagle_output.log"))
#删除填充的input data
system(paste0("rm -rf ",OOOOOOO0OOO0O0O0O0O0))
#解压 .gz 格式			
system(paste0("gzip  -d ",OOOOO0O0OOO0OOOOOOOO,".vcf.gz")) #将.gz文件解压为 .vcf  ,便于进行 vcftools转换
}
 

#Part3:将质控或质控或填充好的数据进行输出, 默认输出的数据为： OOOOO0O0OOO0OOOOOOOO .vcf  
 

 sum_data=geno_format(      
          input_data_name=paste0(output_data_name,".vcf"),  
	     input_data_path=output_data_path,		  
          input_data_type="VCF",  	
		miss_base=miss_base,
		cpu_cores=cpu_cores,
		bigmemory_cal=bigmemory_cal,
		bigmemory_data_type=bigmemory_data_type,
		bigmemory_data_path=bigmemory_data_path,
		bigmemory_data_name=bigmemory_data_name,
		phased_symbol=phased_symbol,
		phased_genotype=phased_genotype,
		haplotype_window_nSNP=haplotype_window_nSNP,
		haplotype_window_kb=haplotype_window_kb,
		haplotype_window_block=haplotype_window_block,		
		output_data_type=output_data_type,
		output_data_name=output_data_name,
		output_data_path=output_data_path,			
		return_result=FALSE)
		
if(file.exists(paste0(beagle_software_name))){system(paste0("rm -rf ",beagle_software_name))}
if(file.exists(paste0(output_data_name,".nosex"))){system(paste0("rm -rf ",output_data_name,".nosex")) }
if(file.exists(paste0(output_data_name,".fam"))){system(paste0("rm -rf ",output_data_name,".fam")) }
if(!"VCF"%in%output_data_type){system(paste0("rm -rf ",output_data_path,"/",output_data_name,".vcf"))}
}
cal_kinship<-function(
		input_data_type=NULL, #"Plink" , "Hapmap" , "VCF" , "BLUPF90","Numeric"  "Haplotype"
		input_data_path=NULL,
		input_data_name=NULL,  #hapmap 为全名，vcf 和 plink 默认不包括后缀，为 name.vcf name.ped name.map 
		input_data_hmp=NULL,
		input_data_plink_ped=NULL,
		input_data_plink_map=NULL,
		input_data_blupf90=NULL,
		input_data_numeric=NULL,
		input_data_vcf=NULL,
		input_data_haplotype_hap=NULL,
		input_data_haplotype_map=NULL,
		input_data_haplotype_sample=NULL,
		input_data_numeric_map=NULL,
		input_data_blupf90_map=NULL,
		bigmemory_cal=FALSE,
		bigmemory_data_type="double",
		bigmemory_data_path=getwd(),
		bigmemory_data_name="blupADC",
		phased_genotype=FALSE,
		haplotype_window_nSNP=NULL,
		haplotype_window_kb=NULL,
		haplotype_window_block=NULL,
		miss_base=NULL,		
		miss_base_num=0,		
	     cpu_cores=1,
	     pedigree_rename=TRUE,
	     dis_progress=FALSE,
	     kinship_type=NULL,         #c("G_A","G_D","P_A","P_D","H_A"),
	     dominance_type="genotypic",     # "genotypic","classical"
	     Metafounder_algorithm=FALSE,
	     SSBLUP_omega=0.05, # ssblup parameter
	     gene_dropping_algorithm=FALSE,
	     gene_dropping_iteration=1000, # gene_dropping_algorithm 的模拟次数
	     inbred_type=NULL, #不同类型的近交系数 "Homozygous","G_Diag","Pedigree","H_diag"
	     kinship_base=FALSE,          # 设置 p = q =0.5, 
	     kinship_trace=FALSE,                     #根据迹计算亲缘矩阵
	     input_pedigree_path=NULL,
	     input_pedigree_name=NULL,
	     input_pedigree=NULL,
	     pedigree_multi_col=FALSE,
	     IND_rename=FALSE, # 根据提供的系谱，得到rename系谱，并将基因型个体进行rename
	     return_result=FALSE,
          output_matrix_type="col_all",  # c("col_three","col_all")
	      matrix_log_det=FALSE,  
	      col3_threshold=0,
	      APY_algorithm=FALSE,
	      APY_eigen_threshold=0.95,
	      APY_n_core=0,
		  col3_bigmemory=FALSE,
		  col3_bigmemory_data_path=NULL,
		  col3_bigmemory_data_name=NULL,
	      output_matrix_name=NULL,
	      output_matrix_path=NULL){ 

if(bigmemory_cal==TRUE){
col3_bigmemory=TRUE
col3_bigmemory_data_path=bigmemory_data_path
col3_bigmemory_data_name=bigmemory_data_name
}
		  
data_block=NULL
library(bigmemory)
options(bigmemory.allow.dimnames=TRUE)
#检查输入的参数
if(sum(is.null(output_matrix_name),is.null(output_matrix_path))==1){
stop("Paramters:output_matrix_name and output_matrix_path couldn't be: the one is non-NULL and the other is NULL")
}

MKL_status=find("getMKLthreads")	
if(length(MKL_status)==0){
cat("In terms of relationship matrix construction, we highly recommend you to switch to R-open which supprots MKL! \n")
}else{

if(cpu_cores==1){

if(getMKLthreads()==1){

max_cores=parallel::detectCores()
setMKLthreads(max_cores)

}}

}

inbred_output_matrix_path=output_matrix_path
inbred_output_matrix_name=output_matrix_name
if(is.null(kinship_type)&!is.null(inbred_type)){
output_matrix_path=NULL
output_matrix_name=NULL
}


#判断输入类型
if("Homozygous"%in%inbred_type){kinship_type=c(kinship_type,"G_A")}
if("G_Diag"%in%inbred_type){kinship_type=c(kinship_type,"G_A")}
if("Pedigree"%in%inbred_type){kinship_type=c(kinship_type,"P_A")}
if("H_diag"%in%inbred_type){kinship_type=c(kinship_type,"H_A")}


#基因型数据
if("G_A"%in%kinship_type | "G_Ainv"%in%kinship_type |"G_D"%in%kinship_type |"G_Dinv"%in%kinship_type| "H_A"%in%kinship_type| "H_Ainv"%in%kinship_type|"H_D"%in%kinship_type|"H_Dinv"%in%kinship_type){ 
#程序自动判断输入类型
input=O0OOO0O0OOO0O0O0OOOO(input_data_type=input_data_type,input_data_hmp=input_data_hmp,input_data_plink_ped=input_data_plink_ped,
                              input_data_plink_map=input_data_plink_map,input_data_blupf90=input_data_blupf90,
						  input_data_numeric=input_data_numeric,input_data_vcf=input_data_vcf,
						  input_data_haplotype_hap=input_data_haplotype_hap,input_data_haplotype_map=input_data_haplotype_map,
						  input_data_haplotype_sample=input_data_haplotype_sample,
						  miss_base=miss_base)

input_data_type=input$input_data_type
miss_base=input$miss_base

if(!is.null(output_matrix_name)){output_matrix_name=paste0(output_matrix_name,"_")}
#读取数据
#获取 input_data_numeric
if(input_data_type=="Numeric"){

if(!is.null(input_data_path)&!is.null(input_data_name)){
input_data_numeric=fread(paste0(input_data_path,"/",input_data_name),header=F,data.table=F)
IND_geno=input_data_numeric[,1]
input_data_numeric=DataFrame_to_arma(input_data_numeric[,-1])
rownames(input_data_numeric)=IND_geno
}

}else{	
sum_data=geno_format(   
          input_data_hmp=input_data_hmp,      
          input_data_name=input_data_name,  
	     input_data_path=input_data_path,		  
          input_data_type=input_data_type,  
	     input_data_plink_map=input_data_plink_map, 
	     input_data_plink_ped=input_data_plink_ped, 
	     input_data_blupf90=input_data_blupf90,
		input_data_numeric=input_data_numeric,
		input_data_vcf=input_data_vcf,
		input_data_haplotype_hap=input_data_haplotype_hap,
		input_data_haplotype_map=input_data_haplotype_map,
		input_data_haplotype_sample=input_data_haplotype_sample,
		input_data_numeric_map=input_data_numeric_map,
		input_data_blupf90_map=input_data_blupf90_map,		
	     phased_genotype=phased_genotype,
		bigmemory_cal=bigmemory_cal,
		bigmemory_data_type="integer",
		bigmemory_data_path=bigmemory_data_path,
		bigmemory_data_name=bigmemory_data_name,		 
	     haplotype_window_nSNP=haplotype_window_nSNP,
	     haplotype_window_kb=haplotype_window_kb,
	     haplotype_window_block=haplotype_window_block,			
		miss_base=miss_base,
		cpu_cores=cpu_cores,
		output_data_type="Numeric",
		return_result=TRUE)	
input_data_numeric=sum_data$numeric
data_block=sum_data$phased_block
rm(sum_data);gc();
}

if(is.null(rownames(input_data_numeric))){stop("Provided numeric data must has rownames!")}

if("big.matrix"%in%class(input_data_numeric)){bigmemory_cal=TRUE}
if(is.data.frame(input_data_numeric)){
IND_geno=rownames(input_data_numeric)
input_data_numeric=DataFrame_to_arma(input_data_numeric)
rownames(input_data_numeric)=IND_geno}

if(typeof(input_data_numeric)!="integer"&(!"big.matrix"%in%class(input_data_numeric))){
input_data_numeric=NumericMatrix_to_arma(input_data_numeric)}
#创建bigmemory-对象
if(bigmemory_cal==TRUE&!("big.matrix"%in%class(input_data_numeric))){
cat("Provided genotype data is not bigmemory-object, convert provided genotype data as bigmemory-object......\n")
IND_geno=rownames(input_data_numeric)
delete_bigmemory_file("convert_numeric",bigmemory_data_name,bigmemory_data_path,TRUE)

input_data_numeric=bigmemory_object_convert(input_data_numeric,bigmemory_data_name,bigmemory_data_path)
gc();
options(bigmemory.allow.dimnames=TRUE)
rownames(input_data_numeric)=IND_geno
colnames(input_data_numeric)=paste0("SNP",1:ncol(input_data_numeric))
cat(paste0("Save converted bigmemory-object as: ",bigmemory_data_path,"/",bigmemory_data_name,"_convert_numeric"," \n"))								   
}

		
if(exists("input_data_hmp")){rm(input_data_hmp)}
if(exists("input_data_plink_map")){rm(input_data_plink_map)}
if(exists("input_data_plink_ped")){rm(input_data_plink_ped)}
if(exists("input_data_blupf90")){rm(input_data_blupf90)}
if(exists("input_data_vcf")){rm(input_data_vcf)}
if(exists("input_data_haplotype_hap")){rm(input_data_haplotype_hap)}
gc();
}

if(bigmemory_cal==TRUE){output_matrix_path=bigmemory_data_path;output_matrix_name=bigmemory_data_name}

#读取系谱数据
if(!is.null(input_pedigree_path)&!is.null(input_pedigree_name)){
input_pedigree=fread(paste0(input_pedigree_path,"/",input_pedigree_name),data.table=F)
}

#rename IND
if(!is.null(input_pedigree)&is.null(input_data_numeric)){

rename_ped=trace_pedigree(input_pedigree,multi_col=pedigree_multi_col)$rename_ped
IND_pedigree=rename_ped[,1]
num_ped=as.matrix(rename_ped[,3:5]) #数字化的3列系谱
num_ped=apply(num_ped,2,as.integer)
if(IND_rename==TRUE){
IND_pedigree=num_ped[,1]
}

if("col_three"%in%output_matrix_type){
if(NA %in% as.numeric(IND_pedigree)){
stop("Error:In terms of constring col_three matrix, provided OOO0OOO0O0OOO0O0OOO0OOOOOO id include character, please set IND_rename=TRUE !")
}
IND_pedigree=as.numeric(IND_pedigree)
}

}else if(!is.null(input_pedigree)&!is.null(input_data_numeric)){
cat("Please make sure the genotype id is accordance with the _pedigree_id! \n")
rename_ped=trace_pedigree(input_pedigree,multi_col=pedigree_multi_col)$rename_ped
IND_pedigree=rename_ped[,1]
num_ped=as.matrix(rename_ped[,3:5]) #数字化的3列系谱
num_ped=apply(num_ped,2,as.integer)
IND_geno=rownames(input_data_numeric)

if(sum(IND_geno%in%IND_pedigree)==0){stop("Provided genotype id is not accordance with the OOO0OOO0O0OOO0O0OOO0OOOOOO id!")}
if(IND_rename==TRUE){
IND_geno=num_ped[match(IND_geno,IND_pedigree),1]
IND_pedigree=num_ped[,1]
}

if("col_three"%in%output_matrix_type){
if(NA %in% as.numeric(IND_geno)){
stop("Error:In terms of constring col_three matrix,, provided genotype id include character, please set the genotype id as integer or provide OOO0OOO0O0OOO0O0OOO0OOOOOO for recoding the genotype id!")
}
IND_geno=as.numeric(IND_geno)
IND_pedigree=as.numeric(IND_pedigree)
}

}else if(is.null(input_pedigree)&!is.null(input_data_numeric)){

IND_geno=rownames(input_data_numeric)
if("col_three"%in%output_matrix_type){
if(NA %in% as.numeric(IND_geno)){
stop("Error:In terms of constring col_three matrix, provided genotype id include character, please set the genotype id as integer or provide OOO0OOO0O0OOO0O0OOO0OOOOOO for recoding the genotype id!")
}
IND_geno=as.numeric(IND_geno)
}
}


G_A=NULL;G_Ainv=NULL;G_D=NULL;G_Dinv=NULL;P_A=NULL;P_Ainv=NULL;P_D=NULL;P_Dinv=NULL;H_A=NULL;H_Ainv=NULL;H_D=NULL;H_Dinv=NULL;
homo_inbred=NULL;diag_inbred=NULL;pedigree_inbred=NULL;O0O0OOOOO0OOOOO0O0OO=NULL

if(bigmemory_cal==FALSE){
#Part1 计算kinship 
########################################## P_A ##########################################
if("P_A"%in%kinship_type&(!"P_Ainv"%in%kinship_type)){
P_A=makeA_cpp(num_ped)
P_Ainv=NULL

}else if("P_A"%in%kinship_type&("P_Ainv"%in%kinship_type)){
P_A=makeA_cpp(num_ped)
P_Ainv=makeAinv_cpp(num_ped)

}else if(!"P_A"%in%kinship_type&("P_Ainv"%in%kinship_type)){
P_Ainv=makeAinv_cpp(num_ped)
P_A=NULL
}


#输出矩阵文件-start
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))|(!is.null(col3_bigmemory_data_name)&!is.null(col3_bigmemory_data_path))){

if(!is.null(P_A)){
if(!is.null(output_matrix_path)){setwd(output_matrix_path)}
if(!is.null(col3_bigmemory_data_path)){setwd(col3_bigmemory_data_path)}
fwrite(data.table(IND_pedigree),"IND_pedigree.txt",quote=F,row.names=F,col.names=F,sep="\t")
if("col_three" %in% output_matrix_type){
if(col3_bigmemory==FALSE){
O0O0OOO0O0OOOOOOOOOO=matrix_col3(P_A,IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
fwrite(data.table(O0O0OOO0O0OOOOOOOOOO),paste0(output_matrix_name,"P_A_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
rm(O0O0OOO0O0OOOOOOOOOO);gc();
}else{
cat(paste0("Save bigmemory-P_A_col_three matrix as ",col3_bigmemory_data_path,"/",col3_bigmemory_data_name,"_P_A_col_three ......\n"))
matrix_col3_memory_alt(P_A,paste0(col3_bigmemory_data_name,"_P_A"),col3_bigmemory_data_path,
						   IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
}
if("col_all" %in% output_matrix_type){
fwrite(data.frame(P_A),paste0(output_matrix_name,"P_A_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(P_A);gc();}
}

if(!is.null(P_Ainv)){
if(!is.null(output_matrix_path)){setwd(output_matrix_path)}
if(!is.null(col3_bigmemory_data_path)){setwd(col3_bigmemory_data_path)}
fwrite(data.table(IND_pedigree),"IND_pedigree.txt",quote=F,row.names=F,col.names=F,sep="\t")
if("col_three" %in% output_matrix_type){

if(col3_bigmemory==FALSE){
O0O0OOO0O0OOOOOOOOOO=matrix_col3(P_Ainv,IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
fwrite(data.table(O0O0OOO0O0OOOOOOOOOO),paste0(output_matrix_name,"P_Ainv_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
rm(O0O0OOO0O0OOOOOOOOOO);gc();
}else{
cat(paste0("Save bigmemory-P_Ainv_col_three matrix as ",col3_bigmemory_data_path,"/",col3_bigmemory_data_name,"_P_Ainv_col_three ......\n"))
matrix_col3_memory_alt(P_Ainv,paste0(col3_bigmemory_data_name,"_P_Ainv"),col3_bigmemory_data_path,
						   IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
}

if("col_all" %in% output_matrix_type){
fwrite(data.frame(P_Ainv),paste0(output_matrix_name,"P_Ainv_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(P_Ainv);gc();}
}
}
#输出矩阵文件-end


########################################## P_D ##########################################
if(gene_dropping_algorithm==FALSE){

if("P_D"%in%kinship_type&(!"P_Dinv"%in%kinship_type)){
P_D=makeD_cpp(num_ped,FALSE)$D
P_Dinv=NULL
}else if("P_D"%in%kinship_type&("P_Dinv"%in%kinship_type)){
P_Dresult=makeD_cpp(num_ped,TRUE)
P_D=P_Dresult$D
P_Dinv=P_Dresult$Dinv
rm(P_Dresult);gc();
}else if(!"P_D"%in%kinship_type&("P_Dinv"%in%kinship_type)){
P_Dinv=makeD_cpp(num_ped,TRUE)$Dinv
P_D=NULL
}

}else{

if("P_D"%in%kinship_type&(!"P_Dinv"%in%kinship_type)){
P_D=gene_dropping_D(num_ped,iteration=gene_dropping_iteration,cpu_cores=cpu_cores,inverse=FALSE)$D
P_Dinv=NULL
}else if("P_D"%in%kinship_type&("P_Dinv"%in%kinship_type)){
P_Dresult=gene_dropping_D(num_ped,iteration=gene_dropping_iteration,cpu_cores=cpu_cores,inverse=TRUE)
P_D=P_Dresult$D
P_Dinv=P_Dresult$Dinv
rm(P_Dresult);gc();
}else if(!"P_D"%in%kinship_type&("P_Dinv"%in%kinship_type)){
P_Dinv=gene_dropping_D(num_ped,iteration=gene_dropping_iteration,cpu_cores=cpu_cores,inverse=TRUE)$Dinv
P_D=NULL
}
}

#输出矩阵文件-start
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))|(!is.null(col3_bigmemory_data_name)&!is.null(col3_bigmemory_data_path))){

if(!is.null(P_D)){
if(!is.null(output_matrix_path)){setwd(output_matrix_path)}
if(!is.null(col3_bigmemory_data_path)){setwd(col3_bigmemory_data_path)}
fwrite(data.table(IND_pedigree),"IND_pedigree.txt",quote=F,row.names=F,col.names=F,sep="\t")

if("col_three" %in% output_matrix_type){
if(col3_bigmemory==FALSE){
O0O0OOO0O0OOOOOOOOOO=matrix_col3(P_D,IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
fwrite(data.table(O0O0OOO0O0OOOOOOOOOO),paste0(output_matrix_name,"P_D_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
rm(O0O0OOO0O0OOOOOOOOOO);gc();
}else{
cat(paste0("Save bigmemory-P_D_col_three matrix as ",col3_bigmemory_data_path,"/",col3_bigmemory_data_name,"_P_D_col_three ......\n"))
matrix_col3_memory_alt(P_D,paste0(col3_bigmemory_data_name,"_P_D"),col3_bigmemory_data_path,
						   IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
}

if("col_all" %in% output_matrix_type){
fwrite(data.frame(P_D),paste0(output_matrix_name,"P_D_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(P_D);gc();}
}

if(!is.null(P_Dinv)){
if(!is.null(output_matrix_path)){setwd(output_matrix_path)}
if(!is.null(col3_bigmemory_data_path)){setwd(col3_bigmemory_data_path)}
fwrite(data.table(IND_pedigree),"IND_pedigree.txt",quote=F,row.names=F,col.names=F,sep="\t")

if("col_three" %in% output_matrix_type){
if(col3_bigmemory==FALSE){
O0O0OOO0O0OOOOOOOOOO=matrix_col3(P_Dinv,IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
fwrite(data.table(O0O0OOO0O0OOOOOOOOOO),paste0(output_matrix_name,"P_Dinv_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
rm(O0O0OOO0O0OOOOOOOOOO);gc();
}else{
cat(paste0("Save bigmemory-P_Dinv_col_three matrix as ",col3_bigmemory_data_path,"/",col3_bigmemory_data_name,"_P_Dinv_col_three ......\n"))
matrix_col3_memory_alt(P_Dinv,paste0(col3_bigmemory_data_name,"_P_Dinv"),col3_bigmemory_data_path,
						   IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
}


if("col_all" %in% output_matrix_type){
fwrite(data.frame(P_Dinv),paste0(output_matrix_name,"P_Dinv_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(P_Dinv);gc();}
}
}
#输出矩阵文件-end


########################################## G_A ##########################################
if(APY_algorithm==FALSE){
if("G_A"%in%kinship_type&(!"G_Ainv"%in%kinship_type)){
G_A=G_matrix_cpp(input_data_numeric,base=kinship_base,trace=kinship_trace,inverse=FALSE)$G
diag_inbred=data.frame(Id=IND_geno,Diag_inbred=diag(G_A)-1,stringsAsFactors=F)
G_Ainv=NULL
}else if("G_A"%in%kinship_type&("G_Ainv"%in%kinship_type)){
G_Aresult=G_matrix_cpp(input_data_numeric,base=kinship_base,trace=kinship_trace,inverse=TRUE)
G_A=G_Aresult$G
diag_inbred=data.frame(Id=IND_geno,Diag_inbred=diag(G_A)-1,stringsAsFactors=F)
G_Ainv=G_Aresult$Ginv
rm(G_Aresult);gc();
}else if(!"G_A"%in%kinship_type&("G_Ainv"%in%kinship_type)){
G_Ainv=G_matrix_cpp(input_data_numeric,base=kinship_base,trace=kinship_trace,inverse=TRUE)$Ginv
G_A=NULL
}

}else{

if("G_A"%in%kinship_type&(!"G_Ainv"%in%kinship_type)){
G_A=APY_inverse_cpp(input_data_numeric,IND_geno,APY_eigen_threshold,APY_n_core,re_inverse=TRUE)$G
G_Ainv=NULL
}else if("G_A"%in%kinship_type&("G_Ainv"%in%kinship_type)){
G_Aresult=APY_inverse_cpp(input_data_numeric,IND_geno,APY_eigen_threshold,APY_n_core,re_inverse=TRUE)
G_A=G_Aresult$G
G_Ainv=G_Aresult$Ginv
rm(G_Aresult);gc();
}else if(!"G_A"%in%kinship_type&("G_Ainv"%in%kinship_type)){
G_Ainv=APY_inverse_cpp(input_data_numeric,IND_geno,APY_eigen_threshold,APY_n_core,re_inverse=FALSE)$Ginv
G_A=NULL
}
}


#输出矩阵文件-start
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))|(!is.null(col3_bigmemory_data_name)&!is.null(col3_bigmemory_data_path))){

if(!is.null(G_A)){
if(!is.null(output_matrix_path)){setwd(output_matrix_path)}
if(!is.null(col3_bigmemory_data_path)){setwd(col3_bigmemory_data_path)}
fwrite(data.table(IND_geno),"IND_geno.txt",quote=F,row.names=F,col.names=F,sep="\t")

if("col_three" %in% output_matrix_type){
if(col3_bigmemory==FALSE){
O0O0OOO0O0OOOOOOOOOO=matrix_col3(G_A,IND_geno=as.numeric(IND_geno),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
fwrite(data.table(O0O0OOO0O0OOOOOOOOOO),paste0(output_matrix_name,"G_A_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
rm(O0O0OOO0O0OOOOOOOOOO);gc();
}else{
cat(paste0("Save bigmemory-G_A_col_three matrix as ",col3_bigmemory_data_path,"/",col3_bigmemory_data_name,"_G_A_col_three ......\n"))
matrix_col3_memory_alt(G_A,paste0(col3_bigmemory_data_name,"_G_A"),col3_bigmemory_data_path,
						   IND_geno=as.numeric(IND_geno),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
}



if("col_all" %in% output_matrix_type){
fwrite(data.frame(G_A),paste0(output_matrix_name,"G_A_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(G_A);gc();}
}

if(!is.null(G_Ainv)){
if(!is.null(output_matrix_path)){setwd(output_matrix_path)}
if(!is.null(col3_bigmemory_data_path)){setwd(col3_bigmemory_data_path)}
fwrite(data.table(IND_geno),"IND_geno.txt",quote=F,row.names=F,col.names=F,sep="\t")

if("col_three" %in% output_matrix_type){
if(col3_bigmemory==FALSE){
O0O0OOO0O0OOOOOOOOOO=matrix_col3(G_Ainv,IND_geno=as.numeric(IND_geno),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
fwrite(data.table(O0O0OOO0O0OOOOOOOOOO),paste0(output_matrix_name,"G_Ainv_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
rm(O0O0OOO0O0OOOOOOOOOO);gc();
}else{
cat(paste0("Save bigmemory-G_Ainv_col_three matrix as ",col3_bigmemory_data_path,"/",col3_bigmemory_data_name,"_G_Ainv_col_three ......\n"))
matrix_col3_memory_alt(G_Ainv,paste0(col3_bigmemory_data_name,"_G_Ainv"),col3_bigmemory_data_path,
						   IND_geno=as.numeric(IND_geno),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
}


if("col_all" %in% output_matrix_type){
fwrite(data.frame(G_Ainv),paste0(output_matrix_name,"G_Ainv_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(G_Ainv);gc();}
}
}
#输出矩阵文件-end

########################################## G_D ##########################################
if("genotypic"%in%dominance_type){
if("G_D"%in%kinship_type&(!"G_Dinv"%in%kinship_type)){
G_D=Dominance_matrix_cpp(input_data_numeric,base=kinship_base,trace=kinship_trace,inverse=FALSE)$D
G_Dinv=NULL
}else if("G_D"%in%kinship_type&("G_Dinv"%in%kinship_type)){
G_Dresult=Dominance_matrix_cpp(input_data_numeric,base=kinship_base,trace=kinship_trace,inverse=TRUE)
G_D=G_Dresult$D
G_Dinv=G_Dresult$Dinv
rm(G_Dresult);gc();
}else if(!"G_D"%in%kinship_type&("G_Dinv"%in%kinship_type)){
G_Dinv=Dominance_matrix_cpp(input_data_numeric,base=kinship_base,trace=kinship_trace,inverse=TRUE)$Dinv
G_D=NULL
}
}else if("classical"%in%dominance_type){

if("G_D"%in%kinship_type&(!"G_Dinv"%in%kinship_type)){
G_D=Deviation_matrix_cpp(input_data_numeric,base=kinship_base,trace=kinship_trace,inverse=FALSE)$D
G_Dinv=NULL
}else if("G_D"%in%kinship_type&("G_Dinv"%in%kinship_type)){
G_Dresult=Deviation_matrix_cpp(input_data_numeric,base=kinship_base,trace=kinship_trace,inverse=TRUE)
G_D=G_Dresult$D
G_Dinv=G_Dresult$Dinv
rm(G_Dresult);gc();
}else if(!"G_D"%in%kinship_type&("G_Dinv"%in%kinship_type)){
G_Dinv=Deviation_matrix_cpp(input_data_numeric,base=kinship_base,trace=kinship_trace,inverse=TRUE)$Dinv
G_D=NULL
}
}


#输出矩阵文件-start
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))|(!is.null(col3_bigmemory_data_name)&!is.null(col3_bigmemory_data_path))){

if(!is.null(G_D)){
if(!is.null(output_matrix_path)){setwd(output_matrix_path)}
if(!is.null(col3_bigmemory_data_path)){setwd(col3_bigmemory_data_path)}
fwrite(data.table(IND_geno),"IND_geno.txt",quote=F,row.names=F,col.names=F,sep="\t")

if("col_three" %in% output_matrix_type){
if(col3_bigmemory==FALSE){
O0O0OOO0O0OOOOOOOOOO=matrix_col3(G_D,IND_geno=as.numeric(IND_geno),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
fwrite(data.table(O0O0OOO0O0OOOOOOOOOO),paste0(output_matrix_name,"G_D_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
rm(O0O0OOO0O0OOOOOOOOOO);gc();
}else{
cat(paste0("Save bigmemory-G_D_col_three matrix as ",col3_bigmemory_data_path,"/",col3_bigmemory_data_name,"_G_D_col_three ......\n"))
matrix_col3_memory_alt(G_D,paste0(col3_bigmemory_data_name,"_G_D"),col3_bigmemory_data_path,
						   IND_geno=as.numeric(IND_geno),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
}


if("col_all" %in% output_matrix_type){
fwrite(data.frame(G_D),paste0(output_matrix_name,"G_D_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(G_D);gc();}
}

if(!is.null(G_Dinv)){
if(!is.null(output_matrix_path)){setwd(output_matrix_path)}
if(!is.null(col3_bigmemory_data_path)){setwd(col3_bigmemory_data_path)}
fwrite(data.table(IND_geno),"IND_geno.txt",quote=F,row.names=F,col.names=F,sep="\t")

if("col_three" %in% output_matrix_type){
if(col3_bigmemory==FALSE){
O0O0OOO0O0OOOOOOOOOO=matrix_col3(G_Dinv,IND_geno=as.numeric(IND_geno),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
fwrite(data.table(O0O0OOO0O0OOOOOOOOOO),paste0(output_matrix_name,"G_Dinv_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
rm(O0O0OOO0O0OOOOOOOOOO);gc();
}else{
cat(paste0("Save bigmemory-G_Dinv_col_three matrix as ",col3_bigmemory_data_path,"/",col3_bigmemory_data_name,"_G_Dinv_col_three ......\n"))
matrix_col3_memory_alt(G_Dinv,paste0(col3_bigmemory_data_name,"_G_Dinv"),col3_bigmemory_data_path,
						   IND_geno=as.numeric(IND_geno),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
}

if("col_all" %in% output_matrix_type){
fwrite(data.frame(G_Dinv),paste0(output_matrix_name,"G_Dinv_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(G_Dinv);gc();}
}
}
#输出矩阵文件-end


########################################## single-step ##########################################
if("H_A"%in%kinship_type|"H_Ainv"%in%kinship_type|"H_D"%in%kinship_type|"H_Dinv"%in%kinship_type){
IND_geno2=setdiff(IND_geno,IND_pedigree)    #有基因型但是无系谱
IND_A22=intersect(IND_pedigree,IND_geno)    #既有基因型又有系谱
IND_A11=setdiff(IND_pedigree,IND_A22)        #有系谱但是无基因型
IND_Additive=c(IND_A11,IND_A22)               #在H矩阵中的个体
pos_geno=match(IND_A22,IND_geno)-1            #既有基因型又有系谱的个体在 基因型数据中的位置
pos_A11=match(IND_A11,IND_pedigree)-1         #有系谱但是无基因型个体在  在系谱中的位置
pos_A22=match(IND_A22,IND_pedigree)-1         #既有基因型又有系谱的个体在  在系谱中的位置
pos_A=match(IND_Additive,IND_pedigree)-1      #H矩阵的个体  在系谱中的位置 
pos_H22=match(IND_A22,IND_Additive)-1         #既有基因型又有系谱的个体在  在H矩阵中的位置
if(length(IND_A22)==0){stop("No aminals both with genotype and_pedigree, Please check you data ! \n")
}else{cat(length(IND_A22),"individuals both have genotype and_pedigree \n")}
if(length(IND_geno2)>0){cat(paste0(length(IND_geno2)," individuals have genotype but not in the_pedigree, these individuals are removed ! \n"))}
}

#return(list(num_ped=num_ped,input_data_numeric=input_data_numeric,pos_geno=pos_geno,
#			pos_A11=pos_A11,pos_A22=pos_A22,pos_A=pos_A,pos_H22=pos_H22))


if(Metafounder_algorithm==FALSE){
if("H_A"%in%kinship_type&(!"H_Ainv"%in%kinship_type)){
H_Aresult=makeHA_cpp(num_ped,input_data_numeric,IND_geno,pos_A11,pos_A22,pos_geno,pos_A,pos_H22,APY_algorithm,APY_eigen_threshold,APY_n_core,direct=TRUE,inverse=FALSE,omega=SSBLUP_omega)
H_A_inbred=H_Aresult$inbred
O0O0OOOOO0OOOOO0O0OO=data.frame(Id=IND_Additive,Diag_inbred=H_A_inbred,stringsAsFactors=F)
H_A=H_Aresult$H
H_Ainv=NULL
rm(H_Aresult);gc();
}else if("H_A"%in%kinship_type&("H_Ainv"%in%kinship_type)){
H_Aresult=makeHA_cpp(num_ped,input_data_numeric,IND_geno,pos_A11,pos_A22,pos_geno,pos_A,pos_H22,APY_algorithm,APY_eigen_threshold,APY_n_core,direct=TRUE,inverse=TRUE,omega=SSBLUP_omega)
H_A_inbred=H_Aresult$inbred
O0O0OOOOO0OOOOO0O0OO=data.frame(Id=IND_Additive,Diag_inbred=H_A_inbred,stringsAsFactors=F)
H_A=H_Aresult$H
H_Ainv=H_Aresult$Hinv
rm(H_Aresult);gc();
}else if(!"H_A"%in%kinship_type&("H_Ainv"%in%kinship_type)){
H_Ainv=makeHA_cpp(num_ped,input_data_numeric,IND_geno,pos_A11,pos_A22,pos_geno,pos_A,pos_H22,APY_algorithm,APY_eigen_threshold,APY_n_core,direct=FALSE,inverse=TRUE,omega=SSBLUP_omega)$Hinv
H_A=NULL
H_A_inbred=NULL
}
}else{

if("H_A"%in%kinship_type&(!"H_Ainv"%in%kinship_type)){
H_Aresult=makeHA_metafounder_cpp(num_ped,input_data_numeric,pos_A11,pos_A22,pos_geno,pos_A,pos_H22,direct=TRUE,inverse=FALSE,omega=SSBLUP_omega)
H_A_inbred=H_Aresult$inbred
O0O0OOOOO0OOOOO0O0OO=data.frame(Id=IND_Additive,Diag_inbred=H_A_inbred,stringsAsFactors=F)
H_A=H_Aresult$H
H_Ainv=NULL
rm(H_Aresult);gc();
}else if("H_A"%in%kinship_type&("H_Ainv"%in%kinship_type)){
H_Aresult=makeHA_metafounder_cpp(num_ped,input_data_numeric,pos_A11,pos_A22,pos_geno,pos_A,pos_H22,direct=TRUE,inverse=TRUE,omega=SSBLUP_omega)
H_A_inbred=H_Aresult$inbred
O0O0OOOOO0OOOOO0O0OO=data.frame(Id=IND_Additive,Diag_inbred=H_A_inbred,stringsAsFactors=F)
H_A=H_Aresult$H
H_Ainv=H_Aresult$Hinv
rm(H_Aresult);gc();
}else if(!"H_A"%in%kinship_type&("H_Ainv"%in%kinship_type)){
H_Ainv=makeHA_metafounder_cpp(num_ped,input_data_numeric,pos_A11,pos_A22,pos_geno,pos_A,pos_H22,direct=FALSE,inverse=TRUE,omega=SSBLUP_omega)$Hinv
H_A=NULL
H_A_inbred=NULL
}

}


#输出矩阵文件-start
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))|(!is.null(col3_bigmemory_data_name)&!is.null(col3_bigmemory_data_path))){

if(!is.null(H_A)){
if(!is.null(output_matrix_path)){setwd(output_matrix_path)}
if(!is.null(col3_bigmemory_data_path)){setwd(col3_bigmemory_data_path)}
fwrite(data.table(IND_Additive),"IND_SSBLUP.txt",quote=F,row.names=F,col.names=F,sep="\t")

if("col_three" %in% output_matrix_type){
if(col3_bigmemory==FALSE){
O0O0OOO0O0OOOOOOOOOO=matrix_col3(H_A,IND_geno=as.numeric(IND_Additive),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
fwrite(data.table(O0O0OOO0O0OOOOOOOOOO),paste0(output_matrix_name,"H_A_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
rm(O0O0OOO0O0OOOOOOOOOO);gc();
}else{
cat(paste0("Save bigmemory-G_A_col_three matrix as ",col3_bigmemory_data_path,"/",col3_bigmemory_data_name,"_H_A_col_three ......\n"))
matrix_col3_memory_alt(H_A,paste0(col3_bigmemory_data_name,"_H_A"),col3_bigmemory_data_path,
						   IND_geno=as.numeric(IND_Additive),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
}

if("col_all" %in% output_matrix_type){
fwrite(data.frame(H_A),paste0(output_matrix_name,"H_A_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(H_A);gc();}
}

if(!is.null(H_Ainv)){
if(!is.null(output_matrix_path)){setwd(output_matrix_path)}
if(!is.null(col3_bigmemory_data_path)){setwd(col3_bigmemory_data_path)}
fwrite(data.table(IND_Additive),"IND_SSBLUP.txt",quote=F,row.names=F,col.names=F,sep="\t")

if("col_three" %in% output_matrix_type){
if(col3_bigmemory==FALSE){
O0O0OOO0O0OOOOOOOOOO=matrix_col3(H_Ainv,IND_geno=as.numeric(IND_Additive),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
fwrite(data.table(O0O0OOO0O0OOOOOOOOOO),paste0(output_matrix_name,"H_Ainv_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
rm(O0O0OOO0O0OOOOOOOOOO);gc();
}else{
cat(paste0("Save bigmemory-H_Ainv_col_three matrix as ",col3_bigmemory_data_path,"/",col3_bigmemory_data_name,"_H_Ainv_col_three ......\n"))
matrix_col3_memory_alt(H_Ainv,paste0(col3_bigmemory_data_name,"_H_Ainv"),col3_bigmemory_data_path,
						   IND_geno=as.numeric(IND_Additive),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
}

if("col_all" %in% output_matrix_type){
fwrite(data.frame(H_Ainv),paste0(output_matrix_name,"H_Ainv_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(H_Ainv);gc();}
}
}


#HD
if(gene_dropping_algorithm==FALSE){
if("H_D"%in%kinship_type&(!"H_Dinv"%in%kinship_type)){
H_Dresult=makeHD_cpp(num_ped,input_data_numeric,IND_geno,pos_A11,pos_A22,pos_geno,pos_A,pos_H22,gene_dropping_algorithm,cpu_cores,gene_dropping_iteration,direct=TRUE,inverse=TRUE,omega=SSBLUP_omega)
H_D=H_Dresult$H
H_Dinv=NULL
rm(H_Dresult);gc();
}else if("H_D"%in%kinship_type&("H_Dinv"%in%kinship_type)){
H_Dresult=makeHD_cpp(num_ped,input_data_numeric,IND_geno,pos_A11,pos_A22,pos_geno,pos_A,pos_H22,gene_dropping_algorithm,cpu_cores,gene_dropping_iteration,direct=TRUE,inverse=TRUE,omega=SSBLUP_omega)
H_D=H_Dresult$H
H_Dinv=H_Dresult$Hinv
rm(H_Dresult);gc();
}else if(!"H_D"%in%kinship_type&("H_Dinv"%in%kinship_type)){

H_Dinv=makeHD_cpp(num_ped,input_data_numeric,IND_geno,pos_A11,pos_A22,pos_geno,pos_A,pos_H22,gene_dropping_algorithm,cpu_cores,gene_dropping_iteration,direct=FALSE,inverse=TRUE,omega=SSBLUP_omega)$Hinv
H_D=NULL
}

}else{
if("H_D"%in%kinship_type&(!"H_Dinv"%in%kinship_type)){
H_Dresult=makeHD_gene_dropping_cpp(num_ped,input_data_numeric,IND_geno,pos_A11,pos_A22,pos_geno,pos_A,pos_H22,gene_dropping_algorithm,cpu_cores,gene_dropping_iteration,direct=TRUE,inverse=TRUE,omega=SSBLUP_omega)
H_D=H_Dresult$H
H_Dinv=NULL
rm(H_Dresult);gc();
}else if("H_D"%in%kinship_type&("H_Dinv"%in%kinship_type)){
H_Dresult=makeHD_gene_dropping_cpp(num_ped,input_data_numeric,IND_geno,pos_A11,pos_A22,pos_geno,pos_A,pos_H22,gene_dropping_algorithm,cpu_cores,gene_dropping_iteration,direct=TRUE,inverse=TRUE,omega=SSBLUP_omega)
H_D=H_Dresult$H
H_Dinv=H_Dresult$Hinv
rm(H_Dresult);gc();
}else if(!"H_D"%in%kinship_type&("H_Dinv"%in%kinship_type)){

H_Dinv=makeHD_gene_dropping_cpp(num_ped,input_data_numeric,IND_geno,pos_A11,pos_A22,pos_geno,pos_A,pos_H22,gene_dropping_algorithm,cpu_cores,gene_dropping_iteration,direct=FALSE,inverse=TRUE,omega=SSBLUP_omega)$Hinv
H_D=NULL
}

}

#输出矩阵文件-start
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))|(!is.null(col3_bigmemory_data_name)&!is.null(col3_bigmemory_data_path))){

if(!is.null(H_D)){
if(!is.null(output_matrix_path)){setwd(output_matrix_path)}
if(!is.null(col3_bigmemory_data_path)){setwd(col3_bigmemory_data_path)}
fwrite(data.table(IND_Additive),"IND_SSBLUP.txt",quote=F,row.names=F,col.names=F,sep="\t")

if("col_three" %in% output_matrix_type){
if(col3_bigmemory==FALSE){
O0O0OOO0O0OOOOOOOOOO=matrix_col3(H_D,IND_geno=as.numeric(IND_Additive),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
fwrite(data.table(O0O0OOO0O0OOOOOOOOOO),paste0(output_matrix_name,"H_D_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
rm(O0O0OOO0O0OOOOOOOOOO);gc();
}else{
cat(paste0("Save bigmemory-H_D_col_three matrix as ",col3_bigmemory_data_path,"/",col3_bigmemory_data_name,"_H_D_col_three ......\n"))
matrix_col3_memory_alt(H_D,paste0(col3_bigmemory_data_name,"_H_D"),col3_bigmemory_data_path,
						   IND_geno=as.numeric(IND_Additive),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
}

if("col_all" %in% output_matrix_type){
fwrite(data.frame(H_D),paste0(output_matrix_name,"H_D_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(H_D);gc();}
}

if(!is.null(H_Dinv)){
if(!is.null(output_matrix_path)){setwd(output_matrix_path)}
if(!is.null(col3_bigmemory_data_path)){setwd(col3_bigmemory_data_path)}
fwrite(data.table(IND_Additive),"IND_SSBLUP.txt",quote=F,row.names=F,col.names=F,sep="\t")

if("col_three" %in% output_matrix_type){
if(col3_bigmemory==FALSE){
O0O0OOO0O0OOOOOOOOOO=matrix_col3(H_Dinv,IND_geno=as.numeric(IND_Additive),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
fwrite(data.table(O0O0OOO0O0OOOOOOOOOO),paste0(output_matrix_name,"H_Dinv_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
rm(O0O0OOO0O0OOOOOOOOOO);gc();
}else{
cat(paste0("Save bigmemory-H_Dinv_col_three matrix as ",col3_bigmemory_data_path,"/",col3_bigmemory_data_name,"_H_Dinv_col_three ......\n"))
matrix_col3_memory_alt(H_Dinv,paste0(col3_bigmemory_data_name,"_H_Dinv"),col3_bigmemory_data_path,
						   IND_geno=as.numeric(IND_Additive),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
}

if("col_all" %in% output_matrix_type){
fwrite(data.frame(H_Dinv),paste0(output_matrix_name,"H_Dinv_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(H_Dinv);gc();}
}
}



#输出矩阵文件-end

#Part2 计算inbreeding coefficients
if("Homozygous"%in%inbred_type){
cat("Start calculate genomic inbreeding coefficients by homozygous site \n")
homo_inbred<-cal_homo_inbred_cpp(input_data_numeric,as.character(IND_geno))
if(!is.null(inbred_output_matrix_path)){
setwd(inbred_output_matrix_path)
cat("Saving inbreeding coefficients......\n")
fwrite(data.table(homo_inbred),"Homozygous_inbreeding_coefficients.txt",quote=F,sep="\t",row.names=F,col.names=T)
}
}else{
homo_inbred=NULL
}

if("G_diag"%in%inbred_type){
cat("Start calculate genomic inbreeding coefficients by diagnoal of genomic additive relationship matrix  \n")
if(!is.null(inbred_output_matrix_path)){
setwd(inbred_output_matrix_path)
cat("Saving inbreeding coefficients......\n")
fwrite(data.table(diag_inbred),"Diagnoal_inbreeding_coefficients.txt",quote=F,sep="\t",row.names=F,col.names=T)
}
}else{
diag_inbred=NULL
}

if("Pedigree"%in%inbred_type){
cat("Start calculate pedigree_inbreeding coefficients  \n")
pedigree_inbred=data.frame(Id=IND_pedigree,Diag_inbred=diag(makeInbreeding_cpp(num_ped))-1,stringsAsFactors=F)
if(!is.null(inbred_output_matrix_path)){
setwd(inbred_output_matrix_path)
cat("Saving inbreeding coefficients......\n")
fwrite(data.table(pedigree_inbred),"Pedigree_inbreeding_coefficients.txt",quote=F,sep="\t",row.names=F,col.names=T)
}
}else{
pedigree_inbred=NULL
}

if("H_diag"%in%inbred_type){
if(Metafounder_algorithm==TRUE){cat("Start calculate  inbreeding coefficients by single-step metafounder method  \n")
}else {cat("Start calculate  inbreeding coefficients by single-step  method  \n")}
#O0O0OOOOO0OOOOO0O0OO=data.frame(Id=c(IND_A11,IND_A22),H_inbreeding=O0O0OOOOO0OOOOO0O0OO,stringsAsFactors=F)
if(!is.null(inbred_output_matrix_path)){
setwd(inbred_output_matrix_path)
cat("Saving inbreeding coefficients......\n")
if(Metafounder_algorithm==TRUE){
fwrite(data.table(O0O0OOOOO0OOOOO0O0OO),"H_metafounder_inbreeding_coefficients.txt",quote=F,sep="\t",row.names=F,col.names=T)
}else {
fwrite(data.table(O0O0OOOOO0OOOOO0O0OO),"H_inbreeding_coefficients.txt",quote=F,sep="\t",row.names=F,col.names=T)}
}
}else{
O0O0OOOOO0OOOOO0O0OO=NULL
}

if(!exists("G_A")){G_A=NULL};if(!exists("G_Ainv")){G_Ainv=NULL};
if(!exists("G_D")){G_D=NULL};if(!exists("G_Dinv")){G_Dinv=NULL};
if(!exists("P_A")){P_A=NULL};if(!exists("P_Ainv")){P_Ainv=NULL};
if(!exists("P_D")){P_D=NULL};if(!exists("P_Dinv")){P_Dinv=NULL};
if(!exists("H_A")){H_A=NULL};if(!exists("H_Ainv")){H_Ainv=NULL};

if(!is.null(G_A)){rownames(G_A)=colnames(G_A)=IND_geno}
if(!is.null(G_Ainv)){rownames(G_Ainv)=colnames(G_Ainv)=IND_geno}
if(!is.null(G_D)){rownames(G_D)=colnames(G_D)=IND_geno}
if(!is.null(G_Dinv)){rownames(G_Dinv)=colnames(G_Dinv)=IND_geno}
if(!is.null(P_A)){rownames(P_A)=colnames(P_A)=IND_pedigree}
if(!is.null(P_Ainv)){rownames(P_Ainv)=colnames(P_Ainv)=IND_pedigree}
if(!is.null(P_D)){rownames(P_D)=colnames(P_D)=IND_pedigree}
if(!is.null(P_Dinv)){rownames(P_Dinv)=colnames(P_Dinv)=IND_pedigree}
if(!is.null(H_A)){rownames(H_A)=colnames(H_A)=IND_Additive}
if(!is.null(H_Ainv)){rownames(H_Ainv)=colnames(H_Ainv)=IND_Additive}


if(return_result==TRUE){
return(list(G_A=list(A=G_A,Ainv=G_Ainv),G_D=list(D=G_D,Dinv=G_Dinv),P_A=list(A=P_A,Ainv=P_Ainv),P_D=list(D=P_D,Dinv=P_Dinv),
		  H_A=list(A=H_A,Ainv=H_Ainv),
		    Inbred=list(Homozygous=homo_inbred,G_diag=diag_inbred,Pedigree=pedigree_inbred,H_diag=O0O0OOOOO0OOOOO0O0OO),
			phased_block=data_block))
}
}else if(bigmemory_cal==TRUE){



###################################big memory##########################
delete_bigmemory_file("double_numeric",bigmemory_data_name,bigmemory_data_path,TRUE);
#Part1 计算kinship 
########################################## P_A ##########################################

#删除bigmemory-文件
if("P_A"%in%kinship_type){
delete_bigmemory_file("P_A",bigmemory_data_name,bigmemory_data_path,TRUE);
}

if("P_Ainv"%in%kinship_type){
delete_bigmemory_file("P_Ainv",bigmemory_data_name,bigmemory_data_path,TRUE);
}	


if("P_A"%in%kinship_type&(!"P_Ainv"%in%kinship_type)){
P_A=K_matrix_cal_memory(numeric_address=NULL,bigmemory_data_name,bigmemory_data_path,
					    IND_geno=c("1","2"),Pedigree=num_ped,pos_A11=c(1,2),pos_A22=c(1,2),
					    pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
				        inverse=FALSE,type=7)$A
P_Ainv=NULL
}else if("P_A"%in%kinship_type&("P_Ainv"%in%kinship_type)){
P_A=K_matrix_cal_memory(numeric_address=NULL,bigmemory_data_name,bigmemory_data_path,
					    IND_geno=c("1","2"),Pedigree=num_ped,pos_A11=c(1,2),pos_A22=c(1,2),
					    pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
				        type=7)$A
P_Ainv=K_matrix_cal_memory(numeric_address=NULL,bigmemory_data_name,bigmemory_data_path,
					       IND_geno=c("1","2"),Pedigree=num_ped,pos_A11=c(1,2),pos_A22=c(1,2),
						  pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
						  type=8)$Ainv

}else if(!"P_A"%in%kinship_type&("P_Ainv"%in%kinship_type)){
P_A=NULL
P_Ainv=K_matrix_cal_memory(numeric_address=NULL,bigmemory_data_name,bigmemory_data_path,
					       IND_geno=c("1","2"),Pedigree=num_ped,pos_A11=c(1,2),pos_A22=c(1,2),
						  pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
						  type=8)$Ainv
}


#输出矩阵文件-start
if(!is.null(output_matrix_path)&!is.null(output_matrix_name)){

if(!is.null(P_A)){
setwd(output_matrix_path)
fwrite(data.table(IND_pedigree),"IND_pedigree.txt",quote=F,row.names=F,col.names=F,sep="\t")
if("col_three" %in% output_matrix_type){
cat(paste0("Save bigmemory-P_A_col_three matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_P_A_col_three ......\n"))
matrix_col3_memory(P_A@address,paste0(bigmemory_data_name,"_P_A"),bigmemory_data_path,
						   IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
if("col_all" %in% output_matrix_type){
cat(paste0("Save bigmemory-P_A_col_all matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_P_A ......\n"))
}
}

if(!is.null(P_Ainv)){
setwd(output_matrix_path)
fwrite(data.table(IND_pedigree),"IND_pedigree.txt",quote=F,row.names=F,col.names=F,sep="\t")
if("col_three" %in% output_matrix_type){
cat(paste0("Save bigmemory-P_Ainv_col_three matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_P_Ainv_col_three ......\n"))
matrix_col3_memory(pBigMat=P_Ainv@address,paste0(bigmemory_data_name,"_P_Ainv"),bigmemory_data_path,
						   IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
if("col_all" %in% output_matrix_type){
cat(paste0("Save bigmemory-P_Ainv_col_all matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_P_Ainv ......\n"))
}
}

}
#输出矩阵文件-end


########################################## P_D ##########################################
#删除bigmemory-文件
if("P_D"%in%kinship_type&(!"P_Dinv"%in%kinship_type)){
delete_bigmemory_file("P_D",bigmemory_data_name,bigmemory_data_path,TRUE);
delete_bigmemory_file("P_A_temp",bigmemory_data_name,bigmemory_data_path,FALSE);
}else if("P_D"%in%kinship_type&("P_Dinv"%in%kinship_type)){
delete_bigmemory_file("P_D",bigmemory_data_name,bigmemory_data_path,TRUE);
delete_bigmemory_file("P_A_temp",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("P_Dinv",bigmemory_data_name,bigmemory_data_path,TRUE);	
}else if(!"P_D"%in%kinship_type&("P_Dinv"%in%kinship_type)){
delete_bigmemory_file("P_D",bigmemory_data_name,bigmemory_data_path,TRUE);
delete_bigmemory_file("P_A_temp",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("P_Dinv",bigmemory_data_name,bigmemory_data_path,TRUE);	
}	


if(gene_dropping_algorithm==FALSE){

if("P_D"%in%kinship_type&(!"P_Dinv"%in%kinship_type)){
P_D=K_matrix_cal_memory(numeric_address=NULL,bigmemory_data_name,bigmemory_data_path,
					       IND_geno=c("1","2"),Pedigree=num_ped,pos_A11=c(1,2),pos_A22=c(1,2),
						  pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
						  type=9,inverse=FALSE)$D
P_Dinv=NULL
}else if("P_D"%in%kinship_type&("P_Dinv"%in%kinship_type)){
P_Dresult=K_matrix_cal_memory(numeric_address=NULL,bigmemory_data_name,bigmemory_data_path,
					       IND_geno=c("1","2"),Pedigree=num_ped,pos_A11=c(1,2),pos_A22=c(1,2),
						  pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
						  type=9,inverse=TRUE)
P_D=P_Dresult$D
P_Dinv=P_Dresult$Dinv
rm(P_Dresult);gc();
}else if(!"P_D"%in%kinship_type&("P_Dinv"%in%kinship_type)){
P_Dinv=K_matrix_cal_memory(numeric_address=NULL,bigmemory_data_name,bigmemory_data_path,
					       IND_geno=c("1","2"),Pedigree=num_ped,pos_A11=c(1,2),pos_A22=c(1,2),
						   pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
						   type=9,inverse=TRUE)$Dinv
P_D=NULL
}

}else{

if("P_D"%in%kinship_type&(!"P_Dinv"%in%kinship_type)){
P_D=K_matrix_cal_memory(numeric_address=NULL,bigmemory_data_name,bigmemory_data_path,
					       IND_geno=c("1","2"),Pedigree=num_ped,pos_A11=c(1,2),pos_A22=c(1,2),
						  pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
						  type=10,inverse=FALSE,cpu_cores=cpu_cores,
						  gene_dropping_iteration=gene_dropping_iteration)$D
P_Dinv=NULL
}else if("P_D"%in%kinship_type&("P_Dinv"%in%kinship_type)){
P_Dresult=K_matrix_cal_memory(numeric_address=NULL,bigmemory_data_name,bigmemory_data_path,
					       IND_geno=c("1","2"),Pedigree=num_ped,pos_A11=c(1,2),pos_A22=c(1,2),
						  pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
						  type=10,inverse=TRUE,cpu_cores=cpu_cores,
						  gene_dropping_iteration=gene_dropping_iteration)
P_D=P_Dresult$D
P_Dinv=P_Dresult$Dinv
rm(P_Dresult);gc();
}else if(!"P_D"%in%kinship_type&("P_Dinv"%in%kinship_type)){
P_Dinv=K_matrix_cal_memory(numeric_address=NULL,bigmemory_data_name,bigmemory_data_path,
					       IND_geno=c("1","2"),Pedigree=num_ped,pos_A11=c(1,2),pos_A22=c(1,2),
						   pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
						   type=10,inverse=TRUE,cpu_cores=cpu_cores,
						  gene_dropping_iteration=gene_dropping_iteration)$Dinv
P_D=NULL
}

}

#删除bigmemory-文件
if("P_D"%in%kinship_type|("P_Dinv"%in%kinship_type)){
delete_bigmemory_file("P_A_temp",bigmemory_data_name,bigmemory_data_path,FALSE);	
}
#输出矩阵文件-start
if(!is.null(output_matrix_path)&!is.null(output_matrix_name)){


if(!is.null(P_D)){
setwd(output_matrix_path)
fwrite(data.table(IND_pedigree),"IND_pedigree.txt",quote=F,row.names=F,col.names=F,sep="\t")
if("col_three" %in% output_matrix_type){
cat(paste0("Save bigmemory-P_D_col_three matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_P_D_col_three ......\n"))
matrix_col3_memory(pBigMat=P_D@address,paste0(bigmemory_data_name,"_P_D"),bigmemory_data_path,
						   IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
if("col_all" %in% output_matrix_type){
cat(paste0("Save bigmemory-P_D_col_all matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_P_D ......\n"))
}
}

if(!is.null(P_Dinv)){
setwd(output_matrix_path)
fwrite(data.table(IND_pedigree),"IND_pedigree.txt",quote=F,row.names=F,col.names=F,sep="\t")
if("col_three" %in% output_matrix_type){
cat(paste0("Save bigmemory-P_Dinv_col_three matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_P_Dinv_col_three ......\n"))
matrix_col3_memory(pBigMat=P_Dinv@address,paste0(bigmemory_data_name,"_P_Dinv"),bigmemory_data_path,
						   IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
if("col_all" %in% output_matrix_type){
cat(paste0("Save bigmemory-P_Dinv_col_all matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_P_Dinv ......\n"))
}
}

}
#输出矩阵文件-end

#删除bigmemory-文件
if("G_A"%in%kinship_type&(!"G_Ainv"%in%kinship_type)){
delete_bigmemory_file("temp_numeric",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_numeric",bigmemory_data_name,bigmemory_data_path,FALSE);
}else if("G_A"%in%kinship_type&("G_Ainv"%in%kinship_type)){
delete_bigmemory_file("temp_numeric",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_numeric",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("G_Ainv",bigmemory_data_name,bigmemory_data_path,TRUE);	
}else if(!"G_A"%in%kinship_type&("G_Ainv"%in%kinship_type)){
delete_bigmemory_file("temp_numeric",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_numeric",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("G_Ainv",bigmemory_data_name,bigmemory_data_path,TRUE);	
}		


########################################## G_A ##########################################
if(APY_algorithm==FALSE){
if("G_A"%in%kinship_type&(!"G_Ainv"%in%kinship_type)){
G_A=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					    IND_geno=IND_geno,Pedigree=matrix(1L,1,1),pos_A11=c(1,2),pos_A22=c(1,2),
					    pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
					    type=1,base=kinship_base,trace=kinship_trace,metafounder=FALSE,inverse=FALSE)$G
	
diag_inbred=data.frame(Id=IND_geno,Diag_inbred=diag(G_A[])-1,stringsAsFactors=F)
G_Ainv=NULL
}else if("G_A"%in%kinship_type&("G_Ainv"%in%kinship_type)){
G_Aresult=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					    IND_geno=IND_geno,Pedigree=matrix(1L,1,1),pos_A11=c(1,2),pos_A22=c(1,2),
					    pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
					    type=1,base=kinship_base,trace=kinship_trace,metafounder=FALSE,inverse=TRUE)


G_A=G_Aresult$G
diag_inbred=data.frame(Id=IND_geno,Diag_inbred=diag(G_A[])-1,stringsAsFactors=F)
G_Ainv=G_Aresult$Ginv
rm(G_Aresult);gc();
}else if(!"G_A"%in%kinship_type&("G_Ainv"%in%kinship_type)){

G_Ainv=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					    IND_geno=IND_geno,Pedigree=matrix(1L,1,1),pos_A11=c(1,2),pos_A22=c(1,2),
					    pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
					    type=1,base=kinship_base,trace=kinship_trace,metafounder=FALSE,inverse=TRUE)$Ginv
G_A=NULL

}

}else{

if("G_A"%in%kinship_type&(!"G_Ainv"%in%kinship_type)){
G_A=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					    IND_geno=IND_geno,Pedigree=matrix(1L,1,1),pos_A11=c(1,2),pos_A22=c(1,2),
					    pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
						APY_eigen_threshold=APY_eigen_threshold,
						APY_n_core=APY_n_core,type=4,re_inverse=TRUE)$G
G_Ainv=NULL
}else if("G_A"%in%kinship_type&("G_Ainv"%in%kinship_type)){
G_Aresult=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					    IND_geno=IND_geno,Pedigree=matrix(1L,1,1),pos_A11=c(1,2),pos_A22=c(1,2),
					    pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
						APY_eigen_threshold=APY_eigen_threshold,
						APY_n_core=APY_n_core,type=4,re_inverse=TRUE)
G_A=G_Aresult$G
G_Ainv=G_Aresult$Ginv
rm(G_Aresult);gc();
}else if(!"G_A"%in%kinship_type&("G_Ainv"%in%kinship_type)){
G_Ainv=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					    IND_geno=IND_geno,Pedigree=matrix(1L,1,1),pos_A11=c(1,2),pos_A22=c(1,2),
					    pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
						APY_eigen_threshold=APY_eigen_threshold,
						APY_n_core=APY_n_core,type=4,re_inverse=FALSE)$Ginv
G_A=NULL
}
}

#删除bigmemory-文件
if("G_A"%in%kinship_type|("G_Ainv"%in%kinship_type)){
delete_bigmemory_file("temp_numeric",bigmemory_data_name,bigmemory_data_path,FALSE);			
}

#输出矩阵文件-start
if(!is.null(output_matrix_path)&!is.null(output_matrix_name)){

if(!is.null(G_A)){
setwd(output_matrix_path)
fwrite(data.table(IND_geno),"IND_geno.txt",quote=F,row.names=F,col.names=F,sep="\t")

if("col_three" %in% output_matrix_type){
cat(paste0("Save bigmemory-G_A_col_three matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_G_A_col_three ......\n"))
matrix_col3_memory(pBigMat=G_A@address,paste0(bigmemory_data_name,"_G_A"),bigmemory_data_path,
						   IND_geno=as.numeric(IND_geno),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
if("col_all" %in% output_matrix_type){
cat(paste0("Save bigmemory-G_A_col_all matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_G_A ......\n"))
}
}

if(!is.null(G_Ainv)){
setwd(output_matrix_path)
fwrite(data.table(IND_geno),"IND_geno.txt",quote=F,row.names=F,col.names=F,sep="\t")

if("col_three" %in% output_matrix_type){
cat(paste0("Save bigmemory-G_Ainv_col_three matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_G_Ainv_col_three ......\n"))
matrix_col3_memory(pBigMat=G_Ainv@address,paste0(bigmemory_data_name,"_G_Ainv"),bigmemory_data_path,
						   IND_geno=as.numeric(IND_geno),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
if("col_all" %in% output_matrix_type){
cat(paste0("Save bigmemory-G_Ainv_col_all matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_G_Ainv ......\n"))
}
}

}
#输出矩阵文件-end

#删除bigmemory-文件
if("G_D"%in%kinship_type&(!"G_Dinv"%in%kinship_type)){
delete_bigmemory_file("temp_numeric",bigmemory_data_name,bigmemory_data_path,FDLSE);
delete_bigmemory_file("temp_numeric",bigmemory_data_name,bigmemory_data_path,FDLSE);
}else if("G_D"%in%kinship_type&("G_Dinv"%in%kinship_type)){
delete_bigmemory_file("temp_numeric",bigmemory_data_name,bigmemory_data_path,FDLSE);
delete_bigmemory_file("temp_numeric",bigmemory_data_name,bigmemory_data_path,FDLSE);
delete_bigmemory_file("G_Dinv",bigmemory_data_name,bigmemory_data_path,TRUE);	
}else if(!"G_D"%in%kinship_type&("G_Dinv"%in%kinship_type)){
delete_bigmemory_file("temp_numeric",bigmemory_data_name,bigmemory_data_path,FDLSE);
delete_bigmemory_file("temp_numeric",bigmemory_data_name,bigmemory_data_path,FDLSE);
delete_bigmemory_file("G_Dinv",bigmemory_data_name,bigmemory_data_path,TRUE);	
}	

########################################## G_D ##########################################
if("genotypic"%in%dominance_type){
if("G_D"%in%kinship_type&(!"G_Dinv"%in%kinship_type)){
G_D=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					    IND_geno=IND_geno,Pedigree=matrix(1L,1,1),pos_A11=c(1,2),pos_A22=c(1,2),
					    pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
					    type=3,base=kinship_base,trace=kinship_trace,inverse=FALSE)$D
G_Dinv=NULL
}else if("G_D"%in%kinship_type&("G_Dinv"%in%kinship_type)){
G_Dresult=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					    IND_geno=IND_geno,Pedigree=matrix(1L,1,1),pos_A11=c(1,2),pos_A22=c(1,2),
					    pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
					    type=3,base=kinship_base,trace=kinship_trace,inverse=TRUE)
G_D=G_Dresult$D
G_Dinv=G_Dresult$Dinv
rm(G_Dresult);gc();
}else if(!"G_D"%in%kinship_type&("G_Dinv"%in%kinship_type)){
G_Dinv=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					    IND_geno=IND_geno,Pedigree=matrix(1L,1,1),pos_A11=c(1,2),pos_A22=c(1,2),
					    pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
					    type=3,base=kinship_base,trace=kinship_trace,inverse=TRUE)$Dinv
G_D=NULL
}
}else if("classical"%in%dominance_type){

if("G_D"%in%kinship_type&(!"G_Dinv"%in%kinship_type)){
G_D=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					    IND_geno=IND_geno,Pedigree=matrix(1L,1,1),pos_A11=c(1,2),pos_A22=c(1,2),
					    pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
					    type=2,base=kinship_base,trace=kinship_trace,inverse=FALSE)$D
G_Dinv=NULL
}else if("G_D"%in%kinship_type&("G_Dinv"%in%kinship_type)){
G_Dresult=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					    IND_geno=IND_geno,Pedigree=matrix(1L,1,1),pos_A11=c(1,2),pos_A22=c(1,2),
					    pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
					    type=2,base=kinship_base,trace=kinship_trace,inverse=TRUE)
G_D=G_Dresult$D
G_Dinv=G_Dresult$Dinv
rm(G_Dresult);gc();
}else if(!"G_D"%in%kinship_type&("G_Dinv"%in%kinship_type)){
G_Dinv=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					    IND_geno=IND_geno,Pedigree=matrix(1L,1,1),pos_A11=c(1,2),pos_A22=c(1,2),
					    pos_geno=c(1,2),pos_A=c(1,2),pos_H22=c(1,2),
					    type=2,base=kinship_base,trace=kinship_trace,inverse=TRUE)$Dinv
G_D=NULL
}
}

#删除bigmemory-文件
if("G_D"%in%kinship_type|("G_Dinv"%in%kinship_type)){
delete_bigmemory_file("temp_numeric",bigmemory_data_name,bigmemory_data_path,FALSE);			
}

#输出矩阵文件-start
if(!is.null(output_matrix_path)&!is.null(output_matrix_name)){


if(!is.null(G_D)){
setwd(output_matrix_path)
fwrite(data.table(IND_geno),"IND_geno.txt",quote=F,row.names=F,col.names=F,sep="\t")

if("col_three" %in% output_matrix_type){
cat(paste0("Save bigmemory-G_D_col_three matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_G_D_col_three ......\n"))
matrix_col3_memory(pBigMat=G_D@address,paste0(bigmemory_data_name,"_G_D"),bigmemory_data_path,
						   IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
if("col_all" %in% output_matrix_type){
cat(paste0("Save bigmemory-G_D_col_all matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_G_D ......\n"))
}
}

if(!is.null(G_Dinv)){
setwd(output_matrix_path)
fwrite(data.table(IND_geno),"IND_geno.txt",quote=F,row.names=F,col.names=F,sep="\t")

if("col_three" %in% output_matrix_type){
cat(paste0("Save bigmemory-G_Dinv_col_three matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_G_Dinv_col_three ......\n"))
matrix_col3_memory(pBigMat=G_Dinv@address,paste0(bigmemory_data_name,"_G_Dinv"),bigmemory_data_path,
						   IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
if("col_all" %in% output_matrix_type){
cat(paste0("Save bigmemory-G_Dinv_col_all matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_G_Dinv ......\n"))
}
}

}
#输出矩阵文件-end


########################################## single-step ##########################################
if("H_A"%in%kinship_type|"H_Ainv"%in%kinship_type){


delete_bigmemory_file("temp_P_A",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_A22",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_O0OOO0O0OOO0OOO0OOO0",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_M_new",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_A11",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_A12",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_G_A",bigmemory_data_name,bigmemory_data_path,FALSE);

if("H_A"%in%kinship_type){
delete_bigmemory_file("H_A",bigmemory_data_name,bigmemory_data_path,TRUE);
}


if("H_Ainv"%in%kinship_type){
delete_bigmemory_file("H_Ainv",bigmemory_data_name,bigmemory_data_path,TRUE);
delete_bigmemory_file("temp_G_Ainv",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_P_Ainv",bigmemory_data_name,bigmemory_data_path,FALSE);
}
}


if("H_A"%in%kinship_type|"H_Ainv"%in%kinship_type){
IND_geno2=setdiff(IND_geno,IND_pedigree)    #有基因型但是无系谱
IND_A22=intersect(IND_pedigree,IND_geno)    #既有基因型又有系谱
IND_A11=setdiff(IND_pedigree,IND_A22)        #有系谱但是无基因型
IND_Additive=c(IND_A11,IND_A22)               #在H矩阵中的个体
pos_geno=match(IND_A22,IND_geno)-1            #既有基因型又有系谱的个体在 基因型数据中的位置
pos_A11=match(IND_A11,IND_pedigree)-1         #有系谱但是无基因型个体在  在系谱中的位置
pos_A22=match(IND_A22,IND_pedigree)-1         #既有基因型又有系谱的个体在  在系谱中的位置
pos_A=match(IND_Additive,IND_pedigree)-1      #H矩阵的个体  在系谱中的位置 
pos_H22=match(IND_A22,IND_Additive)-1         #既有基因型又有系谱的个体在  在H矩阵中的位置
if(length(IND_A22)==0){stop("No aminals both with genotype and_pedigree, Please check you data ! \n")
}else{cat(length(IND_A22),"individuals both have genotype and_pedigree \n")}
if(length(IND_geno2)>0){cat(paste0(length(IND_geno2)," individuals have genotype but not in the_pedigree, these individuals are removed ! \n"))}
}

#return(list(num_ped=num_ped,input_data_numeric=input_data_numeric,pos_geno=pos_geno,
#			pos_A11=pos_A11,pos_A22=pos_A22,pos_A=pos_A,pos_H22=pos_H22))

if(APY_algorithm==TRUE){

if("H_A"%in%kinship_type&(!"H_Ainv"%in%kinship_type)){
re_inverse=TRUE
H_A_direct=TRUE
inverse=FALSE
}else if("H_A"%in%kinship_type&("H_Ainv"%in%kinship_type)){
re_inverse=TRUE
H_A_direct=TRUE
inverse=TRUE
}else if(!"H_A"%in%kinship_type&("H_Ainv"%in%kinship_type)){
re_inverse=FALSE
H_A_direct=FALSE
inverse=TRUE
}

}else{

if("H_A"%in%kinship_type&(!"H_Ainv"%in%kinship_type)){
re_inverse=TRUE
H_A_direct=TRUE
inverse=FALSE
}else if("H_A"%in%kinship_type&("H_Ainv"%in%kinship_type)){
re_inverse=TRUE
H_A_direct=TRUE
inverse=TRUE
}else if(!"H_A"%in%kinship_type&("H_Ainv"%in%kinship_type)){
re_inverse=FALSE
H_A_direct=FALSE
inverse=TRUE
}

}



if(Metafounder_algorithm==FALSE){
if("H_A"%in%kinship_type&(!"H_Ainv"%in%kinship_type)){
H_Aresult=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					IND_geno=IND_geno,Pedigree=num_ped,pos_A11=pos_A11,pos_A22=pos_A22,
					pos_geno=pos_geno,pos_A=pos_A,pos_H22=pos_H22,
					H_A_direct=H_A_direct,omega=SSBLUP_omega,
				     base=kinship_base,trace=kinship_trace,
				     type=5,cpu_cores=cpu_cores,APY_eigen_threshold=APY_eigen_threshold,APY_algorithm=APY_algorithm,
					 APY_n_core=APY_n_core,re_inverse=re_inverse,inverse=inverse)
H_A_inbred=H_Aresult$inbred
O0O0OOOOO0OOOOO0O0OO=data.frame(Id=IND_Additive,Diag_inbred=H_A_inbred,stringsAsFactors=F)
H_A=H_Aresult$H
H_Ainv=NULL
rm(H_Aresult);gc();
}else if("H_A"%in%kinship_type&("H_Ainv"%in%kinship_type)){
H_Aresult=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					IND_geno=IND_geno,Pedigree=num_ped,pos_A11=pos_A11,pos_A22=pos_A22,
					pos_geno=pos_geno,pos_A=pos_A,pos_H22=pos_H22,
					H_A_direct=H_A_direct,omega=SSBLUP_omega,
				     base=kinship_base,trace=kinship_trace,
				     type=5,cpu_cores=cpu_cores,APY_eigen_threshold=APY_eigen_threshold,APY_algorithm=APY_algorithm,
					 APY_n_core=APY_n_core,re_inverse=re_inverse,inverse=inverse)
H_A_inbred=H_Aresult$inbred
O0O0OOOOO0OOOOO0O0OO=data.frame(Id=IND_Additive,Diag_inbred=H_A_inbred,stringsAsFactors=F)
H_A=H_Aresult$H
H_Ainv=H_Aresult$Hinv
rm(H_Aresult);gc();
}else if(!"H_A"%in%kinship_type&("H_Ainv"%in%kinship_type)){
H_Ainv=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					IND_geno=IND_geno,Pedigree=num_ped,pos_A11=pos_A11,pos_A22=pos_A22,
					pos_geno=pos_geno,pos_A=pos_A,pos_H22=pos_H22,
					H_A_direct=H_A_direct,omega=SSBLUP_omega,
				     base=kinship_base,trace=kinship_trace,
				     type=5,cpu_cores=cpu_cores,APY_eigen_threshold=APY_eigen_threshold,APY_algorithm=APY_algorithm,
					 APY_n_core=APY_n_core,re_inverse=re_inverse,inverse=inverse)$Hinv
H_A=NULL
H_A_inbred=NULL
}
}else{

if("H_A"%in%kinship_type&(!"H_Ainv"%in%kinship_type)){
H_Aresult=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					IND_geno=IND_geno,Pedigree=num_ped,pos_A11=pos_A11,pos_A22=pos_A22,
					pos_geno=pos_geno,pos_A=pos_A,pos_H22=pos_H22,
					H_A_direct=TRUE,omega=SSBLUP_omega,type=6,inverse=FALSE)
					
H_A_inbred=H_Aresult$inbred
O0O0OOOOO0OOOOO0O0OO=data.frame(Id=IND_Additive,Diag_inbred=H_A_inbred,stringsAsFactors=F)
H_A=H_Aresult$H
H_Ainv=NULL
rm(H_Aresult);gc();
}else if("H_A"%in%kinship_type&("H_Ainv"%in%kinship_type)){
H_Aresult=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					IND_geno=IND_geno,Pedigree=num_ped,pos_A11=pos_A11,pos_A22=pos_A22,
					pos_geno=pos_geno,pos_A=pos_A,pos_H22=pos_H22,
					H_A_direct=TRUE,omega=SSBLUP_omega,type=6,inverse=TRUE)
H_A_inbred=H_Aresult$inbred
O0O0OOOOO0OOOOO0O0OO=data.frame(Id=IND_Additive,Diag_inbred=H_A_inbred,stringsAsFactors=F)
H_A=H_Aresult$H
H_Ainv=H_Aresult$Hinv
rm(H_Aresult);gc();
}else if(!"H_A"%in%kinship_type&("H_Ainv"%in%kinship_type)){
H_Ainv=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					IND_geno=IND_geno,Pedigree=num_ped,pos_A11=pos_A11,pos_A22=pos_A22,
					pos_geno=pos_geno,pos_A=pos_A,pos_H22=pos_H22,
					H_A_direct=FALSE,omega=SSBLUP_omega,type=6,inverse=TRUE)$Hinv
H_A=NULL
H_A_inbred=NULL
}

}

if("H_A"%in%kinship_type|"H_Ainv"%in%kinship_type){
delete_bigmemory_file("temp_P_A",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_A22",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_O0OOO0O0OOO0OOO0OOO0",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_A11",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_A12",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_M_new",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_G_A",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_G_Ainv",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("temp_P_Ainv",bigmemory_data_name,bigmemory_data_path,FALSE);	
}

#输出矩阵文件-start
if(!is.null(output_matrix_path)&!is.null(output_matrix_name)){


if(!is.null(H_A)){
setwd(output_matrix_path)
fwrite(data.table(IND_Additive),"IND_SSBLUP.txt",quote=F,row.names=F,col.names=F,sep="\t")

if("col_three" %in% output_matrix_type){
cat(paste0("Save bigmemory-H_A_col_three matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_H_A_col_three ......\n"))
matrix_col3_memory(pBigMat=H_A@address,bigmemory_data_name=paste0(bigmemory_data_name,"_H_A"),bigmemory_data_path,
						   IND_geno=as.numeric(IND_Additive),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
if("col_all" %in% output_matrix_type){
cat(paste0("Save bigmemory-H_A_col_all matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_H_A ......\n"))
}
}

if(!is.null(H_Ainv)){
setwd(output_matrix_path)
fwrite(data.table(IND_Additive),"IND_SSBLUP.txt",quote=F,row.names=F,col.names=F,sep="\t")

if("col_three" %in% output_matrix_type){
cat(paste0("Save bigmemory-H_Ainv_col_three matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_H_Ainv_col_three ......\n"))
matrix_col3_memory(pBigMat=H_Ainv@address,bigmemory_data_name=paste0(bigmemory_data_name,"_H_Ainv"),bigmemory_data_path,
						   IND_geno=as.numeric(IND_Additive),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
}
if("col_all" %in% output_matrix_type){
cat(paste0("Save bigmemory-H_Ainv_col_all matrix as ",bigmemory_data_path,"/",bigmemory_data_name,"_H_Ainv ......\n"))
}
}

}
#输出矩阵文件-end


#Part2 计算inbreeding coefficients
if("Homozygous"%in%inbred_type){
cat("Start calculate genomic inbreeding coefficients by homozygous site \n")
homo_inbred<-cal_homo_inbred_memory_cpp(input_data_numeric@address,as.character(IND_geno))
if(!is.null(inbred_output_matrix_path)){
setwd(inbred_output_matrix_path)
cat("Saving inbreeding coefficients......\n")
fwrite(data.table(homo_inbred),"Homozygous_inbreeding_coefficients.txt",quote=F,sep="\t",row.names=F,col.names=T)
}
}else{
homo_inbred=NULL
}

if("G_diag"%in%inbred_type){
cat("Start calculate genomic inbreeding coefficients by diagnoal of genomic additive relationship matrix  \n")
if(!is.null(inbred_output_matrix_path)){
setwd(inbred_output_matrix_path)
cat("Saving inbreeding coefficients......\n")
fwrite(data.table(diag_inbred),"Diagnoal_inbreeding_coefficients.txt",quote=F,sep="\t",row.names=F,col.names=T)
}
}else{
diag_inbred=NULL
}

if("Pedigree"%in%inbred_type){
cat("Start calculate pedigree_inbreeding coefficients  \n")
pedigree_inbred=data.frame(Id=IND_pedigree,Diag_inbred=diag(makeInbreeding_memory_cpp(num_ped))-1,stringsAsFactors=F)
if(!is.null(inbred_output_matrix_path)){
setwd(inbred_output_matrix_path)
cat("Saving inbreeding coefficients......\n")
fwrite(data.table(pedigree_inbred),"Pedigree_inbreeding_coefficients.txt",quote=F,sep="\t",row.names=F,col.names=T)
}
}else{
pedigree_inbred=NULL
}

if("H_diag"%in%inbred_type){
if(Metafounder_algorithm==TRUE){cat("Start calculate  inbreeding coefficients by single-step metafounder method  \n")
}else {cat("Start calculate  inbreeding coefficients by single-step  method  \n")}
#O0O0OOOOO0OOOOO0O0OO=data.frame(Id=c(IND_A11,IND_A22),H_inbreeding=O0O0OOOOO0OOOOO0O0OO,stringsAsFactors=F)
if(!is.null(inbred_output_matrix_path)){
setwd(inbred_output_matrix_path)
cat("Saving inbreeding coefficients......\n")
if(Metafounder_algorithm==TRUE){
fwrite(data.table(O0O0OOOOO0OOOOO0O0OO),"H_metafounder_inbreeding_coefficients.txt",quote=F,sep="\t",row.names=F,col.names=T)
}else {
fwrite(data.table(O0O0OOOOO0OOOOO0O0OO),"H_inbreeding_coefficients.txt",quote=F,sep="\t",row.names=F,col.names=T)}
}
}else{
O0O0OOOOO0OOOOO0O0OO=NULL
}

if(!exists("G_A")){G_A=NULL};if(!exists("G_Ainv")){G_Ainv=NULL};
if(!exists("G_D")){G_D=NULL};if(!exists("G_Dinv")){G_Dinv=NULL};
if(!exists("P_A")){P_A=NULL};if(!exists("P_Ainv")){P_Ainv=NULL};
if(!exists("P_D")){P_D=NULL};if(!exists("P_Dinv")){P_Dinv=NULL};
if(!exists("H_A")){H_A=NULL};if(!exists("H_Ainv")){H_Ainv=NULL};
if(!exists("H_D")){H_D=NULL};if(!exists("H_Dinv")){H_Dinv=NULL};

if(!is.null(G_A)){rownames(G_A)=colnames(G_A)=as.character(IND_geno)}
if(!is.null(G_Ainv)){rownames(G_Ainv)=colnames(G_Ainv)=as.character(IND_geno)}
if(!is.null(G_D)){rownames(G_D)=colnames(G_D)=as.character(IND_geno)}
if(!is.null(G_Dinv)){rownames(G_Dinv)=colnames(G_Dinv)=as.character(IND_geno)}
if(!is.null(P_A)){rownames(P_A)=colnames(P_A)=as.character(IND_pedigree)}
if(!is.null(P_Ainv)){rownames(P_Ainv)=colnames(P_Ainv)=as.character(IND_pedigree)}
if(!is.null(P_D)){rownames(P_D)=colnames(P_D)=as.character(IND_pedigree)}
if(!is.null(P_Dinv)){rownames(P_Dinv)=colnames(P_Dinv)=as.character(IND_pedigree)}
if(!is.null(H_A)){rownames(H_A)=colnames(H_A)=as.character(IND_Additive)}
if(!is.null(H_Ainv)){rownames(H_Ainv)=colnames(H_Ainv)=as.character(IND_Additive)}
if(!is.null(H_D)){rownames(H_D)=colnames(H_D)=as.character(IND_Additive)}
if(!is.null(H_Dinv)){rownames(H_Dinv)=colnames(H_Dinv)=as.character(IND_Additive)}


if(return_result==TRUE){
return(list(G_A=list(A=G_A,Ainv=G_Ainv),G_D=list(D=G_D,Dinv=G_Dinv),P_A=list(A=P_A,Ainv=P_Ainv),P_D=list(D=P_D,Dinv=P_Dinv),
		  H_A=list(A=H_A,Ainv=H_Ainv),H_D=list(D=H_D,Dinv=H_Dinv),
		    Inbred=list(Homozygous=homo_inbred,G_diag=diag_inbred,Pedigree=pedigree_inbred,H_diag=O0O0OOOOO0OOOOO0O0OO),
			phased_block=data_block))
}

}

}



#solve要比 ginv快，如果solve报错，再调用 MASS::ginv求广义逆
OOO0OOO0O0O0OOOOOOO0<-function(Matrix){

	Matrix_inv<- try(chol2inv(chol(Matrix)), silent = TRUE)
	if(class(Matrix_inv) == "try-error"){
	cat("Added small value to diagnoal to ensure matrix invertible \n")
	Matrix_inv <- try(chol2inv(chol(Matrix+diag(0.0001,nrow(Matrix)))),silent=TRUE)
	if(class(Matrix_inv) == "try-error"){
	cat("Using GINV To Calculate Inverse Matrix \n")
	Matrix_inv <- OOO0OOO0O0O0OOOOOOO0_old(Matrix)
	}
	}
	return(Matrix_inv)
}



O0OOO0O0OOO0O0O0OOOO<-function(input_data_type=NULL,input_data_hmp=NULL,input_data_plink_ped=NULL,input_data_plink_map=NULL,
							  input_data_numeric=NULL,input_data_blupf90=NULL,input_data_haplotype_hap=NULL,
							  input_data_haplotype_map=NULL,input_data_haplotype_sample=NULL,input_data_vcf=NULL,miss_base=NULL){

if(is.null(input_data_type)){

if(!is.null(input_data_hmp)){input_data_type=c(input_data_type,"Hapmap")}
if(!is.null(input_data_plink_map)&!is.null(input_data_plink_ped)){input_data_type=c(input_data_type,"Plink")}
if(!is.null(input_data_blupf90)){input_data_type=c(input_data_type,"BLUPF90")}
if(!is.null(input_data_numeric)){input_data_type=c(input_data_type,"Numeric")}
if(!is.null(input_data_vcf)){input_data_type=c(input_data_type,"VCF")}
if(!is.null(input_data_haplotype_hap)&!is.null(input_data_haplotype_sample)&!is.null(input_data_haplotype_map)){input_data_type=c(input_data_type,"Haplotype")}
input_data_type=unique(input_data_type)

if(length(input_data_type)>=2){
stop("Only one kind of input data is allowed in this analysis!Please check your input data carefully! \n")
}else if(length(input_data_type)==0){
stop("No input data avaliable in this analysis \n")
}
}

if(is.null(miss_base)){
miss_base=ifelse(input_data_type=="Plink","0",ifelse(input_data_type=="Hapmap","N","N"))
}
return(list(input_data_type=input_data_type,miss_base=miss_base))
}
LR_evaluation<-function(partial_EBV,  #验证群个体的EBV， 根据partial数据集计算
						whole_EBV,     #完整数据得到的验证群EBV， 根据完整数据计算
						candidate_id=NULL,
						type=c("bias","disperation","accuracy","ratio_of_accuracy","ratio_of_reliability"),
						partial_sigma_a2=NULL, # 利用 partial 数据估计得到的方差组分
						K=NULL,                   # relationship matrix 
						F=NULL,
						f2=NULL
						){
cat("The format of user-provided EBV should contain two columns: id  and EBV! \n")
if(is.matrix(partial_EBV)){partial_EBV=data.frame(partial_EBV,stringsAsFactors=F)}
colnames(partial_EBV)=c("Id","partial_Effect")

result=partial_EBV
result$whole_Effect=whole_EBV[match(result[,1],whole_EBV[,1]),2]

#candidate
result=result[result[,1]%in%candidate_id,]
if(nrow(result)==0){stop("Couldn't find candidate individual in provided EBV data!")}
bias=disperation=accuracy=ratio_of_accuracy=ratio_of_reliability=NULL
# statistics
if("bias"%in%type){bias=mean(result$partial_Effect)-mean(result$whole_Effect)}

if("disperation"%in%type){disperation=as.numeric(lm(result$whole_Effect~result$partial_Effect)[[1]][2])} # cov(w,p)/var(p)

if("accuracy"%in%type){
if(is.null(K)&is.null(F)&is.null(f2)){stop("Please provide Kinship matrix or inbreeding coefficients !")}
if(!is.null(K)){
cat("Please make sure the colnames of Kinship matrix are individual name")
if(is.null(partial_sigma_a2)){stop("Please provide partial_sigma_a2!")}
IND_name=rownames(K)
K=K[match(candidate_id,IND_name),match(candidate_id,IND_name)];gc();
F=mean(diag(K)-1)
f2=mean(K[upper.tri(K)])
rm(K);gc();
}
accuracy=sqrt(cov(result$partial_Effect,result$whole_Effect)/((1+F-f2)*partial_sigma_a2))
}

if("ratio_of_accuracy"%in%type){ratio_of_accuracy=cor(result$partial_Effect,result$whole_Effect)}

if("ratio_of_reliability"%in%type){ratio_of_reliability=cov(result$whole_Effect,result$partial_Effect)/var(result$whole_Effect)}

return(list(bias=bias,disperation=disperation,accuracy=accuracy,
			ratio_of_accuracy=ratio_of_accuracy,ratio_of_reliability=ratio_of_reliability))
}



#相关系数检验
Hotelling_test <- function(gen_value,phe,two_tailed=TRUE){
	
   cat("Please make sure you've been installed package:psych! \n")	
   cat("The format of gen_value is: id,ebv1,ebv2...... \n")   
   cat("The format of phe is: id,phe \n")

   gen_value=data.frame(gen_value,stringsAsFactors=F)
   phe=data.frame(phe,stringsAsFactors=F)
   if(ncol(gen_value)<=2){stop("The column of gen_value should >=3 !")}
   if(!ncol(phe)==2){stop("The column of phe should ==2 !")}
   
   
   intersect_id=intersect(gen_value[,1],phe[,1])
   if(length(intersect_id)==0){stop("There are no identical id between gen_value and phe!")}
   gen_value=gen_value[match(intersect_id,gen_value[,1]),]
   phe=phe[match(intersect_id,phe[,1]),]
   gen_value=gen_value[,-1]
   phe=phe[,-1]
   
cat("Off-digonal is p-value \n")   
   gen_name <- colnames(gen_value)
   results <- matrix(ncol=ncol(gen_value),nrow=ncol(gen_value))
   rownames(results)=colnames(results)=gen_name
   for (OOO0OOO0O0OOO0O0OOO0OOOO in 1:ncol(gen_value)){
     for (OOO0OOO0O0OOO0O0OOO0OOOOO in 1:ncol(gen_value)){
       y <- gen_value[,OOO0OOO0O0OOO0O0OOO0OOOO]
       z <- gen_value[,OOO0OOO0O0OOO0O0OOO0OOOOO]
       if(OOO0OOO0O0OOO0O0OOO0OOOO==OOO0OOO0O0OOO0O0OOO0OOOOO){next}
       xy <- cor(phe,y)
       xz <- cor(phe,z)
       yz <- cor(y,z)
       n <- length(y)
       results [OOO0OOO0O0OOO0O0OOO0OOOO,OOO0OOO0O0OOO0O0OOO0OOOOO] <- psych::paired.r(xy, xz, yz, n,twotailed=two_tailed)$p
     }
   }
   return(results)
 }	#' Print 'Hello world!'
#' @export
#' @param phe_col_names colnames of phenotype
#' @param phe_name phenotype file name
#' @param target_trait_name name of target trait
#' @param fixed_effect_name names of trait's fixed effect
#' @param random_effect_name names of trait's random effect
#' @param covariate_effect_name names of trait's covariate effect
#' @param provided_effect_file_path effect file path
#' @param provided_effect_file_name effect file name
#' @param missing_value missing value in phenotype
#' @param relationship_name relationship file name
#' @param relationship_path relationship file path
#' @param analysis_model BLUPF90 analysis model
#' @param genetic_effect_name name of genetic effect
#' @param included_permanent_effect whehter include permanet effect in analysis
#' @param provided_BLUPF90_prior_file_path prior file path
#' @param provided_BLUPF90_prior_file_name prior file name
#' @param provided_BLUPF90_prior_effect_name prior effect name
#' @param BLUPF90_genumeric_name name
#' @param BLUPF90_map_name map name
generate_renum<-function(
								   phe_col_names=NULL,
								   target_trait_name=NULL,
								   fixed_effect_name=NULL, #列表
								   random_effect_name=NULL, #列表，不包括永久环境效应
								   covariate_effect_name=NULL, #列表
								   provided_effect_file_path=NULL, #各个性状的效应文件
								   provided_effect_file_name="model_define.txt",
								   phe_name=NULL,
								   phe_path=NULL,

								   missing_value="-9999",
								   relationship_name=NULL,
								   relationship_path=NULL,
								   analysis_model="PBLUP_A",
								   genetic_effect_name="Id",
								   user_file_id=NULL,
								   included_permanent_effect=FALSE,
								   provided_BLUPF90_prior_file_path=NULL,
								   provided_BLUPF90_prior_file_name=NULL,
								   provided_BLUPF90_prior_effect_name=NULL, #随机效应的名称, 包括Residual
								   BLUPF90_alt_option=NULL,
								   BLUPF90_genumeric_name=NULL,
								   BLUPF90_map_name=NULL
								   ){

O0O0OOO0O0O0OOOOO0OO=length(target_trait_name)

if(is.null(relationship_name)){
stop("Please provide relationship name!")
}else{
addtive_relationship_name=relationship_name[1]
if(analysis_model%in%c("GBLUP_AD")){dominance_relationship_name=relationship_name[2]}
if(analysis_model%in%c("SSBLUP_A")){SSBLUP_G_matrix_name=relationship_name[2]}
}



if(!is.null(fixed_effect_name)){if(length(fixed_effect_name)!=O0O0OOO0O0O0OOOOO0OO){stop("The number of fixed effect names doesn't match the number of traits")}}
if(!is.null(random_effect_name)){if(length(random_effect_name)!=O0O0OOO0O0O0OOOOO0OO){stop("The number of random effect names doesn't match the number of traits")}}
if(!is.null(covariate_effect_name)){if(length(covariate_effect_name)!=O0O0OOO0O0O0OOOOO0OO){stop("The number of covariate effect names doesn't match the number of traits")}}



#PRIOR
#BLUPF90 PRIOR的结构和DMU不同，和 Additive (CO)VARIANCES 结构一致
if(!is.null(provided_BLUPF90_prior_file_name)&!is.null(provided_BLUPF90_prior_file_path)&!is.null(provided_BLUPF90_prior_effect_name)){
OOOOOOO0O0O0OOOOOOO0=data.table::fread(paste0(provided_BLUPF90_prior_file_path,"/",provided_BLUPF90_prior_file_name),data.table=F,header=F)
OOOOOOO0O0O0OOOOOOO0=as.matrix(OOOOOOO0O0O0OOOOOOO0)
}else{

prior=OOO0OOOOO0OOOOO0OOO0(target_trait_name=target_trait_name,
						      random_effect_name=random_effect_name,
							  included_permanent_effect=included_permanent_effect)

OOOOOOO0O0O0OOOOOOO0=as.matrix(prior)

provided_BLUPF90_prior_effect_name=c(do.call(c,random_effect_name),"Residual")

}

#获取固定效应、随机效应、协变量记录
if(!is.null(provided_effect_file_path)&!is.null(provided_effect_file_name)){

model=data.table::fread(paste0(provided_effect_file_path,"/",provided_effect_file_name),data.table=F,fill=TRUE,header = F)

fixed_effect_name=NULL;random_effect_name=NULL;covariate_effect_name=NULL
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(target_trait_name)){

n=match(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],model[,1]) # which line  is  the trait

if(is.na(n)){stop(paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO]," couldn't find in provided effect_file"))}

OOO0OOOOOOO0O0OOOOOO=grep(pattern="*",as.character(model[n,]),fixed=T)
fixed_effect=as.character(model[n,c((OOO0OOOOOOO0O0OOOOOO[1]+1):(OOO0OOOOOOO0O0OOOOOO[2]-1))])
random_effect=as.character(model[n,c((OOO0OOOOOOO0O0OOOOOO[2]+1):(OOO0OOOOOOO0O0OOOOOO[3]-1))])
covariate_effect=as.character(model[n,c((OOO0OOOOOOO0O0OOOOOO[3]+1):(OOO0OOOOOOO0O0OOOOOO[4]-1))])

if("*"%in%fixed_effect){fixed_effect=NULL}
if("*"%in%random_effect){random_effect=NULL}
if("*"%in%covariate_effect){covariate_effect=NULL}

fixed_effect_name=c(fixed_effect_name,list(fixed_effect))
random_effect_name=c(random_effect_name,list(random_effect))
covariate_effect_name=c(covariate_effect_name,list(covariate_effect))
}
}


####################       #################################
########construct  Header file Part                #########
O0O0O0OOO0O0O0OOO0OO=match(target_trait_name,phe_col_names)


OOOOOOO0OOOOO0O0OOOO=rbind("DATAFILE",phe_name,
							  "TRAITS",paste0(O0O0O0OOO0O0O0OOO0OO,collapse=" "),
							  "FIELDS_PASSED TO OUTPUT","","WEIGHT(S)","")
							  

####################       #################################
########construct  Residual Part                   #########
O0OOOOO0O0O0OOOOOOOO="Residual"
OOOOO0O0O0OOO0O0OOOO=match(O0OOOOO0O0O0OOOOOOOO,provided_BLUPF90_prior_effect_name)
O0O0OOO0OOO0O0OOO0OO_pos=(O0O0OOO0O0O0OOOOO0OO*(OOOOO0O0O0OOO0O0OOOO-1)+1):(O0O0OOO0O0O0OOOOO0OO*OOOOO0O0O0OOO0O0OOOO)
O0O0OOO0OOO0O0OOO0OO=OOOOOOO0O0O0OOOOOOO0[O0O0OOO0OOO0O0OOO0OO_pos,O0O0OOO0OOO0O0OOO0OO_pos]
O0O0OOO0OOO0O0OOO0OO=OOO0OOO0OOOOO0OOOOOO(O0O0OOO0OOO0O0OOO0OO)

O0OOOOO0O0OOOOO0O0O0=rbind("RESIDUAL_VARIANCE",O0O0OOO0OOO0O0OOO0OO) # 残差

#获取各性状下效应对应的位置
####################       #################################
########construct  Fixed-Effect Part               #########

if(!is.null(fixed_effect_name)&(!identical(fixed_effect_name,rep(list(NULL),O0O0OOO0O0O0OOOOO0OO)))){
O0O0OOOOOOOOO0O0OOO0=unique(do.call(c,fixed_effect_name))
O0O0OOOOOOOOO0O0OOO0=sort(O0O0OOOOOOOOO0O0OOO0)
total_O0O0OOO0OOO0OOO0OOO0=match(O0O0OOOOOOOOO0O0OOO0,phe_col_names)

O0OOOOO0OOO0OOO0O0O0=NULL
O0O0OOO0OOO0OOO0OOO0=NULL

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:O0O0OOO0O0O0OOOOO0OO){  #每个性状下固定效应的 pos
	OOO0OOO0OOOOOOOOO0OO=total_O0O0OOO0OOO0OOO0OOO0
	OOO0OOO0OOOOOOOOO0OO[!O0O0OOOOOOOOO0O0OOO0%in%fixed_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]]=0
	O0O0OOO0OOO0OOO0OOO0=c(O0O0OOO0OOO0OOO0OOO0,list(OOO0OOO0OOOOOOOOO0OO))
}
O0O0OOO0OOO0OOO0OOO0_matrix=do.call(rbind,O0O0OOO0OOO0OOO0OOO0) #pos的矩阵格式

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(O0O0OOOOOOOOO0O0OOO0)){
	O0OOOOO0OOO0OOO0O0O0=rbind(O0OOOOO0OOO0OOO0O0O0,"EFFECT",paste0(paste(O0O0OOO0OOO0OOO0OOO0_matrix[,OOO0OOO0O0OOO0O0OOO0OOOO],collapse=" ")," cross alpha "))
}
}else{
O0OOOOO0OOO0OOO0O0O0=data.frame(NULL)}

####################       #################################
########construct  Covariate-Effect Part           #########
if(!is.null(covariate_effect_name)&(!identical(covariate_effect_name,rep(list(NULL),O0O0OOO0O0O0OOOOO0OO)))){
OOOOOOO0O0O0OOOOOOOO=unique(do.call(c,covariate_effect_name))
OOOOOOO0O0O0OOOOOOOO=sort(OOOOOOO0O0O0OOOOOOOO)
O0OOOOO0O0OOO0O0O0OO=match(OOOOOOO0O0O0OOOOOOOO,phe_col_names)
OOO0OOO0O0O0O0OOOOO0=NULL
covariate_effect_pos=NULL

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:O0O0OOO0O0O0OOOOO0OO){  #每个性状下固定效应的 pos
	O0OOO0OOOOO0O0OOOOO0=O0OOOOO0O0OOO0O0O0OO
	O0OOO0OOOOO0O0OOOOO0[!OOOOOOO0O0O0OOOOOOOO%in%covariate_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]]=0
	covariate_effect_pos=c(covariate_effect_pos,list(O0OOO0OOOOO0O0OOOOO0))
}
covariate_effect_pos_matrix=do.call(rbind,covariate_effect_pos) #pos的矩阵格式


for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(OOOOOOO0O0O0OOOOOOOO)){
	OOO0OOO0O0O0O0OOOOO0=rbind(OOO0OOO0O0O0O0OOOOO0,"EFFECT",paste0(paste(covariate_effect_pos_matrix[,OOO0OOO0O0OOO0O0OOO0OOOO],collapse=" ")," cov "))
}
}else{
OOO0OOO0O0O0O0OOOOO0=data.frame(NULL)}


####################       #################################
########construct  Non_genetic-random effect Part  #########

#for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(random_effect_name)){   #remove genetic effect name
#temp_random=random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]
#temp_random=temp_random[!temp_random%in%genetic_effect_name]
#if(length(temp_random)==0){temp_random=NULL}
#random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]=temp_random
#}

O0OOOOOOOOO0O0O0O0O0=unique(do.call(c,random_effect_name))
O0OOOOOOOOO0O0O0O0O0=sort(O0OOOOOOOOO0O0O0O0O0)
OOOOO0OOO0OOOOOOO0OO=match(O0OOOOOOOOO0O0O0O0O0,phe_col_names)

#remove genetic effect
OOOOO0OOO0OOOOOOO0OO=OOOOO0OOO0OOOOOOO0OO[!O0OOOOOOOOO0O0O0O0O0%in%genetic_effect_name]
O0OOOOOOOOO0O0O0O0O0=O0OOOOOOOOO0O0O0O0O0[!O0OOOOOOOOO0O0O0O0O0%in%genetic_effect_name]

if(length(O0OOOOOOOOO0O0O0O0O0)>0&(!identical(random_effect_name,rep(list(NULL),O0O0OOO0O0O0OOOOO0OO)))){
OOOOOOO0OOO0O0O0OOOO=NULL
random_effect_pos=NULL

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:O0O0OOO0O0O0OOOOO0OO){  #每个性状下固定效应的 pos
	temp_random_pos=OOOOO0OOO0OOOOOOO0OO
	temp_random_pos[!O0OOOOOOOOO0O0O0O0O0%in%random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]]=0
	random_effect_pos=c(random_effect_pos,list(temp_random_pos))
}

random_effect_pos_matrix=do.call(rbind,random_effect_pos) #pos的矩阵格式

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(O0OOOOOOOOO0O0O0O0O0)){

	O0OOOOO0O0O0OOOOOOOO=O0OOOOOOOOO0O0O0O0O0[OOO0OOO0O0OOO0O0OOO0OOOO]
	OOOOO0O0O0OOO0O0OOOO=match(O0OOOOO0O0O0OOOOOOOO,provided_BLUPF90_prior_effect_name)
	O0O0OOO0OOO0O0OOO0OO_pos=(O0O0OOO0O0O0OOOOO0OO*(OOOOO0O0O0OOO0O0OOOO-1)+1):(O0O0OOO0O0O0OOOOO0OO*OOOOO0O0O0OOO0O0OOOO)
	O0O0OOO0OOO0O0OOO0OO=OOOOOOO0O0O0OOOOOOO0[O0O0OOO0OOO0O0OOO0OO_pos,O0O0OOO0OOO0O0OOO0OO_pos]
	O0O0OOO0OOO0O0OOO0OO=OOO0OOO0OOOOO0OOOOOO(O0O0OOO0OOO0O0OOO0OO)

	OOOOOOO0OOO0O0O0OOOO=rbind(OOOOOOO0OOO0O0O0OOOO,"EFFECT",paste0(paste(random_effect_pos_matrix[,OOO0OOO0O0OOO0O0OOO0OOOO],collapse=" ")," cross alpha "),
	                              "RANDOM","diagonal","(CO)VARIANCES",O0O0OOO0OOO0O0OOO0OO)
}
}else{
OOOOOOO0OOO0O0O0OOOO=data.frame(NULL)}



OOOOO0OOOOO0O0OOO0OO=NULL
random_effect_pos=NULL
O0O0O0O0OOOOOOO0OOOO=NULL
O0O0OOO0OOOOO0OOOOOO=NULL

####################       #################################
########construct  Permanent Effect Part   		   #########

if(included_permanent_effect==TRUE){
O0O0O0O0OOOOOOO0OOOO=rbind("OPTIONAL","pe")
	O0OOOOO0O0O0OOOOOOOO="Permanent"
	OOOOO0O0O0OOO0O0OOOO=match(O0OOOOO0O0O0OOOOOOOO,provided_BLUPF90_prior_effect_name)
	O0O0OOO0OOO0O0OOO0OO_pos=(O0O0OOO0O0O0OOOOO0OO*(OOOOO0O0O0OOO0O0OOOO-1)+1):(O0O0OOO0O0O0OOOOO0OO*OOOOO0O0O0OOO0O0OOOO)
	O0O0OOO0OOO0O0OOO0OO=OOOOOOO0O0O0OOOOOOO0[O0O0OOO0OOO0O0OOO0OO_pos,O0O0OOO0OOO0O0OOO0OO_pos]
	O0O0OOO0OOO0O0OOO0OO=OOO0OOO0OOOOO0OOOOOO(O0O0OOO0OOO0O0OOO0OO)
	O0O0OOO0OOOOO0OOOOOO=rbind("(CO)VARIANCES_PE",O0O0OOO0OOO0O0OOO0OO)
}

####################       #################################
########construct  Genetic-random effect Part      #########

if(analysis_model=="GBLUP_A"){
	O0OOOOO0O0O0OOOOOOOO=genetic_effect_name
	OOOOO0O0O0OOO0O0OOOO=match(O0OOOOO0O0O0OOOOOOOO,provided_BLUPF90_prior_effect_name)
	O0O0OOO0OOO0O0OOO0OO_pos=(O0O0OOO0O0O0OOOOO0OO*(OOOOO0O0O0OOO0O0OOOO-1)+1):(O0O0OOO0O0O0OOOOO0OO*OOOOO0O0O0OOO0O0OOOO)
	O0O0OOO0OOO0O0OOO0OO=OOOOOOO0O0O0OOOOOOO0[O0O0OOO0OOO0O0OOO0OO_pos,O0O0OOO0OOO0O0OOO0OO_pos]
	O0O0OOO0OOO0O0OOO0OO=OOO0OOO0OOOOO0OOOOOO(O0O0OOO0OOO0O0OOO0OO)

	OOOOO0OOOOO0O0OOO0OO=rbind(OOOOO0OOOOO0O0OOO0OO,"EFFECT",paste0(paste(match(rep(genetic_effect_name,O0O0OOO0O0O0OOOOO0OO),phe_col_names),collapse=" ")," cross alpha "),
	                              "RANDOM","animal","SNP_FILE",addtive_relationship_name,
						     "(CO)VARIANCES",O0O0OOO0OOO0O0OOO0OO,O0O0OOO0OOOOO0OOOOOO)
}else if(analysis_model=="PBLUP_A"){

	O0OOOOO0O0O0OOOOOOOO=genetic_effect_name
	OOOOO0O0O0OOO0O0OOOO=match(O0OOOOO0O0O0OOOOOOOO,provided_BLUPF90_prior_effect_name)
	O0O0OOO0OOO0O0OOO0OO_pos=(O0O0OOO0O0O0OOOOO0OO*(OOOOO0O0O0OOO0O0OOOO-1)+1):(O0O0OOO0O0O0OOOOO0OO*OOOOO0O0O0OOO0O0OOOO)
	O0O0OOO0OOO0O0OOO0OO=OOOOOOO0O0O0OOOOOOO0[O0O0OOO0OOO0O0OOO0OO_pos,O0O0OOO0OOO0O0OOO0OO_pos]
	O0O0OOO0OOO0O0OOO0OO=OOO0OOO0OOOOO0OOOOOO(O0O0OOO0OOO0O0OOO0OO)

	OOOOO0OOOOO0O0OOO0OO=rbind(OOOOO0OOOOO0O0OOO0OO,"EFFECT",paste0(paste(match(rep(genetic_effect_name,O0O0OOO0O0O0OOOOO0OO),phe_col_names),collapse=" ")," cross alpha "),
	                              "RANDOM","animal",O0O0O0O0OOOOOOO0OOOO,"FILE",addtive_relationship_name,"FILE_POS","1 2 3 0 0",
							 "PED_DEPTH","0",
						     "(CO)VARIANCES",O0O0OOO0OOO0O0OOO0OO,O0O0OOO0OOOOO0OOOOOO)

}else if(analysis_model=="User_define"){
    addtive_relationship_name="dummy_pedigree.txt"
	if(is.null(user_file_id)){stop("Please provide user_file_id!")} 
	#generate dummy_pedigree for generating renf90.par file 
	phe=read.table(paste0(phe_path,"/",phe_name),header=F) 
	colnames(phe)=phe_col_names
	Offspring=phe[,genetic_effect_name]
	write.table(data.frame(Offspring=user_file_id,Sire=0,Dam=0),"dummy_pedigree.txt",sep=" ",quote=F,row.names=F,col.names=F)
	#write.table(user_define_ped,"dummy_pedigree.txt",sep=" ",quote=F,row.names=F,col.names=F) 
	effect_pos=(1:length(provided_BLUPF90_prior_effect_name))
	for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(relationship_name)){ 
	O0OOOOO0O0O0OOOOOOOO=genetic_effect_name
	#OOOOO0O0O0OOO0O0OOOO=match(O0OOOOO0O0O0OOOOOOOO,provided_BLUPF90_prior_effect_name)
	OOOOO0O0O0OOO0O0OOOO=effect_pos[provided_BLUPF90_prior_effect_name%in%O0OOOOO0O0O0OOOOOOOO][OOO0OOO0O0OOO0O0OOO0OOOO]
	O0O0OOO0OOO0O0OOO0OO_pos=(O0O0OOO0O0O0OOOOO0OO*(OOOOO0O0O0OOO0O0OOOO-1)+1):(O0O0OOO0O0O0OOOOO0OO*OOOOO0O0O0OOO0O0OOOO)
	O0O0OOO0OOO0O0OOO0OO=OOOOOOO0O0O0OOOOOOO0[O0O0OOO0OOO0O0OOO0OO_pos,O0O0OOO0OOO0O0OOO0OO_pos]
	O0O0OOO0OOO0O0OOO0OO=OOO0OOO0OOOOO0OOOOOO(O0O0OOO0OOO0O0OOO0OO)

	OOOOO0OOOOO0O0OOO0OO=rbind(OOOOO0OOOOO0O0OOO0OO,"EFFECT",paste0(paste(match(rep(genetic_effect_name,O0O0OOO0O0O0OOOOO0OO),phe_col_names),collapse=" ")," cross alpha "),
	                              "RANDOM","animal",O0O0O0O0OOOOOOO0OOOO,"FILE",addtive_relationship_name,"FILE_POS","1 2 3 0 0",
							 "PED_DEPTH","0",
						     "(CO)VARIANCES",O0O0OOO0OOO0O0OOO0OO,O0O0OOO0OOOOO0OOOOOO)
	}
	
}else if(analysis_model=="SSBLUP_A"){

	O0OOOOO0O0O0OOOOOOOO=genetic_effect_name
	OOOOO0O0O0OOO0O0OOOO=match(O0OOOOO0O0O0OOOOOOOO,provided_BLUPF90_prior_effect_name)
	O0O0OOO0OOO0O0OOO0OO_pos=(O0O0OOO0O0O0OOOOO0OO*(OOOOO0O0O0OOO0O0OOOO-1)+1):(O0O0OOO0O0O0OOOOO0OO*OOOOO0O0O0OOO0O0OOOO)
	O0O0OOO0OOO0O0OOO0OO=OOOOOOO0O0O0OOOOOOO0[O0O0OOO0OOO0O0OOO0OO_pos,O0O0OOO0OOO0O0OOO0OO_pos]
	O0O0OOO0OOO0O0OOO0OO=OOO0OOO0OOOOO0OOOOOO(O0O0OOO0OOO0O0OOO0OO)

	OOOOO0OOOOO0O0OOO0OO=rbind(OOOOO0OOOOO0O0OOO0OO,"EFFECT",paste0(paste(match(rep(genetic_effect_name,O0O0OOO0O0O0OOOOO0OO),phe_col_names),collapse=" ")," cross alpha "),
	                              "RANDOM","animal",O0O0O0O0OOOOOOO0OOOO,"FILE",addtive_relationship_name,"FILE_POS","1 2 3 0 0",
							"SNP_FILE",relationship_name[2],
							 "PED_DEPTH","0",
						     "(CO)VARIANCES",O0O0OOO0OOO0O0OOO0OO,O0O0OOO0OOOOO0OOOOOO)


}




################## ##################################
           ###################### ####################
file_genotype_part=data.frame(NULL)
if(!is.null(BLUPF90_genumeric_name)){
SNP_genumeric_file_part=rbind("SNP_FILE",BLUPF90_genumeric_name)
}


####################       #################################
########construct  Option Part                     #########
SNP_map_file_part_option=data.frame(NULL)
if(!is.null(BLUPF90_map_name)){
SNP_map_file_part_option=paste0("OPTION ","chrinfo ",BLUPF90_map_name)
}

OOO0OOO0OOOOO0OOO0OO=data.frame(NULL)
if(!is.null(BLUPF90_genumeric_name)){
OOO0OOO0OOOOO0OOO0OO=paste0("OPTION ","SNP_file ",BLUPF90_genumeric_name)
}

OOO0OOO0O0O0O0OOO0OO=rbind(paste0("OPTION missing ",missing_value),
						 "OPTION alpha_size 30",
						 "OPTION use_yams",
						 "OPTION sol se",
						 #"OPTION se_covar_function H2d G_2_2_1_1/(G_2_2_1_1+G_2_3_1_1+G_3_3_1_1+G_4_4_1_1+R_1_1)",
						 "OPTION residual")

if(!is.null(BLUPF90_alt_option)){
for(option in BLUPF90_alt_option){
OOO0OOO0O0O0O0OOO0OO=rbind(OOO0OOO0O0O0O0OOO0OO,option)
}
}

renum_par_file=plyr::rbind.fill.matrix(OOOOOOO0OOOOO0O0OOOO,
								 O0OOOOO0O0OOOOO0O0O0,
								 O0OOOOO0OOO0OOO0O0O0,
								 OOO0OOO0O0O0O0OOOOO0,
								 OOOOOOO0OOO0O0O0OOOO,
								 OOOOO0OOOOO0O0OOO0OO,
								 SNP_map_file_part_option,
								 OOO0OOO0OOOOO0OOO0OO,
								 OOO0OOO0O0O0O0OOO0OO)
utils::write.table(renum_par_file,"renum.par",quote=F,row.names=F,col.names=F)
}



#将n列矩阵变成一列矩阵
OOO0OOO0OOOOO0OOOOOO<-function(data){
data=as.matrix(data)
OOOOOOO0OOO0O0O0OOO0=matrix(NA,nrow=nrow(data),ncol=1)
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:nrow(OOOOOOO0OOO0O0O0OOO0)){
OOOOOOO0OOO0O0O0OOO0[OOO0OOO0O0OOO0O0OOO0OOOO,]=paste(data[OOO0OOO0O0OOO0O0OOO0OOOO,],collapse=" ")
}
return(OOOOOOO0OOO0O0O0OOO0)
}


OOO0OOOOO0OOOOO0OOO0<-function(target_trait_name=NULL,
								 random_effect_name=NULL,
								 included_permanent_effect=FALSE){


prior=diag(10+length(target_trait_name)+length(do.call(c,random_effect_name))+included_permanent_effect)


return(prior)
}

#*************generate DIR. file **********
generate_DIR<-function(phe_col_names=NULL,
					   target_trait_name=NULL,
					   fixed_effect_name=NULL, #列表
					   random_effect_name=NULL, #列表，不包括永久环境效应
					   covariate_effect_name=NULL, #列表
					   random_regression_effect_name=NULL, #列表，随机回归效应   list(c("",""),c("",""))
					   maternal_effect_name=NULL, #列表
					   
					   include_social_effect=FALSE,    #是否评估social effect
					   group_effect_name=NULL,  #评估social effect时, group的名称,自动生成group_phe,												
					   integer_group_names=NULL,        #整型group的名称
					   real_group_names=NULL,            #实型group的名称
		   
					   provided_effect_file_path=NULL, #各个性状的效应文件
					   provided_effect_file_name="model_define.txt",
					   phe_path=NULL,
					   phe_name=NULL,
					   analysis_model=NULL,
					   ped_inbred_number=1, #构建系谱的时候，用近交的方法
					   genetic_effect_name="Id",
					   included_permanent_effect=FALSE, 
					   included_dominance_effect=FALSE,
					   missing_value= -9999,
					   iteration_criteria= 1.0e-7,
					   relationship_path=NULL,
					   relationship_name=NULL,
					   dmu_module=NULL,
					   dmu_algorithm_code=NULL,
					   provided_prior_file_path=NULL,
					   provided_prior_file_name=NULL,
					   integer_n=NULL,  #整型数目
					   residual_cov_trait=NULL,  #限定残差协方差的性质
					   output_DIR_path=NULL,
					   output_DIR_name=NULL,
					   IND_geno_file_name=NULL, #
					   IND_geno=NULL,
					   SSBLUP_omega=0.05				   
					   ){
#生成 social_effect的表型
#social effect其实可以看做是特殊情况的随机回归，需要将随机效应那一行设置为 1+1+1+1类似的形式即可

if(TRUE%in%include_social_effect){

cat("Performing social genetic effect analysis...... \n")
cat("Group effect should be included as random effect in the social genetic effect analysis \n")
cat("User must specify the integer_group_names and real_group_names or group_effect_name! \n")
if(!is.null(group_effect_name)){

cat("Generating phenotype of social genetic effect automatically......\n")

phe=fread(paste0(phe_path,"/",phe_name),data.table=F,header=F)
colnames(phe)=phe_col_names

social_result=generate_social_phe(phe=phe,
								  integer_n=integer_n,
					                social_effect_group_name=group_effect_name,
					                genetic_effect_name=genetic_effect_name)

phe=social_result$phe
integer_n=social_result$integer_n
integer_group_names=social_result$integer_group_names
real_group_names=social_result$real_group_names

phe_col_names=colnames(phe)
phe_name=paste0("Social_phe_",phe_name)
phe_path=output_DIR_path
fwrite(phe,paste0(phe_path,"/",phe_name),quote=F,sep=" ",row.names=F,col.names=F)

random_regression_effect_name=rep(list(paste0(real_group_names,"&",integer_group_names)),
                                          length(target_trait_name)) #单性状

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(target_trait_name)){
random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]=c(random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]],integer_group_names)  # 单性状
}

}else{

cat("Using user-provided social genetic effect phenotype......\n")

if(is.null(integer_group_names)|is.null(real_group_names)){stop("User must specify the integer_group_names or real_group_names when using provided phenotype!")}
integer_group_names=integer_group_names
real_group_names=real_group_names


random_regression_effect_name=rep(list(paste0(real_group_names,"&",integer_group_names)),
                                          length(target_trait_name)) #单性状

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(target_trait_name)){
random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]=c(random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]],integer_group_names)  # 单性状
}

}
}


if(!is.null(IND_geno)){
write.table(IND_geno,paste0(relationship_path[1],"/","IND_geno.txt"),quote=F,row.names=F,col.names=F,sep=" ")
IND_geno_file_name="IND_geno.txt"
}

					   
O0O0OOO0O0O0OOOOO0OO=length(target_trait_name)
O0OOOOOOOOO0OOO0OOO0=length(phe_col_names)
O0O0OOO0O0OOOOOOOOO0=O0OOOOOOOOO0OOO0OOO0-integer_n
OOOOO0O0O0O0O0OOO0O0=phe_col_names[1:integer_n]
O0OOO0O0OOO0O0OOO0OO=phe_col_names[(integer_n+1):O0OOOOOOOOO0OOO0OOO0]

library(data.table)
#读取各性状的固定效应、随机效应、协变量记录
if(!is.null(random_effect_name)&!is.null(provided_effect_file_path)){stop("Effect_file and  effect_names couldn't provided simultaneously! ")}
if(!is.null(provided_effect_file_path)&!is.null(provided_effect_file_name)){

model=fread(paste0(provided_effect_file_path,"/",provided_effect_file_name),data.table=F,fill=TRUE,header = F)

fixed_effect_name=NULL;random_effect_name=NULL;covariate_effect_name=NULL
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(target_trait_name)){

n=match(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],model[,1]) # which line  is  the trait

if(is.na(n)){stop(paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO]," couldn't find in provided effect_file; "))}

OOO0OOOOOOO0O0OOOOOO=grep(pattern="*",as.character(model[n,]),fixed=T)
fixed_effect=as.character(model[n,c((OOO0OOOOOOO0O0OOOOOO[1]+1):(OOO0OOOOOOO0O0OOOOOO[2]-1))])
random_effect=as.character(model[n,c((OOO0OOOOOOO0O0OOOOOO[2]+1):(OOO0OOOOOOO0O0OOOOOO[3]-1))])
covariate_effect=as.character(model[n,c((OOO0OOOOOOO0O0OOOOOO[3]+1):(OOO0OOOOOOO0O0OOOOOO[4]-1))])

if("*"%in%fixed_effect){fixed_effect=NULL}
if("*"%in%random_effect){random_effect=NULL}
if("*"%in%covariate_effect){covariate_effect=NULL}

fixed_effect_name=c(fixed_effect_name,list(fixed_effect))
random_effect_name=c(random_effect_name,list(random_effect))
covariate_effect_name=c(covariate_effect_name,list(covariate_effect))
}
}


#software create DIR file 
COMMENT=matrix("",nrow=2,ncol=length(phe_col_names))

ANALYSE=matrix("",nrow=2,ncol=length(phe_col_names))

DATA=matrix("",nrow=2,ncol=length(phe_col_names))

VARIABLE=matrix("",nrow=10,ncol=length(phe_col_names))

VAR_STR=matrix("",nrow=4+length(relationship_name),ncol=length(phe_col_names))

SOLUTION=matrix("",nrow=2,ncol=length(phe_col_names))

MODEL=matrix("",nrow=1+O0O0OOO0O0O0OOOOO0OO*(O0O0OOO0O0O0OOOOO0OO+6),ncol=length(phe_col_names))

PRIOR=matrix("",nrow=(O0O0OOO0O0O0OOOOO0OO*(O0O0OOO0O0O0OOOOO0OO+1)/2)*6,ncol=length(phe_col_names))

OOO0OOO0OOOOO0O0OOOO=matrix("",nrow=10,ncol=length(phe_col_names))




#$COMMENT
COMMENT[1,1]="$COMMENT"

#$ANALYSE 
if(is.null(dmu_module)){dmu_module="dmuai"}
OOOOO0OOO0OOOOO0O0OO=ifelse(dmu_module=="dmuai",1,
                   ifelse(dmu_module=="rjmc",2,
		          ifelse(dmu_module=="dmu4",11,
			      ifelse(dmu_module=="dmu5",12,NA))))
	
if(is.null(dmu_algorithm_code)){	
dmu_algorithm_code=ifelse(dmu_module=="dmuai",1,
                   ifelse(dmu_module=="rjmc",0,
		          ifelse(dmu_module=="dmu4",14,
			      ifelse(dmu_module=="dmu5",2,NA))))}		

		  
ANALYSE[1,1]=paste0("$ANALYSE ",OOOOO0OOO0OOOOO0O0OO," ",dmu_algorithm_code," 0 0 ")

#$DATA 
DATA[1,1]=paste0("$DATA  ASCII (",integer_n,",",O0O0OOO0O0OOOOOOOOO0," ,",missing_value,") ",phe_path,"/",phe_name)


#$VARIABLE
VARIABLE[1,1]="$VARIABLE"
VARIABLE[3,1:integer_n]=paste0("I",1:integer_n)
VARIABLE[4,1:integer_n]=OOOOO0O0O0O0O0OOO0O0
VARIABLE[6,1:O0O0OOO0O0OOOOOOOOO0]=paste0("R",1:O0O0OOO0O0OOOOOOOOO0)
VARIABLE[7,1:O0O0OOO0O0OOOOOOOOO0]=O0OOO0O0OOO0O0OOO0OO
VARIABLE[3,1]="#I1";VARIABLE[6,1]="#R1";


#将 I1.... R1.... 和变量名称对齐
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:ncol(VARIABLE)){
integer_space=nchar(VARIABLE[4,OOO0OOO0O0OOO0O0OOO0OOOO])
VARIABLE[3,OOO0OOO0O0OOO0O0OOO0OOOO]=sprintf(paste0("%-",integer_space,"s"),VARIABLE[3,OOO0OOO0O0OOO0O0OOO0OOOO])
real_space=nchar(VARIABLE[7,OOO0OOO0O0OOO0O0OOO0OOOO])
VARIABLE[6,OOO0OOO0O0OOO0O0OOO0OOOO]=sprintf(paste0("%-",real_space,"s"),VARIABLE[6,OOO0OOO0O0OOO0O0OOO0OOOO])
}


#$MODEL
MODEL[1,1]="$MODEL"
MODEL[2,1]=paste0(O0O0OOO0O0O0OOOOO0OO," ",O0O0OOO0O0O0OOOOO0OO," 0 0 0")
MODEL[3:(2+O0O0OOO0O0O0OOOOO0OO),1]="0"

if(!is.null(fixed_effect_name)){if(length(fixed_effect_name)!=O0O0OOO0O0O0OOOOO0OO){stop("The number of fixed effect names doesn't match the number of traits")}}
if(!is.null(random_effect_name)){if(length(random_effect_name)!=O0O0OOO0O0O0OOOOO0OO){stop("The number of random effect names doesn't match the number of traits")}}
if(!is.null(covariate_effect_name)){if(length(covariate_effect_name)!=O0O0OOO0O0O0OOOOO0OO){stop("The number of covariate effect names doesn't match the number of traits")}}



#User_define模型下，添加多个遗传效应 放入 random_effect_name, 后续只需要在 OOO0O0O0OOOOO0O0O0O0 中，把遗传效应的位置替换成 Id的位置，其余不变
if(("User_define"%in%analysis_model)&length(relationship_name)>=2){
User_define_effect=paste0("User_define",1:(length(relationship_name)-1))
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(random_effect_name)){
random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]=c(random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]],User_define_effect)
random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]=c(genetic_effect_name,User_define_effect,setdiff(random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]],c(genetic_effect_name,User_define_effect)))
}
OOOOO0O0O0O0O0OOO0O0=c(OOOOO0O0O0O0O0OOO0O0,User_define_effect)
}


#添加显性效应 放入 random_effect_name, 后续只需要在 OOO0O0O0OOOOO0O0O0O0 中，把dominance效应的位置替换成 Id的位置，其余不变
if(included_dominance_effect==TRUE){
dominance_effect="Added_dominance_effect"
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(random_effect_name)){
random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]=c(random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]],dominance_effect)
random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]=c(genetic_effect_name,dominance_effect,setdiff(random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]],c(genetic_effect_name,dominance_effect)))
}
OOOOO0O0O0O0O0OOO0O0=c(OOOOO0O0O0O0O0OOO0O0,dominance_effect)
}


#添加永久环境效应 放入 random_effect_name, 后续只需要在 OOO0O0O0OOOOO0O0O0O0 中，把permanent的位置替换成 Id的位置，其余不变
if(identical(included_permanent_effect,TRUE)){included_permanent_effect=rep(list(included_permanent_effect),length(target_trait_name))}
if(identical(included_permanent_effect,FALSE)){included_permanent_effect=rep(list(included_permanent_effect),length(target_trait_name))}

permanent_effect="pe_effect"
if(permanent_effect%in%OOOOO0O0O0O0O0OOO0O0){stop("pe_effect already stands for permant effect,please change this effect name in the phe_col_name!!! ")}
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(random_effect_name)){

if(included_permanent_effect[[OOO0OOO0O0OOO0O0OOO0OOOO]]==TRUE){

random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]=c(random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]],permanent_effect)
random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]=c(genetic_effect_name,permanent_effect,setdiff(random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]],c(genetic_effect_name,permanent_effect)))

}
}

if(TRUE %in% included_permanent_effect)OOOOO0O0O0O0O0OOO0O0=c(OOOOO0O0O0O0O0OOO0O0,permanent_effect)


#添加母性效应 放入 random_effect_name, 后续只需要将母性效应对应的随机效应编号 设置与 遗传效应的编号一样即可

if(!is.null(maternal_effect_name)){

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(random_effect_name)){

maternal=maternal_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]] #列表

if(!is.null(maternal)){

if(!maternal%in%OOOOO0O0O0O0O0OOO0O0){stop("Maternal effect name could not find in provided phe_col_names; ")}

random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]=c(random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]],maternal)
random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]=c(genetic_effect_name,maternal,setdiff(random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]],c(genetic_effect_name,maternal)))


}
}}



#记录所有性状的名称，并自定义随机效应序号

random_lev=as.factor(unique(do.call(c,random_effect_name))) #统计所有的随机效应名称，并将随机效应赋值
levels(random_lev)=c(genetic_effect_name,setdiff(as.character(random_lev),genetic_effect_name))
OOO0OOOOOOOOOOOOO0O0=as.numeric(random_lev)
OOOOOOOOOOOOO0O0OOOO=as.character(random_lev)
OOO0OOO0O0OOOOOOO0O0=OOO0OOOOOOOOOOOOO0O0[match(do.call(c,random_effect_name),OOOOOOOOOOOOO0O0OOOO)]



#添加残差效应序号：
O0OOO0OOOOOOOOOOO0OO=length(OOO0OOOOOOOOOOOOO0O0)+1
OOO0OOO0O0OOOOOOO0O0=c(OOO0OOO0O0OOOOOOO0O0,rep(O0OOO0OOOOOOOOOOO0OO,O0O0OOO0O0O0OOOOO0OO))


#多个性状循环
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:O0O0OOO0O0O0OOOOO0OO){
trait=target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO]
fixed=fixed_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]
random=random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]
if(!genetic_effect_name%in%random){stop("Genetic effect name could not find in provided random effect name! ")}
random=c(genetic_effect_name,setdiff(random,c(genetic_effect_name)))#确保 genetic effect name 为随机效应第一位
covariate=covariate_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]
regression=random_regression_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]
maternal=maternal_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]] #列表

#判断性状名称、效应名称是否都在 提供的 phe_col_names 中
if(is.na(match(trait,O0OOO0O0OOO0O0OOO0OO))){
stop(paste0("Trait:",trait," could not find in provided phe_col_names; "))
}else{trait_pos=match(trait,O0OOO0O0OOO0O0OOO0OO)}



if(!is.null(fixed)){
OOO0O0O0O0OOOOO0O0OO=match(fixed,OOOOO0O0O0O0O0OOO0O0)
if(NA%in%OOO0O0O0O0OOOOO0O0OO){stop(paste0("Fixed effect:",fixed[is.na(OOO0O0O0O0OOOOO0O0OO)]," could not find in provided phe_col_names; "))}
}else{OOO0O0O0O0OOOOO0O0OO=NULL}

if(!is.null(random)){

random_pos=match(random,OOOOO0O0O0O0O0OOO0O0)

if(NA%in%random_pos){stop(paste0("Random effect:",random[is.na(random_pos)]," could not find in provided phe_col_names; "))}
}else{random_pos=NULL}


if(!is.null(covariate)){
OOOOOOO0O0OOOOOOOOO0=match(covariate,O0OOO0O0OOO0O0OOO0OO)
if(NA%in%OOOOOOO0O0OOOOOOOOO0){stop(paste0("Covariate effect:",covariate[is.na(OOOOOOO0O0OOOOOOOOO0)]," could not find in provided phe_col_names; "))}
}else {OOOOOOO0O0OOOOOOOOO0=NULL}


#随机回归效应-在covariate生成后进行修改
#嵌套的随机效应必须要出现在随机效应中
#嵌套的随机效应同时需要考虑在PRIOR文件中，每个随机效应下，有多少个对角的效应 eg.随机效应代码为，嵌套两个，那么就需要变为：
																				# 1 1 1 1
																				# 1 2 2 1
																				# 1 3 3 1
linked_number=0
if(!is.null(regression)){

if(TRUE %in% included_permanent_effect){cat("Permanent effect name in random regression model is:pe_effect! \n")}

##确定回归系数的种类
#regression_name_coef_set=sapply(regression,function(x)unlist(strsplit(x,split = "&"))[1])
#regression_name_coef_set=unique(regression_name_coef_set)
#n_regre=length(regression_name_coef_set)
#
#regression_split=strsplit(regression,split = "&")
#regression_name_effect_set=as.list(rep(NA,n_regre))
##确定每个回归系数匹配的效应
#for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:n_regre){
#
#for(OOO0OOO0O0OOO0O0OOO0OOOOO in 1:length(regression_split)){
#regression_name_coef=regression_split[[OOO0OOO0O0OOO0O0OOO0OOOOO]][1]       # 第一个名称是 多项式回归系数的名称
#regression_name_effect=regression_split[[OOO0OOO0O0OOO0O0OOO0OOOOO]][-1]    # 非第一个，剩余的是 该多项式回归系数嵌套的所有随机效应
#
#if(regression_name_coef==regression_name_coef_set[OOO0OOO0O0OOO0O0OOO0OOOO]){
#regression_name_effect_set[[OOO0OOO0O0OOO0O0OOO0OOOO]]=c(regression_name_effect_set[[OOO0OOO0O0OOO0O0OOO0OOOO]],regression_name_effect)
#}
#}
#tmp_reg=unique(as.character(na.omit(regression_name_effect_set[[OOO0OOO0O0OOO0O0OOO0OOOO]])))
#linked_number=linked_number+length(tmp_reg)
#pos_In_blanket=match(match(tmp_reg,OOOOO0O0O0O0O0OOO0O0),
#                                   c(OOO0O0O0O0OOOOO0O0OO,random_pos)) # 嵌套的变量在括号中的位置
#pos_In_blanket=paste(pos_In_blanket,collapse = " ")    # 合并成一个字符串
#
#regression_name_coef_effect_pos=paste0(match(regression_name_coef,O0OOO0O0OOO0O0OOO0OO),
#                                                   "(",
#										   pos_In_blanket,
#										   ") "
#										   )
#print(regression_name_coef_effect_pos)										   
#OOOOOOO0O0OOOOOOOOO0=c(OOOOOOO0O0OOOOOOOOO0,regression_name_coef_effect_pos)
#}

reg_k=length(regression) #随机回归效应的个数

for(OOO0OOO0O0OOO0O0OOO0OOOOO in 1:reg_k){

regression_name=unlist(strsplit(regression[OOO0OOO0O0OOO0O0OOO0OOOOO],split = "&"))
  
regression_name_coef=regression_name[1]       # 第一个名称是 多项式回归系数的名称
regression_name_effect=regression_name[-1]             #非第一个，剩余的是 该多项式回归系数嵌套的所有随机效应

linked_number=linked_number+length(regression_name_effect)-1 #针对括号内有2个及以上嵌套效应的情况，需要将协变量的总数目进行修改

#if(FALSE%in%(regression_name_effect%in%random)){stop("Random effect in random_regression_effect_name don't exist in the random_effect_name!")}

if(FALSE%in%(regression_name_effect%in%c(fixed,random))){
stop("Regression effect in random_regression_effect_name don't exist in the random_effect_name or fixed_effect_name!")
}


pos_In_blanket=match(match(regression_name_effect,OOOOO0O0O0O0O0OOO0O0),
                                   c(OOO0O0O0O0OOOOO0O0OO,random_pos)) # 嵌套的变量在括号中的位置
pos_In_blanket=paste(pos_In_blanket,collapse = " ")    # 合并成一个字符串

								   
regression_name_coef_effect_pos=paste0(match(regression_name_coef,O0OOO0O0OOO0O0OOO0OO),
                                                   "(",
										   pos_In_blanket,
										   ") "
										   )
OOOOOOO0O0OOOOOOOOO0=c(OOOOOOO0O0OOOOOOOOO0,regression_name_coef_effect_pos)
}

}


OOO0O0O0OOOOO0O0O0O0=paste0(trait_pos," 0 ",length(OOO0O0O0O0OOOOO0O0OO)+length(random_pos)," ",paste(c(OOO0O0O0O0OOOOO0O0OO,random_pos),collapse=" "))
random_line=paste(c(length(random),OOO0OOOOOOOOOOOOO0O0[match(random,OOOOOOOOOOOOO0O0OOOO)]),collapse=" ")
OOO0O0O0O0O0OOO0O0OO=paste(c(length(OOOOOOO0O0OOOOOOOOO0)+linked_number,OOOOOOO0O0OOOOOOOOO0),collapse=" ")


#考虑永久环境效应，将  OOO0O0O0OOOOO0O0O0O0 中的 Added_permanent_effect 的位置替换成 Id 位置
if(TRUE%in%included_permanent_effect[[OOO0OOO0O0OOO0O0OOO0OOOO]]){
O0OOO0O0OOO0O0OOOOOO=match(permanent_effect,OOOOO0O0O0O0O0OOO0O0)
random_pos[random_pos%in%O0OOO0O0OOO0O0OOOOOO]=match(genetic_effect_name,OOOOO0O0O0O0O0OOO0O0)
OOO0O0O0OOOOO0O0O0O0=paste0(trait_pos," 0 ",length(OOO0O0O0O0OOOOO0O0OO)+length(random_pos)," ",paste(c(OOO0O0O0O0OOOOO0O0OO,random_pos),collapse=" "))
}
 
#考虑显性效应，将  OOO0O0O0OOOOO0O0O0O0 中的 Added_dominance_effect 的位置替换成 Id 位置
if(included_dominance_effect==TRUE){
dominance_effect_pos=match(dominance_effect,OOOOO0O0O0O0O0OOO0O0)
random_pos[random_pos%in%dominance_effect_pos]=match(genetic_effect_name,OOOOO0O0O0O0O0OOO0O0)
OOO0O0O0OOOOO0O0O0O0=paste0(trait_pos," 0 ",length(OOO0O0O0O0OOOOO0O0OO)+length(random_pos)," ",paste(c(OOO0O0O0O0OOOOO0O0OO,random_pos),collapse=" "))
}

 
#考虑User_define模型，将  OOO0O0O0OOOOO0O0O0O0 中的 User_define_effect 的位置替换成 Id 位置
if(("User_define"%in%analysis_model)&length(relationship_name)>=2){
User_define_effect_pos=match(User_define_effect,OOOOO0O0O0O0O0OOO0O0)
random_pos[random_pos%in%User_define_effect_pos]=match(genetic_effect_name,OOOOO0O0O0O0O0OOO0O0)
OOO0O0O0OOOOO0O0O0O0=paste0(trait_pos," 0 ",length(OOO0O0O0O0OOOOO0O0OO)+length(random_pos)," ",paste(c(OOO0O0O0O0OOOOO0O0OO,random_pos),collapse=" "))
}


#考虑母性效应，将  random_line 中的 maternal effect 的随机效应编号替换为遗传效应的随机效应编号
if(!is.null(maternal)){

#将maternal effect的水平修改成 和  遗传效应水平一致
#其他随机效应的水平均减去1

tmp_OOO0OOOOOOOOOOOOO0O0=OOO0OOOOOOOOOOOOO0O0

tmp_OOO0OOOOOOOOOOOOO0O0[match(maternal,OOOOOOOOOOOOO0O0OOOO)]=tmp_OOO0OOOOOOOOOOOOO0O0[match(genetic_effect_name,OOOOOOOOOOOOO0O0OOOO)]

tmp_OOO0OOOOOOOOOOOOO0O0[tmp_OOO0OOOOOOOOOOOOO0O0!=1]=tmp_OOO0OOOOOOOOOOOOO0O0[tmp_OOO0OOOOOOOOOOOOO0O0!=1]-1

random_line=paste(c(length(random),tmp_OOO0OOOOOOOOOOOOO0O0[match(random,OOOOOOOOOOOOO0O0OOOO)]),collapse=" ")

}
 

#考虑social genetic effect, 需要额外对DIR文件中，固定效应-随机效应行 和 随机效应编号 行进行修改
if(TRUE%in%include_social_effect){
if(TRUE%in%include_social_effect[[OOO0OOO0O0OOO0O0OOO0OOOO]]){

#修改固定效应-随机效应行， 在integer_group_name后添加(0)

random_pos[match(integer_group_names,random)]=paste0(random_pos[match(integer_group_names,random)],"(","0",")")
OOO0O0O0OOOOO0O0O0O0=paste0(trait_pos," 0 ",length(OOO0O0O0O0OOOOO0O0OO)+length(random_pos)," ",paste(c(OOO0O0O0O0OOOOO0O0OO,random_pos),collapse=" "))

#修改随机效应编号 行，形如 +1+1+1

temp=OOO0OOOOOOOOOOOOO0O0[match(random,OOOOOOOOOOOOO0O0OOOO)]

temp1=OOO0OOOOOOOOOOOOO0O0[match(setdiff(random,integer_group_names),OOOOOOOOOOOOO0O0OOOO)] #从随机效应名称中排除 integer_group_names
temp2=paste(rep(1,length(random)-length(temp1)),collapse = "+")

random_line=paste(c(length(random),temp1,temp2),collapse=" ")

}}



MODEL[OOO0OOO0O0OOO0O0OOO0OOOO+2+O0O0OOO0O0O0OOOOO0OO,1]=OOO0O0O0OOOOO0O0O0O0
MODEL[(OOO0OOO0O0OOO0O0OOO0OOOO+2+O0O0OOO0O0O0OOOOO0OO)+O0O0OOO0O0O0OOOOO0OO,1]=random_line
MODEL[(OOO0OOO0O0OOO0O0OOO0OOOO+2+O0O0OOO0O0O0OOOOO0OO)+2*O0O0OOO0O0O0OOOOO0OO,1]=OOO0O0O0O0O0OOO0O0OO

}

#残差限定
if(is.null(residual_cov_trait)){
MODEL[(O0O0OOO0O0O0OOOOO0OO+2+O0O0OOO0O0O0OOOOO0OO)+2*O0O0OOO0O0O0OOOOO0OO+1,1]="0"
}else{
if(length(residual_cov_trait)>O0O0OOO0O0O0OOOOO0OO){stop("The number of residual covariance traits is larger than the number of analyse traits \n")}
residual_cov_trait_n=0

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(residual_cov_trait)){

O0O0O0OOO0OOO0O0O0O0=residual_cov_trait[[OOO0OOO0O0OOO0O0OOO0OOOO]]

if(!is.null(O0O0O0OOO0OOO0O0O0O0)){
residual_cov_trait_n=residual_cov_trait_n+1
trait1=O0O0O0OOO0OOO0O0O0O0[1]
trait2=O0O0O0OOO0OOO0O0O0O0[2]
MODEL[(O0O0OOO0O0O0OOOOO0OO+2+O0O0OOO0O0O0OOOOO0OO)+2*O0O0OOO0O0O0OOOOO0OO+1+OOO0OOO0O0OOO0O0OOO0OOOO,1]=paste0(match(trait1,target_trait_name)," ",match(trait2,target_trait_name))
}
}

MODEL[(O0O0OOO0O0O0OOOOO0OO+2+O0O0OOO0O0O0OOOOO0OO)+2*O0O0OOO0O0O0OOOOO0OO+1,1]=residual_cov_trait_n
}



#$VAR_STR
#定义 $VAR_STR 中的文件名称
if(is.null(relationship_name)){
stop("Please provide relationship name!")
}else{

addtive_relationship_name=relationship_name[1]
if(analysis_model%in%c("GBLUP_AD")){dominance_relationship_name=relationship_name[2]}
if(analysis_model%in%c("SSBLUP_A")){SSBLUP_G_matrix_name=relationship_name[2]}
}

if(analysis_model=="PBLUP_A"){
VAR_STR[1,1]=paste0("$VAR_STR 1 PED ",ped_inbred_number," ASCII   ",relationship_path,"/",addtive_relationship_name)
}else if (analysis_model=="PBLUP_AD"){
VAR_STR[1,1]=paste0("$VAR_STR 1 PED ",ped_inbred_number," ASCII   ",relationship_path,"/",addtive_relationship_name)
VAR_STR[3,1]=paste0("$VAR_STR 2 COR ASCII     ",relationship_path,"/",dominance_relationship_name)
}else if (analysis_model=="GBLUP_A"){
VAR_STR[1,1]=paste0("$VAR_STR 1 COR ASCII     ",relationship_path,"/",addtive_relationship_name)
}else if (analysis_model=="GBLUP_AD"){
VAR_STR[1,1]=paste0("$VAR_STR 1 COR ASCII     ",relationship_path,"/",addtive_relationship_name)
VAR_STR[3,1]=paste0("$VAR_STR 2 COR ASCII     ",relationship_path,"/",dominance_relationship_name)
}else if (analysis_model=="SSBLUP_A"){
VAR_STR[1,1]=paste0("$VAR_STR 1 PGMIX 1 ASCII  ",relationship_path,"/",addtive_relationship_name,"    ",
                      relationship_path,"/",IND_geno_file_name,"   ",
				   relationship_path,"/",SSBLUP_G_matrix_name,"  ",
				   SSBLUP_omega,"  G-ADJUST")
}else if (analysis_model=="User_define"){

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(relationship_name)){
VAR_STR[OOO0OOO0O0OOO0O0OOO0OOOO,1]=paste0("$VAR_STR ",OOO0OOO0O0OOO0O0OOO0OOOO," COR ASCII     ",relationship_path,"/",relationship_name[OOO0OOO0O0OOO0O0OOO0OOOO])
}}


#$SOLUTION
SOLUTION[1,1]="$RESIDUALS ASCII"
SOLUTION[2,1]="$SOLUTION"

#$PRIOR

if(include_social_effect==FALSE){  #social effect 不知道该如何构建PRIOR文件，

PRIOR[1,1]="$PRIOR"
if(!is.null(provided_prior_file_name)&!is.null(provided_prior_file_path)){
OOOOOOO0O0O0OOOOOOO0=fread(paste0(provided_prior_file_path,"/",provided_prior_file_name),data.table=F)
OOOOOOO0O0O0OOOOOOO0=as.matrix(OOOOOOO0O0O0OOOOOOO0)
PRIOR[2:(nrow(OOOOOOO0O0O0OOOOOOO0)+1),1:ncol(OOOOOOO0O0O0OOOOOOO0)]=OOOOOOO0O0O0OOOOOOO0
}else{

iteration_pos=1
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:(length(OOO0OOOOOOOOOOOOO0O0)+1)){  #随机效应的数量 

	for(OOO0OOO0O0OOO0O0OOO0OOOOO in 1:sum(OOO0OOO0O0OOOOOOO0O0%in%OOO0OOO0O0OOO0O0OOO0OOOO)){
	
		for(k in 1:OOO0OOO0O0OOO0O0OOO0OOOOO){
			iteration_pos=iteration_pos+1
			PRIOR[iteration_pos,1]=OOO0OOO0O0OOO0O0OOO0OOOO
			PRIOR[iteration_pos,2]=OOO0OOO0O0OOO0O0OOO0OOOOO
			PRIOR[iteration_pos,3]=k
			if(OOO0OOO0O0OOO0O0OOO0OOOOO!=k){
			PRIOR[iteration_pos,4]=0.5
			}else {PRIOR[iteration_pos,4]=1}
	}}}
}	
}

#残差限定

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(residual_cov_trait)){

O0O0O0OOO0OOO0O0O0O0=residual_cov_trait[[OOO0OOO0O0OOO0O0OOO0OOOO]]
O0OOO0OOOOOOOOOOO0OO=O0OOO0OOOOOOOOOOO0OO

if(!is.null(O0O0O0OOO0OOO0O0O0O0)){

O0OOO0OOOOOOOOOOOOO0=match(O0O0O0OOO0OOO0O0O0O0,target_trait_name)

trait1_pos=O0OOO0OOOOOOOOOOOOO0[1]
trait2_pos=O0OOO0OOOOOOOOOOOOO0[2]
PRIOR=PRIOR[!(PRIOR[,1]%in%O0OOO0OOOOOOOOOOO0OO&PRIOR[,2]%in%trait1_pos&PRIOR[,3]%in%trait2_pos),]
PRIOR=PRIOR[!(PRIOR[,1]%in%O0OOO0OOOOOOOOOOO0OO&PRIOR[,2]%in%trait2_pos&PRIOR[,3]%in%trait1_pos),]		
}}


#不显示 PRIOR结果， 因为针对随机回归模型，PRIOR结果存在些许问题，使用DMU默认值即可
#如果不是人为提供PRIOR文件的话，PRIOR设置为NULL
if(!is.null(provided_prior_file_name)&!is.null(provided_prior_file_path)){

}else{
PRIOR=NULL
}


# module_parameter

if(dmu_module=="dmuai"){

OOO0OOO0OOOOO0O0OOOO[1,1]="$DMUAI"
OOO0OOO0OOOOO0O0OOOO[2,1]="10"
OOO0OOO0OOOOO0O0OOOO[3,1]=iteration_criteria
OOO0OOO0OOOOO0O0OOOO[4,1]="1.0d-6"
OOO0OOO0OOOOO0O0OOOO[5,1]="1"
OOO0OOO0OOOOO0O0OOOO[6,1]="0"

}else if(dmu_module=="dmu5"){
OOO0OOO0OOOOO0O0OOOO[1,1]="$DMU5"
OOO0OOO0OOOOO0O0OOOO[2,2]="1000000 1e-9"
OOO0OOO0OOOOO0O0OOOO[3,3]="512"
}






#输出结果
#DIR=rbind(COMMENT,ANALYSE,DATA,VARIABLE,VAR_STR,MODEL,SOLUTION,PRIOR,OOO0OOO0OOOOO0O0OOOO)
#DIR=DIR[!DIR[,1]%in%"",]

COMMENT=COMMENT[!COMMENT[,1]%in%"",]
ANALYSE=ANALYSE[!ANALYSE[,1]%in%"",]
DATA=DATA[!DATA[,1]%in%"",]
VARIABLE=VARIABLE[!VARIABLE[,1]%in%"",]
VAR_STR=VAR_STR[!VAR_STR[,1]%in%"",]
MODEL=MODEL[!MODEL[,1]%in%"",]
SOLUTION=SOLUTION[!SOLUTION[,1]%in%"",]
PRIOR=PRIOR[!PRIOR[,1]%in%"",]
OOO0OOO0OOOOO0O0OOOO=OOO0OOO0OOOOO0O0OOOO[!OOO0OOO0OOOOO0O0OOOO[,1]%in%"",]

DIR=rbind(COMMENT,"",ANALYSE,"",DATA,"",VARIABLE,"",VAR_STR,"",MODEL,"",SOLUTION,"",PRIOR,"",OOO0OOO0OOOOO0O0OOOO)

setwd(output_DIR_path)
if(is.null(output_DIR_name)){output_DIR_name="Trait.DIR"}
write.table(DIR,output_DIR_name,quote=F,row.names=F,col.names=F,sep=" ")
return(list(DIR=OOO0O0OOO0O0O0OOOOO0(DIR),random_effect_name=random_effect_name))
}


OOO0O0OOO0O0O0OOOOO0<-function(data){
union_data=matrix(NA,nrow=nrow(data),ncol=1)
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:nrow(data)){
union_data[OOO0OOO0O0OOO0O0OOO0OOOO,1]=paste(data[OOO0OOO0O0OOO0O0OOO0OOOO,],collapse=" ")
}
return(union_data)
}

#

generate_social_phe<-function(phe=NULL,
							  integer_n=3,
							  social_effect_group_name=NULL,
							  genetic_effect_name=NULL
							  ){

if(!social_effect_group_name%in%colnames(phe)){stop("Group name is not included in provided phe_col_names")}
if(!genetic_effect_name%in%colnames(phe)){stop("Genetic effect name is not included in provided phe_col_names")}
phe=as.data.frame(phe,stringsAsFactors=F)

phe_col_names=colnames(phe)
O0OOOOOOOOO0OOO0OOO0=length(phe_col_names)
O0O0OOO0O0OOOOOOOOO0=O0OOOOOOOOO0OOO0OOO0-integer_n
OOOOO0O0O0O0O0OOO0O0=phe_col_names[1:integer_n]
O0OOO0O0OOO0O0OOO0OO=phe_col_names[(integer_n+1):O0OOOOOOOOO0OOO0OOO0]

res_OOOOO0O0O0O0O0OOO0O0=setdiff(OOOOO0O0O0O0O0OOO0O0,c(genetic_effect_name,social_effect_group_name)) #除遗传和group之外的整型变量名称
res_O0OOO0O0OOO0O0OOO0OO=O0OOO0O0OOO0O0OOO0OO #剩余的实型变量名称


#统计最大的group_size
group_size=as.data.frame(table(phe[,social_effect_group_name]),stringsAsFactors=F)
group_size[,1]=as.numeric(group_size[,1])
max_size=max(group_size[,2])-1
colnames(group_size)=c("Group","Size")

if(max_size==0){stop("The number of indivuals in each group is equal to 1, please check your group record!!!")}


#构建group_ind_set和group的对应列表
group_ind_set <- vector(mode = "list", length = nrow(group_size))
group_i_set=as.numeric(group_size$Group)
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:nrow(group_size)){
group_ind_set[[OOO0OOO0O0OOO0O0OOO0OOOO]]=phe[phe[,social_effect_group_name]%in%group_i_set[OOO0OOO0O0OOO0O0OOO0OOOO],genetic_effect_name]
}


#group信息和 id名称所在的列
group_column=match(social_effect_group_name,colnames(phe))-1
id_column=match(genetic_effect_name,colnames(phe))-1


group_result=get_rest_id(ref=group_ind_set,ref_index=group_i_set,
                 phe=as.matrix(phe),column_n=max_size,
			  group_column=group_column ,id_column=id_column)

integer_group=group_result$integer_group
real_group=group_result$real_group
rm(group_result);gc();

integer_group_names=paste0("Gr_id",1:ncol(integer_group))
real_group_names=paste0("Status_Gr_id",1:ncol(real_group))

colnames(integer_group)=integer_group_names
colnames(real_group)=real_group_names

phe=cbind(phe,integer_group,real_group)


#按整型和实型将表型数据的列重新排列
phe=phe[,c(genetic_effect_name,res_OOOOO0O0O0O0O0OOO0O0,social_effect_group_name,integer_group_names,
            res_O0OOO0O0OOO0O0OOO0OO,real_group_names)]

return(list(phe=phe,integer_n=integer_n+length(integer_group_names),
			integer_group_names=integer_group_names,real_group_names=real_group_names))			  
}
geno_format<-function(             
		input_data_type=NULL, #"Plink" , "Hapmap" , "VCF" , "BLUPF90","Numeric"  "Haplotype"
		input_data_path=NULL,
		input_data_name=NULL,  #hapmap 为全名，vcf 和 plink 默认不包括后缀，为 name.vcf name.ped name.map 
		input_data_hmp=NULL,
		input_data_plink_ped=NULL,
		input_data_plink_map=NULL,
		input_data_blupf90=NULL,
		input_data_blupf90_map=NULL, #格式和 Haplotype-map格式一样
		input_data_numeric=NULL,	 
		input_data_numeric_map=NULL, #格式和 Haplotype-map格式一样
		input_data_vcf=NULL,
		input_data_haplotype_hap=NULL,
		input_data_haplotype_map=NULL,
		input_data_haplotype_sample=NULL,
		bigmemory_cal=FALSE,
		bigmemory_data_type="integer",
		bigmemory_data_path=getwd(),
		bigmemory_data_name="blupADC",
		phased_symbol="/",
		phased_genotype=FALSE,
		haplotype_window_nSNP=NULL,
		haplotype_window_kb=NULL,
		haplotype_window_block=NULL,	
		miss_base=NULL,		
		miss_base_num=0,
		output_data_name="blupADC_geno_convertion",
		output_data_type=NULL,  # "Plink" , "Hapmap" , "VCF"  "BLUPF90"  "Numeric"  "Haplotype"
		output_data_path=NULL,
		return_result=FALSE,
		cpu_cores=1  #调用的cpu数目，用于加速计算
		){	
options (warn =-1)	
library(data.table)
library(bigmemory)
options(bigmemory.allow.dimnames=TRUE)
options(bigmemory.typecast.warning=FALSE)
if(!is.null(output_data_path)){
if(!file.exists(output_data_path)){
dir.create(output_data_path,recursive=T)
cat("Output data path is not exist,software will constructs this path  automatically \n")
}
setwd(output_data_path)
}

output_data_type=unique(output_data_type)
if(bigmemory_cal==TRUE){return_result=TRUE}
#程序自动判断输入类型
	
input=O0OOO0O0OOO0O0O0OOOO(input_data_type=input_data_type,input_data_hmp=input_data_hmp,input_data_plink_ped=input_data_plink_ped,
                                          input_data_plink_map=input_data_plink_map,input_data_blupf90=input_data_blupf90,
								    input_data_numeric=input_data_numeric,input_data_vcf=input_data_vcf,
									input_data_haplotype_hap=input_data_haplotype_hap,input_data_haplotype_map=input_data_haplotype_map,
									input_data_haplotype_sample=input_data_haplotype_sample,
								     miss_base=miss_base)

input_data_type=input$input_data_type
miss_base=input$miss_base

if(identical(input_data_type,output_data_type)){return(NULL)}
##############################################初始化

data_hmp=NULL;data_ped=NULL;data_map=NULL;data_blupf90=NULL;data_numeric=NULL;data_vcf=NULL
data_haplotype_hap=NULL;data_haplotype_map=NULL;data_haplotype_sample=NULL;output_block=NULL
haplotype_allele=NULL;

#判断是否存在本地文件
if(!is.null(input_data_path)&!is.null(input_data_name)){
if(input_data_type=="Plink"){
if(!file.exists(paste0(input_data_path,"/",input_data_name,".ped"))|!file.exists(paste0(input_data_path,"/",input_data_name,".map"))){
stop("Couldn't found provided input data!")}

}else if(input_data_type=="Haplotype"){
if(!file.exists(paste0(input_data_path,"/",input_data_name,".hap"))|!file.exists(paste0(input_data_path,"/",input_data_name,".sample"))|!file.exists(paste0(input_data_path,"/",input_data_name,".map"))){
stop("Couldn't found provided input data!")}
}else{
if(!file.exists(paste0(input_data_path,"/",input_data_name))){stop("Couldn't found provided input data!")}
}
}

#判断输入的bigmeory类型
if(!bigmemory_data_type%in%(c("double","float","integer","short","char"))){

stop("Unknow input bigmemory_data_type, bigmemory_data_type should be double or float or integer or short or char ")

}


####################################input data type is Hapmap#################################### 
if(input_data_type=="Hapmap"){
if(!is.null(input_data_path)&!is.null(input_data_name)){
cat("Start read the Hapmap format genotype data \n")
input_data_hmp=fread(paste0(input_data_path,"/",input_data_name),data.table=F)
cat("Complete read the Hapmap format genotype data  \n")}

if(!is.matrix(input_data_hmp)){input_data_hmp=as.matrix(input_data_hmp)}
input_data_hmp[,3]=as.numeric(input_data_hmp[,3])
input_data_hmp[,4]=as.numeric(input_data_hmp[,4])

IND_name=colnames(input_data_hmp)[-c(1:11)]
n_snp=nrow(input_data_hmp)
n_ind=ncol(input_data_hmp)-11

output_type_number=get_output_type_number(input_data_type=input_data_type,
										  output_data_type=output_data_type,
										  bigmemory_cal=bigmemory_cal,
										  phased_genotype=phased_genotype)

if(length(output_type_number)==0){stop("Please provided standard output data type!")}

for(type in output_type_number){
if(exists(".Random.seed")){fixed_openMP(.Random.seed);set.seed(19980204)}

if(type==1){
data_ped=hapmap_convertion(input_data_hmp,type,miss_base,phased_symbol,miss_base_num,cpu_cores)$ped;
data_map=input_data_hmp[,c(3,1,4,4)];data_map[,3]=as.numeric(data_map[,3])/1000000;
data_map[,3]=as.numeric(data_map[,3])
data_map[,4]=as.numeric(data_map[,4])
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_ped),paste0(output_data_name,".ped"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(data_map),paste0(output_data_name,".map"),quote=F,row.names=F,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_map,data_ped);gc();}
}

if(type==2){data_vcf=hapmap_convertion(input_data_hmp,type,miss_base,phased_symbol,miss_base_num,cpu_cores)$vcf;
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(get_vcf_header()),paste0(output_data_name,".vcf"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(data_vcf),paste0(output_data_name,".vcf"),append=TRUE,quote=F,row.names=F,col.names=T,sep="\t")}
if(return_result==FALSE){rm(data_vcf);gc();}
}

if(type==3){data_numeric=hapmap_convertion(input_data_hmp,type,miss_base,phased_symbol,miss_base_num,cpu_cores)$numeric
rownames(data_numeric)=IND_name
if(!is.null(output_data_path)&!is.null(output_data_name)){
fwrite(data.frame(data_numeric),paste0(output_data_name,".numeric"),quote=F,row.names=T,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_numeric);gc();}
}

if(type==4){data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=hapmap_convertion(input_data_hmp,type,miss_base,phased_symbol,miss_base_num,cpu_cores)$blupf90,stringsAsFactors=F)
if(!is.null(output_data_path)&!is.null(output_data_name)){
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=F,sep=" ")}
if(return_result==FALSE){rm(data_blupf90);gc();}
}

#构建bigmemory对象 
if(type==13){

if(exists("pBigMat_num_address")){rm(pBigMat_num_address)}
gc()
delete_bigmemory_file("numeric",bigmemory_data_name,bigmemory_data_path,FALSE);

pBigMat_num_address=make_bigmemory_object_address_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
																
hapmap_to_numeric_memory_cpp(input_data_hmp,pBigMat_num_address,miss_base,miss_base_num,cpu_cores,type=bigmemory_data_type)									 
data_numeric=make_bigmemory_object_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)	
options(bigmemory.allow.dimnames=TRUE)
rownames(data_numeric)=IND_name
colnames(data_numeric)=paste0("SNP",1:ncol(data_numeric))																
cat(paste0("bigmemory-Numeric data is save as ",bigmemory_data_path,"/",bigmemory_data_name,"_numeric \n"))
																
}


}
rm(input_data_hmp);gc();
}


####################################input data type is Plink#################################### 
if(input_data_type=="Plink"){
if(!is.null(input_data_path)&!is.null(input_data_name)){
cat("Start read the Plink format genotype data \n")
input_data_plink_ped=fread(paste0(input_data_path,"/",input_data_name,".ped"),data.table=F)
input_data_plink_map=fread(paste0(input_data_path,"/",input_data_name,".map"),data.table=F)
cat("Complete read the Plink format genotype data \n")}

IND_name=input_data_plink_ped[,1]
n_ind=length(IND_name)
n_snp=nrow(input_data_plink_map)

if(!is.matrix(input_data_plink_map)){input_data_plink_map=as.matrix(input_data_plink_map)}
if(!is.matrix(input_data_plink_ped)){input_data_plink_ped=as.matrix(input_data_plink_ped)}
input_data_plink_map[,1]=as.numeric(input_data_plink_map[,1])
input_data_plink_map[,3]=as.numeric(input_data_plink_map[,3])
input_data_plink_map[,4]=as.numeric(input_data_plink_map[,4])

output_type_number=get_output_type_number(input_data_type=input_data_type,
										  output_data_type=output_data_type,
										  bigmemory_cal=bigmemory_cal,
										  phased_genotype=phased_genotype)

if(length(output_type_number)==0){stop("Please provided standard output data type!")}

for(type in output_type_number){
if(exists(".Random.seed")){fixed_openMP(.Random.seed);set.seed(19980204)}
if(type==1){data_hmp=plink_convertion(input_data_plink_ped,input_data_plink_map,	type, miss_base,phased_symbol,miss_base_num,cpu_cores)$hmp
data_hmp[,1]=input_data_plink_map[,2]
data_hmp[,3]=input_data_plink_map[,1]
data_hmp[,4]=input_data_plink_map[,4]
data_hmp[,3]=as.numeric(data_hmp[,3])
data_hmp[,4]=as.numeric(data_hmp[,4])
data_hmp[,c(2,5:11)]="NA"
colnames(data_hmp)=c("rs#","alleles","chrom","pos","strand","assembly","center","protLSID","assayLSID","panelLSID","QCcode",IND_name)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.table(data_hmp),paste0(output_data_name,".hmp.txt"),quote=F,row.names=F,col.names=T,sep="\t")}
if(return_result==FALSE){rm(data_hmp);gc();}
}

if(type==2){data_vcf=plink_convertion(input_data_plink_ped,input_data_plink_map,	type, miss_base,phased_symbol,miss_base_num,cpu_cores)$vcf
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(get_vcf_header()),paste0(output_data_name,".vcf"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(data_vcf),paste0(output_data_name,".vcf"),append=TRUE,quote=F,row.names=F,col.names=T,sep="\t")}
if(return_result==FALSE){rm(data_vcf);gc();}
}

if(type==3){data_numeric=plink_convertion(input_data_plink_ped,input_data_plink_map,	type, miss_base,phased_symbol,miss_base_num,cpu_cores)$numeric
rownames(data_numeric)=IND_name
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_numeric),paste0(output_data_name,".numeric"),quote=F,row.names=T,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_numeric);gc();}
}

if(type==4){data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=plink_convertion(input_data_plink_ped,input_data_plink_map,	type, miss_base,phased_symbol,miss_base_num,cpu_cores)$blupf90,stringsAsFactors=F)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=F,sep=" ")}
if(return_result==FALSE){rm(data_blupf90);gc();}
}


#构建bigmemory对象 
if(type==13){

if(exists("pBigMat_num_address")){rm(pBigMat_num_address)}
gc()
delete_bigmemory_file("numeric",bigmemory_data_name,bigmemory_data_path,FALSE);

pBigMat_num_address=make_bigmemory_object_address_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
																
plink_to_numeric_memory_cpp(input_data_plink_ped,input_data_plink_map,pBigMat_num_address,cpu_cores,miss_base,miss_base_num,type=bigmemory_data_type)								 
data_numeric=make_bigmemory_object_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)	
options(bigmemory.allow.dimnames=TRUE)
rownames(data_numeric)=IND_name
colnames(data_numeric)=paste0("SNP",1:ncol(data_numeric))																
cat(paste0("bigmemory-Numeric data is save as ",bigmemory_data_path,"/",bigmemory_data_name,"_numeric \n"))
																
}


}
rm(input_data_plink_ped,input_data_plink_map);gc();
}


####################################input data type is VCF#################################### 
if(input_data_type=="VCF"){
if(!is.null(input_data_path)&!is.null(input_data_name)){
cat("Start read the VCF format genotype data \n")
input_data_vcf=fread(paste0(input_data_path,"/",input_data_name),data.table=F)
cat("Complete read the VCF format genotype data \n")}

if(!is.matrix(input_data_vcf)){input_data_vcf=as.matrix(input_data_vcf)}

input_data_vcf[,1]=as.numeric(input_data_vcf[,1])
input_data_vcf[,2]=as.numeric(input_data_vcf[,2])

IND_name=colnames(input_data_vcf)[-c(1:9)]
input_data_haplotype_map=input_data_vcf[,c(1,3,2,4,5)]


#构建block
if(!is.null(haplotype_window_nSNP)|!is.null(haplotype_window_kb)|!is.null(haplotype_window_block)){
block_result=get_haplotype_block(input_data_haplotype_map,haplotype_window_nSNP,
				                    haplotype_window_kb,haplotype_window_block)
block_hap=block_result$block_rcpp
block_start=block_hap[,1];block_end=block_hap[,2]
output_block=block_result$block
}else{
block_start=1:5
block_end=1:5
output_block=NULL
}

output_type_number=get_output_type_number(input_data_type=input_data_type,
										  output_data_type=output_data_type,
										  bigmemory_cal=bigmemory_cal,
										  phased_genotype=phased_genotype)
n_ind=ncol(input_data_vcf)-9
n_snp=nrow(input_data_vcf)
n_window=length(block_start)
data_haplotype_map=input_data_vcf[,c(1,3,2,4,5)];IND_name=colnames(input_data_vcf)[-c(1:9)]

if(length(output_type_number)==0){stop("Please provided standard output data type!")}
for(type in output_type_number){
if(exists(".Random.seed")){fixed_openMP(.Random.seed);set.seed(19980204)}
if(type==1){data_hmp=vcf_convertion(input_data_vcf,block_start,block_end,type,miss_base,miss_base_num,cpu_cores)$hmp
data_hmp[,1]=input_data_vcf[,3]
data_hmp[,3]=input_data_vcf[,1]
data_hmp[,4]=input_data_vcf[,2]
data_hmp[,c(2,5:11)]="NA"
data_hmp[,3]=as.numeric(data_hmp[,3])
data_hmp[,4]=as.numeric(data_hmp[,4])
colnames(data_hmp)=c("rs#","alleles","chrom","pos","strand","assembly","center","protLSID","assayLSID","panelLSID","QCcode",IND_name)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.table(data_hmp),paste0(output_data_name,".hmp.txt"),quote=F,row.names=F,col.names=T,sep="\t")}
if(return_result==FALSE){rm(data_hmp);gc();}
}

if(type==2){data_ped=vcf_convertion(input_data_vcf,block_start,block_end,type,miss_base,miss_base_num,cpu_cores)$ped
data_map=input_data_vcf[,c(1,3,2,2)];data_map[,3]=(as.numeric(data_map[,3])/1000000)
data_map[,3]=as.numeric(data_map[,3])
data_map[,4]=as.numeric(data_map[,4])
data_ped[,3:6]=0
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_ped),paste0(output_data_name,".ped"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(data_map),paste0(output_data_name,".map"),quote=F,row.names=F,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_map,data_ped);gc();}
}


if(type==8){
tmp=vcf_convertion(input_data_vcf,block_start,block_end,type,miss_base,miss_base_num,cpu_cores)
data_ped=tmp$ped
haplotype_allele=tmp$SNP
origin_haplotype_allele=haplotype_allele
rm(tmp);gc();

haplotype_allele_num=sapply(haplotype_allele,length)
block_pos_allele=output_block$block_pos
block_pos_allele$chr=input_data_haplotype_map[block_pos_allele[,2],1]
block_pos_allele$pos=paste0(input_data_haplotype_map[block_pos_allele[,2],3],"_",
						     input_data_haplotype_map[block_pos_allele[,3],3])

haplotype_allele_chr=rep(block_pos_allele$chr,times=haplotype_allele_num)
haplotype_allele_pos=rep(block_pos_allele$pos,times=haplotype_allele_num)

haplotype_allele=do.call(c,haplotype_allele)
haplotype_allele=data.frame(Chr=haplotype_allele_chr,
							SNP=paste0(1:length(haplotype_allele)),
							cM=0,
							Pos=haplotype_allele_pos,
							Allele=haplotype_allele)

data_map=haplotype_allele[,1:4]
haplotype_allele=origin_haplotype_allele
data_ped[,3:6]=0
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_ped),paste0(output_data_name,".ped"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(data_map),paste0(output_data_name,".map"),quote=F,row.names=F,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_map,data_ped);gc();}
}

if(type==3){data_haplotype_hap=vcf_convertion(input_data_vcf,block_start,block_end,type,miss_base,miss_base_num,cpu_cores)$hap
data_haplotype_map=input_data_vcf[,c(1,3,2,4,5)];IND_name=colnames(input_data_vcf)[-c(1:9)]
data_haplotype_sample=data.frame(IND_name)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_haplotype_hap),paste0(output_data_name,"_haplotype.hap"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(data_haplotype_map),paste0(output_data_name,"_haplotype.map"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(data_haplotype_sample),paste0(output_data_name,"_haplotype.sample"),quote=F,row.names=F,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_haplotype_hap,data_haplotype_map,data_haplotype_sample);gc();}
}

if(type==4){

tmp=vcf_convertion(input_data_vcf,block_start,block_end,type,miss_base,miss_base_num,cpu_cores)
haplotype_allele=tmp$SNP
data_numeric=tmp$numeric
rm(tmp);gc();
rownames(data_numeric)=IND_name

if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_numeric),paste0(output_data_name,".numeric"),quote=F,row.names=T,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_numeric);gc();}
}

if(type==5){
tmp=vcf_convertion(input_data_vcf,block_start,block_end,type,miss_base,miss_base_num,cpu_cores)
haplotype_allele=tmp$SNP
data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=tmp$blupf90,stringsAsFactors=F)
rm(tmp);gc();
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=F,sep=" ")}
if(return_result==FALSE){rm(data_blupf90);gc();}
}

if(type==6){data_numeric=vcf_convertion(input_data_vcf,block_start,block_end,type,miss_base,miss_base_num,cpu_cores)$numeric
rownames(data_numeric)=IND_name
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_numeric),paste0(output_data_name,".numeric"),quote=F,row.names=T,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_numeric);gc();}
}

if(type==7){data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=vcf_convertion(input_data_vcf,block_start,block_end,type,miss_base,miss_base_num,cpu_cores)$blupf90,stringsAsFactors=F)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=F,sep=" ")}
if(return_result==FALSE){rm(data_blupf90);gc();}
}
#构建bigmemory对象 
if(type%in%c(13,14,15)){
delete_bigmemory_file("haplotype_hap",bigmemory_data_name,bigmemory_data_path,FALSE);
pBigMat_hap_address=make_bigmemory_object_address_cpp(nrow=n_snp,ncol=n_ind*2,file_name=paste0(bigmemory_data_name,"_haplotype_hap"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
}

if(type%in%c(14,15)){
delete_bigmemory_file("haplotype_list",bigmemory_data_name,bigmemory_data_path,FALSE);

if(type==14){
delete_bigmemory_file("numeric",bigmemory_data_name,bigmemory_data_path,TRUE)
pBigMat_list_address=make_bigmemory_object_address_cpp(nrow=n_window,ncol=n_ind*2,file_name=paste0(bigmemory_data_name,"_haplotype_list"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
}
}

if(type%in%c(16)){
delete_bigmemory_file("numeric",bigmemory_data_name,bigmemory_data_path,TRUE);

pBigMat_num_address=make_bigmemory_object_address_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
}


if(type==13){
phased_vcf_to_haplotype_memory_cpp(input_data_vcf, pBigMat_hap_address, cpu_cores)
rm(pBigMat_hap_address);gc();
}


if(type==14){
haplotype_allele=phased_vcf_to_numeric_memory_cpp(pBigMat_hap_address,pBigMat_list_address,block_start,block_end,
									 numeric_file_name=paste0(bigmemory_data_name,"_numeric"),numeric_file_path=bigmemory_data_path,
									 input_data_vcf,cpu_cores)
rm(pBigMat_hap_address,pBigMat_list_address);gc();									 
}

if(type==15){	
tmp=phased_vcf_blupf90_memory_cpp(pBigMat_hap_address,pBigMat_list_address,block_start,block_end,input_data_vcf,cpu_cores)
haplotype_allele=tmp$haplotype_allele
data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=tmp$allele_string_vector,stringsAsFactors=F)									 
rm(tmp,pBigMat_hap_address,pBigMat_list_address);gc();	
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=T,sep=" ")}
}

if(type==16){unphased_vcf_to_numeric_memory_cpp(pBigMat_num_address,input_data_vcf,cpu_cores,miss_base,miss_base_num)								 
rm(pBigMat_num_address);gc();	
}
}
rm(input_data_vcf);gc();
#删除bigmemory中间文件

delete_bigmemory_file("haplotype_list",bigmemory_data_name,bigmemory_data_path,FALSE);

if(!13%in%output_type_number){

delete_bigmemory_file("haplotype_hap",bigmemory_data_name,bigmemory_data_path,TRUE);

}


if(13%in%output_type_number){
data_haplotype_hap=make_bigmemory_object_cpp(nrow=n_snp,ncol=n_ind*2,file_name=paste0(bigmemory_data_name,"_haplotype_hap"),
                                                               file_path=bigmemory_data_path,type=bigmemory_data_type)
data_haplotype_sample=data.frame(IND_name)
cat(paste0("bigmemory-Haplotype-hap data is save as ",bigmemory_data_path,"/",bigmemory_data_name,"_haplotype_hap \n"))
}

if(14%in%output_type_number|16%in%output_type_number){									 
data_numeric=make_bigmemory_object_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)	
options(bigmemory.allow.dimnames=TRUE)
rownames(data_numeric)=IND_name
colnames(data_numeric)=paste0("SNP",1:ncol(data_numeric))																
cat(paste0("bigmemory-Numeric data is save as ",bigmemory_data_path,"/",bigmemory_data_name,"_numeric \n"))
}

}
####################################input data type is Haplotype#################################### 
if(input_data_type=="Haplotype"){
if(!is.null(input_data_path)&!is.null(input_data_name)){
cat("Start read the Haplotype format genotype data \n")
input_data_haplotype_hap=fread(paste0(input_data_path,"/",input_data_name,".hap"),header=F,data.table=F)
input_data_haplotype_sample=fread(paste0(input_data_path,"/",input_data_name,".sample"),header=F,data.table=F)
input_data_haplotype_map=fread(paste0(input_data_path,"/",input_data_name,".map"),header=F,data.table=F)
input_data_haplotype_map[,1]=as.numeric(input_data_haplotype_map[,1])
input_data_haplotype_map[,3]=as.numeric(input_data_haplotype_map[,3])
cat("Complete read the Haplotype format genotype data \n")}

if(is.data.frame(input_data_haplotype_hap)){input_data_haplotype_hap=DataFrame_to_arma(input_data_haplotype_hap)}
if(typeof(input_data_haplotype_hap)!="integer"){input_data_haplotype_hap=NumericMatrix_to_arma(input_data_haplotype_hap)}

IND_name=as.character(input_data_haplotype_sample[,1])
#构建block
if(!is.null(phased_genotype==TRUE)&all(is.null(haplotype_window_nSNP),is.null(haplotype_window_kb),is.null(haplotype_window_block))){
stop("Once set phased_genotype==TRUE,  haplotype_window_nSNP or haplotype_window_kb or haplotype_window_block must be provided!")
}

if(!is.null(haplotype_window_nSNP)|!is.null(haplotype_window_kb)|!is.null(haplotype_window_block)){
block_result=get_haplotype_block(input_data_haplotype_map,haplotype_window_nSNP,
				                    haplotype_window_kb,haplotype_window_block)
block_hap=block_result$block_rcpp
block_start=block_hap[,1];block_end=block_hap[,2]
output_block=block_result$block
}else{
block_start=1:5
block_end=1:5
output_block=NULL
}

n_ind=ncol(input_data_haplotype_hap)/2
n_snp=nrow(input_data_haplotype_hap)
n_window=length(block_start)


output_type_number=get_output_type_number(input_data_type=input_data_type,
										  output_data_type=output_data_type,
										  bigmemory_cal=bigmemory_cal,
										  phased_genotype=phased_genotype)

if(length(output_type_number)==0){stop("Please provided standard output data type!")}

for(type in output_type_number){
if(exists(".Random.seed")){fixed_openMP(.Random.seed);set.seed(19980204)}
if(type==1){
tmp=haplotype_convertion(input_data_haplotype_hap,block_start,block_end,IND_name,type,cpu_cores)

haplotype_allele=tmp$SNP
data_numeric=tmp$numeric
rm(tmp);gc();
rownames(data_numeric)=IND_name
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_numeric),paste0(output_data_name,".numeric"),quote=F,row.names=T,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_numeric);gc();}
}

if(type==2){
tmp=haplotype_convertion(input_data_haplotype_hap,block_start,block_end,IND_name,type,cpu_cores)
haplotype_allele=tmp$SNP
data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=tmp$blupf90,stringsAsFactors=F)
rm(tmp);gc();
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=F,sep=" ")}
if(return_result==FALSE){rm(data_blupf90);gc();}
}

if(type==3){

tmp=haplotype_convertion(input_data_haplotype_hap,block_start,block_end,IND_name,type,cpu_cores)
haplotype_allele=tmp$SNP
data_numeric=tmp$numeric
rm(tmp);gc();
rownames(data_numeric)=IND_name
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_numeric),paste0(output_data_name,".numeric"),quote=F,row.names=T,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_numeric);gc();}
}

if(type==4){
tmp=haplotype_convertion(input_data_haplotype_hap,block_start,block_end,IND_name,type,cpu_cores)
haplotype_allele=tmp$SNP
data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=tmp$blupf90,stringsAsFactors=F)
rm(tmp);gc();
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=F,sep=" ")}
if(return_result==FALSE){rm(data_blupf90);gc();}
}

#构建bigmemory对象 
if(type%in%c(11,12,13,14)){

if(exists("pBigMat_hap")){rm(pBigMat_hap)}
gc()
delete_bigmemory_file("haplotype_hap",bigmemory_data_name,bigmemory_data_path,FALSE);
pBigMat_hap=make_bigmemory_object_cpp(nrow=n_snp,ncol=n_ind*2,file_name=paste0(bigmemory_data_name,"_haplotype_hap"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
pBigMat_hap[]=input_data_haplotype_hap
rm(input_data_haplotype_hap,pBigMat_hap);gc();
pBigMat_hap_address=make_bigmemory_object_address_cpp(nrow=n_snp,ncol=n_ind*2,file_name=paste0(bigmemory_data_name,"_haplotype_hap"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)																
}

if(type%in%c(11,12)){

if(exists("pBigMat_list_address")){rm(pBigMat_list_address)}
gc()
delete_bigmemory_file("haplotype_list",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("numeric",bigmemory_data_name,bigmemory_data_path,FALSE);

pBigMat_list_address=make_bigmemory_object_address_cpp(nrow=n_window,ncol=n_ind*2,file_name=paste0(bigmemory_data_name,"_haplotype_list"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
}

if(type%in%c(13)){

if(exists("pBigMat_num_address")){rm(pBigMat_num_address)}
gc()
delete_bigmemory_file("numeric",bigmemory_data_name,bigmemory_data_path,FALSE);

pBigMat_num_address=make_bigmemory_object_address_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
}


if(type==11){
haplotype_allele=phased_haplotype_to_numeric_memory_cpp(pBigMat_hap_address,pBigMat_list_address,
									   numeric_file_name=paste0(bigmemory_data_name,"_numeric"),
									   numeric_file_path=bigmemory_data_path,
                                             block_start,block_end,cpu_cores)

data_numeric=make_bigmemory_object_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
options(bigmemory.allow.dimnames=TRUE)
rownames(data_numeric)=IND_name
colnames(data_numeric)=paste0("SNP",1:ncol(data_numeric))																		
cat(paste0("bigmemory-Numeric data is save as ",bigmemory_data_path,"/",bigmemory_data_name,"_numeric \n"))
}


if(type==12){

tmp=phased_haplotype_to_blupf90_memory_cpp(block_start,block_end,pBigMat_hap_address,pBigMat_list_address,IND_name,cpu_cores)
haplotype_allele=tmp$SNP

data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=tmp$allele_string_vector,stringsAsFactors=F)								 
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=T,sep=" ")}
}

if(type==13){
unphased_haplotype_to_numeric_memory_cpp(pBigMat_hap_address,pBigMat_num_address,cpu_cores,miss_base,miss_base_num)

data_numeric=make_bigmemory_object_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
options(bigmemory.allow.dimnames=TRUE)
rownames(data_numeric)=IND_name
colnames(data_numeric)=paste0("SNP",1:ncol(data_numeric))																		
cat(paste0("bigmemory-Numeric data is save as ",bigmemory_data_path,"/",bigmemory_data_name,"_numeric \n"))
}

if(type==14){
data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=unphased_haplotype_to_blupf90_memory_cpp(pBigMat_hap_address,IND_name,cpu_cores),stringsAsFactors=F)								 
if(!is.null(output_data_path)&!is.null(output_data_name)){
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=T,sep=" ")}
}

}
}



####################################input data type is Numeric#################################### 
if(input_data_type=="Numeric"){
if(!is.null(input_data_path)&!is.null(input_data_name)){
cat("Start read the Numeric format genotype data \n")
input_data_numeric=fread(paste0(input_data_path,"/",input_data_name),data.table=F)
IND_name=as.character(input_data_numeric[,1])
input_data_numeric=as.matrix(input_data_numeric[,-1])
if(file.exists(paste0(input_data_path,"/",input_data_name,".map"))){
cat("Start read the Numeric_map data \n")
input_data_numeric_map=fread(paste0(input_data_path,"/",input_data_name,".map"),data.table=F)
}
cat("Complete read the Numeric format genotype data \n")}

output_type_number=get_output_type_number(input_data_type=input_data_type,
										  output_data_type=output_data_type,
										  bigmemory_cal=bigmemory_cal,
										  phased_genotype=phased_genotype)

if(length(output_type_number)==0){stop("Please provided standard output data type!")}

if(!is.null(input_data_numeric_map)){
input_data_numeric_map=as.matrix(input_data_numeric_map)
input_data_numeric_map[,1]=as.numeric(input_data_numeric_map[,1])
input_data_numeric_map[,3]=as.numeric(input_data_numeric_map[,3])
}


for(type in output_type_number){

if(type==1){data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=numeric_to_blupf90_cpp(input_data_numeric,IND_name,cpu_cores),stringsAsFactors=F)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=F,sep=" ")}
if(return_result==FALSE){rm(data_blupf90);gc();}
}

if(type==2){
if(is.null(input_data_numeric_map)){stop("Please provide Numeric map information!")}
data_ped=numeric_to_ped_cpp(IND_name,input_data_numeric_map,input_data_numeric,cpu_cores,miss_base)
data_map=input_data_numeric_map[,c(1,2,3,3)];data_map[,3]=as.numeric(data_map[,3])/1000000;
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_ped),paste0(output_data_name,".ped"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(data_map),paste0(output_data_name,".map"),quote=F,row.names=F,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_map,data_ped);gc();}

}


if(type==3){
if(is.null(input_data_numeric_map)){stop("Please provide Numeric map information!")}
input_data_numeric_map=as.matrix(input_data_numeric_map)
data_hmp=numeric_to_hapmap_cpp(IND_name,input_data_numeric_map,input_data_numeric,cpu_cores,miss_base)
data_hmp[,c(2,5:11)]="NA"
colnames(data_hmp)=c("rs#","alleles","chrom","pos","strand","assembly","center","protLSID","assayLSID","panelLSID","QCcode",IND_name)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.table(data_hmp),paste0(output_data_name,".hmp.txt"),quote=F,row.names=F,col.names=T,sep="\t")}
if(return_result==FALSE){rm(data_hmp);gc();}
}


if(type==4){
if(is.null(input_data_numeric_map)){stop("Please provide Numeric map information!")}
data_vcf=numeric_to_vcf_cpp(IND_name,input_data_numeric_map,input_data_numeric,phased_symbol,cpu_cores)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(get_vcf_header()),paste0(output_data_name,".vcf"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(data_vcf),paste0(output_data_name,".vcf"),append=TRUE,quote=F,row.names=F,col.names=T,sep="\t")}
if(return_result==FALSE){rm(data_vcf);gc();}
}

}
rm(input_data_numeric);gc();
}


####################################input data type is BLUPF90#################################### 
if(input_data_type=="BLUPF90"){
if(!is.null(input_data_path)&!is.null(input_data_name)){
cat("Start read the BLUPF90 format genotype data \n")
input_data_blupf90=fread(paste0(input_data_path,"/",input_data_name),data.table=F)
if(file.exists(paste0(input_data_path,"/",input_data_name,".map"))){
cat("Start read the BLUPF90_map data \n")
input_data_blupf90_map=fread(paste0(input_data_path,"/",input_data_name,".map"),data.table=F)
}
cat("Complete read the BLUPF90 format genotype data \n")}
IND_name=as.character(input_data_blupf90[,1]);input_data_blupf90=input_data_blupf90[,2];gc();
n_ind=length(IND_name)
n_snp=nchar(input_data_blupf90[1])
output_type_number=get_output_type_number(input_data_type=input_data_type,
										  output_data_type=output_data_type,
										  bigmemory_cal=bigmemory_cal,
										  phased_genotype=phased_genotype)

if(length(output_type_number)==0){stop("Please provided standard output data type!")}

if(!is.null(input_data_blupf90_map)){
input_data_blupf90_map=as.matrix(input_data_blupf90_map)
input_data_blupf90_map[,1]=as.numeric(input_data_blupf90_map[,1])
input_data_blupf90_map[,3]=as.numeric(input_data_blupf90_map[,3])
}


for(type in output_type_number){
if(exists(".Random.seed")){fixed_openMP(.Random.seed);set.seed(19980204)}

if(type==1){data_numeric=blupf90_to_numeric_cpp(input_data_blupf90,cpu_cores)
rownames(data_numeric)=IND_name
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_numeric),paste0(output_data_name,".numeric"),quote=F,row.names=T,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_numeric);gc();}
}


if(type==2){
if(is.null(input_data_blupf90_map)){stop("Please provide BLUPF90 map information!")}
data_ped=blupf90_to_ped_cpp(IND_name,input_data_blupf90_map,input_data_blupf90,cpu_cores,miss_base)
data_map=input_data_blupf90_map[,c(1,2,3,3)];data_map[,3]=as.numeric(data_map[,3])/1000000;
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_ped),paste0(output_data_name,".ped"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(data_map),paste0(output_data_name,".map"),quote=F,row.names=F,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_map,data_ped);gc();}

}


if(type==3){
if(is.null(input_data_blupf90_map)){stop("Please provide BLUPF90 map information!")}
data_hmp=blupf90_to_hapmap_cpp(IND_name,input_data_blupf90_map,input_data_blupf90,cpu_cores,miss_base)
data_hmp[,c(2,5:11)]="NA"
colnames(data_hmp)=c("rs#","alleles","chrom","pos","strand","assembly","center","protLSID","assayLSID","panelLSID","QCcode",IND_name)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.table(data_hmp),paste0(output_data_name,".hmp.txt"),quote=F,row.names=F,col.names=T,sep="\t")}
if(return_result==FALSE){rm(data_hmp);gc();}
}


if(type==4){
if(is.null(input_data_blupf90_map)){stop("Please provide BLUPF90 map information!")}
data_vcf=blupf90_to_vcf_cpp(IND_name,input_data_blupf90_map,input_data_blupf90,phased_symbol,cpu_cores)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(get_vcf_header()),paste0(output_data_name,".vcf"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(data_vcf),paste0(output_data_name,".vcf"),append=TRUE,quote=F,row.names=F,col.names=T,sep="\t")}
if(return_result==FALSE){rm(data_vcf);gc();}
}


#构建bigmemory对象 
if(type==11){

if(exists("pBigMat_num_address")){rm(pBigMat_num_address)}
gc()
delete_bigmemory_file("numeric",bigmemory_data_name,bigmemory_data_path,FALSE);

pBigMat_num_address=make_bigmemory_object_address_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
																
blupf90_to_numeric_memory_cpp(input_data_blupf90,pBigMat_num_address,cpu_cores,type=bigmemory_data_type)								 
data_numeric=make_bigmemory_object_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)	
options(bigmemory.allow.dimnames=TRUE)
rownames(data_numeric)=IND_name
colnames(data_numeric)=paste0("SNP",1:ncol(data_numeric))																
cat(paste0("bigmemory-Numeric data is save as ",bigmemory_data_path,"/",bigmemory_data_name,"_numeric \n"))
																
}


}
rm(input_data_blupf90);gc();
}

options (warn =1)
##############################################

if(!exists("data_hmp")){data_hmp=NULL}
if(!exists("data_ped")){data_ped=NULL}
if(!exists("data_map")){data_map=NULL}
if(!exists("data_blupf90")){data_blupf90=NULL}
if(!exists("data_numeric")){data_numeric=NULL}
if(!exists("data_vcf")){data_vcf=NULL}
if(!exists("data_haplotype_hap")){data_haplotype_hap=NULL}
if(!exists("data_haplotype_map")){data_haplotype_map=NULL}
if(!exists("data_haplotype_sample")){data_haplotype_sample=NULL}
if(!exists("output_block")){output_block=NULL}


if(!is.null(haplotype_allele)&!is.null(output_block)){

haplotype_allele_num=sapply(haplotype_allele,length)
block_pos_allele=output_block$block_pos
block_pos_allele$chr=input_data_haplotype_map[block_pos_allele[,2],1]
block_pos_allele$pos=paste0(input_data_haplotype_map[block_pos_allele[,2],3],"_",
						     input_data_haplotype_map[block_pos_allele[,3],3])

haplotype_allele_chr=rep(block_pos_allele$chr,times=haplotype_allele_num)
haplotype_allele_pos=rep(block_pos_allele$pos,times=haplotype_allele_num)

haplotype_allele=do.call(c,haplotype_allele)
haplotype_allele=data.frame(Chr=haplotype_allele_chr,
							SNP=paste0(1:length(haplotype_allele)),
							cM=0,
							Pos=haplotype_allele_pos,
							Allele=haplotype_allele)
output_block=c(output_block,list(allele=haplotype_allele))							
}

if(return_result==TRUE){return(list(hmp=data_hmp,
                                            ped=data_ped,
								     map=data_map,
									 blupf90=data_blupf90,
									 numeric=data_numeric,
									 vcf=data_vcf,
									 haplotype_hap=data_haplotype_hap,
									 haplotype_map=data_haplotype_map,
									 haplotype_sample=data_haplotype_sample,
									 phased_block=output_block))}
} 


O0OOO0O0OOO0O0O0OOOO<-function(input_data_type=NULL,input_data_hmp=NULL,input_data_plink_ped=NULL,input_data_plink_map=NULL,
							  input_data_numeric=NULL,input_data_blupf90=NULL,input_data_haplotype_hap=NULL,
							  input_data_haplotype_map=NULL,input_data_haplotype_sample=NULL,input_data_vcf=NULL,miss_base=NULL){

if(is.null(input_data_type)){

if(!is.null(input_data_hmp)){input_data_type=c(input_data_type,"Hapmap")}
if(!is.null(input_data_plink_map)&!is.null(input_data_plink_ped)){input_data_type=c(input_data_type,"Plink")}
if(!is.null(input_data_blupf90)){input_data_type=c(input_data_type,"BLUPF90")}
if(!is.null(input_data_numeric)){input_data_type=c(input_data_type,"Numeric")}
if(!is.null(input_data_vcf)){input_data_type=c(input_data_type,"VCF")}
if(!is.null(input_data_haplotype_hap)&!is.null(input_data_haplotype_sample)&!is.null(input_data_haplotype_map)){input_data_type=c(input_data_type,"Haplotype")}
input_data_type=unique(input_data_type)

if(length(input_data_type)>=2){
stop("Only one kind of input data is allowed in this analysis!Please check your input data carefully! \n")
}else if(length(input_data_type)==0){
stop("No input data avaliable in this analysis \n")
}
}

if(is.null(miss_base)){
miss_base=ifelse(input_data_type=="Plink","0",ifelse(input_data_type=="Hapmap","N","N"))
}
return(list(input_data_type=input_data_type,miss_base=miss_base))
}


get_vcf_header<-function(date=Sys.time()){
H1="##fileformat=VCFv4.2"
H2=paste0("##filedate=",date)
H3="##source=R-package:blupADC"
H4="##INFO=<ID=AF,Number=A,Type=Float"
H5="##INFO=<ID=DR2,Number=A,Type=Float"
H6="##INFO=<ID=IMP,Number=0,Type=Flag"
H7="##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype"
H8="##FORMAT=<ID=DS,Number=A,Type=Float"
return(rbind(H1,H2,H3,H4,H5,H6,H7,H8))
}


get_output_type_number<-function(input_data_type=NULL,
								 output_data_type=NULL,
								 bigmemory_cal=FALSE,
								 phased_genotype=FALSE)
								 {
data_type_set=c("Hapmap","Plink","VCF","Haplotype","Numeric","BLUPF90")
output_data_type=unique(output_data_type)
input_data_type=unique(input_data_type)

if(length(unique(input_data_type))!=1){stop("Provided input_data_type  more than two or equal zero!")}
if(sum(!input_data_type%in%data_type_set)>=1){stop("Provided input_data_type contain non-standard type!")}
if(sum(!output_data_type%in%data_type_set)>=1){stop("Provided output_data_type contain non-standard type!")}

output_data_type=output_data_type[!output_data_type%in%input_data_type]

type_number=NULL

if(input_data_type=="Plink"){

   if(bigmemory_cal==FALSE){
		if("Hapmap"%in%output_data_type){type_number=c(type_number,1)}
		if("VCF"%in%output_data_type){type_number=c(type_number,2)}	
		if("Numeric"%in%output_data_type){type_number=c(type_number,3)}
		if("BLUPF90"%in%output_data_type){type_number=c(type_number,4)}   
   }else{
		
		if("Hapmap"%in%output_data_type){type_number=c(type_number,1)}
		if("VCF"%in%output_data_type){type_number=c(type_number,2)}	
		if("Numeric"%in%output_data_type){type_number=c(type_number,13)}
		if("BLUPF90"%in%output_data_type){type_number=c(type_number,4)}    		
   }

}else if(input_data_type=="Hapmap"){

   if(bigmemory_cal==FALSE){
		if("Plink"%in%output_data_type){type_number=c(type_number,1)}
		if("VCF"%in%output_data_type){type_number=c(type_number,2)}	
		if("Numeric"%in%output_data_type){type_number=c(type_number,3)}
		if("BLUPF90"%in%output_data_type){type_number=c(type_number,4)}   
   }else{
		if("Plink"%in%output_data_type){type_number=c(type_number,1)}
		if("VCF"%in%output_data_type){type_number=c(type_number,2)}	
		if("Numeric"%in%output_data_type){type_number=c(type_number,13)}
		if("BLUPF90"%in%output_data_type){type_number=c(type_number,4)}   		
   }

}else if(input_data_type=="VCF"){

   if(bigmemory_cal==FALSE){
		if("Hapmap"%in%output_data_type){type_number=c(type_number,1)}
		if("Haplotype"%in%output_data_type){type_number=c(type_number,3)}	
		if("Numeric"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,4)}
		if("BLUPF90"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,5)}
		if("Numeric"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,6)}
		if("BLUPF90"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,7)} 
		if("Plink"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,8)} 
		if("Plink"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,2)} 
   }else{		
		if("Hapmap"%in%output_data_type){type_number=c(type_number,1)}
		if("Plink"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,8)} 
		if("Plink"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,2)} 
		if("Haplotype"%in%output_data_type){type_number=c(type_number,13)}	
		if("Numeric"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,14)}
		if("BLUPF90"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,15)}
		if("Numeric"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,16)}
		if("BLUPF90"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,7)}  
   }

}else if(input_data_type=="Haplotype"){

   if(bigmemory_cal==FALSE){
		if("Numeric"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,1)}  
		if("BLUPF90"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,2)} 			
		if("Numeric"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,3)}		
		if("BLUPF90"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,4)} 
	
   }else{ 
            
		if("Numeric"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,11)}  
		if("BLUPF90"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,12)} 			
		if("Numeric"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,13)}		
		if("BLUPF90"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,14)} 
   
   }
   
}else if(input_data_type=="Numeric"){

   if(bigmemory_cal==FALSE){
		if("BLUPF90"%in%output_data_type){type_number=c(type_number,1)}
		if("Plink"%in%output_data_type){type_number=c(type_number,2)}
		if("Hapmap"%in%output_data_type){type_number=c(type_number,3)}
		if("VCF"%in%output_data_type){type_number=c(type_number,4)}		
	
   }else{ 
         if("BLUPF90"%in%output_data_type){type_number=c(type_number,1)}
		if("Plink"%in%output_data_type){type_number=c(type_number,2)}
		if("Hapmap"%in%output_data_type){type_number=c(type_number,3)}
		if("VCF"%in%output_data_type){type_number=c(type_number,4)}		  
   
   }
   
}else if(input_data_type=="BLUPF90"){

   if(bigmemory_cal==FALSE){
		if("Numeric"%in%output_data_type){type_number=c(type_number,1)}
		if("Plink"%in%output_data_type){type_number=c(type_number,2)}
		if("Hapmap"%in%output_data_type){type_number=c(type_number,3)}
		if("VCF"%in%output_data_type){type_number=c(type_number,4)}
	
   }else{ 
         if("Numeric"%in%output_data_type){type_number=c(type_number,11)} 
		if("Plink"%in%output_data_type){type_number=c(type_number,2)}
		if("Hapmap"%in%output_data_type){type_number=c(type_number,3)}
		if("VCF"%in%output_data_type){type_number=c(type_number,4)}		 
   
   }
}
return(type_number)
}








#获取单倍型 block
get_haplotype_block<-function(data_map=NULL,
							  haplotype_window_nSNP=NULL,
							  haplotype_window_kb=NULL,
							  haplotype_window_block=NULL){

#同时提供三种参数中至少两种会报错
if(sum(c(!is.null(haplotype_window_nSNP),!is.null(haplotype_window_kb),!is.null(haplotype_window_block)))>=2){stop("Provided too much parateters: haplotype_window_nSNP, haplotype_window_kb, haplotype_window_block!")}
if(sum(c(!is.null(haplotype_window_nSNP),!is.null(haplotype_window_kb),!is.null(haplotype_window_block)))==0){stop("Provided no parateters: haplotype_window_nSNP, haplotype_window_kb, haplotype_window_block!")}

#获取SNP位置-考虑过了染色体
data_map=data.frame(data_map,stringsAsFactors=F)
data_map[,1]=as.numeric(data_map[,1])
data_map[,3]=as.numeric(data_map[,3])
chr_set=unique(data_map[,1])
chr_set_snp_n=0
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(chr_set)){
tmp_data_map=data_map[data_map[,1]==chr_set[OOO0OOO0O0OOO0O0OOO0OOOO],]
chr_set_snp_n=c(chr_set_snp_n,nrow(tmp_data_map))
}
chr_set_snp_n=cumsum(chr_set_snp_n)
chr_set_snp_n=chr_set_snp_n[-length(chr_set_snp_n)]

if(!is.null(haplotype_window_block)){

block_start=haplotype_window_block[,1]-1  #window_block为两列数据，均为位置信息
block_end=haplotype_window_block[,2]-1

}else if(!is.null(haplotype_window_nSNP)){ #根据SNP数目划分block

block_start=NULL
block_end=NULL
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(chr_set)){
tmp_data_map=data_map[data_map[,1]==chr_set[OOO0OOO0O0OOO0O0OOO0OOOO],]
block_start_tmp=seq(1,nrow(tmp_data_map),haplotype_window_nSNP)-1
block_end_tmp=seq(1,nrow(tmp_data_map),haplotype_window_nSNP)-2
block_end_tmp=c(block_end_tmp,nrow(tmp_data_map)-1)
block_end_tmp=block_end_tmp[-1]

#添加位置信息
block_start_tmp=block_start_tmp+chr_set_snp_n[OOO0OOO0O0OOO0O0OOO0OOOO]
block_end_tmp=block_end_tmp+chr_set_snp_n[OOO0OOO0O0OOO0O0OOO0OOOO]
block_start=c(block_start,block_start_tmp)
block_end=c(block_end,block_end_tmp)
}

}else if(!is.null(haplotype_window_kb)){  #根据物理距离划分block

block_start=NULL
block_end=NULL
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(chr_set)){
tmp_data_map=data_map[data_map[,1]==chr_set[OOO0OOO0O0OOO0O0OOO0OOOO],]
block=seq(min(tmp_data_map[,3]),max(tmp_data_map[,3]),haplotype_window_kb*1000)

block_1=block
block_2=block-1;
block_2=block_2[-1]; #去除最后一列
block_2=c(block_2,max(tmp_data_map[,3]))

#R-function too slow
#block_start_tmp=NULL
#block_end_tmp=NULL
#for(OOO0OOO0O0OOO0O0OOO0OOOOO in 1:length(block)){
#block_start_tmp=c(block_start_tmp,min(which(tmp_data_map[,3]>=block_1[OOO0OOO0O0OOO0O0OOO0OOOOO])))
#block_end_tmp=c(block_end_tmp,max(which(tmp_data_map[,3]<=block_2[OOO0OOO0O0OOO0O0OOO0OOOOO])))
#}

#Rcpp function 
block_result=define_block_window_kb_cpp(block_1,block_2,tmp_data_map[,3])
block_start_tmp=block_result[[1]]
block_end_tmp=block_result[[2]]

pos_status=block_start_tmp<=block_end_tmp
block_start_tmp=block_start_tmp[pos_status]
block_end_tmp=block_end_tmp[pos_status]
block_start_tmp=block_start_tmp-1
block_end_tmp=block_end_tmp-1
#添加位置信息
block_start_tmp=block_start_tmp+chr_set_snp_n[OOO0OOO0O0OOO0O0OOO0OOOO]
block_end_tmp=block_end_tmp+chr_set_snp_n[OOO0OOO0O0OOO0O0OOO0OOOO]
block_start=c(block_start,block_start_tmp)
block_end=c(block_end,block_end_tmp)
}
}
#进行单倍型构建分析
#写出block的信息
snp_map=as.character(data_map[,2])
block_pos=data.frame(block_num=paste0("block",1:length(block_start)),window_start_pos=block_start+1,window_end_pos=block_end+1)
block_snp=data.frame(block_num=paste0("block",1:length(block_start)),window_start_snp=snp_map[block_start+1],window_end_snp=snp_map[block_end+1])
return(list(block=list(block_pos=block_pos,block_snp=block_snp),
			block_rcpp=cbind(block_start,block_end) #c++程序的输入文件
			))
}

fixed_openMP<-function(x){
  Sx <- deparse(substitute(x))
  rm(list=Sx,envir=sys.frame(-1))
}

get_bigmemory_address<-function(data){

	return(data@address)
}

delete_bigmemory_file<-function(matrix_type,
								    bigmemory_data_name,
						             bigmemory_data_path,
							         message=FALSE){
								   
		
		file1=paste0(bigmemory_data_path,"/",bigmemory_data_name,"_",matrix_type,".desc")
		file2=paste0(bigmemory_data_path,"/",bigmemory_data_name,"_",matrix_type,".bk")

	    if(file.exists(file1)&file.exists(file2)){
			if(message==TRUE){
			cat(paste0("Found ",matrix_type,".desc & ",matrix_type,".bk files in path:", bigmemory_data_path,"\n"))
			cat("Software will delete these files automatically!\n")
			}
			unlink(file2)
			unlink(file1)
			
		}
	    if(file.exists(file1)){unlink(file1)}
		 if(file.exists(file2)){unlink(file2)}

		
}
#
trace_pedigree<-function(input_pedigree=NULL,  #第一列为个体号,第二列为父亲,第三列为母亲
						 input_pedigree_path=NULL,
						 input_pedigree_name=NULL,
						 multi_col=FALSE,#不对系谱数据进行重新整理
						 dup_error_check=TRUE,
						 sex_error_check=TRUE,
						 breed_error_check=FALSE,
						 birth_date_error_check=FALSE,
						 trace_id=NULL,   #默认追溯所有系谱中的个体
						 trace_sibs=FALSE,#是否追溯同胞
						 trace_fullsibs=FALSE, #追溯全同胞
						 max_fullsibs=NULL,  #追溯全同胞时，全同胞最大数量
						 trace_direction="forward", #reverse 为反向
						 #trace_offspring=FALSE, #追溯个体的子代
						 trace_generation=NULL, #追溯的代数	
						 trace_birth_date=NULL, #仅追溯出生日期在给定日期之后的个体
						 output_rename_pedigree=TRUE, #输出rename后的系谱
						 output_pedigree_path=NULL,
						 output_pedigree_name=NULL,
						 output_pedigree_tree=FALSE,
						 pedigree_tree_depth=3,
						 summary_sibs=FALSE
						 ){

options(scipen=200)

library(data.table)


#检查输入的参数
if(sum(is.null(output_pedigree_path),is.null(output_pedigree_name))==1){
stop("Paramters:output_pedigree_path and output_pedigree_name couldn't be: the one is non-NULL and the other is NULL")
}


#检查输入的系谱
if(is.null(input_pedigree)){
OOO0OOO0O0OOO0O0OOO0OOOOOO=data.table::fread(paste0(input_pedigree_path,"/",input_pedigree_name),data.table=F)
}else{
OOO0OOO0O0OOO0O0OOO0OOOOOO=input_pedigree
}

if(ncol(OOO0OOO0O0OOO0O0OOO0OOOOOO)==3){	
cat("Peidgree provided has three columns,please make sure the format of pedigree_data has three columns: Offspring Sire Dam  \n")
}else if(ncol(OOO0OOO0O0OOO0O0OOO0OOOOOO)==4){
cat("Peidgree provided has four columns,please make sure the format of pedigree_data has four columns: Offspring Sire Dam  Birth_Date \n")
}else if(multi_col==TRUE){
cat("Peidgree provided has multiple columns,please make sure the format of pedigree_data similar to: Offspring Sire Dam  SireSire SireDam ......  \n")
OOO0OOO0O0OOO0O0OOO0OOOOOO=O0OOOOO0OOO0O0OOOOO0(OOO0OOO0O0OOO0O0OOO0OOOOOO)
}else{
stop("Error:peidgree_provided is not standard format!")
}

#将缺失值变为NA
OOO0OOO0O0OOO0O0OOO0OOOOOO[OOO0OOO0O0OOO0O0OOO0OOOOOO=="-9999"]=NA
OOO0OOO0O0OOO0O0OOO0OOOOOO[OOO0OOO0O0OOO0O0OOO0OOOOOO=="0"]=NA
OOO0OOO0O0OOO0O0OOO0OOOOOO[OOO0OOO0O0OOO0O0OOO0OOOOOO==""]=NA

##将系谱变成长系谱，确保第一列包含所有个体
#id_F=na.omit(unique(OOO0OOO0O0OOO0O0OOO0OOOOOO[!OOO0OOO0O0OOO0O0OOO0OOOOOO[,2]%in%OOO0OOO0O0OOO0O0OOO0OOOOOO[,1],2]))
#id_M=na.omit(unique(OOO0OOO0O0OOO0O0OOO0OOOOOO[!OOO0OOO0O0OOO0O0OOO0OOOOOO[,3]%in%OOO0OOO0O0OOO0O0OOO0OOOOOO[,1],3]))
#if(length(id_F)==0){id_F=NA}
#if(length(id_M)==0){id_M=NA}
#
#OOO0OOO0O0OOO0O0OOO0OOOOOO=rbind(as.matrix(OOO0OOO0O0OOO0O0OOO0OOOOOO[,1:3]),cbind(id_F,NA,NA),cbind(id_M,NA,NA))
#
#if(ncol(input_pedigree)==4){
#OOO0OOO0O0OOO0O0OOO0OOOOOO=cbind(OOO0OOO0O0OOO0O0OOO0OOOOOO,input_pedigree[match(OOO0OOO0O0OOO0O0OOO0OOOOOO[,1],input_pedigree[,1]),4])
#colnames(OOO0OOO0O0OOO0O0OOO0OOOOOO)[4]="Birth_Date"}

ped=as.matrix(OOO0OOO0O0OOO0O0OOO0OOOOOO)
if(!mode(ped)=="character"){ped=apply(ped,2,as.character)}
colnames(ped)[1:3]=c("Offspring","Sire","Dam")
ped=ped[!is.na(ped[,1]),]
#检查系谱错误
O0O0OOOOOOOOOOO0OOO0=OOOOOOO0O0O0O0OOOOO0(ped=ped, dup_error_check=dup_error_check, sex_error_check=sex_error_check,
						                        breed_error_check=breed_error_check, birth_date_error_check=birth_date_error_check)

ped=O0O0OOOOOOOOOOO0OOO0$ped

error_duplicated_id=O0O0OOOOOOOOOOO0OOO0$error_duplicated_id
error_sex_id=O0O0OOOOOOOOOOO0OOO0$error_sex_id
error_breed_id=O0O0OOOOOOOOOOO0OOO0$error_breed_id
error_birth_date_id=O0O0OOOOOOOOOOO0OOO0$error_birth_date_id

rename_ped=NULL;rename_phenotype=NULL
#系谱排序, 数字化
if(!is.null(trace_birth_date)){
cat("Please make sure_pedigree data has four columns \n")
trace_id=na.omit(ped[as.numeric(ped[,4])>=trace_birth_date,1])
}


#trace 选择的个体
#之前的trace_id思路，感觉有点错误
#if(!is.null(trace_id)){                     #根据给定的个体id, trace所有和这些个体id有联系的个体的系谱(默认trace这些个体以及它们的父母)
#
#trace_id=as.character(trace_id)       
#trace_id_parent=as.character(match_parents(trace_id,ped)[[1]])
#
#if(trace_offspring==TRUE){                 #根据给定的个体id, trace这些个体的后代
#O0OOOOOOOOOOO0O0O0O0=ped[ped[,2]%in%trace_id|ped[,3]%in%trace_id,1]
#trace_id=unique(c(trace_id,O0OOOOOOOOOOO0O0O0O0))
#}
#
#trace_id=unique(na.omit(c(trace_id,trace_id_parent))) #追溯给定个体的父母
#ped=ped[ped[,1]%in%trace_id,]
#}
ped[ped==0]=NA
#指定个体追溯系谱
if(!is.null(trace_id)){                     #根据给定的个体id, trace所有和这些个体id有联系的个体的系谱(默认trace这些个体以及它们的父母)
cat("Tracing porvided id......\n")
OOO0OOO0O0OOO0O0OOO0OOOO=0
trace_type="Match"
trace_id_set=trace_id

trace_id_sibs=NULL
#找寻 trace_id 这一代的同胞
if(trace_fullsibs==TRUE){
trace_id_sibs=match_fullsibs(trace_id,ped,max_fullsibs)
}


while(trace_type=="Match"){
OOO0OOO0O0OOO0O0OOO0OOOO=OOO0OOO0O0OOO0O0OOO0OOOO+1
if(!is.null(trace_generation)){
if(OOO0OOO0O0OOO0O0OOO0OOOO>=trace_generation)break
}

if(trace_fullsibs==TRUE){
sibs=match_fullsibs(trace_id,ped,max_fullsibs)
trace_id=c(trace_id,sibs)
}

trace_result=match_parents(trace_id,ped)
trace_id=na.omit(unique(as.character(trace_result[[1]])))

trace_type=trace_result[[2]]
trace_id_set=unique(c(trace_id_set,trace_id))
}
trace_id_set=unique(na.omit(trace_id_set))
trace_id_set=unique(c(trace_id_set,trace_id_sibs))
ped=ped[ped[,1]%in%trace_id_set,]
}

#通过for 循环判断哪些个体的代数

rename_pedigree=ped

base_result=single_pedigree_cpp(ped)

IND_base=base_result[[1]]
IND_middle=base_result[[2]]
IND_offspring=base_result[[3]]

if(trace_direction=="reverse"){

OOOOO0OOOOOOOOO0OOOO=data.frame(Offspring=c(IND_base,IND_middle,IND_offspring),Generation="NA",stringsAsFactors=F)
OOOOO0OOOOOOOOO0OOOO[OOOOO0OOOOOOOOO0OOOO[,1]%in%IND_base,2]=paste0("A",sprintf("%5s",0))
OOOOO0OOOOOOOOO0OOOO[OOOOO0OOOOOOOOO0OOOO[,1]%in%IND_offspring,2]=paste0("B",sprintf("%5s",1001))

#接下来的目的是将 IND_middle 中的个体划分出  父母与子代
tmp_ped=ped
OOO0OOO0O0OOO0O0OOO0OOOO=0
IND_base_set=IND_base
while(length(IND_middle)>0){

if(OOO0OOO0O0OOO0O0OOO0OOOO>1000){stop("Found generations of provided OOO0OOO0O0OOO0O0OOO0OOOOOO larger than 1000, please check your OOO0OOO0O0OOO0O0OOO0OOOOOO carefully!")}
OOO0OOO0O0OOO0O0OOO0OOOO=OOO0OOO0O0OOO0O0OOO0OOOO+1
tmp_ped=tmp_ped[tmp_ped[,1]%in%IND_middle,]
tmp_base_result=single_pedigree(tmp_ped)
tmp_IND_base=setdiff(tmp_base_result[[1]],IND_base_set);IND_base_set=c(IND_base_set,tmp_IND_base)
tmp_IND_middle=tmp_base_result[[2]]
tmp_IND_offspring=tmp_base_result[[3]]

OOOOO0OOOOOOOOO0OOOO[OOOOO0OOOOOOOOO0OOOO[,1]%in%tmp_IND_base,2]=paste0("A",sprintf("%5s",OOO0OOO0O0OOO0O0OOO0OOOO))
OOOOO0OOOOOOOOO0OOOO[OOOOO0OOOOOOOOO0OOOO[,1]%in%tmp_IND_offspring,2]=paste0("B",sprintf("%5s",1000-OOO0OOO0O0OOO0O0OOO0OOOO))

IND_middle=tmp_IND_middle
}

OOOOO0OOOOOOOOO0OOOO$Generation=as.numeric(as.factor(OOOOO0OOOOOOOOO0OOOO$Generation))

OOOOO0OOOOOOOOO0OOOO$Generation=OOOOO0OOOOOOOOO0OOOO$Generation-1
}else if(trace_direction=="forward"){

#OOOOO0OOOOOOOOO0OOOO=data.frame(Offspring=c(IND_base,IND_middle,IND_offspring),Generation="NA",stringsAsFactors=F)
#tmp_ped=ped
#OOO0OOO0O0OOO0O0OOO0OOOO=0
#IND_base_set=IND_base
#OOOOO0OOOOOOOOO0OOOO[OOOOO0OOOOOOOOO0OOOO[,1]%in%IND_base_set,2]=OOO0OOO0O0OOO0O0OOO0OOOO
#while(length(IND_base_set)>0){
#OOO0OOO0O0OOO0O0OOO0OOOO=OOO0OOO0O0OOO0O0OOO0OOOO+1
#IND_base_set=unique(na.omit(tmp_ped[tmp_ped[,2]%in%IND_base_set|tmp_ped[,3]%in%IND_base_set,1]))
#OOOOO0OOOOOOOOO0OOOO[OOOOO0OOOOOOOOO0OOOO[,1]%in%IND_base_set,2]=OOO0OOO0O0OOO0O0OOO0OOOO
#if(OOO0OOO0O0OOO0O0OOO0OOOO>1000){stop("Found generations of provided OOO0OOO0O0OOO0O0OOO0OOOOOO larger than 1000, please check your OOO0OOO0O0OOO0O0OOO0OOOOOO carefully!")}
#}

OOOOO0OOOOOOOOO0OOOO=data.frame(Offspring=c(IND_base,IND_middle,IND_offspring),stringsAsFactors=F)
OOOOO0OOOOOOOOO0OOOO$Sire=ped[match(OOOOO0OOOOOOOOO0OOOO[,1],ped[,1]),2]
OOOOO0OOOOOOOOO0OOOO$Dam=ped[match(OOOOO0OOOOOOOOO0OOOO[,1],ped[,1]),3]
IND_base_set=IND_base
pos=get_offspring_generation_cpp(OOOOO0OOOOOOOOO0OOOO,IND_base_set)
OOOOO0OOOOOOOOO0OOOO$Generation=pos
OOOOO0OOOOOOOOO0OOOO=OOOOO0OOOOOOOOO0OOOO[,c(1,4)]


}else{
stop("Non-standard input of trace_direction parameter")
}


OOOOO0OOOOOOOOO0OOOO=OOOOO0OOOOOOOOO0OOOO[!duplicated(OOOOO0OOOOOOOOO0OOOO[,1]),]
OOOOO0OOOOOOOOO0OOOO=OOOOO0OOOOOOOOO0OOOO[!OOOOO0OOOOOOOOO0OOOO[,1]%in%NA,]
options (warn =-1)	
if(NA%in%as.numeric(OOOOO0OOOOOOOOO0OOOO[,1])){
OOOOO0OOOOOOOOO0OOOO=OOOOO0OOOOOOOOO0OOOO[order(OOOOO0OOOOOOOOO0OOOO[,2],OOOOO0OOOOOOOOO0OOOO[,1]),]
}else{
OOOOO0OOOOOOOOO0OOOO=OOOOO0OOOOOOOOO0OOOO[order(OOOOO0OOOOOOOOO0OOOO[,2],as.numeric(OOOOO0OOOOOOOOO0OOOO[,1])),]
}
options (warn =1)	

if(!is.null(trace_generation)&is.null(trace_id)){
trace_generation=as.numeric(trace_generation)
O0O0OOO0O0O0O0OOO0OO=max(as.numeric(OOOOO0OOOOOOOOO0OOOO[,2]))
OOO0O0OOOOO0OOOOOOOO=O0O0OOO0O0O0O0OOO0OO-trace_generation+1

if(OOO0O0OOOOO0OOOOOOOO< 0){
message(paste0("Attention:the max generation of this_pedigree is ",O0O0OOO0O0O0O0OOO0OO))
OOO0O0OOOOO0OOOOOOOO=0
}

OOOOO0OOOOOOOOO0OOOO=OOOOO0OOOOOOOOO0OOOO[OOOOO0OOOOOOOOO0OOOO[,2]%in%c(OOO0O0OOOOO0OOOOOOOO:O0O0OOO0O0O0O0OOO0OO),]
OOOOO0OOOOOOOOO0OOOO[,2]=OOOOO0OOOOOOOOO0OOOO[,2]-OOO0O0OOOOO0OOOOOOOO
}

OOOOO0OOOOOOOOO0OOOO=cbind(OOOOO0OOOOOOOOO0OOOO,1:nrow(OOOOO0OOOOOOOOO0OOOO));
colnames(OOOOO0OOOOOOOOO0OOOO)[3]="Offspring_Id"


for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:ncol(rename_pedigree)){

rename_pedigree[,OOO0OOO0O0OOO0O0OOO0OOOO]=OOOOO0OOOOOOOOO0OOOO[match(rename_pedigree[,OOO0OOO0O0OOO0O0OOO0OOOO],OOOOO0OOOOOOOOO0OOOO[,"Offspring"]),"Offspring_Id"]

}



OOOOO0OOOOOOOOO0OOOO=cbind(OOOOO0OOOOOOOOO0OOOO,match_parents(as.character(as.numeric(OOOOO0OOOOOOOOO0OOOO[,"Offspring_Id"])),rename_pedigree)[[1]])
OOOOO0OOOOOOOOO0OOOO=cbind(OOOOO0OOOOOOOOO0OOOO,1:nrow(OOOOO0OOOOOOOOO0OOOO));
colnames(OOOOO0OOOOOOOOO0OOOO)[4:6]=c("Sire_Id","Dam_Id","Order")
OOOOO0OOOOOOOOO0OOOO[,4]=as.numeric(as.character(OOOOO0OOOOOOOOO0OOOO[,4]));
OOOOO0OOOOOOOOO0OOOO[,5]=as.numeric(as.character(OOOOO0OOOOOOOOO0OOOO[,5]))


O0OOOOO0O0OOO0OOO0O0=OOOOO0OOOOOOOOO0OOOO[,c(1,4,5,2)];O0OOOOO0O0OOO0OOO0O0[,2]=OOOOO0OOOOOOOOO0OOOO[match(O0OOOOO0O0OOO0OOO0O0[,2],OOOOO0OOOOOOOOO0OOOO[,3]),1]
O0OOOOO0O0OOO0OOO0O0[,3]=OOOOO0OOOOOOOOO0OOOO[match(O0OOOOO0O0OOO0OOO0O0[,3],OOOOO0OOOOOOOOO0OOOO[,3]),1]
O0OOOOO0O0OOO0OOO0O0[is.na(O0OOOOO0O0OOO0OOO0O0)]=0
colnames(O0OOOOO0O0OOO0OOO0O0)=c("Offspring","Sire","Dam","Generation")

if(!is.null(output_pedigree_path)&!is.null(output_pedigree_name)){
setwd(output_pedigree_path)


if(output_rename_pedigree==TRUE){
OOOOO0OOOOOOOOO0OOOO[is.na(OOOOO0OOOOOOOOO0OOOO)]=0
data.table::fwrite(data.table(OOOOO0OOOOOOOOO0OOOO[,-c(1:2)]),paste0(output_pedigree_name,".txt"),col.names=F,quote=F,sep=" ",row.names=F)
tmp=OOOOO0OOOOOOOOO0OOOO[,c(1,3)]
colnames(tmp)=c("Original_Id","Renamed_Id")
data.table::fwrite(data.table(tmp),paste0(output_pedigree_name,"_renamed_key.txt"),col.names=T,quote=F,sep=" ",row.names=F)
}else{
cat("Output ordered but non-renamed OOO0OOO0O0OOO0O0OOO0OOOOOO! \n")
tmp=O0OOOOO0O0OOO0OOO0O0
tmp[is.na(tmp)]=0
tmp[,4]=1:nrow(tmp)
data.table::fwrite(data.table(tmp),paste0(output_pedigree_name,".txt"),col.names=F,quote=F,sep=" ",row.names=F)
}
}
OOOOO0OOOOOOOOO0OOOO[is.na(OOOOO0OOOOOOOOO0OOOO)]=0

if(sum((OOOOO0OOOOOOOOO0OOOO$Sire_Id-OOOOO0OOOOOOOOO0OOOO$Offspring_Id)>0)>=1 |sum((OOOOO0OOOOOOOOO0OOOO$Dam_Id-OOOOO0OOOOOOOOO0OOOO$Offspring_Id)>0)>=1){
stop("Rename error: sire_id or dam_id larger than offspring_id \n")
}


#追溯全系谱
#将系谱 rename排序后，
#依次从上到下，找个体的祖先，
#利用迭代的思想，因为系谱在上面的个体是系谱在下面个体的祖先，如果上面个体的祖先已知，那么下面个体的祖先就也是已知的。
OOOOOOOOO0OOO0OOOOOO=NULL
if(output_pedigree_tree==TRUE){
cat("Constructing pedigree_tree......\n")

#original_ped=OOOOO0OOOOOOOOO0OOOO
#original_ped[,4]=original_ped[match(original_ped[,4],original_ped[,3]),1]
#original_ped[,5]=original_ped[match(original_ped[,5],original_ped[,3]),1]
#original_ped[is.na(original_ped)]="0"
#original_ped=as.matrix(original_ped[,3:5])
#original_ped=apply(original_ped,2,as.character)

generation_names=O0OOO0O0OOO0OOO0O0OO_colnames(generation=pedigree_tree_depth)

O0OOOOO0O0OOO0OOO0O0[is.na(O0OOOOO0O0OOO0OOO0O0)]="0"
O0OOOOO0O0OOO0OOO0O0=as.matrix(O0OOOOO0O0OOO0OOO0O0)
O0OOOOO0O0OOO0OOO0O0=apply(O0OOOOO0O0OOO0OOO0O0,2,as.character)


OOOOOOOOO0OOO0OOOOOO=full_generation_conversion(generation_names=generation_names,ped=O0OOOOO0O0OOO0OOO0O0[,1:3])

colnames(OOOOOOOOO0OOO0OOOOOO)=c("Offspring",generation_names)

}

#统计同胞数目
group_full_sibs=NULL
family_size_full_sibs=NULL
if(summary_sibs==TRUE){

#统计全同胞数目
full_sibs_pedigree=data.frame(O0OOOOO0O0OOO0OOO0O0,stringsAsFactors=F)
full_sibs_pedigree$MF=paste0(full_sibs_pedigree[,2],"_",full_sibs_pedigree[,3])
full_sibs_pedigree=full_sibs_pedigree[!(full_sibs_pedigree[,2]%in%0|full_sibs_pedigree[,3]%in%0),]
full_sibs_com=aggregate(full_sibs_pedigree,by=list(full_sibs_pedigree$MF),length)
group_full_sibs=nrow(full_sibs_com)
family_size_full_sibs=mean(full_sibs_com[,2])
}



return(list(ped=O0OOOOO0O0OOO0OOO0O0,rename_ped=OOOOO0OOOOOOOOO0OOOO,pedigree_tree=OOOOOOOOO0OOO0OOOOOO,
				   error_id_set=list(error_duplicated_id=error_duplicated_id,error_sex_id=error_sex_id,
			         error_breed_id=error_breed_id,error_birth_date_id=error_birth_date_id),
				sibs=list(group_full_sibs=group_full_sibs,family_size_full_sibs=family_size_full_sibs)))

}



#检查系谱错误
OOOOOOO0O0O0O0OOOOO0<-function(ped,  #matrix format
                               dup_error_check=TRUE, 
						   sex_error_check=TRUE,
						   breed_error_check=FALSE,
						   birth_date_error_check=FALSE)
{
#将缺失值变为NA
ped[ped=="-9999"]=NA
ped[ped=="0"]=NA
ped[ped==""]=NA
ped=as.matrix(ped)
error_duplicated_id=NULL;error_sex_id=NULL;error_breed_id=NULL;error_birth_date_id=NULL;

#错误0：个体既出在子代又出现在父母列



#错误1：个体名相同，但是父母至少有一个不相同
if(dup_error_check==TRUE){
duplicated_id=unique(ped[duplicated(ped[,1]),1])

error_duplicated_id=NULL

if(length(duplicated_id)>0){

for(id in duplicated_id){
duplicated_set=ped[ped[,1]%in%id,]
if(length(unique(duplicated_set[,2]))>=2| length(unique(duplicated_set[,3]))>=2){
error_duplicated_id=c(error_duplicated_id,id)
}
}

if(length(error_duplicated_id)>0){
cat(paste0("Found ",length(error_duplicated_id)," duplicated id error records: offsprings with  same id but have different sire or dam records, records of sire and dam would be treated as missing value \n"))
ped[ped[,1]%in%error_duplicated_id,2:3]=NA
}
ped=ped[!duplicated(ped[,1]),]  #选择非重复的行
}
}


#错误2  性别有误，个体既出现在父亲这一列，又出现在母亲这一列
if(sex_error_check==TRUE){
error_sex_id=unique(na.omit(c(ped[,2][ped[,2]%in%ped[,3]],ped[ped[,1]==ped[,2]|ped[,1]==ped[,3],1])))
if(length(error_sex_id)>0){

cat(paste0("Found ",length(error_sex_id)," sex error records: ids in the sire column also appear in the dam column or offspring column, these ids would be treated as missing value \n"))


ped[ped[,2]%in%error_sex_id,2:3]=NA
ped[ped[,3]%in%error_sex_id,2:3]=NA
}}

#错误3  品种错误，对于某些系谱记录来说，个体ID 和 父母ID 记录有品种信息，因此可根据品种记录找出父母和子女品种记录不一致的个体，仅针对纯种
#默认品种信息记录在Id的开始两位，eg. YY201201
if(breed_error_check==TRUE){

O0OOO0OOOOO0OOO0OOOO=substr(ped[,1],1,2)
breed_sire=substr(ped[,2],1,2)
breed_dam=substr(ped[,3],1,2)

O0O0O0OOOOO0OOOOOOOO=(O0OOO0OOOOO0OOO0OOOO==breed_sire) & (O0OOO0OOOOO0OOO0OOOO==breed_dam) #逻辑值：父母的品种记录均和后代一致

if(sum(O0O0O0OOOOO0OOOOOOOO==FALSE)>=1){
cat("Found breed error: the breed of offspring, sire and dam are not equal  \n")
error_breed_id=ped[!O0O0O0OOOOO0OOOOOOOO,1]
ped[!O0O0O0OOOOO0OOOOOOOO,2:3]=NA
}
}

#错误4 根据出生日期的记录，判断后代的出生日期是否早于父母的出生日期，日期需要是数字格式，如：20180204
if(birth_date_error_check==TRUE){  

if(ncol(ped)<4){
stop("Error: peidgree provided doesn't have birth date records, please turn off birth_date_check function")
}

O0OOO0O0O0OOOOO0O0O0=ped[!(is.na(ped[,2])&is.na(ped[,3])),]  #去除父母为缺失的个体

ordered_offspring=O0OOO0O0O0OOOOO0O0O0[,4][order(as.numeric(O0OOO0O0O0OOOOO0O0O0[,4]))]

O0OOO0O0OOOOO0O0O0OO=(1:length(ordered_offspring))[match(O0OOO0O0O0OOOOO0O0O0[,4],ordered_offspring)]

sire_order=O0OOO0O0OOOOO0O0O0OO[match(O0OOO0O0O0OOOOO0O0O0[,2],O0OOO0O0O0OOOOO0O0O0[,1])]

dam_order=O0OOO0O0OOOOO0O0O0OO[match(O0OOO0O0O0OOOOO0O0O0[,3],O0OOO0O0O0OOOOO0O0O0[,1])]

OOOOO0OOOOOOO0OOOOO0=(O0OOO0O0OOOOO0O0O0OO>sire_order)|(O0OOO0O0OOOOO0O0O0OO>dam_order) #逻辑值

if(sum(na.omit(OOOOO0OOOOOOO0OOOOO0==FALSE))>=1){
error_birth_date_id=O0OOO0O0O0OOOOO0O0O0[!OOOOO0OOOOOOO0OOOOO0,1]
error_birth_date_id=na.omit(error_birth_date_id)
ped[ped[,1]%in%error_birth_date_id,2:3]=NA
cat(paste0("Found ",length(error_birth_date_id)," birth date error records: the birth date of offspring is before than the birth date of sire and dam \n"))
}
}



return(list(ped=ped,error_duplicated_id=error_duplicated_id,error_sex_id=error_sex_id,
			         error_breed_id=error_breed_id,error_birth_date_id=error_birth_date_id))
}


match_parents<-function(id,ped){ # ped 为三列向量（个体 父亲 母亲）， id 为一列向量 

parents=matrix(NA,ncol=2,nrow=length(id))

parents[,1]=ped[match(id,ped[,1]),2]
parents[,2]=ped[match(id,ped[,1]),3]
match_state=ifelse(sum(is.na(c(parents[,1],parents[,2])))==c(length(id)*2),"Complete_Unmatch","Match")
return(c(list(parents),list(match_state)))
}


single_pedigree<-function(ped){   #针对单次的排序，挑选出初代,中代,终代， 
if(is.data.frame(ped)){ped=as.matrix(ped)}
if(!is.matrix(ped)){ped=matrix(ped,ncol=length(ped))}
f0=ped[is.na(ped[,2])&is.na(ped[,3]),1]

ped=ped[!is.na(ped[,2])|!is.na(ped[,3]),]
if(!is.matrix(ped)){ped=matrix(ped,ncol=length(ped))}
n=unique(na.omit(ped[,1]))
P=unique(na.omit(c(ped[,2],ped[,3])))   #假设缺失值是NA, 这些个体是父母，它们有后代
f1=P[!P %in% n] # 这些个体没有再往上的系谱了， 所以假定它们是初代
f1=unique(c(f0,f1))
f2=P[P %in% n]  #这些个体可以往上追溯的父母，我们主要是将这些个体进行排序， 看谁是谁的曾祖代,祖代。。。。。。。
fn=n[!n %in% P] #这些个体不是父母，在父母那一栏里找不到这些个体， 所以这些个体是系谱最底层的那些个体
return(list(sort(f1),sort(f2),sort(fn)))
}

#get pedigree_tree_depth names
O0OOO0O0OOO0OOO0O0OO_colnames<-function(generation=1){
generation_colnames=c("Sire","Dam")
total=generation_colnames

if(generation>=2){
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:(generation-1)){
temp=NULL
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(generation_colnames)){
temp=c(temp,paste0(generation_colnames[OOO0OOO0O0OOO0O0OOO0OOOO],c("Sire","Dam")))
}
generation_colnames=temp
total=c(total,generation_colnames)
}}

return(total)
}


 
O0OOOOO0OOO0O0OOOOO0<-function(OOO0OOO0O0OOO0O0OOO0OOOOOO){

#check OOO0OOO0O0OOO0O0OOO0OOOOOO names
col=colnames(OOO0OOO0O0OOO0O0OOO0OOOOOO)
O0O0OOO0O0O0O0OOO0OO=round(nchar(col[1])/3)
O0O0OOO0O0O0O0OOO0OO_names=c("Offspring",O0OOO0O0OOO0OOO0O0OO_colnames(O0O0OOO0O0O0O0OOO0OO))


if(sum(col%in%O0O0OOO0O0O0O0OOO0OO_names==FALSE)>0){
error_colnames=col[!col%in%O0O0OOO0O0O0O0OOO0OO_names]
stop(paste0("Colnames of provided pedigree_data, including << ",error_colnames, ">> doesn't meet the requirement, please modify these columns names!"))
}

Total_result=list(OOO0OOO0O0OOO0O0OOO0OOOOOO[,1:3])
O0O0OOO0OOO0O0O0O0OO<-function(name,OOO0OOO0O0OOO0O0OOO0OOOOOO){
part=OOO0OOO0O0OOO0O0OOO0OOOOOO[,colnames(OOO0OOO0O0OOO0O0OOO0OOOOOO)%in%name]

if(paste0(name,"Sire")%in%colnames(OOO0OOO0O0OOO0O0OOO0OOOOOO)){
part=cbind(part,OOO0OOO0O0OOO0O0OOO0OOOOOO[,colnames(OOO0OOO0O0OOO0O0OOO0OOOOOO)%in%paste0(name,"Sire")])
}else{part=cbind(part,NA)}

if(paste0(name,"Dam")%in%colnames(OOO0OOO0O0OOO0O0OOO0OOOOOO)){
part=cbind(part,OOO0OOO0O0OOO0O0OOO0OOOOOO[,colnames(OOO0OOO0O0OOO0O0OOO0OOOOOO)%in%paste0(name,"Dam")])
}else{part=cbind(part,NA)}
colnames(part)=c("Offspring","Sire","Dam")
return(part)
}

for(OOO0OOO0O0OOO0O0OOO0OOOO in colnames(OOO0OOO0O0OOO0O0OOO0OOOOOO)){
if(paste0(OOO0OOO0O0OOO0O0OOO0OOOO,"Sire")%in%colnames(OOO0OOO0O0OOO0O0OOO0OOOOOO) | paste0(OOO0OOO0O0OOO0O0OOO0OOOO,"Dam")%in%colnames(OOO0OOO0O0OOO0O0OOO0OOOOOO)){
result=list(O0O0OOO0OOO0O0O0O0OO(OOO0OOO0O0OOO0O0OOO0OOOO,OOO0OOO0O0OOO0O0OOO0OOOOOO))
Total_result=c(Total_result,result)
}
}
OOO0OOO0O0OOO0O0OOO0OOOOOO=do.call(rbind,Total_result)
OOO0OOO0O0OOO0O0OOO0OOOOOO=OOO0OOO0O0OOO0O0OOO0OOOOOO[!duplicated(paste0(OOO0OOO0O0OOO0O0OOO0OOOOOO[,1],OOO0OOO0O0OOO0O0OOO0OOOOOO[,2],OOO0OOO0O0OOO0O0OOO0OOOOOO[,3])),]
message("Complete multi-coulums pedigree_data convert into standard 3 columns pedigree_data!")
return(OOO0OOO0O0OOO0O0OOO0OOOOOO)
}


#追溯全同胞
match_fullsibs<-function(trace_id=NULL,ped=NULL,max_fullsibs=NULL){
ped=ped[!is.na(ped[,2])|!is.na(ped[,3]),]
id_F=ped[ped[,1]%in%trace_id,3]
id_M=ped[ped[,1]%in%trace_id,2]
MF=paste0(id_M,"_",id_F)
sibs_ped=ped[paste0(ped[,2],"_",ped[,3])%in%MF,]

if(!is.null(max_fullsibs)){

if(max_fullsibs<1){stop("Please provide valid max_fullsibs!")}
sibs_MF=paste0(sibs_ped[,2],"_",sibs_ped[,3])
pos=rep(FALSE,length(sibs_MF))
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:max_fullsibs){
pos_tmp=!duplicated(sibs_MF)
sibs_MF[pos_tmp]=NA
pos=pos|pos_tmp
}
sibs=unique(sibs_ped[pos,1])
}else{
sibs=unique(sibs_ped[,1])
}
return(setdiff(sibs,trace_id))
}
#绘制方差组分的图
plot_dmu_blupf90_prior<-function(target_trait_name=NULL,random_effect_name=NULL,output_path=NULL,genetic_effect_name="Id"){
library(ggplot2)
library(RColorBrewer)
color_set=brewer.pal(9, "Set1")
color_set=color_set[c(1,3,2,4:9)]
unique_random=c(unique(do.call(c,random_effect_name)),"Residual")
group_color=color_set[1:length(unique_random)]
names(group_color)=unique_random

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(target_trait_name)){

a=data.table::fread(paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],"_heritability_result.txt"),data.table=F)
a[,-1]=round(a[,-1],2)

a$prior_se_blanket=paste0(a$prior,"(",a$prior_se,")")
a$h2_se_blanket=paste0(a$h2,"(",a$h2_se,")")
a$h2=a$prior

colnames(a)[2]="Proportion of variance"
colnames(a)[4]="Variance"

b=reshape2::melt(a)
b=b[b[,"variable"]%in%c("Proportion of variance","Variance"),]
#b$variable=factor(b$variable,levels = c("Variance","Proportion of variance"))
b$prior_se_blanket=ifelse(b$variable=="Variance",b$prior_se_blanket,b$h2_se_blanket)

random_effect_order=c(genetic_effect_name,sort(setdiff(unique(b$Random_effect_name),genetic_effect_name)))
b$Random_effect_name=factor(b$Random_effect_name,levels =rev(random_effect_order))

p=ggplot2::ggplot(data = b, aes(x=variable, y=value, fill = Random_effect_name)) +
  geom_col(aes(fill = Random_effect_name), position = "dodge") +
  geom_text(aes(label =prior_se_blanket), position = position_dodge(0.9),size=2*length(target_trait_name),vjust = 0.5,hjust=-0.25,family = "serif")+
  labs(fill='Random effect')+  #修改图例
  scale_fill_manual(values =group_color)+  #设置颜色
  theme_bw() +ylab("Estimates(SE)")+
  coord_flip()+
  theme(axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(), #去除x轴名称
        axis.text.y = element_text(hjust =0.5,face="bold",size=6.5*length(target_trait_name), family = "serif"),#调整x轴标签的内容显示
        axis.title.x=element_text(face="bold",size=6.5*length(target_trait_name), family = "serif"),
		axis.ticks.x=element_blank(), #去除y轴刻度线
	    axis.text.x=element_blank(),   #去除y轴刻度内容
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
	   legend.text=element_text(size=5*length(target_trait_name)), #调整图例文本字体大小
	   legend.title=element_text(size=5*length(target_trait_name)), #调整图例标题文本字体大小	   
        panel.background = element_blank())+ 
        scale_y_continuous(limits = c(0,max(b$value)*1.75))+
		theme(legend.position='bottom')+
		guides(fill = guide_legend(reverse = TRUE))
		

if(OOO0OOO0O0OOO0O0OOO0OOOO==1){total=p}
if(OOO0OOO0O0OOO0O0OOO0OOOO>=2){library(patchwork);total=total+p}
}

if(OOO0OOO0O0OOO0O0OOO0OOOO>=2){
final=total+ patchwork::plot_annotation(tag_levels = 'A')&theme(plot.tag = element_text(size = 10*length(target_trait_name),face="bold"))
}else{
final=total
}


if(!is.null(output_path)){
if(length(target_trait_name)==1){
ggsave(final,filename = paste0(output_path,"/Variance_plot.png"),width=5*length(target_trait_name),height=2.7*length(target_trait_name))
}

if(length(target_trait_name)>=2){
ggsave(final,filename = paste0(output_path,"/Variance_plot.png"),width=6*length(target_trait_name),height=2*length(target_trait_name))
}
}else{
return(final)
}

}



detect_pca_plot<-function(input_data_numeric=NULL,Breed=NULL){
library(ggplot2)
pca<-data.frame(rownames(input_data_numeric),prcomp(input_data_numeric)$x[,1:2],stringsAsFactors=F);colnames(pca)=c("Id",paste0("PC",1:2))
PCA_Breed=data.frame(Breed,stringsAsFactors=F);colnames(PCA_Breed)=c("Id","Breed")
pca=merge(pca,PCA_Breed,by="Id")

if(nrow(pca)==0){stop("Can't match the records  between genotype data and provided breed information!!!")}
message("Start detect the breed of population by PCA......")

#detect outlier breed by kmeans
target_Breed_number=length(unique(PCA_Breed[,2]))
km=kmeans(pca[,2:3],target_Breed_number)
aaa <- data.frame(pca, km$cluster)

Data=NULL
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:target_Breed_number){

tmp_Data <- pca[which(aaa$km.cluster == OOO0OOO0O0OOO0O0OOO0OOOO),]
tmp_Data$Expeced_Breed=names(table(tmp_Data[,"Breed"]))[table(tmp_Data[,"Breed"]) == max(table(tmp_Data[,"Breed"]))]
Data=rbind(Data,tmp_Data)
}

outlier_Data=Data[Data[,"Breed"]!=Data[,"Expeced_Breed"],]

library(RColorBrewer)
color_set=brewer.pal(8, "Set1")
color_set=color_set[1:target_Breed_number]

origin_pca=ggplot2::ggplot(pca)+geom_point(aes(x=PC1,y=PC2,shape=Breed,color=Breed),size=5)+
				theme(legend.title = element_text(size=20),
				      legend.text = element_text(size = 15),
					  axis.title = element_text(size=20, color="black", face= "bold", vjust=0.5,hjust=0.5),
					  axis.text = element_text(size=15, color="black", face= "bold", vjust=0.5,hjust=0.5),
				     panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border = element_blank(),
	                    panel.background = element_blank(),
                        axis.line = element_line(colour = "black")						
					 )+
				   scale_shape_manual(values=c(18, 17,15))+ 
				    scale_colour_manual(values=color_set)  #设置颜色

kmeans_pca=ggplot2::ggplot(Data)+geom_point(aes(x=PC1,y=PC2,shape=Expeced_Breed,color=Expeced_Breed),size=5)+
				theme(legend.title = element_text(size=20),
				      legend.text = element_text(size = 15),
					  axis.title = element_text(size=20, color="black", face= "bold", vjust=0.5,hjust=0.5),
					  axis.text = element_text(size=15, color="black", face= "bold", vjust=0.5,hjust=0.5),
				     panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border = element_blank(),
	                    panel.background = element_blank(),
                        axis.line = element_line(colour = "black")						
					 )+ 
					 labs(shape='Expect Breed',color='Expect Breed')+			
					 scale_shape_manual(values=c(18, 17,15))+ 
				    scale_colour_manual(values=color_set)  #设置颜色


library(patchwork)
p=origin_pca+kmeans_pca+patchwork::plot_annotation(tag_levels = 'A')&theme(plot.tag = element_text(size = 20,face="bold"))
#ggsave(origin_pca+kmeans_pca+patchwork::plot_annotation(tag_levels = 'A')&theme(plot.tag = element_text(size = 20,face="bold")),
#  filename = paste0(output_plot_path,"/PCA_detect.png"),width=20,height=10)
return(list(outlier=outlier_Data[,-c(2:3)],plot=p))
}




test123<-function(){

#从这个例子就可以看出来 trace_pedigree 追溯的系谱代数计算方法了。DD18922905仍然被视作初代，尽管其和第一代是夫妻，但是由于其无父母记录，因此将其视为初代。
ped=data.frame(Offspring=c("DD20967904","DD18923611"),
                  Sire=c("DD18922905","DD17773903"),
			    Dam=c("DD18923611","DD17516512"))
				
				
ped2=data.frame(
 Offspring=c("DD20967904","DD18923611","DD17516512","DD16546520","DD12546527","DD10546538","DD10246533"),
 Sire=c("DD18922905","DD17773903","DD16526514","DD12546527","DD10546538","DD10246533","DD09246524"),
 Dam=c("DD18923611","DD17516512","DD16546520","DD13546529","DD11546519","DD10546579","DD08246517"))				


}



#计算谱系数据
OOOOO0O0O0OOO0O0O0OO<-function(ped=NULL){



OOO0O0OOO0O0O0O0OOOO=plot_ped_single_pedigree(ped)
Ind=na.omit(unique(c(ped[,1],ped[,2],ped[,3])))
data_set=data.frame(Id=Ind,Family=NA)

Fn=OOO0O0OOO0O0O0O0OOOO[[3]]

get_related<-function(id,ped){
parents=as.vector(match_parents(id,ped)[[1]])
sibs=do.call(c,OOOOOOOOO0O0O0OOOOO0(id,ped))
offsprings=O0O0O0OOO0OOO0O0O0OO(id,ped)
related_ids=na.omit(unique(c(id,parents,sibs,offsprings)))
ped=as.matrix(ped)
related_ids=na.omit(unique(as.vector(ped[ped[,1]%in%related_ids|ped[,2]%in%related_ids|ped[,3]%in%related_ids,])))
return(related_ids)
}

k=0
for(OOO0OOO0O0OOO0O0OOO0OOOOO in 1:length(Fn)){
id=Fn[OOO0OOO0O0OOO0O0OOO0OOOOO]
if(all(!is.na(data_set[,2]))){break} #当所有个体都有family记录时，退出for循环
if(!is.na(data_set[data_set[,1]%in%id,2])){next} #当前个体有family记录时，退出当前循环
k=k+1
related_ids=get_related(id,ped)
n=1
while(length(related_ids)>n){
n=length(related_ids)
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(related_ids)){
related_ids=unique(na.omit(c(related_ids,get_related(related_ids[OOO0OOO0O0OOO0O0OOO0OOOO],ped))))
}
}
data_set[data_set[,1]%in%related_ids,2]=k
}
return(data_set)
}


plot_ped_single_pedigree<-function(ped){   #针对单次的排序，挑选出初代,中代,终代， 

n=unique(na.omit(ped[,1]))

P=unique(na.omit(c(ped[,2],ped[,3])))   #假设缺失值是NA, 这些个体是父母，它们有后代

f1=P[!P %in% n] # 这些个体没有再往上的系谱了， 所以假定它们是初代

f2=P[P %in% n]  #这些个体可以往上追溯的父母，我们主要是将这些个体进行排序， 看谁是谁的曾祖代,祖代。。。。。。。

fn=n[!n %in% P] #这些个体不是父母，在父母那一栏里找不到这些个体， 所以这些个体是系谱最底层的那些个体

return(list(f1,f2,fn))

}

#获得系谱中所有个体的代数
O0OOO0O0OOO0OOO0O0OO<-function(ped){
ped[ped==0]=NA
OOO0O0OOO0O0O0O0OOOO=plot_ped_single_pedigree(ped)
Whole_ped=data.frame(Id=unique(do.call(c,OOO0O0OOO0O0O0O0OOOO)))
Whole_ped$Generation=NA
Whole_ped$X_score=NA  #父亲-1，母亲+1
Whole_ped$Y_score=NA  #亲本+1
Fn=OOO0O0OOO0O0O0O0OOOO[[3]] #无父母,  存在 个体无后代，但其全同胞有后代的情况，这些个体的代数有问题，后面将解决代数问题

#去除Fn中-全同胞有后代的个体 
fullsibs=NULL
for(OOO0OOO0O0OOO0O0OOO0OOOO in Fn){
#将全同胞为父母的个体，从基础群中去除
if(length(OOOOOOOOO0O0O0OOOOO0(OOO0OOO0O0OOO0O0OOO0OOOO,ped)$full_sibs)>0){
if(sum(OOOOOOOOO0O0O0OOOOO0(OOO0OOO0O0OOO0O0OOO0OOOO,ped)$full_sibs%in%ped[,2])>0|sum(OOOOOOOOO0O0O0OOOOO0(OOO0OOO0O0OOO0O0OOO0OOOO,ped)$full_sibs%in%ped[,3])>0){
Fn=setdiff(Fn,OOO0OOO0O0OOO0O0OOO0OOOO)
}
}
}
OOO0OOO0O0OOO0O0OOO0OOOOO=0
data=data.frame(Id=Fn,Generation=OOO0OOO0O0OOO0O0OOO0OOOOO,stringsAsFactors=F)

while(length(Fn)>=1&match_parents(Fn,ped)[[2]]=="Match"){

  OOO0OOO0O0OOO0O0OOO0OOOOO=OOO0OOO0O0OOO0O0OOO0OOOOO+1
  Fn=setdiff(na.omit(unique(as.vector(match_parents(Fn,ped)[[1]]))),Fn)
 
  for( tmp in Fn){
  Fn=unique(c(Fn,do.call(c,OOOOOOOOO0O0O0OOOOO0(tmp,ped))))
  }
  if(length(Fn)==0){break}
  data=rbind(data,data.frame(Id=Fn,Generation=OOO0OOO0O0OOO0O0OOO0OOOOO))
}
return(data)
}



#获得个体及父母的 X轴，Y轴位置

O0OOOOO0O0O0O0OOO0O0_old<-function(X_score=0,Y_score=0,Fn=NULL,ped=NULL){
final_data=NULL
#根据nodes将 Fn排序
Fn=Fn[order(OOOOOOOOO0O0O0O0OOO0(Fn,ped),decreasing = T)]

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(Fn)){
Full_sibs_M=NULL
Full_sibs_F=NULL

id=Fn[OOO0OOO0O0OOO0O0OOO0OOOO]
individual_data=data.frame(Id=id,
		                       X_score=X_score,
				             Y_score=Y_score,
						   X_score_end=X_score,
						   Y_score_end=Y_score+1,
						   Type="Individual")

if(OOO0OOO0O0OOO0O0OOO0OOOO>1){
if(length(final_data[final_data[,1]%in%id,"Type"])>0){
if(final_data[final_data[,1]%in%id,"Type"]=="Sibs"){
next 
#已经在final_data里的个体，代表其位置信息已经确定，不需要再进行计算
}	
}
}					   
						   

if(match_parents(id,ped)[[2]]=="Complete_Unmatch"){
final_data=rbind(final_data,individual_data)
}else{

F=match_parents(id,ped)[[1]][1,1]  #父亲个体号
M=match_parents(id,ped)[[1]][1,2]  #母亲个体号

F_X_score=X_score-1
F_Y_score=Y_score+1

M_X_score=X_score+1
M_Y_score=Y_score+1

#判断个体有无全同胞，有全同胞的话，母亲X轴坐标向右移动(移动的长度等于全同胞的个数)
sibs_data=NULL
if(length(OOOOOOOOO0O0O0OOOOO0(id,ped)$full_sibs)>0){
sibs_Id=OOOOOOOOO0O0O0OOOOO0(id,ped)$full_sibs#全同胞的个体号
M_X_score=M_X_score+1*length(sibs_Id)

sibs_X_score=X_score+1:length(sibs_Id)
sibs_Y_score=Y_score

sibs_data=cbind(sibs_Id,sibs_X_score,sibs_Y_score)

sibs_data=data.frame(Id=sibs_Id,
		                X_score=sibs_X_score,
				       Y_score=sibs_Y_score,
					  X_score_end=sibs_X_score,
					  Y_score_end=sibs_Y_score+1,
					  Type="Sibs")

}

#从母亲出发，计算母亲上面分支的节点数
M_X_score=M_X_score+OOOOOOOOO0O0O0O0OOO0(M,ped)*1 #母亲X轴向右平移
F_X_score=F_X_score-OOOOOOOOO0O0O0O0OOO0(M,ped)*1 #父亲X轴向左平移


#从父亲出发，计算父亲上面分支的节点数
M_X_score=M_X_score+OOOOOOOOO0O0O0O0OOO0(F,ped)*1 #母亲X轴向右平移
F_X_score=F_X_score-OOOOOOOOO0O0O0O0OOO0(F,ped)*1 #父亲X轴向左平移

#统计父母的全同胞数目
Full_sibs_F=OOOOOOOOO0O0O0OOOOO0(F,ped)$full_sibs
Full_sibs_M=OOOOOOOOO0O0O0OOOOO0(M,ped)$full_sibs

parents_data=data.frame(Id=c(F,M),
		                   X_score=c(F_X_score,M_X_score),
				         Y_score=c(F_Y_score,M_Y_score),
					    X_score_end=c(M_X_score,M_X_score),
					    Y_score_end=c(M_Y_score,M_Y_score),
						Type=c("Father","Mother"))


final_data=rbind(final_data,sibs_data,parents_data,individual_data)
final_data=final_data[!is.na(final_data[,1]),]
#return(final_data)
}

X_score=max(final_data$X_score,final_data$X_score_end)+2*OOOOOOOOO0O0O0O0OOO0(Fn[OOO0OOO0O0OOO0O0OOO0OOOO+1],ped)
Y_score=0

#考虑母亲有多个全同胞，需要将下一个个体的X轴往右偏移
X_score=X_score+length(Full_sibs_M)

}

return(final_data)
}


O0OOOOO0O0O0O0OOO0O0<-function(X_score=0,Y_score=0,Fn=NULL,ped=NULL){
final_data=NULL
#根据nodes将 Fn排序

Fn=Fn[order(OOOOOOOOO0O0O0O0OOO0(Fn,ped),decreasing = T)]

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(Fn)){
Full_sibs_M=NULL
Full_sibs_F=NULL

id=Fn[OOO0OOO0O0OOO0O0OOO0OOOO]
individual_data=data.frame(Id=id,
		                       X_score=X_score,
				             Y_score=Y_score,
						   X_score_end=X_score,
						   Y_score_end=Y_score+1,
						   Type="Individual")

if(OOO0OOO0O0OOO0O0OOO0OOOO>1){
if(length(final_data[final_data[,1]%in%id,"Type"])>0){
if(final_data[final_data[,1]%in%id,"Type"]=="Sibs"){
next 
#已经在final_data里的个体，代表其位置信息已经确定，不需要再进行计算
}	
}
}					   
						   

if(match_parents(id,ped)[[2]]=="Complete_Unmatch"){
final_data=rbind(final_data,individual_data)
}else{

F=match_parents(id,ped)[[1]][1,1]  #父亲个体号
M=match_parents(id,ped)[[1]][1,2]  #母亲个体号

F_X_score=X_score-1
F_Y_score=Y_score+1

M_X_score=X_score+1
M_Y_score=Y_score+1

#判断个体有无全同胞，有全同胞的话，母亲X轴坐标向右移动(移动的长度等于全同胞的个数)
sibs_data=NULL
if(length(OOOOOOOOO0O0O0OOOOO0(id,ped)$full_sibs)>0){
sibs_Id=OOOOOOOOO0O0O0OOOOO0(id,ped)$full_sibs#全同胞的个体号


if(!is.na(id)&id%in%ped[,2]){
sibs_X_score=X_score-1:length(sibs_Id)  #当个体性别为公时，全同胞的X轴向左移动
F_X_score=F_X_score-1*length(sibs_Id)
}else {
sibs_X_score=X_score+1:length(sibs_Id)  #当个体性别为母时，全同胞的X轴向右移动
M_X_score=M_X_score+1*length(sibs_Id)
}

sibs_Y_score=Y_score

sibs_data=cbind(sibs_Id,sibs_X_score,sibs_Y_score)

sibs_data=data.frame(Id=sibs_Id,
		                X_score=sibs_X_score,
				       Y_score=sibs_Y_score,
					  X_score_end=sibs_X_score,
					  Y_score_end=sibs_Y_score+1,
					  Type="Sibs")

}

#从母亲出发，计算母亲上面分支的节点数
M_X_score=M_X_score+OOOOOOOOO0O0O0O0OOO0(M,ped)*1 #母亲X轴向右平移
F_X_score=F_X_score-OOOOOOOOO0O0O0O0OOO0(M,ped)*1 #父亲X轴向左平移


#从父亲出发，计算父亲上面分支的节点数
M_X_score=M_X_score+OOOOOOOOO0O0O0O0OOO0(F,ped)*1 #母亲X轴向右平移
F_X_score=F_X_score-OOOOOOOOO0O0O0O0OOO0(F,ped)*1 #父亲X轴向左平移

#统计父母的全同胞数目
Full_sibs_F=OOOOOOOOO0O0O0OOOOO0(F,ped)$full_sibs
Full_sibs_M=OOOOOOOOO0O0O0OOOOO0(M,ped)$full_sibs

parents_data=data.frame(Id=c(rep(F,length(F_X_score)),rep(M,length(M_X_score))),
		                   X_score=c(F_X_score,M_X_score),
				         Y_score=c(F_Y_score,M_Y_score),
					    X_score_end=c(M_X_score,M_X_score),
					    Y_score_end=c(M_Y_score,M_Y_score),
					    Type=c("Father","Mother"))


final_data=rbind(final_data,sibs_data,parents_data,individual_data)
final_data=final_data[!is.na(final_data[,1]),]
#return(final_data)
}

#X_score=max(final_data$X_score,final_data$X_score_end)+2*OOOOOOOOO0O0O0O0OOO0(Fn[OOO0OOO0O0OOO0O0OOO0OOOO+1],ped)
#id_max_X_score=final_data[which.max(final_data$X_score),1][1]

X_score=max(final_data$X_score,final_data$X_score_end)+4*OOOOOOOOO0O0O0O0OOO0(Fn[OOO0OOO0O0OOO0O0OOO0OOOO+1],ped)
#X_score=max(final_data$X_score,final_data$X_score_end)+2*OOOOOOOOO0O0O0O0OOO0(id_max_X_score,ped)+2*OOOOOOOOO0O0O0O0OOO0(Fn[OOO0OOO0O0OOO0O0OOO0OOOO+1],ped)
if(OOOOOOOOO0O0O0O0OOO0(Fn[OOO0OOO0O0OOO0O0OOO0OOOO+1],ped)==0){X_score=X_score+1}
Y_score=0

#考虑母亲有多个全同胞，需要将下一个个体的X轴往右偏移
X_score=X_score+length(Full_sibs_M)

}

return(final_data)
}

#获得个体上面的代数(等价于系谱树的分支数)
OOOOOOOOO0O0O0O0OOO0<-function(id,ped){
total_nodes=NULL
for(OOO0OOO0O0OOO0O0OOO0OOOO in id){
nodes=0
while(match_parents(OOO0OOO0O0OOO0O0OOO0OOOO,ped)[[2]]=="Match"){
nodes=nodes+1
OOO0OOO0O0OOO0O0OOO0OOOO=na.omit(as.vector(match_parents(OOO0OOO0O0OOO0O0OOO0OOOO,ped)[[1]]))
}
total_nodes=c(total_nodes,nodes)
}
return(total_nodes)
}

match_parents<-function(id,ped){ # ped 为三列向量（个体 父亲 母亲）， id 为一列向量 

parents=matrix(NA,ncol=2,nrow=length(id))

parents[,1]=ped[match(id,ped[,1]),2]
parents[,2]=ped[match(id,ped[,1]),3]
match_state=ifelse(sum(is.na(c(parents[,1],parents[,2])))==c(length(id)*2),"Complete_Unmatch","Match")
return(c(list(parents),list(match_state)))
}

#寻找个体的子代
O0O0O0OOO0OOO0O0O0OO<-function(id=NULL,ped=NULL){
offspring=NULL
if(!is.na(id)){
offspring=na.omit(unique(c(ped[ped[,2]%in%id,1],ped[ped[,3]%in%id,1])))
}
return(offspring)
}


#找寻个体的全同胞合半同胞
OOOOOOOOO0O0O0OOOOO0<-function(id=NULL,ped=NULL){
id_F=ped[ped[,1]%in%id,2]
id_M=ped[ped[,1]%in%id,3]
if(length(id_F)==0){id_F=NA}
if(length(id_M)==0){id_M=NA}

sibs_F=NULL
sibs_M=NULL
if(!is.na(id_F)){sibs_F=ped[ped[,2]%in%id_F,1]}
if(!is.na(id_M)){sibs_M=ped[ped[,3]%in%id_M,1]}

full_sibs=NULL
half_sibs=NULL
sibs=unique(c(id,sibs_F,sibs_M))
sibs=setdiff(sibs,id)

if(!is.na(id_F)&!is.na(id_M)){ #当父母其中有一方为缺失值的时候，sibs必定为半同胞
for(OOO0OOO0O0OOO0O0OOO0OOOO in sibs){
if(ped[ped[,1]%in%OOO0OOO0O0OOO0O0OOO0OOOO,2]%in%id_F&ped[ped[,1]%in%OOO0OOO0O0OOO0O0OOO0OOOO,3]%in%id_M){
full_sibs=c(full_sibs,OOO0OOO0O0OOO0O0OOO0OOOO)
}else{
half_sibs=c(half_sibs,OOO0OOO0O0OOO0O0OOO0OOOO)
}
}
}else{
half_sibs=c(half_sibs,sibs)
}
if(length(half_sibs)==0){half_sibs=NULL}
if(length(full_sibs)==0){full_sibs=NULL}

return(list(full_sibs=full_sibs,half_sibs=half_sibs))
}

#性别的形状

get_image<-function(shape_type=NULL,male_url=NULL,female_url=NULL){

if(!is.null(male_url)&!is.null(female_url)){
cat("Using user-provided figure to plot......")
male=male_url
female=female_url
}else{

if(shape_type==2){
female="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210708210603.png"
male="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210708210612.png"
}

if(shape_type==6){
male="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210709155300.png"
female="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210709155303.png"
}

if(shape_type==5){
female="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210708213611.png"
male="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210708213621.png"
}

if(shape_type==3){
male="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210711221548.png"
female="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210711221551.png"
}

if(shape_type==7){
male="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210711220705.png"
female="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210711220709.png"
}

if(shape_type==8){
male="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210711223116.png"
female="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210711223103.png"
}

if(shape_type==9){
male="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210711223704.png"
female="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210711223713.png"
}

}
return(list(male=male,female=female))
}

#ped=tmp_ped2
ggped<-function(input_pedigree=NULL,ind_text_size=4,ind_text_vjust=3,value_text_size=2,shape_size=8,plot_family=FALSE,trace_id=NULL,trace_sibs=FALSE,ind_sex=NULL,plot_inbred=FALSE,
				shape_type=4,male_url=NULL,female_url=NULL,show_curve=FALSE){

ped=input_pedigree

if(plot_inbred==TRUE){
K=cal_kinship(input_pedigree = ped,kinship_type = "P_A",return_result = TRUE)$P_A$A
Inbred=data.frame(Id=rownames(K),Inbred=diag(K)-1,stringsAsFactors=F)
rm(K);gc();
}

if(!is.null(trace_id)){
ped=trace_pedigree(input_pedigree=ped,trace_id=trace_id,trace_sibs=trace_sibs)[[1]][,1:3]
}

ped[ped==0]=NA
ped=ped[!(is.na(ped[,2])&is.na(ped[,3])),]
data_set=O0OOO0O0OOO0OOO0O0OO(ped)

Fn=data_set[data_set$Generation==0,1]

plot_pos=O0OOOOO0O0O0O0OOO0O0(X_score=0,Y_score=0,Fn=Fn,ped=ped) #起始位置

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:nrow(data_set)){

Start_id=data_set[OOO0OOO0O0OOO0O0OOO0OOOO,1]

if(Start_id %in% plot_pos[,1]){
if(match_parents(Start_id,ped)[[2]]=="Match"&all(plot_pos[plot_pos[,1]%in%Start_id,"Type"]!="Sibs")){
i_X_score=unique(plot_pos[plot_pos[,1]%in%Start_id,"X_score"])
i_Y_score=unique(plot_pos[plot_pos[,1]%in%Start_id,"Y_score"])
i_plot_pos=O0OOOOO0O0O0O0OOO0O0(X_score=i_X_score,Y_score=i_Y_score,Fn=Start_id,ped=ped)

plot_pos=rbind(plot_pos,i_plot_pos)
}

}else{
next  #已经在plot_pos里的个体，代表其位置信息已经确定，不需要再进行计算
}
plot_pos=plot_pos[!duplicated(paste0(plot_pos[,1],plot_pos[,4],plot_pos[,5])),]
}


if(!is.null(ind_sex)){

plot_pos$Sex=ind_sex[match(plot_pos[,1],ind_sex[,1]),2]
if(sum(is.na(plot_pos$Sex))>0){
stop("Please provide sex information for all individuals!")
}

}else{
plot_pos$Sex="Male"
plot_pos[plot_pos[,1]%in%ped[,3],"Sex"]="Female"
}

#统计同胞个体
Start_X_Y=paste0(plot_pos$X_score,"_",plot_pos$Y_score)
half_sibs_id=plot_pos$Id[duplicated(plot_pos$Id)&!duplicated(Start_X_Y)] #个体号相同，位置信息却不相同的id
half_sibs_id=unique(half_sibs_id)

#全同胞个体用虚线连接起来
total_half_sibs=NULL
if(length(half_sibs_id)>0){
for(OOO0OOO0O0OOO0O0OOO0OOOO in half_sibs_id){
half_sibs_score=plot_pos[plot_pos$Id%in%OOO0OOO0O0OOO0O0OOO0OOOO,]
half_sibs_score[1:(nrow(half_sibs_score)-1),4:5]=half_sibs_score[-1,2:3]
half_sibs_score=half_sibs_score[-nrow(half_sibs_score),]
total_half_sibs=rbind(total_half_sibs,half_sibs_score)
}

#check是否有重复的连线
status1=rep(NA,nrow(total_half_sibs))
status2=rep(NA,nrow(total_half_sibs))

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:nrow(total_half_sibs)){
status1[OOO0OOO0O0OOO0O0OOO0OOOO]=paste(total_half_sibs[OOO0OOO0O0OOO0O0OOO0OOOO,1:5],collapse = "_")
status2[OOO0OOO0O0OOO0O0OOO0OOOO]=paste(total_half_sibs[OOO0OOO0O0OOO0O0OOO0OOOO,c(1,4,5,2,3)],collapse = "_")
}
total_half_sibs=total_half_sibs[!duplicated(status1%in%status2),]
}




library(ggplot2)
library(RColorBrewer)
color_set=brewer.pal(9, "Set1")

#将近交系数添加到数据集里
if(plot_inbred==TRUE){
plot_pos$Inbred_value=Inbred[match(plot_pos$Id,Inbred$Id),2]
}



if((shape_type!=1&shape_type!=4)|(!is.null(male_url)&!is.null(female_url))){
sex_url=get_image(shape_type,male_url=male_url,female_url=female_url)
male=sex_url$male
female=sex_url$female
plot_pos$Image=male
plot_pos[plot_pos$Sex=="Female","Image"]=female
}



if(plot_family==FALSE){
p=ggplot(plot_pos,aes(x=X_score,y=Y_score))+
        geom_segment(aes(x=X_score,
		              xend=X_score_end,
				    y=Y_score,
				    yend=Y_score_end))

if(show_curve==TRUE){
if(length(half_sibs_id)>=1){

p=p+geom_curve(aes(x =X_score,
                 y =Y_score,
                 xend=X_score_end, 
                 yend =Y_score_end),
				 data=total_half_sibs)

}}

if(shape_type==4){
p=p+geom_point(aes(x=X_score,y=Y_score,shape=Sex,color=Sex),size=shape_size)
}else if(shape_type==1){
p=p+geom_point(aes(x=X_score,y=Y_score,shape=Sex),size=shape_size)
}else{
p=p+ggimage::geom_image(size =shape_size/150,aes(image=Image))
}		
	
	p=p+geom_text(aes(label =Id),size=ind_text_size,vjust = ind_text_vjust,family = "serif")+
   theme(
        axis.title.y=element_blank(), #去除y轴名称
	   axis.title.x=element_blank(), #去除x轴名称
        axis.ticks.x=element_blank(), #去除x轴刻度线
	   axis.ticks.y=element_blank(), #去除y轴刻度线	
	   axis.text.x=element_blank(),   #去除y轴刻度内容
	   axis.text.y=element_blank(),   #去除y轴刻度内容
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
	   legend.position = "none", #去除所有图例内容   
        panel.background = element_blank())+		
		scale_y_continuous(limits = c(min(plot_pos$Y_score)-1/max(plot_pos$Y_score),max(plot_pos$Y_score)))+
	    scale_x_continuous(limits = c(min(plot_pos$X_score),max(plot_pos$X_score)))+
	   scale_shape_manual(values = c(16,15))+
	   scale_colour_manual(values =color_set)  #设置颜色
	   
if(plot_inbred==TRUE){
p=p+geom_text(aes(label =Inbred_value),size=value_text_size,vjust = 0,family = "serif")
}	   
	   
return(p)
}else{
family_set=OOOOO0O0O0OOO0O0O0OO(ped)
plot_pos$Family=family_set[match(plot_pos$Id,family_set$Id),2]
total_half_sibs$Family=family_set[match(total_half_sibs$Id,family_set$Id),2]
Family_value=unique(plot_pos$Family)


tmp_plot_pos=plot_pos[plot_pos$Family%in%1,]
p=ggplot(tmp_plot_pos,aes(x=X_score,y=Y_score))+
        geom_segment(aes(x=X_score,
		              xend=X_score_end,
				    y=Y_score,
				    yend=Y_score_end))

if(show_curve==TRUE){
if(length(half_sibs_id)>=1){
tmp_total_half_sibs=total_half_sibs[total_half_sibs$Family%in%1,]
p=p+geom_curve(aes(x =X_score,
                 y =Y_score,
                 xend=X_score_end, 
                 yend =Y_score_end),
			data=tmp_total_half_sibs)
}}		

if(shape_type==4){
p=p+geom_point(aes(x=X_score,y=Y_score,shape=Sex,color=Sex),size=shape_size)
}else if(shape_type==1){
p=p+geom_point(aes(x=X_score,y=Y_score,shape=Sex),size=shape_size)
}else{
p=p+ggimage::geom_image(size =shape_size/150,aes(image=Image))
}
	p=p+geom_text(aes(label =Id),size=ind_text_size,vjust = ind_text_vjust,family = "serif")+
	  xlab("Family 1")+
   theme(
	   axis.title.x=element_text(face="bold",size=15, family = "serif"),
        axis.title.y=element_blank(), #去除y轴名称
	   #axis.title.x=element_blank(), #去除x轴名称
        axis.ticks.x=element_blank(), #去除x轴刻度线
	   axis.ticks.y=element_blank(), #去除y轴刻度线	
	   axis.text.x=element_blank(),   #去除x轴刻度内容
	   axis.text.y=element_blank(),   #去除y轴刻度内容
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
	   legend.position = "none", #去除所有图例内容   
        panel.background = element_blank())+		
		scale_y_continuous(limits = c(min(tmp_plot_pos$Y_score)-1/max(plot_pos$Y_score),max(tmp_plot_pos$Y_score)))+
	    scale_x_continuous(limits = c(min(tmp_plot_pos$X_score),max(tmp_plot_pos$X_score)))+
	   scale_shape_manual(values = c(16,15))+
	   scale_colour_manual(values =color_set)  #设置颜色

if(plot_inbred==TRUE){
p=p+geom_text(aes(label =Inbred_value),size=value_text_size,vjust = 0,family = "serif")
}
	   
sum_p=p	   
if(length(Family_value)>=2){

for(OOO0OOO0O0OOO0O0OOO0OOOO in 2:length(Family_value)){

tmp_plot_pos=plot_pos[plot_pos$Family%in%OOO0OOO0O0OOO0O0OOO0OOOO,]
p=ggplot(tmp_plot_pos,aes(x=X_score,y=Y_score))+
        geom_segment(aes(x=X_score,
		              xend=X_score_end,
				    y=Y_score,
				    yend=Y_score_end))

if(length(half_sibs_id)>=1){
tmp_total_half_sibs=total_half_sibs[total_half_sibs$Family%in%OOO0OOO0O0OOO0O0OOO0OOOO,]
p=p+geom_curve(aes(x =X_score,
                 y =Y_score,
                 xend=X_score_end, 
                 yend =Y_score_end),
				 data=tmp_total_half_sibs)
}

if(shape_type==4){
p=p+geom_point(aes(x=X_score,y=Y_score,shape=Sex,color=Sex),size=shape_size)
}else if(shape_type==1){
p=p+geom_point(aes(x=X_score,y=Y_score,shape=Sex),size=shape_size)
}else{
p=p+ggimage::geom_image(size =shape_size/150,aes(image=Image))
}	
	p=p+geom_text(aes(label =Id),size=ind_text_size,vjust = ind_text_vjust,family = "serif")+
		xlab(paste0("Family ",OOO0OOO0O0OOO0O0OOO0OOOO))+
   theme(
	   axis.title.x=element_text(face="bold",size=15, family = "serif"),
        axis.title.y=element_blank(), #去除y轴名称
	   #axis.title.x=element_blank(), #去除x轴名称
        axis.ticks.x=element_blank(), #去除x轴刻度线
	   axis.ticks.y=element_blank(), #去除y轴刻度线	
	   axis.text.x=element_blank(),   #去除x轴刻度内容
	   axis.text.y=element_blank(),   #去除y轴刻度内容
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
	   legend.position = "none", #去除所有图例内容   
        panel.background = element_blank())+		
		scale_y_continuous(limits = c(min(tmp_plot_pos$Y_score)-1/max(plot_pos$Y_score),max(tmp_plot_pos$Y_score)))+
	    scale_x_continuous(limits = c(min(tmp_plot_pos$X_score),max(tmp_plot_pos$X_score)))+
	   scale_shape_manual(values = c(16,15))+
	   scale_colour_manual(values =color_set)  #设置颜色
library(patchwork)	 
sum_p=sum_p+p
}}
sum_p=sum_p+ patchwork::plot_annotation(tag_levels = 'A')&theme(plot.tag = element_text(size = 20,face="bold",, family = "serif"))

return(sum_p)

}	   
}


show_ped_type<-function(){
library(ggplot2)
library(patchwork)
Theme=theme(
        axis.title.y=element_blank(), #去除y轴名称
	   axis.title.x=element_blank(), #去除x轴名称
        axis.ticks.x=element_blank(), #去除x轴刻度线
	   axis.ticks.y=element_blank(), #去除y轴刻度线	
	   axis.text.x=element_blank(),   #去除y轴刻度内容
	   axis.text.y=element_blank(),   #去除y轴刻度内容
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
	   legend.position = "none", #去除所有图例内容   
        panel.background = element_blank())	

plot_pos=data.frame(Id=c("Male","Female"),Sex=c("Male","Female"),X_score=c(1,1),Y_score=c(1,2.3))

four=ggplot(plot_pos,aes(x=X_score,y=Y_score))+
  geom_point(aes(x=X_score,y=Y_score,shape=Sex,color=Sex),size=10)+
  geom_text(aes(label =Id),size=4,vjust =3,family = "serif")+
  scale_shape_manual(values = c(16,15))+
	   scale_colour_manual(values =color_set)+Theme+
	   scale_y_continuous(limits = c(-1,3))
  
one=ggplot(plot_pos,aes(x=X_score,y=Y_score))+
  geom_point(aes(x=X_score,y=Y_score,shape=Sex),size=10)+
  geom_text(aes(label =Id),size=4,vjust =3,family = "serif")+
  scale_shape_manual(values = c(16,15))+
	   scale_colour_manual(values =color_set)+Theme+
	   scale_y_continuous(limits = c(-1,3))

for(shape_type in c(2:3)){
sex_url=get_image(shape_type=shape_type)
male=sex_url$male
female=sex_url$female
plot_pos$Image=male
plot_pos[plot_pos$Sex=="Female","Image"]=female  
p=ggplot(plot_pos,aes(x=X_score,y=Y_score))+
  geom_text(aes(label =Id),size=4,vjust =3,family = "serif")+
  ggimage::geom_image(size =.08,aes(image=Image))+Theme+
  scale_y_continuous(limits = c(-1,3))
one=one+p
}

one=one+four

for(shape_type in c(5:6)){
sex_url=get_image(shape_type=shape_type)
male=sex_url$male
female=sex_url$female
plot_pos$Image=male
plot_pos[plot_pos$Sex=="Female","Image"]=female  
p=ggplot(plot_pos,aes(x=X_score,y=Y_score))+
  geom_text(aes(label =Id),size=4,vjust =3,family = "serif")+
  ggimage::geom_image(size =.07,aes(image=Image))+Theme+
  scale_y_continuous(limits = c(-1,3))
one=one+p
} 

for(shape_type in c(7:9)){
sex_url=get_image(shape_type=shape_type)
male=sex_url$male
female=sex_url$female
plot_pos$Image=male
plot_pos[plot_pos$Sex=="Female","Image"]=female  
p=ggplot(plot_pos,aes(x=X_score,y=Y_score))+
  geom_text(aes(label =Id),size=4,vjust =3,family = "serif")+
  ggimage::geom_image(size =0.1,aes(image=Image))+Theme+
  scale_y_continuous(limits = c(-1,3))
one=one+p
} 
one=one+plot_annotation(tag_levels = list(paste0("Shape type=",1:9)))
return(one)
}

#ggsave(show_ped_type(),filename="c:/Users/26564/Desktop/test1.png",width = 10,height =7)		


#ped_plot=ggped(input_pedigree = blupADC::plot_pedigree,trace_id = c("121","125"),trace_sibs = T,shape_type =6,shape_size = 4)
#ggsave(ped_plot,filename="c:/Users/26564/Desktop/test-t.png",width = 10,height =7)



#ped_plot=ggped(input_pedigree = blupADC::plot_pedigree,trace_sibs = T,shape_type =4,
#               shape_size =8,trace_id = c("121","122","210","214","205","206","207"),plot_family=T)
#ggsave(ped_plot,filename="c:/Users/26564/Desktop/test-4.png",width =13,height =8)#' RUN BLUPF90 software in R
#' @export
#' @param phe_col_names colnames of phenotype
#' @param phe_name phenotype file name
#' @param target_trait_name name of target trait
#' @param fixed_effect_name names of trait's fixed effect
#' @param random_effect_name names of trait's random effect
#' @param covariate_effect_name names of trait's covariate effect
#' @param provided_effect_file_path effect file path
#' @param provided_effect_file_name effect file name
#' @param missing_value missing value in phenotype
#' @param relationship_name relationship file name
#' @param relationship_path relationship file path
#' @param analysis_model BLUPF90 analysis model
#' @param genetic_effect_name name of genetic effect
#' @param included_permanent_effect whehter include permanet effect in analysis
#' @param provided_BLUPF90_prior_file_path prior file path
#' @param provided_BLUPF90_prior_file_name prior file name
#' @param provided_BLUPF90_prior_effect_name prior effect name
#' @param BLUPF90_genumeric_name name
#' @param BLUPF90_map_name map name
#' @param output_result_path name
#' @param output_ebv_path map name
#' @param output_ebv_name map name
#' @param BLUPF90_algorithm map name
#' @param BLUPF90_software_path map name
#' @param phe_path phenotype file path
run_BLUPF90<-function(
								   phe_col_names=NULL,
								   target_trait_name=NULL,
								   fixed_effect_name=NULL, #列表
								   random_effect_name=NULL, #列表，不包括永久环境效应
								   covariate_effect_name=NULL, #列表
								   provided_effect_file_path=NULL, #各个性状的效应文件
								   provided_effect_file_name="model_define.txt",
								   phe_path=NULL,
								   phe_name=NULL,
								   user_file_id=NULL, # user_file_id, 
								   output_result_path=NULL,
								   output_ebv_path=NULL,
						            output_ebv_name=NULL,
								   missing_value="-9999",
								   relationship_name=NULL,
								   relationship_path=NULL,
								   analysis_model="PBLUP_A",
								   BLUPF90_algorithm="AI_REML",
								   genetic_effect_name="Id",
								   included_permanent_effect=FALSE,
								   provided_BLUPF90_prior_file_path=NULL,
								   provided_BLUPF90_prior_file_name=NULL,
								   provided_BLUPF90_prior_effect_name=NULL, #随机效应的名称, 包括Residual
					                 provided_renf90_par_file_path=NULL,
					                 provided_renf90_par_file_name=NULL,
								   BLUPF90_genumeric_name=NULL,
								   BLUPF90_map_name=NULL,
								   plot_variance=FALSE,
								   BLUPF90_software_path=ifelse(as.character(Sys.info()["sysname"])=="Linux",system.file("extdata/bin", package = "blupADC"),system.file("extdata/bin_windows", package = "blupADC"))

){

cat("R package:blupADC is only the wrapper of blupf90-package in the field of academic research! \n")
cat("For some very complicated models, a direct use of the blupf90-package is needed ! \n")
cat("For commercial use of the blupf90-package please contact Ignacy Misztal. Email: ignacy@uga.edu  !\n")

cat(paste0("Start the ",analysis_model," analyse of ",length(target_trait_name)," trait model:",paste(target_trait_name,collapse = " & ")," \n"))

if(length(grep(" ",output_result_path))>0){
stop("output_result_path shouldn't contain blank, please change into another path!")
}

output_result_path=paste0(output_result_path,"/",analysis_model,"_",paste(target_trait_name,collapse = "_"))
if(!file.exists(output_result_path)){dir.create(output_result_path,recursive=TRUE)}

setwd(output_result_path)
output_result_path=getwd()
cat(paste0("Results are saved in path: ",output_result_path,"\n"))
file.copy(from=paste0(phe_path,"/",phe_name),to=output_result_path,overwrite=TRUE)

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(relationship_name)){
file.copy(from=paste0(relationship_path,"/",relationship_name[OOO0OOO0O0OOO0O0OOO0OOOO]),to=output_result_path,overwrite=TRUE)
}

if(!is.null(provided_renf90_par_file_path)&!is.null(provided_renf90_par_file_name)){ #用户提供 blupf90的参数文件

file.copy(from=paste0(provided_renf90_par_file_path,"/",provided_renf90_par_file_name),to=paste0(output_result_path,"renf90.par"),overwrite=TRUE)

}else{#自动生成 blupf90的参数文件

generate_renum(
			    phe_col_names=phe_col_names,
			    phe_path=phe_path,
			    phe_name=phe_name,
				user_file_id=user_file_id,
			    target_trait_name=target_trait_name,
			    fixed_effect_name=fixed_effect_name, #列表
			    random_effect_name=random_effect_name, #列表，不包括永久环境效应
			    covariate_effect_name=covariate_effect_name, #列表
			    provided_effect_file_path=provided_effect_file_path, #各个性状的效应文件
			    provided_effect_file_name=provided_effect_file_name,
			    analysis_model=analysis_model,
			    relationship_name=relationship_name,
			    relationship_path=relationship_path,
			    genetic_effect_name=genetic_effect_name,
				missing_value=missing_value,
			    included_permanent_effect=included_permanent_effect,
			    provided_BLUPF90_prior_file_path=provided_BLUPF90_prior_file_path,
			    provided_BLUPF90_prior_file_name=provided_BLUPF90_prior_file_name,
			    provided_BLUPF90_prior_effect_name=provided_BLUPF90_prior_effect_name, #随机效应的名称, 包括Residual
			    BLUPF90_genumeric_name=BLUPF90_genumeric_name,
			    BLUPF90_map_name=BLUPF90_map_name)

cat("Start running renumf90 module of BLUPF90......\n")
system2(paste0(BLUPF90_software_path,"/renumf90"),"renum.par",stdout="renumf90.log")

if(analysis_model=="User_define"){
#rename renf90.dat 
#phe=read.table(phe_name,header=F);colnames(phe)=phe_col_names
#renf90.phe=read.table("renf90.dat",header=F)
#renf90.phe[,ncol(renf90.phe)]=phe[,genetic_effect_name]
#write.table(renf90.phe,"renf90.dat",quote=F,row.names=F,sep=" ",col.names=F)

#modify renf90.par

renf90=data.table::fread("renf90.par",data.table=F,fill=T)
pos=match("add_animal",renf90[,1])
renf90[pos,1]="user_file"
renf90[pos+2,1]=relationship_name[1]
renf90[is.na(renf90)]=""

#读取 renadd.ped 将表型数据中的新id

write.table(renf90,"renf90.par",quote=F,row.names=F,col.names=T)
}

}
system("ulimit -s unlimited")
#system2(paste0(BLUPF90_software_path,"/renumf90"),"renum.par",stdout="renumf90.log")
if(BLUPF90_algorithm=="AI_REML"){
cat("Start running AI REML module of BLUPF90......\n")
system2(paste0(BLUPF90_software_path,"/airemlf90"),"renf90.par",stdout="ai_remlf90.log")
cat("Complete running AI REML module of BLUPF90\n")
}else if(BLUPF90_algorithm=="EM_REML"){
cat("Start running EM REML module of BLUPF90......\n")
system2(paste0(BLUPF90_software_path,"/remlf90"),"renf90.par",stdout="em_remlf90.log")
cat("Complete running AI REML module of BLUPF90\n")
}else if(BLUPF90_algorithm=="BLUP"){
cat("Start running BLUP module of BLUPF90......\n")
cat("No need to estimate variace components......\n")
system2(paste0(BLUPF90_software_path,"/blupf90"),"renf90.par",stdout="blup_remlf90.log")
cat("Complete running BLUP module of BLUPF90\n")
}

#计算校正表型
O0O0OOO0O0OOOOO0O0OO=cal_corrected_phe_BLUPF90(target_trait_name=target_trait_name)
O0O0OOO0O0OOOOO0O0OO=O0O0OOO0O0OOOOO0O0OO[order(O0O0OOO0O0OOOOO0O0OO[,1]),]
if(is.null(output_ebv_path)){output_ebv_path=output_result_path}
if(is.null(output_ebv_name)){output_ebv_name="corrected_phe_BLUPF90"}

cat("Output the corrected phenotype......\n")
setwd(output_ebv_path)
utils::write.table(O0O0OOO0O0OOOOO0O0OO,paste0("colnames_",output_ebv_name,".txt"),quote=F,row.names=F,col.names=T,sep=" ")

O0O0OOO0O0OOOOO0O0OO[is.na(O0O0OOO0O0OOOOO0O0OO)]=missing_value
setwd(output_ebv_path)
utils::write.table(O0O0OOO0O0OOOOO0O0OO,paste0(output_ebv_name,".txt"),quote=F,row.names=F,col.names=F,sep=" ")


#统计方差组分结果-遗传力及标准误

if(BLUPF90_algorithm!="BLUP"){
cal_blupf90_se_reml(target_trait_name=target_trait_name,random_effect_name=random_effect_name,genetic_effect_name=genetic_effect_name,BLUPF90_algorithm=BLUPF90_algorithm)
}

if(plot_variance==TRUE){
message("Plot the result of variance components......")
plot_dmu_blupf90_prior(genetic_effect_name=genetic_effect_name,target_trait_name=target_trait_name,random_effect_name=random_effect_name,output_path=output_result_path)
}

}


cal_corrected_phe_BLUPF90<-function(target_trait_name=NULL){

ebv=data.table::fread("solutions",data.table=F,header=F)
colnames(ebv)=c("OOOOO0O0O0O0OOOOOOO0","Effect_number","Id","EBV","SE_EBV")
max_effect_number=max(ebv[,"Effect_number"])
ebv=ebv[ebv[,"Effect_number"]%in%max_effect_number,]


residual=data.table::fread("yhat_residual",data.table=F,header=F)
O0O0O0O0O0O0OOOOOOO0=data.table::fread("renf90.dat",data.table=F,header=F)
residual$Id=O0O0O0O0O0O0OOOOOOO0[,ncol(O0O0O0O0O0O0OOOOOOO0)]


O0O0OOO0O0OOOOO0O0OO=matrix(NA,nrow=length(unique(ebv[,"Id"])),ncol=1+4*length(target_trait_name))
O0O0OOO0O0OOOOO0O0OO=as.data.frame(O0O0OOO0O0OOOOO0O0OO,stringsAsFactors=F)
O0O0OOO0O0OOOOO0O0OO[,1]=sort(unique(ebv[,"Id"]));colnames(O0O0OOO0O0OOOOO0O0OO)[1]="Id"

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(target_trait_name)){

trait_ebv=ebv[ebv$OOOOO0O0O0O0OOOOOOO0%in%OOO0OOO0O0OOO0O0OOO0OOOO,]
O0O0OOO0O0OOOOO0O0OO[,1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*4+1]= trait_ebv[match(O0O0OOO0O0OOOOO0O0OO[,1],trait_ebv[,"Id"]),"EBV"]
O0O0OOO0O0OOOOO0O0OO[,1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*4+2]= residual[match(O0O0OOO0O0OOOOO0O0OO[,1],residual[,"Id"]),length(target_trait_name)+OOO0OOO0O0OOO0O0OOO0OOOO]
O0O0OOO0O0OOOOO0O0OO[,1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*4+3]=O0O0OOO0O0OOOOO0O0OO[,1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*4+1]+O0O0OOO0O0OOOOO0O0OO[,1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*4+2]

colnames(O0O0OOO0O0OOOOO0O0OO)[1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*4+1]=paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],"_EBV")
colnames(O0O0OOO0O0OOOOO0O0OO)[1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*4+2]=paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],"_Res")
colnames(O0O0OOO0O0OOOOO0O0OO)[1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*4+3]=paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],"_EBV_Plus_Res")
colnames(O0O0OOO0O0OOOOO0O0OO)[1+(OOO0OOO0O0OOO0O0OOO0OOOO-1)*4+4]=paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],"_dEBV")

}


#将rename 后的名称映射回来
OOOOO0OOOOOOO0O0O0OO=list.files(getwd(),pattern = '.ped$')
if(length(OOOOO0OOOOOOO0O0O0OO)>1){
OOOOO0OOOOOOO0O0O0OO=paste0("renadd0",max_effect_number,".ped")
}
BLUPF90_key=data.table::fread(OOOOO0OOOOOOO0O0O0OO,data.table=F,header=F)
O0O0OOO0O0OOOOO0O0OO[,"Id"]=BLUPF90_key[match(O0O0OOO0O0OOOOO0O0OO[,"Id"],BLUPF90_key[,1]),10]

return(O0O0OOO0O0OOOOO0O0OO)


}




cal_blupf90_se_reml<-function(target_trait_name=NULL,random_effect_name=NULL,genetic_effect_name="Id",BLUPF90_algorithm=NULL){

O0OOO0O0O0O0O0OOOOOO<-function(infor_matrix=NULL,OOOOOOOOOOOOOOO0OOOO=NULL,Prior=NULL){  #根据信息矩阵计算标准误

if(is.null(OOOOOOOOOOOOOOO0OOOO)){
Var=OOO0OOO0O0O0OOOOOOO0(infor_matrix)  #信息矩阵的逆矩阵为 方差-协方差矩阵
}else{
Var=OOOOOOOOOOOOOOO0OOOO
}

heritablilty=Prior/sum(Prior[,1])
SE_h2=matrix(NA,nrow=nrow(Prior),ncol=1)
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:nrow(Prior)){

derivate=NULL
for(OOO0OOO0O0OOO0O0OOO0OOOOO in 1:nrow(Prior)){

if(OOO0OOO0O0OOO0O0OOO0OOOOO==OOO0OOO0O0OOO0O0OOO0OOOO){derivate=cbind(derivate,(sum(Prior[,1])-Prior[OOO0OOO0O0OOO0O0OOO0OOOO,1])/(sum(Prior[,1])^2) )}
if(OOO0OOO0O0OOO0O0OOO0OOOOO!=OOO0OOO0O0OOO0O0OOO0OOOO){derivate=cbind(derivate,(-Prior[OOO0OOO0O0OOO0O0OOO0OOOO,1])/(sum(Prior[,1])^2 ))}

}
SE_h2[OOO0OOO0O0OOO0O0OOO0OOOO,1]=sqrt(derivate%*%Var%*%t(derivate)) # 开平方根才是最终结果
}
return(SE_h2)
}

temp=fread("AI_inv.txt",data.table=F)
O0OOO0OOO0O0O0OOOOOO=diag(max(temp[,1]))
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:nrow(temp)){
O0OOO0OOO0O0O0OOOOOO[temp[OOO0OOO0O0OOO0O0OOO0OOOO,1],temp[OOO0OOO0O0OOO0O0OOO0OOOO,2]]=temp[OOO0OOO0O0OOO0O0OOO0OOOO,7]
O0OOO0OOO0O0O0OOOOOO[temp[OOO0OOO0O0OOO0O0OOO0OOOO,2],temp[OOO0OOO0O0OOO0O0OOO0OOOO,1]]=temp[OOO0OOO0O0OOO0O0OOO0OOOO,7]
}

if(BLUPF90_algorithm=="AI_REML"){
log=fread("airemlf90.log",data.table=F,fill = T)
}else if(BLUPF90_algorithm=="EM_REML"){
log=fread("emremlf90.log",data.table=F,fill = T)
}else if(BLUPF90_algorithm=="BLUP"){
log=fread("blupf90.log",data.table=F,fill = T)
}

#分性状计算
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(target_trait_name)){

OOOOOOOOOOOOOOO0OOOO=O0OOO0OOO0O0O0OOOOOO[seq(OOO0OOO0O0OOO0O0OOO0OOOO,nrow(O0OOO0OOO0O0O0OOOOOO),length(target_trait_name)),seq(OOO0OOO0O0OOO0O0OOO0OOOO,nrow(O0OOO0OOO0O0O0OOOOOO),length(target_trait_name))]
OOOOOOOOOOOOOOO0OOOO=OOOOOOOOOOOOOOO0OOOO[!diag(OOOOOOOOOOOOOOO0OOOO)==0,!diag(OOOOOOOOOOOOOOO0OOOO)==0]

Prior_pos=(1:nrow(log))[log[,1]%in%c("Genetic","Residual")]+OOO0OOO0O0OOO0O0OOO0OOOO
Prior=matrix(as.numeric(log[Prior_pos,OOO0OOO0O0OOO0O0OOO0OOOO]))

OOO0O0OOO0O0O0OOO0O0=(1:nrow(log))[log[,1]%in%c("SE")]+OOO0OOO0O0OOO0O0OOO0OOOO
Prior_se=matrix(as.numeric(log[OOO0O0OOO0O0O0OOO0O0,OOO0OOO0O0OOO0O0OOO0OOOO]))

Prior_se=Prior_se[!Prior[,1]==1,]
Prior=Prior[!Prior[,1]==1,] #去除冗余的随机效应(e.g. 性状没有该随机效应，但是BLUPF90会自动赋予该效应的方差为1)
Prior=as.matrix(Prior)

h2_se=O0OOO0O0O0O0O0OOOOOO(OOOOOOOOOOOOOOO0OOOO=OOOOOOOOOOOOOOO0OOOO,Prior=Prior)
h2=Prior/sum(Prior[,1])

temp_random_effect_name=random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]]
temp_random_effect_name=sort(temp_random_effect_name[!temp_random_effect_name%in%genetic_effect_name])
temp_random_effect_name=c(temp_random_effect_name,"Id","Residual")

data_h2=cbind(temp_random_effect_name,Prior,Prior_se,h2,h2_se)
colnames(data_h2)=c("Random_effect_name","prior","prior_se","h2","h2_se")
write.table(data_h2,paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],"_heritability_result.txt"),quote=F,row.names=F,col.names=T,sep="\t")

}

}





run_DMU<-function(
	                   phe_col_names=NULL,
					   target_trait_name=NULL,
					   fixed_effect_name=NULL, #列表
					   random_effect_name=NULL, #列表，不包括永久环境效应			   
					   covariate_effect_name=NULL, #列表
					   random_regression_effect_name=NULL, #列表，随机回归效应   list(c("",""),c("",""))
					   maternal_effect_name=NULL,    #列表		   
					   include_social_effect=FALSE,    #是否评估social effect
					   group_effect_name=NULL,  #评估social effect时, group的名称,自动生成group_phe,												
					   integer_group_names=NULL,        #整型group的名称
					   real_group_names=NULL,            #实型group的名称					   
					   ped_inbred_number=2, #是否考虑近交构建A
					   provided_effect_file_path=NULL, #各个性状的效应文件
					   provided_effect_file_name=NULL,
					   provided_DIR_file_path=NULL,
					   provided_DIR_file_name=NULL,
					   phe_path=NULL,
					   phe_name=NULL,
					   analysis_model=NULL,
					   genetic_effect_name="Id",
					   included_permanent_effect=FALSE, 
					   included_dominance_effect=FALSE,
					   missing_value= -9999,
					   iteration_criteria=1.0e-7,
					   relationship_name=NULL,
					   relationship_path=NULL,
					   dmu_module="dmuai",
					   dmu_algorithm_code=NULL,
					   provided_prior_file_path=NULL,
					   provided_prior_file_name=NULL,
					   integer_n=NULL,  #整型数目
					   residual_cov_trait=NULL,  #限定残差协方差的性状
					   output_result_path=NULL,
					   #O0O0OOO0O0OOOOO0O0OO
					   genetic_effect_number=NULL,
					   residual_average=TRUE,
					   debv_pedigree_path=NULL, #计算debv需要提供系谱
					   debv_pedigree_name=NULL, #计算debv需要提供系谱						
					   cal_debv=FALSE,      #是否计算debv
					   cal_reliability=FALSE,
					   debv_id=NULL,    #计算这些个体的校正表型
					   output_ebv_path=NULL,
					   output_ebv_name=NULL,			   
					   DMU_software_path=ifelse(as.character(Sys.info()["sysname"])=="Linux",system.file("extdata/bin", package = "blupADC"),system.file("extdata/bin_windows", package = "blupADC")),
					   IND_geno_file_name="IND_geno.txt", #
					   IND_geno=NULL,
					   SSBLUP_omega=0.05,
					   plot_variance=FALSE,
					   return_result=FALSE
					   ){
cat("R package:blupADC is only the wrapper of dmu-package in the field of academic research! \n")
cat("The release of pre-installed dmu-package is 5.2 ! \n")
cat("For commercial use of the dmu-package please contact QGG. Email: Per.Madsen@agrsci.dk  !\n")


#检查输入的参数
if(sum(is.null(provided_prior_file_path),is.null(provided_prior_file_name))==1){
stop("Paramters:provided_prior_file_path and provided_prior_file_name couldn't be: the one is non-NULL and the other is NULL")
}


if(is.null(genetic_effect_number)){

if(analysis_model%in%c("User_define","GBLUP_A","GBLUP_AD")){genetic_effect_number=3}
if(analysis_model=="SSBLUP_A"|analysis_model=="PBLUP_A"){genetic_effect_number=4}
}

library(data.table)

if(length(grep(" ",output_result_path))>0){
stop("output_result_path shouldn't contain blank, please change into another path!")
}

output_result_path=paste0(output_result_path,"/",analysis_model,"_",paste(target_trait_name,collapse = "_"))	
if(!file.exists(output_result_path)){
dir.create(output_result_path,recursive=TRUE)
}
	
if(is.null(output_ebv_path)){output_ebv_path=output_result_path}
if(is.null(output_ebv_name)){output_ebv_name="corrected_phe_dmu"}
if(!file.exists(output_ebv_path)){dir.create(output_ebv_path,recursive=TRUE)}			
		

setwd(output_result_path)
output_result_path=getwd()
cat(paste0("Results are saved in path: ",output_result_path,"\n"))
#针对路径含有空格的情况，如果路径存在空格，将表型数据copy 到 输出路径下
cat(paste0("Start the ",analysis_model," analyse of ",length(target_trait_name)," trait model:",paste(target_trait_name,collapse = " & ")," \n"))

if(length(grep(" ",phe_path))>0){
#cat("Attention: Phe path include space!")
file.copy(paste0(phe_path,"/",phe_name),output_result_path)
phe_path=output_result_path
}

if(length(grep(" ",relationship_path))>0){
#cat("Attention: Relationship_path path include space!")

if(file.exists(paste0(relationship_path,"/",IND_geno_file_name))){file.copy(paste0(relationship_path,"/",IND_geno_file_name),output_result_path,recursive=TRUE)}

tmp=relationship_path
for(OOO0OOO0O0OOO0O0OOO0OOOO in relationship_name){
file.copy(paste0(tmp,"/",OOO0OOO0O0OOO0O0OOO0OOOO),output_result_path,recursive=TRUE)
relationship_path=output_result_path
}
}

#################
#*************generate DIR. file **********

if(!is.null(provided_DIR_file_path)&!is.null(provided_DIR_file_name)){

system(paste0("cp -r ",provided_DIR_file_path,"/",provided_DIR_file_name,"  ",output_result_path,"/","Trait.DIR"))

}else{

       DIR=generate_DIR(phe_col_names=phe_col_names,
					    target_trait_name=target_trait_name,
					   fixed_effect_name=fixed_effect_name, #列表
					   random_effect_name=random_effect_name, #列表，不包括永久环境效应
					   covariate_effect_name=covariate_effect_name, #列表
					   random_regression_effect_name=random_regression_effect_name, #列表
					   maternal_effect_name=maternal_effect_name, #列表
					   include_social_effect=include_social_effect, #列表
					   integer_group_names=integer_group_names,
					   real_group_names=real_group_names,
					   group_effect_name=group_effect_name,           
					   ped_inbred_number=ped_inbred_number, #是否考虑近交构建A
					   provided_effect_file_path=provided_effect_file_path, #各个性状的效应文件
					   provided_effect_file_name=provided_effect_file_name,
					   phe_path=phe_path,
					   phe_name=phe_name,
					   analysis_model=analysis_model,
					   genetic_effect_name=genetic_effect_name,
					   included_permanent_effect=included_permanent_effect, 
					   included_dominance_effect=included_dominance_effect,
					   missing_value=missing_value,
					   iteration_criteria=iteration_criteria,
					   relationship_name=relationship_name,
					   relationship_path=relationship_path,
					   dmu_module=dmu_module,
					   dmu_algorithm_code=dmu_algorithm_code,
					   provided_prior_file_path=provided_prior_file_path,
					   provided_prior_file_name=provided_prior_file_name,
					   integer_n=integer_n,  #整型数目
					   residual_cov_trait=residual_cov_trait,  #限定残差协方差的性质
					   output_DIR_path=output_result_path,
					   output_DIR_name="Trait.DIR",
					   IND_geno_file_name=IND_geno_file_name, #
					   IND_geno=IND_geno,
					   SSBLUP_omega=SSBLUP_omega
					   )
}

##################

#################
#***********run DMU1 and DMUai***************
setwd(output_result_path)

#判断操作系统的类型

OOO0OOO0O0OOOOO0OOO0=as.character(Sys.info()["sysname"])

if(OOO0OOO0O0OOOOO0OOO0=="Linux"){
system(paste0(DMU_software_path,"/dmu1  release=5.3 < ",output_result_path,"/Trait.DIR"," >  ",paste(target_trait_name,collapse = "_"),".dmu1.lst"))
system(paste0(DMU_software_path,"/",dmu_module," release=5.3 < ",output_result_path,"/Trait.DIR"," >  ",paste(target_trait_name,collapse = "_"),".",dmu_module,".lst"))
}else{
file.copy(from=paste0(DMU_software_path,"/dmu1.exe"),to=output_result_path)
file.copy(from=paste0(DMU_software_path,"/",dmu_module,".exe"),to=output_result_path)
#生成 .bat文件运行在windows下运行dmu
OOO0OOOOO0OOOOOOO0O0("dmu1")
OOO0OOOOO0OOOOOOO0O0(dmu_module)
system("run_dmu1.bat Trait.DIR")
system(paste0("run_",dmu_module,".bat"," Trait.DIR"))
file.remove(list.files(pattern = "*.exe"))
file.remove(list.files(pattern = "*.bat"))
}
if(file.exists("HINV1"))file.remove("HINV1")
if(file.exists("COR2"))file.remove("COR2")
if(file.exists("COR1"))file.remove("COR1")
if(file.exists("COR3"))file.remove("COR3")
#************finished run DMU1 and DMUai*****************
#################

#将PAROUT 拷贝到输出的路径下
if(!is.null(provided_prior_file_path)&!is.null(provided_prior_file_name)){
file.copy(paste0(provided_prior_file_path,"/",provided_prior_file_name),
            paste0(output_result_path,"/","PAROUT"))
}

#calculate model_reliability  and  corrected_phe 
if("PAROUT_STD"%in%list.files()&(is.null(maternal_effect_name))&(FALSE%in%include_social_effect)){
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(target_trait_name)){

prior=round(cal_dmu_se_reml(target_trait_name=target_trait_name,"PAROUT_STD")$prior[[OOO0OOO0O0OOO0O0OOO0OOOO]],5)
prior_se=round(cal_dmu_se_reml(target_trait_name=target_trait_name,"PAROUT_STD")$prior_se[[OOO0OOO0O0OOO0O0OOO0OOOO]],5)
h2=round(cal_dmu_se_reml(target_trait_name=target_trait_name,"PAROUT_STD")$h2[[OOO0OOO0O0OOO0O0OOO0OOOO]],5)
h2_se=round(cal_dmu_se_reml(target_trait_name=target_trait_name,"PAROUT_STD")$h2_se[[OOO0OOO0O0OOO0O0OOO0OOOO]],5)

data=cbind(c(DIR$random_effect_name[[OOO0OOO0O0OOO0O0OOO0OOOO]],"Residual"),prior,prior_se,h2,h2_se)

colnames(data)=c("Random_effect_name","prior","prior_se","h2","h2_se")
write.table(data,paste0(target_trait_name[OOO0OOO0O0OOO0O0OOO0OOOO],"_heritability_result.txt"),quote=F,row.names=F,col.names=T,sep="\t")
}
}

if(!is.null(provided_prior_file_path) & dmu_module%in%c("dmu4","dmu5")){
system(paste0("cp -r ",provided_prior_file_path,"/",provided_prior_file_name,"  ",output_result_path))
}
#calcualte corrected phenotype
dmu_result=cal_corrected_phe(target_trait_name=target_trait_name,
				   phe_col_names=phe_col_names,
				   genetic_effect_number=genetic_effect_number,
				   genetic_effect_name=genetic_effect_name,
				   dmu_result_path=output_result_path,
				   residual_average=TRUE,
				   dmu_module=dmu_module,
				   cal_reliability=cal_reliability,
				   phe_path=phe_path,  #没有列名 
				   phe_name=phe_name,
				   pedigree_path=debv_pedigree_path, #计算debv需要提供系谱
				   pedigree_name=debv_pedigree_name, #计算debv需要提供系谱
				   provided_prior_file_path=provided_prior_file_path,
				   provided_prior_file_name=provided_prior_file_name,							
				   cal_debv=cal_debv,
				   analysis_model=analysis_model,
				   relationship_name=relationship_name,
				   debv_id=debv_id, #计算有基因型个体的debv
				   output_ebv_path=output_ebv_path,
				   output_ebv_name=output_ebv_name,
				   return_result=return_result
							)

if(plot_variance==TRUE){
message("Plot the result of variance components......")
plot_dmu_blupf90_prior(genetic_effect_name=genetic_effect_name,target_trait_name=target_trait_name,random_effect_name=random_effect_name,output_path=output_result_path)
}

# calculate heritablilty and SE
cat(paste0("Completed the ",analysis_model," analyse of ",length(target_trait_name)," trait model:",paste(target_trait_name,collapse = " & ")," \n"))

if(return_result==TRUE){
return(dmu_result)
}

}


#在 windows下运行DMU
OOO0OOOOO0OOOOOOO0O0<-function(dmu_module="dmu1"){
file=rbind("ECHO OFF",
		   "set CWD=%CD%",
		   "set FN=%1%",
		   "IF EXIST %FN% (",
           paste0("echo Start running ",dmu_module),
		   paste0(dmu_module," < %1 > ",dmu_module,".lst"),
		   ")")
write.table(file,paste0("run_",dmu_module,".bat"),quote=F,row.names=F,col.names=F,sep="")
}



#计算标准误
cal_dmu_se_reml<-function(target_trait_name=target_trait_name,file="PAROUT_STD",file_path=NULL){  #根据方差-协方差矩阵计算标准误

if(!is.null(file_path)){
setwd(file_path)
}

Prior=as.matrix(read.table("PAROUT",header=F))
npar<-as.numeric(read.table(file,nrows=1))

#prior-se
prior_se=read.table(file,skip=1,nrows=npar)

#covariance 
cov=read.table(file,skip=1+npar)
OOO0OOOOO0O0OOO0OOOO=diag(1,npar)    #第一列和第二列对应的是 PRIOR文件的方差组分的行号

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:nrow(cov)){
 OOO0OOOOO0O0OOO0OOOO[cov[OOO0OOO0O0OOO0O0OOO0OOOO,1],cov[OOO0OOO0O0OOO0O0OOO0OOOO,2]]=cov[OOO0OOO0O0OOO0O0OOO0OOOO,4]
 OOO0OOOOO0O0OOO0OOOO[cov[OOO0OOO0O0OOO0O0OOO0OOOO,2],cov[OOO0OOO0O0OOO0O0OOO0OOOO,1]]=cov[OOO0OOO0O0OOO0O0OOO0OOOO,4]
}

h2=NULL;h2_se=NULL;O0O0OOOOOOOOO0OOOOOO=NULL;O0OOOOO0O0OOO0O0OOO0=NULL;

for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:length(target_trait_name)){

OOO0OOOOOOOOOOOOO0OO_pos=(1:nrow(Prior))[Prior[,2]%in%OOO0OOO0O0OOO0O0OOO0OOOO & Prior[,3]%in%OOO0OOO0O0OOO0O0OOO0OOOO]
OOO0OOOOOOOOOOOOO0OO=as.matrix(Prior[OOO0OOOOOOOOOOOOO0OO_pos,4])
O0O0O0OOOOOOO0O0O0OO=OOO0OOOOO0O0OOO0OOOO[OOO0OOOOOOOOOOOOO0OO_pos,OOO0OOOOOOOOOOOOO0OO_pos]

Var=O0O0O0OOOOOOO0O0O0OO #信息矩阵的逆矩阵为 方差-协方差矩阵
heritablilty=OOO0OOOOOOOOOOOOO0OO/sum(OOO0OOOOOOOOOOOOO0OO[,1])
SE_h2=matrix(NA,nrow=nrow(OOO0OOOOOOOOOOOOO0OO),ncol=1)
for(OOO0OOO0O0OOO0O0OOO0OOOO in 1:nrow(OOO0OOOOOOOOOOOOO0OO)){
derivate=NULL
for(OOO0OOO0O0OOO0O0OOO0OOOOO in 1:nrow(OOO0OOOOOOOOOOOOO0OO)){
if(OOO0OOO0O0OOO0O0OOO0OOOOO==OOO0OOO0O0OOO0O0OOO0OOOO){derivate=cbind(derivate,(sum(OOO0OOOOOOOOOOOOO0OO[,1])-OOO0OOOOOOOOOOOOO0OO[OOO0OOO0O0OOO0O0OOO0OOOO,1])/(sum(OOO0OOOOOOOOOOOOO0OO[,1])^2) )}
if(OOO0OOO0O0OOO0O0OOO0OOOOO!=OOO0OOO0O0OOO0O0OOO0OOOO){derivate=cbind(derivate,(-OOO0OOOOOOOOOOOOO0OO[OOO0OOO0O0OOO0O0OOO0OOOO,1])/(sum(OOO0OOOOOOOOOOOOO0OO[,1])^2 ))}
}
SE_h2[OOO0OOO0O0OOO0O0OOO0OOOO,1]=sqrt(derivate%*%Var%*%t(derivate)) # 开平方根才是最终结果
}

O0O0OOOOOOOOO0OOOOOO=c(O0O0OOOOOOOOO0OOOOOO,list(OOO0OOOOOOOOOOOOO0OO))
O0OOOOO0O0OOO0O0OOO0=c(O0OOOOO0O0OOO0O0OOO0,list(prior_se[OOO0OOOOOOOOOOOOO0OO_pos,3]))
h2=c(h2,list(heritablilty))
h2_se=c(h2_se,list(SE_h2))
}

return(list(prior=O0O0OOOOOOOOO0OOOOOO,prior_se=O0OOOOO0O0OOO0O0OOO0,h2=h2,h2_se=h2_se))
}
