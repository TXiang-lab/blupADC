
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
Trait_n=length(target_trait_name)
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
colnames(SOL)=c("Genetic_number","Random_number","Trait","Id","EBV","SE")
tmp_SOL=SOL
additive_id=unique(SOL[(SOL$Genetic_number%in%genetic_effect_number)&(SOL$Random_number%in%(1:length(relationship_name))),c("Id")])
cor_phe=matrix(NA,nrow=length(additive_id),ncol=1+(2+length(relationship_name))*Trait_n)
cor_phe=as.data.frame(cor_phe,stringsAsFactors=F)
cor_phe$V1=additive_id
colnames(cor_phe)[1]="Id"
RESIDUAL=fread(paste0(dmu_result_path,"/RESIDUAL"),heade=F,data.table=F)

for(i in 1:length(target_trait_name)){
residual=RESIDUAL[,c(1,1+i+2*Trait_n)];colnames(residual)=c("Pos","RESIDUAL")
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
RES_pos=match(cor_phe[,1],residual[,1])

for(j in 1:length(relationship_name)){
additive=SOL[(SOL$Genetic_number%in%genetic_effect_number)&(SOL$Trait%in%i)&(SOL$Random_number%in%j),c("Id","EBV")]
EBV_pos=match(cor_phe[,1],additive[,1])
cor_phe[,(i-1)*(length(relationship_name)+2)+j+1]=additive[EBV_pos,"EBV"]
cor_phe[,(i-1)*(length(relationship_name)+2)+j+2]=residual[RES_pos,"RESIDUAL"]
cor_phe[,(i-1)*(length(relationship_name)+2)+j+3]=cor_phe[,(i-1)*(length(relationship_name)+2)+j+2]+cor_phe[,(i-1)*(length(relationship_name)+2)+j+1]
colnames(cor_phe)[(i-1)*(length(relationship_name)+2)+j+1]=paste0(target_trait_name[i],"_R",j,"_Value")
colnames(cor_phe)[(i-1)*(length(relationship_name)+2)+j+2]=paste0(target_trait_name[i],"_Residual")
colnames(cor_phe)[(i-1)*(length(relationship_name)+2)+j+3]=paste0(target_trait_name[i],"_R",j,"_Value_Plus_Residual")

}
#cor_phe[,i*(length(relationship_name)+2)+1]= residual[RES_pos,"RESIDUAL"]
#cor_phe[,i*(length(relationship_name)+2)+2]= cor_phe[,i*(length(relationship_name)+2)+0]+cor_phe[,i*(length(relationship_name)+2)+1]
#colnames(cor_phe)[i*(length(relationship_name)+2)+1]=paste0(target_trait_name[i],"_Residual")
#colnames(cor_phe)[(i-1)*(length(relationship_name)+2)+j+2]=paste0(target_trait_name[i],"_R",j,"_Value_Plus_Residual")
}


# 常规模型-计算EBV
}else{
SOL=fread(paste0(dmu_result_path,"/SOL"),heade=F,data.table=F)
SOL=SOL[,c(1,2,5,8,9)];
colnames(SOL)=c("Genetic_number","Trait","Id","EBV","SE")
tmp_SOL=SOL
additive_id=unique(SOL[(SOL$Genetic_number%in%genetic_effect_number),c("Id")])
cor_phe=matrix(NA,nrow=length(additive_id),ncol=1+4*Trait_n)
cor_phe=as.data.frame(cor_phe,stringsAsFactors=F)
cor_phe$V1=additive_id
colnames(cor_phe)[1]="Id"



RESIDUAL=fread(paste0(dmu_result_path,"/RESIDUAL"),heade=F,data.table=F)



for(i in 1:length(target_trait_name)){

additive=SOL[(SOL$Genetic_number%in%genetic_effect_number)&(SOL$Trait%in%i),c("Id","EBV")]
residual=RESIDUAL[,c(1,1+i+2*Trait_n)];colnames(residual)=c("Pos","RESIDUAL")
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
EBV_pos=match(cor_phe[,1],additive[,1])
RES_pos=match(cor_phe[,1],residual[,1])

cor_phe[,1+(i-1)*4+1]=   additive[EBV_pos,"EBV"]
cor_phe[,1+(i-1)*4+2]= residual[RES_pos,"RESIDUAL"]
cor_phe[,1+(i-1)*4+3]=cor_phe[,1+(i-1)*3+i]+cor_phe[,1+(i-1)*3+i+1]

colnames(cor_phe)[1+(i-1)*4+1]=paste0(target_trait_name[i],"_EBV")
colnames(cor_phe)[1+(i-1)*4+2]=paste0(target_trait_name[i],"_Res")
colnames(cor_phe)[1+(i-1)*4+3]=paste0(target_trait_name[i],"_EBV_Plus_Residual")
colnames(cor_phe)[1+(i-1)*4+4]=paste0(target_trait_name[i],"_dEBV")
}


#计算debv 
if(cal_debv==TRUE){
if(is.null(pedigree_path)&is.null(pedigree_name)){
stop("Please provide pedigree path and name for calculating dEBV!")
}
cat("Output dEBV......\n")
debv_result=calculate_debv(target_trait_name=target_trait_name,
				   genetic_effect_number=genetic_effect_number,
				   genetic_effect_name=genetic_effect_name,
				   dmu_result_path=dmu_result_path,
				   pedigree_path=pedigree_path,
				   pedigree_name=pedigree_name,							
				   provided_prior_file_path=provided_prior_file_path,
				   provided_prior_file_name=provided_prior_file_name,
				   debv_id=debv_id #计算有基因型个体的debv
				   )

for(i in 1:length(target_trait_name)){
cor_phe[,1+(i-1)*4+i+3]=as.numeric(debv_result[match(cor_phe[,1],debv_result[,1]),i+1])
}
original_cor_phe=cor_phe #带有NA的
}
}

if(cal_reliability==TRUE){

Sum_reliability=NULL

if(dmu_module=="dmu5"){stop("dmu5 module coundn't be used to estimate reliability,please use dmuai module!")}
if(!file.exists("PAROUT")){stop("PAROUT file is needed to estimate reliability,please use dmuai module!")}
cat("Output trait reliability result! \n")
vc_value<-data.frame(read.table("PAROUT",header=F))
se=tmp_SOL[tmp_SOL[,1]%in%genetic_effect_number,]
se=se[,c("Id","Trait","EBV","SE")]
se=se[!se[,4]%in%0,]  #去除为0的SE

for(i in 1:length(target_trait_name)){
trait_se=se[se[,2]%in%i,]
trait_var_Additive=vc_value[vc_value[,1]%in%1&vc_value[,2]%in%i&vc_value[,3]%in%i,4]
reliability=round(mean(1-(trait_se[,"SE"]^2)/trait_var_Additive),4)
Sum_reliability=rbind(Sum_reliability,c(target_trait_name[[i]],reliability))
}

colnames(Sum_reliability)=c("Trait","Reliability")
write.table(Sum_reliability,paste("Trait_reliability.txt",sep=""),row.names=F,col.names=T,quote=F)
}


if(!is.null(output_ebv_path)&!is.null(output_ebv_name)){
setwd(output_ebv_path)
write.table(cor_phe,paste0("colnames_",output_ebv_name,".txt"),quote=F,row.names=F,sep=" ")

cor_phe[is.na(cor_phe)]="-9999"
write.table(cor_phe,paste0(output_ebv_name,".txt"),quote=F,row.names=F,col.names=F,sep=" ")
}

if(return_result==TRUE){return(original_cor_phe)}
}


#calcualte corrected phenotype
calculate_debv_old<-function(target_trait_name=NULL,
				   genetic_effect_number=NULL,
				   genetic_effect_name="Id",
				   dmu_result_path=NULL,
				   pedigree_path=NULL,
				   pedigree_name=NULL,							
				   provided_prior_file_path=NULL,
				   provided_prior_file_name=NULL,
				   debv_id=NULL #计算有基因型个体的debv
				   ){

ebv_reliability=calculate_reliability(return_individual_result=TRUE,
					    target_trait_name=target_trait_name,
					    genetic_effect_number=genetic_effect_number,
						provided_prior_file_path=provided_prior_file_path,
						provided_prior_file_name=provided_prior_file_name,
						dmu_result_path=dmu_result_path,
						debv_id=NULL)

deregressed_ebv=matrix(nrow=length(debv_id),ncol=length(target_trait_name)+1)
deregressed_ebv=as.data.frame(deregressed_ebv,stringsAsFactors=F)
deregressed_ebv[,1]=debv_id

for(i in 1:length(target_trait_name)){

ebv_r2=ebv_reliability$ebv_r2[[i]]
h2=ebv_reliability$h2[[i]]
pedigree=fread(paste0(pedigree_path,"/",pedigree_name),header=F,data.table=F)

pedigree=pedigree[,1:3]
colnames(pedigree)=c("ID","Sire","Dam")
p.varSNP=0.5

debv_id=data.frame(debv_id,stringsAsFactors=F);

colnames(debv_id)[1]="Id"

PedTraits <- merge(x=ebv_r2,y=pedigree,by.x=1,by.y=1)
  PedTraitsgenoIDs <- merge(x=PedTraits,y=debv_id,by.x=1,by.y=1)
  colnames(x=PedTraitsgenoIDs) <- c("ID","ID_EBV","ID_Rel","SireID","DamID") 
  Rel.sire <- merge(x=PedTraits,y=PedTraitsgenoIDs,by.x=1,by.y=4)
  Rel.sire <- Rel.sire[,-4:-5]  
  colnames(Rel.sire)[1:3] <- c('SireID','Sire_EBV','Sire_R2')
  Rel.dam <- merge(x=PedTraits,y=Rel.sire,by.x=1,by.y=7)
  Rel.dam <- Rel.dam[,-4:-5]  
  colnames(Rel.dam)[1:3] <- c('DamID','Dam_EBV','Dam_R2')
  data <- Rel.dam[,c('ID','ID_EBV','ID_Rel','SireID','Sire_EBV','Sire_R2','DamID','Dam_EBV','Dam_R2')]
  dEBV <- data[,-c(1,4,7)]
  dEBV$h2 <- h2
  dEBV$p.varSNP <- 1-p.varSNP
  dEBV=as.matrix(dEBV)
  dEBV=apply(dEBV,2,as.numeric)
  Deregress <- t(apply(X=dEBV,MARGIN=1,FUN=Debv_Garrick))
  Deregress <- cbind.data.frame(Id=data[,1],Deregress)

  deregressed_ebv[na.omit(match(deregressed_ebv[,1],Deregress[,1])),i+1]=Deregress[na.omit(match(deregressed_ebv[,1],Deregress[,1])),"dEBV"]
}

colnames(deregressed_ebv)=c("Id",paste0(target_trait_name,"_dEBV"))
return(deregressed_ebv)

}

#calcualte corrected phenotype
calculate_debv<-function(target_trait_name=NULL,
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

ebv_reliability=calculate_reliability(return_individual_result=TRUE,
					    target_trait_name=target_trait_name,
					    genetic_effect_number=genetic_effect_number,
						provided_prior_file_path=provided_prior_file_path,
						provided_prior_file_name=provided_prior_file_name,
						dmu_result_path=dmu_result_path,
						debv_id=NULL)

pedigree=fread(paste0(pedigree_path,"/",pedigree_name),header=F,data.table=F)
pedigree=pedigree[,1:3]
colnames(pedigree)=c("ID","Sire","Dam")
p_var=p_var


ebv_r2=data.frame(ebv_reliability$ebv_r2[[1]],stringsAsFactors=F)

h2=ebv_reliability$h2[[1]]

debv=debv_method1(ebv=ebv_r2,ped=pedigree,h2=h2,p_var=p_var,debv_id=debv_id)[,c("Id","dEBV")]

if(length(target_trait_name)>=2){
for(i in 2:length(target_trait_name)){

ebv_r2=data.frame(ebv_reliability$ebv_r2[[i]],stringsAsFactors=F)
h2=ebv_reliability$h2[[i]]

debv_trait=debv_method1(ebv=ebv_r2,ped=pedigree,h2=h2,p_var=p_var,debv_id=debv_id)[,c("Id","dEBV")]

debv=cbind(debv,debv_trait[match(debv[,1],debv_trait[,1]),2])

}
}


return(debv)

}

debv_method1<-function(ebv,ped,h2=1,p_var=0.5,debv_id=NULL){
 Debv_Garrick <- function (ebv_mat){
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

Deregress <- t(apply(X=as.matrix(ped[,-c(1:3)]),MARGIN=1,FUN=Debv_Garrick))
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
colnames(SOL)=c("Genetic_number","Trait","Id","EBV","SE")
SOL=SOL[SOL$Genetic_number%in%genetic_effect_number,]

if(is.null(provided_prior_file_path)&is.null(provided_prior_file_name)){

vc_value<-data.frame(read.table("PAROUT",header=F),stringsAsFactors=F)
}else{

vc_value<-data.frame(read.table(paste0(provided_prior_file_path,"/",provided_prior_file_name),header=F),stringsAsFactors=F)}

for(i in 1:length(target_trait_name)){

trait_se=SOL[SOL$Trait%in%i,]
trait_var_Additive=as.numeric(vc_value[(vc_value[,1]%in%1)&(vc_value[,2]%in%i)&(vc_value[,3]%in%i),4])

reliability=mean(1-(trait_se[,"SE"]^2)/trait_var_Additive)
r2=c(r2,list(reliability))
write.table(reliability,paste0("The_reliability_of_",target_trait_name[i],".txt"),row.names=F,col.names=F,quote=F)

ebv_r2=cbind(SOL[,"Id"],SOL[,"EBV"],(1-(trait_se[,"SE"]^2)/trait_var_Additive))
h2=c(h2,list(trait_var_Additive/sum(vc_value[,4])))
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

		