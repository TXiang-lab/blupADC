#' RUN BLUPF90 software in R
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
								   maternal_effect_option=NULL, # pe mat mpe
								   mat_cov=NULL,
								   mat_pe_cov=NULL,
								   mat_mpe_cov=NULL,
								   phe_path=NULL,
								   phe_name=NULL,
								   random_regression_effect_name=NULL,
								   reg_pe_cov=NULL,
								   reg_gen_cov=NULL,
								   user_file_id=NULL, # user_file_id, 
								   output_result_path=NULL,
								   output_ebv_path=NULL,
						            output_ebv_name=NULL,
								   missing_value="-9999",
								   relationship_name=NULL,
								   relationship_path=NULL,
								   analysis_model="PBLUP_A",
								   BLUPF90_algorithm="AI_REML",
								   gibbs_sampler="thrgibbs1f90",
								   genetic_effect_name="Id",
								   included_permanent_effect=FALSE,
								   provided_BLUPF90_prior_file_path=NULL,
								   provided_BLUPF90_prior_file_name=NULL,
								   provided_BLUPF90_prior_effect_name=NULL, #随机效应的名称, 包括Residual
					                 provided_renf90_par_file_path=NULL,
					                 provided_renf90_par_file_name=NULL,
								  cal_se_reml=FALSE, #if output the se of heritablilty
								   BLUPF90_alt_option=NULL,
								   BLUPF90_genumeric_name=NULL,
								   BLUPF90_map_name=NULL,
								   plot_variance=FALSE,
								   BLUPF90_software_path=ifelse(as.character(Sys.info()["sysname"])=="Linux",
																system.file("extdata/bin_linux", package = "blupADC"),
																ifelse(as.character(Sys.info()["sysname"])=="Windows",
																system.file("extdata/bin_windows", package = "blupADC"),
																system.file("extdata/bin_mac", package = "blupADC")
																))

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

for(i in 1:length(relationship_name)){
file.copy(from=paste0(relationship_path,"/",relationship_name[i]),to=output_result_path,overwrite=TRUE)
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
				maternal_effect_option=maternal_effect_option, # pe mat mpe
				mat_cov=mat_cov,
				mat_pe_cov=mat_pe_cov,
				mat_mpe_cov=mat_mpe_cov,
				random_regression_effect_name=random_regression_effect_name,
				reg_pe_cov=reg_pe_cov,
				reg_gen_cov=reg_gen_cov,				
			    included_permanent_effect=included_permanent_effect,
			    provided_BLUPF90_prior_file_path=provided_BLUPF90_prior_file_path,
			    provided_BLUPF90_prior_file_name=provided_BLUPF90_prior_file_name,
			    provided_BLUPF90_prior_effect_name=provided_BLUPF90_prior_effect_name, #随机效应的名称, 包括Residual
				BLUPF90_alt_option=BLUPF90_alt_option,
			    BLUPF90_genumeric_name=BLUPF90_genumeric_name,
			    BLUPF90_map_name=BLUPF90_map_name)

cat("Start running renumf90 module of BLUPF90......\n")
system2(paste0(BLUPF90_software_path,"/renumf90"),"renum.par",stdout="renumf90.log")

if(analysis_model=="User_define"){
#rename renf90.dat 
phe=read.table(phe_name,header=F);colnames(phe)=phe_col_names
renf90.phe=read.table("renf90.dat",header=F)
renf90.phe[,(ncol(renf90.phe)-length(relationship_name)+1):ncol(renf90.phe)]=phe[,genetic_effect_name]
write.table(renf90.phe,"renf90.dat",quote=F,row.names=F,sep=" ",col.names=F)

#modify renf90.par

renf90=data.table::fread("renf90.par",data.table=F,fill=T)
pos=match("add_animal",renf90[,1])
renf90[pos,1]="user_file"
renf90[pos+2,1]=relationship_name[1]
renf90[is.na(renf90)]=""

renf90=data.table::fread("renf90.par",data.table=F,fill=T)
pos=(1:nrow(renf90))[renf90[,1]%in%"add_animal"]
renf90[pos,1]="user_file"
renf90[pos+2,1]=relationship_name
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
}else if(BLUPF90_algorithm=="Gibbs"){
cat("Start running Gibbs module of BLUPF90......\n")
cat(paste0("Using gibbs_sampler: ",gibbs_sampler," to estimate variace components......\n"))
cat("Please enter the following three papamters in each line! \n")
cat("1. Number of samples \n")
cat("2. Length of burn-in \n")
cat("3. Give n to store every n-th sample? (1 means store all samples) \n")
system2(paste0(BLUPF90_software_path,"/",gibbs_sampler),"renf90.par",stdout=paste0(gibbs_sampler,".log"))
cat(paste0("Using gibbs_sampler: ",gibbs_sampler," to estimate variace components......\n"))
cat("Complete running Gibbs module of BLUPF90 \n")
cat("Please enter the following three papamters in each line! \n")
cat("1. Additional number of Burn-in \n")
cat("2. Give n to store every n-th sample? (1 means store all samples) \n")
cat("3. Enter 0 for exiting \n")
system2(paste0(BLUPF90_software_path,"/postgibbsf90"),"renf90.par",stdout=paste0("postgibbsf90.log"))
cat("Complete postgibbsf90 analysis! \n")
}

if(!"Gibbs"%in%BLUPF90_algorithm){
#计算校正表型
cor_phe=cal_corrected_phe_BLUPF90(target_trait_name=target_trait_name)
cor_phe=cor_phe[order(cor_phe[,1]),]
if(is.null(output_ebv_path)){output_ebv_path=output_result_path}
if(is.null(output_ebv_name)){output_ebv_name="corrected_phe_BLUPF90"}

cat("Output the corrected phenotype......\n")
setwd(output_ebv_path)
utils::write.table(cor_phe,paste0("colnames_",output_ebv_name,".txt"),quote=F,row.names=F,col.names=T,sep=" ")

cor_phe[is.na(cor_phe)]=missing_value
setwd(output_ebv_path)
utils::write.table(cor_phe,paste0(output_ebv_name,".txt"),quote=F,row.names=F,col.names=F,sep=" ")


#统计方差组分结果-遗传力及标准误

if(cal_se_reml==TRUE){
cal_blupf90_se_reml(target_trait_name=target_trait_name,random_effect_name=random_effect_name,genetic_effect_name=genetic_effect_name,BLUPF90_algorithm=BLUPF90_algorithm)
}

if(plot_variance==TRUE){
message("Plot the result of variance components......")
plot_dmu_blupf90_prior(genetic_effect_name=genetic_effect_name,target_trait_name=target_trait_name,random_effect_name=random_effect_name,output_path=output_result_path)
}
}
}


cal_corrected_phe_BLUPF90<-function(target_trait_name=NULL){

ebv=data.table::fread("solutions",data.table=F,header=F)
colnames(ebv)=c("Trait_number","Effect_number","Id","EBV","SE_EBV")
max_effect_number=max(ebv[,"Effect_number"])
ebv=ebv[ebv[,"Effect_number"]%in%max_effect_number,]


residual=data.table::fread("yhat_residual",data.table=F,header=F)
data_original_phe=data.table::fread("renf90.dat",data.table=F,header=F)
residual$Id=data_original_phe[,ncol(data_original_phe)]


cor_phe=matrix(NA,nrow=length(unique(ebv[,"Id"])),ncol=1+4*length(target_trait_name))
cor_phe=as.data.frame(cor_phe,stringsAsFactors=F)
cor_phe[,1]=sort(unique(ebv[,"Id"]));colnames(cor_phe)[1]="Id"

for(i in 1:length(target_trait_name)){

trait_ebv=ebv[ebv$Trait_number%in%i,]
cor_phe[,1+(i-1)*4+1]= trait_ebv[match(cor_phe[,1],trait_ebv[,"Id"]),"EBV"]
cor_phe[,1+(i-1)*4+2]= residual[match(cor_phe[,1],residual[,"Id"]),length(target_trait_name)+i]
cor_phe[,1+(i-1)*4+3]=cor_phe[,1+(i-1)*4+1]+cor_phe[,1+(i-1)*4+2]

colnames(cor_phe)[1+(i-1)*4+1]=paste0(target_trait_name[i],"_EBV")
colnames(cor_phe)[1+(i-1)*4+2]=paste0(target_trait_name[i],"_Res")
colnames(cor_phe)[1+(i-1)*4+3]=paste0(target_trait_name[i],"_EBV_Plus_Res")
colnames(cor_phe)[1+(i-1)*4+4]=paste0(target_trait_name[i],"_dEBV")

}


#将rename 后的名称映射回来
rename_ped_name=list.files(getwd(),pattern = '.ped$')
if(length(rename_ped_name)>1){
rename_ped_name=paste0("renadd0",max_effect_number,".ped")
}
BLUPF90_key=data.table::fread(rename_ped_name,data.table=F,header=F)
cor_phe[,"Id"]=BLUPF90_key[match(cor_phe[,"Id"],BLUPF90_key[,1]),10]

return(cor_phe)


}




cal_blupf90_se_reml<-function(target_trait_name=NULL,random_effect_name=NULL,genetic_effect_name="Id",BLUPF90_algorithm=NULL){

cal_se_reml<-function(infor_matrix=NULL,infor_matrix_inv=NULL,Prior=NULL){  #根据信息矩阵计算标准误

if(is.null(infor_matrix_inv)){
Var=my.solve(infor_matrix)  #信息矩阵的逆矩阵为 方差-协方差矩阵
}else{
Var=infor_matrix_inv
}

heritablilty=Prior/sum(Prior[,1])
SE_h2=matrix(NA,nrow=nrow(Prior),ncol=1)
for(i in 1:nrow(Prior)){

derivate=NULL
for(j in 1:nrow(Prior)){

if(j==i){derivate=cbind(derivate,(sum(Prior[,1])-Prior[i,1])/(sum(Prior[,1])^2) )}
if(j!=i){derivate=cbind(derivate,(-Prior[i,1])/(sum(Prior[,1])^2 ))}

}
SE_h2[i,1]=sqrt(derivate%*%Var%*%t(derivate)) # 开平方根才是最终结果
}
return(SE_h2)
}

temp=fread("AI_inv.txt",data.table=F)
total_infor_matrix_inv=diag(max(temp[,1]))
for(i in 1:nrow(temp)){
total_infor_matrix_inv[temp[i,1],temp[i,2]]=temp[i,7]
total_infor_matrix_inv[temp[i,2],temp[i,1]]=temp[i,7]
}

if(BLUPF90_algorithm=="AI_REML"){
log=fread("airemlf90.log",data.table=F,fill = T)
}else if(BLUPF90_algorithm=="EM_REML"){
log=fread("emremlf90.log",data.table=F,fill = T)
}else if(BLUPF90_algorithm=="BLUP"){
log=fread("blupf90.log",data.table=F,fill = T)
}

#分性状计算
for(i in 1:length(target_trait_name)){

infor_matrix_inv=total_infor_matrix_inv[seq(i,nrow(total_infor_matrix_inv),length(target_trait_name)),seq(i,nrow(total_infor_matrix_inv),length(target_trait_name))]
infor_matrix_inv=infor_matrix_inv[!diag(infor_matrix_inv)==0,!diag(infor_matrix_inv)==0]

Prior_pos=(1:nrow(log))[log[,1]%in%c("Genetic","Residual")]+i
Prior=matrix(as.numeric(log[Prior_pos,i]))

Prior_se_pos=(1:nrow(log))[log[,1]%in%c("SE")]+i
Prior_se=matrix(as.numeric(log[Prior_se_pos,i]))

Prior_se=Prior_se[!Prior[,1]==1,]
Prior=Prior[!Prior[,1]==1,] #去除冗余的随机效应(e.g. 性状没有该随机效应，但是BLUPF90会自动赋予该效应的方差为1)
Prior=as.matrix(Prior)

h2_se=cal_se_reml(infor_matrix_inv=infor_matrix_inv,Prior=Prior)
h2=Prior/sum(Prior[,1])

temp_random_effect_name=random_effect_name[[i]]
temp_random_effect_name=sort(temp_random_effect_name[!temp_random_effect_name%in%genetic_effect_name])
temp_random_effect_name=c(temp_random_effect_name,"Id","Residual")

data_h2=cbind(temp_random_effect_name,Prior,Prior_se,h2,h2_se)
colnames(data_h2)=c("Random_effect_name","prior","prior_se","h2","h2_se")
write.table(data_h2,paste0(target_trait_name[i],"_heritability_result.txt"),quote=F,row.names=F,col.names=T,sep="\t")

}

}




