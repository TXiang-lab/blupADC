
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
					   #cor_phe
					   genetic_effect_number=NULL,
					   residual_average=TRUE,
					   debv_pedigree_path=NULL, #计算debv需要提供系谱
					   debv_pedigree_name=NULL, #计算debv需要提供系谱						
					   cal_debv=FALSE,      #是否计算debv
					   cal_reliability=FALSE,
					   debv_id=NULL,    #计算这些个体的校正表型
					   output_ebv_path=NULL,
					   output_ebv_name=NULL,			   
					   DMU_software_path=ifelse(as.character(Sys.info()["sysname"])=="Linux",
																system.file("extdata/bin_linux", package = "blupADC"),
																ifelse(as.character(Sys.info()["sysname"])=="Windows",
																system.file("extdata/bin_windows", package = "blupADC"),
																system.file("extdata/bin_mac", package = "blupADC")
																)),
					   IND_geno_file_name="IND_geno.txt", #
					   IND_geno=NULL,
					   SSBLUP_omega=0.05,
					   plot_variance=FALSE,
					   return_result=FALSE,
					   EM_iter=10,
					   dmu5_iter=1000000,
					   gibbs_n_rounds=10000,
					   gibbs_n_burnin=1000,
					   gibbs_n_gaps=5,
					   gibbs_seeds=c(1994,1998)
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
for(i in relationship_name){
file.copy(paste0(tmp,"/",i),output_result_path,recursive=TRUE)
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
					   SSBLUP_omega=SSBLUP_omega,
					   EM_iter=EM_iter,
					   dmu5_iter=dmu5_iter,
					   gibbs_n_rounds=gibbs_n_rounds,
					   gibbs_n_burnin=gibbs_n_burnin,
					   gibbs_n_gaps=gibbs_n_gaps,
					   gibbs_seeds=gibbs_seeds
					   )
}

##################

#################
#***********run DMU1 and DMUai***************
setwd(output_result_path)

#判断操作系统的类型

os_type=as.character(Sys.info()["sysname"])

if(os_type!="Windows"){
system(paste0(DMU_software_path,"/dmu1  release=5.3 < ",output_result_path,"/Trait.DIR"," >  ",paste(target_trait_name,collapse = "_"),".dmu1.lst"))
system(paste0(DMU_software_path,"/",dmu_module," release=5.3 < ",output_result_path,"/Trait.DIR"," >  ",paste(target_trait_name,collapse = "_"),".",dmu_module,".lst"))
}else{
file.copy(from=paste0(DMU_software_path,"/dmu1.exe"),to=output_result_path)
file.copy(from=paste0(DMU_software_path,"/",dmu_module,".exe"),to=output_result_path)
#生成 .bat文件运行在windows下运行dmu
generate_dmu_bat("dmu1")
generate_dmu_bat(dmu_module)
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

h2_data=NULL
#将PAROUT 拷贝到输出的路径下
if(!is.null(provided_prior_file_path)&!is.null(provided_prior_file_name)){
file.copy(paste0(provided_prior_file_path,"/",provided_prior_file_name),
            paste0(output_result_path,"/","PAROUT"))
}

#calculate model_reliability  and  corrected_phe 
if("PAROUT_STD"%in%list.files()&(is.null(maternal_effect_name))&(FALSE%in%include_social_effect)){
for(i in 1:length(target_trait_name)){

prior=round(cal_dmu_se_reml(target_trait_name=target_trait_name,"PAROUT_STD")$prior[[i]],5)
prior_se=round(cal_dmu_se_reml(target_trait_name=target_trait_name,"PAROUT_STD")$prior_se[[i]],5)
h2=round(cal_dmu_se_reml(target_trait_name=target_trait_name,"PAROUT_STD")$h2[[i]],5)
h2_se=round(cal_dmu_se_reml(target_trait_name=target_trait_name,"PAROUT_STD")$h2_se[[i]],5)

h2_data=cbind(c(DIR$random_effect_name[[i]],"Residual"),prior,prior_se,h2,h2_se)

colnames(h2_data)=c("Random_effect_name","prior","prior_se","h2","h2_se")
write.table(h2_data,paste0(target_trait_name[i],"_heritability_result.txt"),quote=F,row.names=F,col.names=T,sep="\t")
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
return(list(dmu_result,h2_data,DIR))
}

}


#在 windows下运行DMU
generate_dmu_bat<-function(dmu_module="dmu1"){
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
cov_matrix=diag(1,npar)    #第一列和第二列对应的是 PRIOR文件的方差组分的行号

for(i in 1:nrow(cov)){
 cov_matrix[cov[i,1],cov[i,2]]=cov[i,4]
 cov_matrix[cov[i,2],cov[i,1]]=cov[i,4]
}

h2=NULL;h2_se=NULL;total_trait_prior=NULL;trait_prior_se=NULL;

for(i in 1:length(target_trait_name)){

Trait_prior_pos=(1:nrow(Prior))[Prior[,2]%in%i & Prior[,3]%in%i]
Trait_prior=as.matrix(Prior[Trait_prior_pos,4])
Trait_cov_matrix=cov_matrix[Trait_prior_pos,Trait_prior_pos]

Var=Trait_cov_matrix #信息矩阵的逆矩阵为 方差-协方差矩阵
heritablilty=Trait_prior/sum(Trait_prior[,1])
SE_h2=matrix(NA,nrow=nrow(Trait_prior),ncol=1)
for(i in 1:nrow(Trait_prior)){
derivate=NULL
for(j in 1:nrow(Trait_prior)){
if(j==i){derivate=cbind(derivate,(sum(Trait_prior[,1])-Trait_prior[i,1])/(sum(Trait_prior[,1])^2) )}
if(j!=i){derivate=cbind(derivate,(-Trait_prior[i,1])/(sum(Trait_prior[,1])^2 ))}
}
SE_h2[i,1]=sqrt(derivate%*%Var%*%t(derivate)) # 开平方根才是最终结果
}

total_trait_prior=c(total_trait_prior,list(Trait_prior))
trait_prior_se=c(trait_prior_se,list(prior_se[Trait_prior_pos,3]))
h2=c(h2,list(heritablilty))
h2_se=c(h2_se,list(SE_h2))
}

return(list(prior=total_trait_prior,prior_se=trait_prior_se,h2=h2,h2_se=h2_se))
}


Multitasks_run_DMU_BLUPF90<-function(tasks){
cat("Perform multitasks of DMU simultaneously! \n")
library(future)
plan(cluster, workers = length(tasks))

final_tasks=NULL
for(i in 1:length(tasks)){
final_tasks=c(final_tasks,future({eval(tasks[i])} ))
}
future::value(final_tasks)
}