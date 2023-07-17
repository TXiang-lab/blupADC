
run_DMU<-function(
	                   phe_col_names=NULL,
					target_trait_name=NULL,
					fixed_effect_name=NULL,  #list 
					random_effect_name=NULL, #list		   
					covariate_effect_name=NULL, #list
					random_regression_effect_name=NULL, #list, e.g. list(c("",""),c("",""))
					maternal_effect_name=NULL,    #list		   
					include_social_effect=FALSE, #whether including social effect in genetic evaluation
					group_effect_name=NULL,       #for social genetic effect evaluation 											
					integer_group_names=NULL,    #integer name for genetic group
					real_group_names=NULL,       #real name for genetic group	   
					ped_inbred_number=2,         #whether consider inbreeding 
					provided_effect_file_path=NULL, #
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
					integer_n=NULL,  #number of integer varable in phenotype
					residual_cov_trait=NULL,  #resitricted-residual trait 
					output_result_path=NULL,
					#cor_phe
					genetic_effect_number=NULL,
					residual_average=TRUE,
					debv_pedigree_path=NULL, 
					debv_pedigree_name=NULL, 			
					cal_debv=FALSE,     
					cal_reliability=FALSE,
					debv_id=NULL,   
					output_ebv_path=NULL,
					output_ebv_name=NULL,			   
					DMU_software_path=blupADC_software_path(),
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
					gibbs_seeds=c(1994,1998),
					#newly added parameters for the new feature in blupADC
					phe_file=NULL,
					kinship_file=NULL,
					prior_file=NULL,
					show_message=TRUE
					   ){
if(show_message)cat("R package:blupADC is only the wrapper of dmu-package in the field of academic research! \n")
if(show_message)cat("The release of pre-installed dmu-package is 5.2 ! \n")
if(show_message)cat("For commercial use of the dmu-package please contact QGG. Email: Per.Madsen@agrsci.dk  !\n")


if(!is.null(phe_file)){


last_slash <- sub(".*/([^/]+)$", "\\1", phe_file)
before_last_slash <- sub("/[^/]*$", "", phe_file)

phe_path=before_last_slash
phe_name=last_slash

}


if(!is.null(kinship_file)){

	
	for(i in 1:length(kinship_file)){

		last_slash <- sub(".*/([^/]+)$", "\\1", kinship_file[i])
		before_last_slash <- sub("/[^/]*$", "", kinship_file[i])

		relationship_name=c(relationship_name,last_slash)
		relationship_path=c(relationship_path,before_last_slash)
	}
	relationship_path=unique(relationship_path)
	
}

if(!is.null(prior_file)){


last_slash <- sub(".*/([^/]+)$", "\\1", prior_file)
before_last_slash <- sub("/[^/]*$", "", prior_file)

provided_prior_file_path=before_last_slash
provided_prior_file_name=last_slash

}

#checking input parameter 
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
if(show_message)cat(paste0("Results are saved in path: ",output_result_path,"\n"))

if(show_message)cat(paste0("Start the ",analysis_model," analyse of ",length(target_trait_name)," trait model:",paste(target_trait_name,collapse = " & ")," \n"))


#be patient to space that appears in provided path
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

 if(length(grep(" ",phe_file))>0){
	# last_slash <- sub(".*/([^/]+)$", "\\1", my_string)
	# before_last_slash <- sub("/[^/]*$", "", my_string)
	stop("phe_file shouldn't contain blank!")
 }

 # if(length(grep(" ",kinship_file))>0){
	# stop("kinship_file shouldn't contain blank!")
 # }

#################
#*************generate DIR. file **********

if(!is.null(provided_DIR_file_path)&!is.null(provided_DIR_file_name)){

system(paste0("cp -r ",provided_DIR_file_path,"/",provided_DIR_file_name,"  ",output_result_path,"/","Trait.DIR"))

}else{

       DIR=generate_DIR(phe_col_names=phe_col_names,
					   target_trait_name=target_trait_name,
					   fixed_effect_name=fixed_effect_name, 
					   random_effect_name=random_effect_name, 
					   covariate_effect_name=covariate_effect_name, 
					   random_regression_effect_name=random_regression_effect_name, 
					   maternal_effect_name=maternal_effect_name, 
					   include_social_effect=include_social_effect, 
					   integer_group_names=integer_group_names,
					   real_group_names=real_group_names,
					   group_effect_name=group_effect_name,           
					   ped_inbred_number=ped_inbred_number, 
					   provided_effect_file_path=provided_effect_file_path,
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
					   integer_n=integer_n,  
					   residual_cov_trait=residual_cov_trait, 
					   output_DIR_path=output_result_path,
					   output_DIR_name="Trait.DIR",
					   IND_geno_file_name=IND_geno_file_name, 
					   IND_geno=IND_geno,
					   SSBLUP_omega=SSBLUP_omega,
					   EM_iter=EM_iter,
					   dmu5_iter=dmu5_iter,
					   gibbs_n_rounds=gibbs_n_rounds,
					   gibbs_n_burnin=gibbs_n_burnin,
					   gibbs_n_gaps=gibbs_n_gaps,
					   gibbs_seeds=gibbs_seeds,
					   phe_file=phe_file,
					   kinship_file=kinship_file
					   )
}


##################

#################
#***********run DMU1 and DMUai***************
setwd(output_result_path)


#get the system used 
os_type=as.character(Sys.info()["sysname"])

if(os_type!="Windows"){

dmu1= try(system(paste0(DMU_software_path,"/dmu1  release=5.3 < ",output_result_path,"/Trait.DIR"," >  ",paste(target_trait_name,collapse = "_"),".dmu1.lst")), silent = TRUE)
dmuai=try(system(paste0(DMU_software_path,"/",dmu_module," release=5.3 < ",output_result_path,"/Trait.DIR"," >  ",paste(target_trait_name,collapse = "_"),".",dmu_module,".lst")), silent = TRUE)

}else{
file.copy(from=paste0(DMU_software_path,"/dmu1.exe"),to=output_result_path)
file.copy(from=paste0(DMU_software_path,"/",dmu_module,".exe"),to=output_result_path)
#generate .bat file for running dmu under windows system
generate_dmu_bat("dmu1")
generate_dmu_bat(dmu_module)
system("run_dmu1.bat Trait.DIR")
system(paste0("run_",dmu_module,".bat"," Trait.DIR"))
file.remove(list.files(pattern = "*.exe"))
file.remove(list.files(pattern = "*.bat"))
}


if(dmu1 != 0){
	system(paste0("tail ",paste(target_trait_name,collapse = "_"),".dmu1.lst"))
	stop(paste0("Found errors in running dmu1, please check the log file:",paste(target_trait_name,collapse = "_"),".dmu1.lst"))

}	

if(dmuai != 0 & dmuai != 232){ #for GBLUP errors, returning code 1000
	system(paste0("tail ",paste(target_trait_name,collapse = "_"),".dmuai.lst"))
	stop(paste0("Found errors in running dmu1 or dmuai, please check the log file:",paste(target_trait_name,collapse = "_"),".dmu1.lst & ",paste(target_trait_name,collapse = "_"),".dmuai.lst"))
}	




if(file.exists("HINV1"))file.remove("HINV1")
if(file.exists("COR2"))file.remove("COR2")
if(file.exists("COR1"))file.remove("COR1")
if(file.exists("COR3"))file.remove("COR3")
#************finished run DMU1 and DMUai*****************
#################

h2_data=NULL

if(!is.null(provided_prior_file_path)&!is.null(provided_prior_file_name)){
file.copy(paste0(provided_prior_file_path,"/",provided_prior_file_name),
            paste0(output_result_path,"/","PAROUT"))
}

vars=NULL


#calculate model_reliability  and  corrected_phe 
if("PAROUT_STD"%in%list.files()&(is.null(maternal_effect_name))&(FALSE%in%include_social_effect)){

vars$inv_ai_mat=cal_dmu_se_reml(target_trait_name=target_trait_name,"PAROUT_STD")$inv_ai_mat

prior_mat=cal_dmu_se_reml(target_trait_name=target_trait_name,"PAROUT_STD")$mat_prior
prior_mat_se=cal_dmu_se_reml(target_trait_name=target_trait_name,"PAROUT_STD")$mat_prior_se


vars$vars_mat=cbind(prior_mat[,4],prior_mat_se[,3])

#print(inv_ai_mat)
#print(prior_mat)

for(i in 1:length(target_trait_name)){

prior=round(cal_dmu_se_reml(target_trait_name=target_trait_name,"PAROUT_STD")$prior[[i]],5)
prior_se=round(cal_dmu_se_reml(target_trait_name=target_trait_name,"PAROUT_STD")$prior_se[[i]],5)
h2=round(cal_dmu_se_reml(target_trait_name=target_trait_name,"PAROUT_STD")$h2[[i]],5)
h2_se=round(cal_dmu_se_reml(target_trait_name=target_trait_name,"PAROUT_STD")$h2_se[[i]],5)

h2_data=cbind(c(DIR$random_effect_name[[i]],"Residual"),prior,prior_se,h2,h2_se)

colnames(h2_data)=c("Random_effect_name","prior","prior_se","h2","h2_se")
write.table(h2_data,paste0(target_trait_name[i],"_heritability_result.txt"),quote=F,row.names=F,col.names=T,sep="\t")
}


#assgin name for AI-inv matrix 
#calculate the se of variance components
# my_cal_se<-function(expr,vars_value,inv_ai=NULL){  #expr( R expresson ), vars_value(variance components)
			# ai_name=colnames(inv_ai)  #random effect names of ai matrix
			# vars_value=as.numeric(vars_value[1,])
			# for (i in 1:length(ai_name))
				# assign(ai_name[i], vars_value[i])
				# grad <- t(as.numeric(attr(eval(deriv(expr, ai_name)), "gradient")))
				# h2<-eval(parse(text=as.character(expr)[-1]))
				# return(list(value=h2,value_se=ifelse(diag(grad %*% inv_ai %*% t(grad))<=0,0,sqrt(diag(grad%*% inv_ai %*% t(grad))))))
# }


ai_lev=DIR$ai_name$lev
ai_lev_name=DIR$ai_name$lev_name;
#ai_lev_name[1]="Gen" #set Gen as genetic effect name
ai_prior_lev=prior_mat[,1:3] #the first three columns of PRIOR
ai_trait_effect=DIR$random_effect_name

ai_mat_names=NULL
for(i in 1:nrow(ai_prior_lev)){

	if(ai_prior_lev[i,2]==ai_prior_lev[i,3]){
	
		i_ai_names=paste0(ai_lev_name[ai_prior_lev[i,1]],ai_prior_lev[i,2])
	
	}else{
	
		i_ai_names=paste0(ai_lev_name[ai_prior_lev[i,1]],paste(sort(ai_prior_lev[i,2:3]),collapse="_"))		
		
	}

	ai_mat_names=c(ai_mat_names,i_ai_names)

}


rownames(vars$inv_ai_mat)=colnames(vars$inv_ai_mat)=ai_mat_names
rownames(vars$vars_mat)=ai_mat_names
colnames(vars$vars_mat)=c("Esimated Variance","SE")

vars$gen_cor=NULL
vars$gen_cor_se=NULL

#calculate h2 and SE for each trait
trait_h2=NULL
k=0
for(i_vars_name in ai_trait_effect){

	k=k+1#which trait 
	
	i_vars_name=c(i_vars_name,"Res")
	i_vars_name=paste0(i_vars_name,k)	

	
	h2=NULL 
	h2_se=NULL 
	
	for (j in i_vars_name){
		expr=eval(parse(text=paste0(paste0("~",j,"/","(",paste(i_vars_name,collapse="+"),")"))))	
		result=my_cal_se(expr,vars$vars_mat[,1],vars$inv_ai_mat)

		h2=c(h2,result[[1]])
		h2_se=c(h2_se,result[[2]])
	}
	i_trait_h2=data.frame(Heritability=h2,SE=h2_se,row.names=i_vars_name,stringsAsFactors=F)

trait_h2=c(trait_h2,list(i_trait_h2))	

}


vars$h2=trait_h2
names(vars$h2)=target_trait_name

#calculate genetic correlation and SE for each trait 
#genetic correlation
if(length(target_trait_name)>=2){
	n_trait=length(target_trait_name)
	gen_cor=diag(n_trait)
	gen_cor_se=diag(0,n_trait) #diagonal is 0
		r2=NULL
		r2_se=NULL
		for(i in 1:(n_trait-1)){
			if(i!=n_trait){
			
				for(j in (i+1):n_trait){
				
					g_lev=paste0(genetic_effect_name,i,"_",j)
					expr=eval(parse(text=paste0(paste0("~",g_lev,"/","sqrt(",genetic_effect_name,i,"*",genetic_effect_name,j,")"))))
					result=my_cal_se(expr,vars$vars_mat[,1],vars$inv_ai_mat)
					r2=c(r2,result[[1]]) #genetic correlation
					r2_se=c(r2_se,result[[2]]) #se of genetic correlation
				}
			}
		}
	gen_cor[upper.tri(gen_cor,diag=F)]<-r2;
	gen_cor[lower.tri(gen_cor)]=t(gen_cor)[lower.tri(t(gen_cor))]
	gen_cor_se[upper.tri(gen_cor_se,diag=F)]<-r2_se;
	gen_cor_se[lower.tri(gen_cor_se)]=t(gen_cor_se)[lower.tri(t(gen_cor_se))]

	colnames(gen_cor)=rownames(gen_cor)=target_trait_name
	colnames(gen_cor_se)=rownames(gen_cor_se)=target_trait_name
	
	vars$gen_cor=gen_cor
	vars$gen_cor_se=gen_cor_se
	
}



if(!is.null(provided_prior_file_path) & dmu_module%in%c("dmu4","dmu5")){
system(paste0("cp -r ",provided_prior_file_path,"/",provided_prior_file_name,"  ",output_result_path))
}


}


if(is.null(phe_file)){
phe_file=paste0(phe_path,"/",phe_name)
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
				   phe_path=phe_path,  
				   phe_name=phe_name,
				   pedigree_path=debv_pedigree_path, 
				   pedigree_name=debv_pedigree_name, 
				   provided_prior_file_path=provided_prior_file_path,
				   provided_prior_file_name=provided_prior_file_name,							
				   cal_debv=cal_debv,
				   analysis_model=analysis_model,
				   relationship_name=relationship_name,
				   debv_id=debv_id, 
				   output_ebv_path=output_ebv_path,
				   output_ebv_name=output_ebv_name,
				   return_result=return_result
							)

if(plot_variance==TRUE){
message("Plot the result of variance components......")
plot_dmu_blupf90_prior(genetic_effect_name=genetic_effect_name,target_trait_name=target_trait_name,random_effect_name=random_effect_name,output_path=output_result_path)
}

# calculate heritablilty and SE
if(show_message)cat(paste0("Completed the ",analysis_model," analyse of ",length(target_trait_name)," trait model:",paste(target_trait_name,collapse = " & ")," \n"))

if(return_result==TRUE){
return(list(dmu_result,h2_data,DIR,vars=vars))
}

}


#ruing DMU under windows system
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



#calculte standared error 
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
cov_matrix=diag(1,npar)   

for(i in 1:nrow(cov)){
 cov_matrix[cov[i,1],cov[i,2]]=cov[i,4]
 cov_matrix[cov[i,2],cov[i,1]]=cov[i,4]
}

h2=NULL;h2_se=NULL;total_trait_prior=NULL;trait_prior_se=NULL;

for(i in 1:length(target_trait_name)){

Trait_prior_pos=(1:nrow(Prior))[Prior[,2]%in%i & Prior[,3]%in%i]
Trait_prior=as.matrix(Prior[Trait_prior_pos,4])
Trait_cov_matrix=cov_matrix[Trait_prior_pos,Trait_prior_pos]

Var=Trait_cov_matrix 
heritablilty=Trait_prior/sum(Trait_prior[,1])
SE_h2=matrix(NA,nrow=nrow(Trait_prior),ncol=1)
for(i in 1:nrow(Trait_prior)){
derivate=NULL
for(j in 1:nrow(Trait_prior)){
if(j==i){derivate=cbind(derivate,(sum(Trait_prior[,1])-Trait_prior[i,1])/(sum(Trait_prior[,1])^2) )}
if(j!=i){derivate=cbind(derivate,(-Trait_prior[i,1])/(sum(Trait_prior[,1])^2 ))}
}
SE_h2[i,1]=sqrt(derivate%*%Var%*%t(derivate)) 
}

total_trait_prior=c(total_trait_prior,list(Trait_prior))
trait_prior_se=c(trait_prior_se,list(prior_se[Trait_prior_pos,3]))
h2=c(h2,list(heritablilty))
h2_se=c(h2_se,list(SE_h2))
}

return(list(prior=total_trait_prior,prior_se=trait_prior_se,h2=h2,h2_se=h2_se,inv_ai_mat=cov_matrix,mat_prior=Prior,mat_prior_se=prior_se))
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