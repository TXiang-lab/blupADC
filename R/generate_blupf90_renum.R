#' Print 'Hello world!'
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

Trait_n=length(target_trait_name)

if(is.null(relationship_name)){
stop("Please provide relationship name!")
}else{
addtive_relationship_name=relationship_name[1]
if(analysis_model%in%c("GBLUP_AD")){dominance_relationship_name=relationship_name[2]}
if(analysis_model%in%c("SSBLUP_A")){SSBLUP_G_matrix_name=relationship_name[2]}
}



if(!is.null(fixed_effect_name)){if(length(fixed_effect_name)!=Trait_n){stop("The number of fixed effect names doesn't match the number of traits")}}
if(!is.null(random_effect_name)){if(length(random_effect_name)!=Trait_n){stop("The number of random effect names doesn't match the number of traits")}}
if(!is.null(covariate_effect_name)){if(length(covariate_effect_name)!=Trait_n){stop("The number of covariate effect names doesn't match the number of traits")}}



#PRIOR
#BLUPF90 PRIOR的结构和DMU不同，和 Additive (CO)VARIANCES 结构一致
if(!is.null(provided_BLUPF90_prior_file_name)&!is.null(provided_BLUPF90_prior_file_path)&!is.null(provided_BLUPF90_prior_effect_name)){
given_prior=data.table::fread(paste0(provided_BLUPF90_prior_file_path,"/",provided_BLUPF90_prior_file_name),data.table=F,header=F)
given_prior=as.matrix(given_prior)
}else{

prior=BLUPF90_generate_prior(target_trait_name=target_trait_name,
						      random_effect_name=random_effect_name,
							  included_permanent_effect=included_permanent_effect)

given_prior=as.matrix(prior)

provided_BLUPF90_prior_effect_name=c(do.call(c,random_effect_name),rep("Residual",length(target_trait_name)))

}

#获取固定效应、随机效应、协变量记录
if(!is.null(provided_effect_file_path)&!is.null(provided_effect_file_name)){

model=data.table::fread(paste0(provided_effect_file_path,"/",provided_effect_file_name),data.table=F,fill=TRUE,header = F)

fixed_effect_name=NULL;random_effect_name=NULL;covariate_effect_name=NULL
for(i in 1:length(target_trait_name)){

n=match(target_trait_name[i],model[,1]) # which line  is  the trait

if(is.na(n)){stop(paste0(target_trait_name[i]," couldn't find in provided effect_file"))}

effect_col=grep(pattern="*",as.character(model[n,]),fixed=T)
fixed_effect=as.character(model[n,c((effect_col[1]+1):(effect_col[2]-1))])
random_effect=as.character(model[n,c((effect_col[2]+1):(effect_col[3]-1))])
covariate_effect=as.character(model[n,c((effect_col[3]+1):(effect_col[4]-1))])

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
trait_position=match(target_trait_name,phe_col_names)


header_information_part=rbind("DATAFILE",phe_name,
							  "TRAITS",paste0(trait_position,collapse=" "),
							  "FIELDS_PASSED TO OUTPUT","","WEIGHT(S)","")
							  

####################       #################################
########construct  Residual Part                   #########
temp_effect="Residual"
#temp_effect_pos=match(temp_effect,provided_BLUPF90_prior_effect_name)
#co_variance_pos=(Trait_n*(temp_effect_pos-1)+1):(Trait_n*temp_effect_pos)
temp_effect_pos=which(provided_BLUPF90_prior_effect_name%in%temp_effect)+length(relationship_name)-1
co_variance_pos=temp_effect_pos
co_variance=given_prior[co_variance_pos,co_variance_pos]
co_variance=inner_matrix(co_variance)

residual_effect_part=rbind("RESIDUAL_VARIANCE",co_variance) # 残差

#获取各性状下效应对应的位置
####################       #################################
########construct  Fixed-Effect Part               #########

if(!is.null(fixed_effect_name)&(!identical(fixed_effect_name,rep(list(NULL),Trait_n)))){
total_fixed_effect=unique(do.call(c,fixed_effect_name))
total_fixed_effect=sort(total_fixed_effect)
total_fixed_effect_pos=match(total_fixed_effect,phe_col_names)

fixed_effect_part=NULL
fixed_effect_pos=NULL

for(i in 1:Trait_n){  #每个性状下固定效应的 pos
	temp_fixed_pos=total_fixed_effect_pos
	temp_fixed_pos[!total_fixed_effect%in%fixed_effect_name[[i]]]=0
	fixed_effect_pos=c(fixed_effect_pos,list(temp_fixed_pos))
}
fixed_effect_pos_matrix=do.call(rbind,fixed_effect_pos) #pos的矩阵格式

for(i in 1:length(total_fixed_effect)){
	fixed_effect_part=rbind(fixed_effect_part,"EFFECT",paste0(paste(fixed_effect_pos_matrix[,i],collapse=" ")," cross alpha "))
}
}else{
fixed_effect_part=data.frame(NULL)}

####################       #################################
########construct  Covariate-Effect Part           #########
if(!is.null(covariate_effect_name)&(!identical(covariate_effect_name,rep(list(NULL),Trait_n)))){
total_covariate_effect=unique(do.call(c,covariate_effect_name))
total_covariate_effect=sort(total_covariate_effect)
total_covariate_effect_pos=match(total_covariate_effect,phe_col_names)
covariate_effect_part=NULL
covariate_effect_pos=NULL

for(i in 1:Trait_n){  #每个性状下固定效应的 pos
	temp_covariate_pos=total_covariate_effect_pos
	temp_covariate_pos[!total_covariate_effect%in%covariate_effect_name[[i]]]=0
	covariate_effect_pos=c(covariate_effect_pos,list(temp_covariate_pos))
}
covariate_effect_pos_matrix=do.call(rbind,covariate_effect_pos) #pos的矩阵格式


for(i in 1:length(total_covariate_effect)){
	covariate_effect_part=rbind(covariate_effect_part,"EFFECT",paste0(paste(covariate_effect_pos_matrix[,i],collapse=" ")," cov "))
}
}else{
covariate_effect_part=data.frame(NULL)}


####################       #################################
########construct  Non_genetic-random effect Part  #########

#for(i in 1:length(random_effect_name)){   #remove genetic effect name
#temp_random=random_effect_name[[i]]
#temp_random=temp_random[!temp_random%in%genetic_effect_name]
#if(length(temp_random)==0){temp_random=NULL}
#random_effect_name[[i]]=temp_random
#}

total_random_effect=unique(do.call(c,random_effect_name))
total_random_effect=sort(total_random_effect)
total_random_effect_pos=match(total_random_effect,phe_col_names)

#remove genetic effect
total_random_effect_pos=total_random_effect_pos[!total_random_effect%in%genetic_effect_name]
total_random_effect=total_random_effect[!total_random_effect%in%genetic_effect_name]

if(length(total_random_effect)>0&(!identical(random_effect_name,rep(list(NULL),Trait_n)))){
non_genetic_random_effect_part=NULL
random_effect_pos=NULL

for(i in 1:Trait_n){  #每个性状下固定效应的 pos
	temp_random_pos=total_random_effect_pos
	temp_random_pos[!total_random_effect%in%random_effect_name[[i]]]=0
	random_effect_pos=c(random_effect_pos,list(temp_random_pos))
}

random_effect_pos_matrix=do.call(rbind,random_effect_pos) #pos的矩阵格式

for(i in 1:length(total_random_effect)){

	temp_effect=total_random_effect[i]
	#temp_effect_pos=match(temp_effect,provided_BLUPF90_prior_effect_name)
	#co_variance_pos=(Trait_n*(temp_effect_pos-1)+1):(Trait_n*temp_effect_pos)
	temp_effect_pos=which(provided_BLUPF90_prior_effect_name%in%temp_effect)
	co_variance_pos=temp_effect_pos
	co_variance=given_prior[co_variance_pos,co_variance_pos]
	co_variance=inner_matrix(co_variance)

	non_genetic_random_effect_part=rbind(non_genetic_random_effect_part,"EFFECT",paste0(paste(random_effect_pos_matrix[,i],collapse=" ")," cross alpha "),
	                              "RANDOM","diagonal","(CO)VARIANCES",co_variance)
}
}else{
non_genetic_random_effect_part=data.frame(NULL)}



genetic_random_effect_part=NULL
random_effect_pos=NULL
Optional_genetic_effect=NULL
permanent_effect_part=NULL

####################       #################################
########construct  Permanent Effect Part   		   #########

if(included_permanent_effect==TRUE){
Optional_genetic_effect=rbind("OPTIONAL","pe")
	temp_effect="Permanent"
	temp_effect_pos=match(temp_effect,provided_BLUPF90_prior_effect_name)
	co_variance_pos=(Trait_n*(temp_effect_pos-1)+1):(Trait_n*temp_effect_pos)
	co_variance=given_prior[co_variance_pos,co_variance_pos]
	co_variance=inner_matrix(co_variance)
	permanent_effect_part=rbind("(CO)VARIANCES_PE",co_variance)
}

####################       #################################
########construct  Genetic-random effect Part      #########

if(analysis_model=="GBLUP_A"){
	temp_effect=genetic_effect_name
	temp_effect_pos=match(temp_effect,provided_BLUPF90_prior_effect_name)
	co_variance_pos=(Trait_n*(temp_effect_pos-1)+1):(Trait_n*temp_effect_pos)
	co_variance=given_prior[co_variance_pos,co_variance_pos]
	co_variance=inner_matrix(co_variance)

	genetic_random_effect_part=rbind(genetic_random_effect_part,"EFFECT",paste0(paste(match(rep(genetic_effect_name,Trait_n),phe_col_names),collapse=" ")," cross alpha "),
	                              "RANDOM","animal","SNP_FILE",addtive_relationship_name,
						     "(CO)VARIANCES",co_variance,permanent_effect_part)
}else if(analysis_model=="PBLUP_A"){

	temp_effect=genetic_effect_name
	temp_effect_pos=match(temp_effect,provided_BLUPF90_prior_effect_name)
	co_variance_pos=(Trait_n*(temp_effect_pos-1)+1):(Trait_n*temp_effect_pos)
	co_variance=given_prior[co_variance_pos,co_variance_pos]
	co_variance=inner_matrix(co_variance)

	genetic_random_effect_part=rbind(genetic_random_effect_part,"EFFECT",paste0(paste(match(rep(genetic_effect_name,Trait_n),phe_col_names),collapse=" ")," cross alpha "),
	                              "RANDOM","animal",Optional_genetic_effect,"FILE",addtive_relationship_name,"FILE_POS","1 2 3 0 0",
							 "PED_DEPTH","0",
						     "(CO)VARIANCES",co_variance,permanent_effect_part)

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

	for(i in 1:length(relationship_name)){ 
	temp_effect=genetic_effect_name
	#temp_effect_pos=match(temp_effect,provided_BLUPF90_prior_effect_name)
	#temp_effect_pos=effect_pos[provided_BLUPF90_prior_effect_name%in%temp_effect][i]
	temp_effect_pos=which(provided_BLUPF90_prior_effect_name%in%temp_effect)
	co_variance_pos=temp_effect_pos+i-1
	#co_variance_pos=(Trait_n*(temp_effect_pos-1)+1):(Trait_n*temp_effect_pos)
	co_variance=given_prior[co_variance_pos,co_variance_pos]
	co_variance=inner_matrix(co_variance)

	genetic_random_effect_part=rbind(genetic_random_effect_part,"EFFECT",paste0(paste(match(rep(genetic_effect_name,Trait_n),phe_col_names),collapse=" ")," cross alpha "),
	                              "RANDOM","animal",Optional_genetic_effect,"FILE",addtive_relationship_name,"FILE_POS","1 2 3 0 0",
							 "PED_DEPTH","0",
						     "(CO)VARIANCES",co_variance,permanent_effect_part)
	}
	
	
}else if(analysis_model=="SSBLUP_A"){

	temp_effect=genetic_effect_name
	temp_effect_pos=match(temp_effect,provided_BLUPF90_prior_effect_name)
	co_variance_pos=(Trait_n*(temp_effect_pos-1)+1):(Trait_n*temp_effect_pos)
	co_variance=given_prior[co_variance_pos,co_variance_pos]
	co_variance=inner_matrix(co_variance)

	genetic_random_effect_part=rbind(genetic_random_effect_part,"EFFECT",paste0(paste(match(rep(genetic_effect_name,Trait_n),phe_col_names),collapse=" ")," cross alpha "),
	                              "RANDOM","animal",Optional_genetic_effect,"FILE",addtive_relationship_name,"FILE_POS","1 2 3 0 0",
							"SNP_FILE",relationship_name[2],
							 "PED_DEPTH","0",
						     "(CO)VARIANCES",co_variance,permanent_effect_part)


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

SNP_genumeric_file_part_option=data.frame(NULL)
if(!is.null(BLUPF90_genumeric_name)){
SNP_genumeric_file_part_option=paste0("OPTION ","SNP_file ",BLUPF90_genumeric_name)
}

file_OPTION_part=rbind(paste0("OPTION missing ",missing_value),
						 "OPTION alpha_size 30",
						 "OPTION use_yams",
						 "OPTION sol se",
						 #"OPTION se_covar_function H2d G_2_2_1_1/(G_2_2_1_1+G_2_3_1_1+G_3_3_1_1+G_4_4_1_1+R_1_1)",
						 "OPTION residual")

if(!is.null(BLUPF90_alt_option)){
for(option in BLUPF90_alt_option){
file_OPTION_part=rbind(file_OPTION_part,option)
}
}

renum_par_file=plyr::rbind.fill.matrix(header_information_part,
								 residual_effect_part,
								 fixed_effect_part,
								 covariate_effect_part,
								 non_genetic_random_effect_part,
								 genetic_random_effect_part,
								 SNP_map_file_part_option,
								 SNP_genumeric_file_part_option,
								 file_OPTION_part)
utils::write.table(renum_par_file,"renum.par",quote=F,row.names=F,col.names=F)
}



#将n列矩阵变成一列矩阵
inner_matrix<-function(data){
data=as.matrix(data)
temp_matrix=matrix(NA,nrow=nrow(data),ncol=1)
for(i in 1:nrow(temp_matrix)){
temp_matrix[i,]=paste(data[i,],collapse=" ")
}
return(temp_matrix)
}


BLUPF90_generate_prior<-function(target_trait_name=NULL,
								 random_effect_name=NULL,
								 included_permanent_effect=FALSE){


prior=diag(10+length(target_trait_name)+length(do.call(c,random_effect_name))+included_permanent_effect)
prior[prior==0]=0.5

return(prior)
}

