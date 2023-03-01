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
					   SSBLUP_omega=0.05,
					   EM_iter=10,
					   dmu5_iter=1000000,
					   gibbs_n_rounds=10000,
					   gibbs_n_burnin=1000,
					   gibbs_n_gaps=5,
					   gibbs_seeds=c(1994,1998)
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

for(i in 1:length(target_trait_name)){
random_effect_name[[i]]=c(random_effect_name[[i]],integer_group_names)  # 单性状
}

}else{

cat("Using user-provided social genetic effect phenotype......\n")

if(is.null(integer_group_names)|is.null(real_group_names)){stop("User must specify the integer_group_names or real_group_names when using provided phenotype!")}
integer_group_names=integer_group_names
real_group_names=real_group_names


random_regression_effect_name=rep(list(paste0(real_group_names,"&",integer_group_names)),
                                          length(target_trait_name)) #单性状

for(i in 1:length(target_trait_name)){
random_effect_name[[i]]=c(random_effect_name[[i]],integer_group_names)  # 单性状
}

}
}


if(!is.null(IND_geno)){
write.table(IND_geno,paste0(relationship_path[1],"/","IND_geno.txt"),quote=F,row.names=F,col.names=F,sep=" ")
IND_geno_file_name="IND_geno.txt"
}

					   
Trait_n=length(target_trait_name)
Variable_n=length(phe_col_names)
real_n=Variable_n-integer_n
Integer_col_names=phe_col_names[1:integer_n]
Real_col_names=phe_col_names[(integer_n+1):Variable_n]

library(data.table)
#读取各性状的固定效应、随机效应、协变量记录
if(!is.null(random_effect_name)&!is.null(provided_effect_file_path)){stop("Effect_file and  effect_names couldn't provided simultaneously! ")}
if(!is.null(provided_effect_file_path)&!is.null(provided_effect_file_name)){

model=fread(paste0(provided_effect_file_path,"/",provided_effect_file_name),data.table=F,fill=TRUE,header = F)

fixed_effect_name=NULL;random_effect_name=NULL;covariate_effect_name=NULL
for(i in 1:length(target_trait_name)){

n=match(target_trait_name[i],model[,1]) # which line  is  the trait

if(is.na(n)){stop(paste0(target_trait_name[i]," couldn't find in provided effect_file; "))}

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


#software create DIR file 
COMMENT=matrix("",nrow=2,ncol=length(phe_col_names)+10)

ANALYSE=matrix("",nrow=2,ncol=length(phe_col_names)+10)

DATA=matrix("",nrow=2,ncol=length(phe_col_names)+10)

VARIABLE=matrix("",nrow=10,ncol=length(phe_col_names)+10)

VAR_STR=matrix("",nrow=4+length(relationship_name),ncol=length(phe_col_names)+10)

SOLUTION=matrix("",nrow=2,ncol=length(phe_col_names)+10)

MODEL=matrix("",nrow=1+Trait_n*(Trait_n+6),ncol=length(phe_col_names)+10)

PRIOR=matrix("",nrow=(Trait_n*(Trait_n+1)/2)*6,ncol=length(phe_col_names)+10)

PARAMETER=matrix("",nrow=10*length(target_trait_name),ncol=length(phe_col_names)+10)




#$COMMENT
COMMENT[1,1]="$COMMENT"

#$ANALYSE 
if(is.null(dmu_module)){dmu_module="dmuai"}
dmu_module_code=ifelse(dmu_module=="dmuai",1,
                   ifelse(dmu_module=="rjmc",2,
		          ifelse(dmu_module=="dmu4",11,
			      ifelse(dmu_module=="dmu5",12,NA))))
	
if(is.null(dmu_algorithm_code)){	
dmu_algorithm_code=ifelse(dmu_module=="dmuai",1,
                   ifelse(dmu_module=="rjmc",0,
		          ifelse(dmu_module=="dmu4",14,
			      ifelse(dmu_module=="dmu5",2,NA))))}		

		  
ANALYSE[1,1]=paste0("$ANALYSE ",dmu_module_code," ",dmu_algorithm_code," 0 0 ")

#$DATA 
DATA[1,1]=paste0("$DATA  ASCII (",integer_n,",",real_n," ,",missing_value,") ",phe_path,"/",phe_name)


#$VARIABLE
VARIABLE[1,1]="$VARIABLE"
VARIABLE[3,1:integer_n]=paste0("I",1:integer_n)
VARIABLE[4,1:integer_n]=Integer_col_names
VARIABLE[6,1:real_n]=paste0("R",1:real_n)
VARIABLE[7,1:real_n]=Real_col_names
VARIABLE[3,1]="#I1";VARIABLE[6,1]="#R1";


#将 I1.... R1.... 和变量名称对齐
for(i in 1:ncol(VARIABLE)){
integer_space=nchar(VARIABLE[4,i])
VARIABLE[3,i]=sprintf(paste0("%-",integer_space,"s"),VARIABLE[3,i])
real_space=nchar(VARIABLE[7,i])
VARIABLE[6,i]=sprintf(paste0("%-",real_space,"s"),VARIABLE[6,i])
}


#$MODEL
MODEL[1,1]="$MODEL"
MODEL[2,1]=paste0(Trait_n," ",Trait_n," 0 0 0")
MODEL[3:(2+Trait_n),1]="0"

if(!is.null(fixed_effect_name)){if(length(fixed_effect_name)!=Trait_n){stop("The number of fixed effect names doesn't match the number of traits")}}
if(!is.null(random_effect_name)){if(length(random_effect_name)!=Trait_n){stop("The number of random effect names doesn't match the number of traits")}}
if(!is.null(covariate_effect_name)){if(length(covariate_effect_name)!=Trait_n){stop("The number of covariate effect names doesn't match the number of traits")}}



#User_define模型下，添加多个遗传效应 放入 random_effect_name, 后续只需要在 trait_fixed_random_line 中，把遗传效应的位置替换成 Id的位置，其余不变
if(("User_define"%in%analysis_model)&length(relationship_name)>=2){
User_define_effect=paste0("User_define",1:(length(relationship_name)-1))
for(i in 1:length(random_effect_name)){
random_effect_name[[i]]=c(random_effect_name[[i]],User_define_effect)
random_effect_name[[i]]=c(genetic_effect_name,User_define_effect,setdiff(random_effect_name[[i]],c(genetic_effect_name,User_define_effect)))
}
Integer_col_names=c(Integer_col_names,User_define_effect)
}


#添加显性效应 放入 random_effect_name, 后续只需要在 trait_fixed_random_line 中，把dominance效应的位置替换成 Id的位置，其余不变
if(included_dominance_effect==TRUE){
dominance_effect="Added_dominance_effect"
for(i in 1:length(random_effect_name)){
random_effect_name[[i]]=c(random_effect_name[[i]],dominance_effect)
random_effect_name[[i]]=c(genetic_effect_name,dominance_effect,setdiff(random_effect_name[[i]],c(genetic_effect_name,dominance_effect)))
}
Integer_col_names=c(Integer_col_names,dominance_effect)
}


#添加永久环境效应 放入 random_effect_name, 后续只需要在 trait_fixed_random_line 中，把permanent的位置替换成 Id的位置，其余不变
if(identical(included_permanent_effect,TRUE)){included_permanent_effect=rep(list(included_permanent_effect),length(target_trait_name))}
if(identical(included_permanent_effect,FALSE)){included_permanent_effect=rep(list(included_permanent_effect),length(target_trait_name))}

permanent_effect="pe_effect"
if(permanent_effect%in%Integer_col_names){stop("pe_effect already stands for permant effect,please change this effect name in the phe_col_name!!! ")}
for(i in 1:length(random_effect_name)){

if(included_permanent_effect[[i]]==TRUE){

random_effect_name[[i]]=c(random_effect_name[[i]],permanent_effect)
random_effect_name[[i]]=c(genetic_effect_name,permanent_effect,setdiff(random_effect_name[[i]],c(genetic_effect_name,permanent_effect)))

}
}

if(TRUE %in% included_permanent_effect)Integer_col_names=c(Integer_col_names,permanent_effect)


#添加母性效应 放入 random_effect_name, 后续只需要将母性效应对应的随机效应编号 设置与 遗传效应的编号一样即可

if(!is.null(maternal_effect_name)){

for(i in 1:length(random_effect_name)){

maternal=maternal_effect_name[[i]] #列表

if(!is.null(maternal)){

if(!maternal%in%Integer_col_names){stop("Maternal effect name could not find in provided phe_col_names; ")}

random_effect_name[[i]]=c(random_effect_name[[i]],maternal)
random_effect_name[[i]]=c(genetic_effect_name,maternal,setdiff(random_effect_name[[i]],c(genetic_effect_name,maternal)))


}
}}



#记录所有性状的名称，并自定义随机效应序号

random_lev=as.factor(unique(do.call(c,random_effect_name))) #统计所有的随机效应名称，并将随机效应赋值
levels(random_lev)=c(genetic_effect_name,setdiff(as.character(random_lev),genetic_effect_name))
random_lev_number=as.numeric(random_lev)
random_lev_names=as.character(random_lev)
random_total_lev_number=random_lev_number[match(do.call(c,random_effect_name),random_lev_names)]



#添加残差效应序号：
residual_number=length(random_lev_number)+1
random_total_lev_number=c(random_total_lev_number,rep(residual_number,Trait_n))


#多个性状循环
for(i in 1:Trait_n){
trait=target_trait_name[i]
fixed=fixed_effect_name[[i]]
random=random_effect_name[[i]]
if(!genetic_effect_name%in%random){stop("Genetic effect name could not find in provided random effect name! ")}
random=c(genetic_effect_name,setdiff(random,c(genetic_effect_name)))#确保 genetic effect name 为随机效应第一位
covariate=covariate_effect_name[[i]]
regression=random_regression_effect_name[[i]]
maternal=maternal_effect_name[[i]] #列表

#判断性状名称、效应名称是否都在 提供的 phe_col_names 中
if(is.na(match(trait,Real_col_names))){
stop(paste0("Trait:",trait," could not find in provided phe_col_names; "))
}else{trait_pos=match(trait,Real_col_names)}



if(!is.null(fixed)){
fixed_pos=match(fixed,Integer_col_names)
if(NA%in%fixed_pos){stop(paste0("Fixed effect:",fixed[is.na(fixed_pos)]," could not find in provided phe_col_names; "))}
}else{fixed_pos=NULL}

if(!is.null(random)){

random_pos=match(random,Integer_col_names)

if(NA%in%random_pos){stop(paste0("Random effect:",random[is.na(random_pos)]," could not find in provided phe_col_names; "))}
}else{random_pos=NULL}


if(!is.null(covariate)){
covariate_pos=match(covariate,Real_col_names)
if(NA%in%covariate_pos){stop(paste0("Covariate effect:",covariate[is.na(covariate_pos)]," could not find in provided phe_col_names; "))}
}else {covariate_pos=NULL}


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
#for(i in 1:n_regre){
#
#for(j in 1:length(regression_split)){
#regression_name_coef=regression_split[[j]][1]       # 第一个名称是 多项式回归系数的名称
#regression_name_effect=regression_split[[j]][-1]    # 非第一个，剩余的是 该多项式回归系数嵌套的所有随机效应
#
#if(regression_name_coef==regression_name_coef_set[i]){
#regression_name_effect_set[[i]]=c(regression_name_effect_set[[i]],regression_name_effect)
#}
#}
#tmp_reg=unique(as.character(na.omit(regression_name_effect_set[[i]])))
#linked_number=linked_number+length(tmp_reg)
#pos_In_blanket=match(match(tmp_reg,Integer_col_names),
#                                   c(fixed_pos,random_pos)) # 嵌套的变量在括号中的位置
#pos_In_blanket=paste(pos_In_blanket,collapse = " ")    # 合并成一个字符串
#
#regression_name_coef_effect_pos=paste0(match(regression_name_coef,Real_col_names),
#                                                   "(",
#										   pos_In_blanket,
#										   ") "
#										   )
#print(regression_name_coef_effect_pos)										   
#covariate_pos=c(covariate_pos,regression_name_coef_effect_pos)
#}

reg_k=length(regression) #随机回归效应的个数

for(j in 1:reg_k){

regression_name=unlist(strsplit(regression[j],split = "&"))
  
regression_name_coef=regression_name[1]       # 第一个名称是 多项式回归系数的名称
regression_name_effect=regression_name[-1]             #非第一个，剩余的是 该多项式回归系数嵌套的所有随机效应

linked_number=linked_number+length(regression_name_effect)-1 #针对括号内有2个及以上嵌套效应的情况，需要将协变量的总数目进行修改

#if(FALSE%in%(regression_name_effect%in%random)){stop("Random effect in random_regression_effect_name don't exist in the random_effect_name!")}

if(FALSE%in%(regression_name_effect%in%c(fixed,random))){
stop("Regression effect in random_regression_effect_name don't exist in the random_effect_name or fixed_effect_name!")
}


pos_In_blanket=match(match(regression_name_effect,Integer_col_names),
                                   c(fixed_pos,random_pos)) # 嵌套的变量在括号中的位置
pos_In_blanket=paste(pos_In_blanket,collapse = " ")    # 合并成一个字符串

								   
regression_name_coef_effect_pos=paste0(match(regression_name_coef,Real_col_names),
                                                   "(",
										   pos_In_blanket,
										   ") "
										   )
covariate_pos=c(covariate_pos,regression_name_coef_effect_pos)
}

}


trait_fixed_random_line=paste0(trait_pos," 0 ",length(fixed_pos)+length(random_pos)," ",paste(c(fixed_pos,random_pos),collapse=" "))
random_line=paste(c(length(random),random_lev_number[match(random,random_lev_names)]),collapse=" ")
covariate_line=paste(c(length(covariate_pos)+linked_number,covariate_pos),collapse=" ")


#考虑永久环境效应，将  trait_fixed_random_line 中的 Added_permanent_effect 的位置替换成 Id 位置
if(TRUE%in%included_permanent_effect[[i]]){
permanent_effect_pos=match(permanent_effect,Integer_col_names)
random_pos[random_pos%in%permanent_effect_pos]=match(genetic_effect_name,Integer_col_names)
trait_fixed_random_line=paste0(trait_pos," 0 ",length(fixed_pos)+length(random_pos)," ",paste(c(fixed_pos,random_pos),collapse=" "))
}
 
#考虑显性效应，将  trait_fixed_random_line 中的 Added_dominance_effect 的位置替换成 Id 位置
if(included_dominance_effect==TRUE){
dominance_effect_pos=match(dominance_effect,Integer_col_names)
random_pos[random_pos%in%dominance_effect_pos]=match(genetic_effect_name,Integer_col_names)
trait_fixed_random_line=paste0(trait_pos," 0 ",length(fixed_pos)+length(random_pos)," ",paste(c(fixed_pos,random_pos),collapse=" "))
}

 
#考虑User_define模型，将  trait_fixed_random_line 中的 User_define_effect 的位置替换成 Id 位置
if(("User_define"%in%analysis_model)&length(relationship_name)>=2){
User_define_effect_pos=match(User_define_effect,Integer_col_names)
random_pos[random_pos%in%User_define_effect_pos]=match(genetic_effect_name,Integer_col_names)
trait_fixed_random_line=paste0(trait_pos," 0 ",length(fixed_pos)+length(random_pos)," ",paste(c(fixed_pos,random_pos),collapse=" "))
}


#考虑母性效应，将  random_line 中的 maternal effect 的随机效应编号替换为遗传效应的随机效应编号
if(!is.null(maternal)){

#将maternal effect的水平修改成 和  遗传效应水平一致
#其他随机效应的水平均减去1

tmp_random_lev_number=random_lev_number

tmp_random_lev_number[match(maternal,random_lev_names)]=tmp_random_lev_number[match(genetic_effect_name,random_lev_names)]

tmp_random_lev_number[tmp_random_lev_number!=1]=tmp_random_lev_number[tmp_random_lev_number!=1]-1

random_line=paste(c(length(random),tmp_random_lev_number[match(random,random_lev_names)]),collapse=" ")

}
 

#考虑social genetic effect, 需要额外对DIR文件中，固定效应-随机效应行 和 随机效应编号 行进行修改
if(TRUE%in%include_social_effect){
if(TRUE%in%include_social_effect[[i]]){

#修改固定效应-随机效应行， 在integer_group_name后添加(0)

random_pos[match(integer_group_names,random)]=paste0(random_pos[match(integer_group_names,random)],"(","0",")")
trait_fixed_random_line=paste0(trait_pos," 0 ",length(fixed_pos)+length(random_pos)," ",paste(c(fixed_pos,random_pos),collapse=" "))

#修改随机效应编号 行，形如 +1+1+1

temp=random_lev_number[match(random,random_lev_names)]

temp1=random_lev_number[match(setdiff(random,integer_group_names),random_lev_names)] #从随机效应名称中排除 integer_group_names
temp2=paste(rep(1,length(random)-length(temp1)),collapse = "+")

random_line=paste(c(length(random),temp1,temp2),collapse=" ")

}}



MODEL[i+2+Trait_n,1]=trait_fixed_random_line
MODEL[(i+2+Trait_n)+Trait_n,1]=random_line
MODEL[(i+2+Trait_n)+2*Trait_n,1]=covariate_line

}

#残差限定
if(is.null(residual_cov_trait)){
MODEL[(Trait_n+2+Trait_n)+2*Trait_n+1,1]="0"
}else{
if(length(residual_cov_trait)>Trait_n){stop("The number of residual covariance traits is larger than the number of analyse traits \n")}
residual_cov_trait_n=0

for(i in 1:length(residual_cov_trait)){

residual_cov_combination=residual_cov_trait[[i]]

if(!is.null(residual_cov_combination)){
residual_cov_trait_n=residual_cov_trait_n+1
trait1=residual_cov_combination[1]
trait2=residual_cov_combination[2]
MODEL[(Trait_n+2+Trait_n)+2*Trait_n+1+i,1]=paste0(match(trait1,target_trait_name)," ",match(trait2,target_trait_name))
}
}

MODEL[(Trait_n+2+Trait_n)+2*Trait_n+1,1]=residual_cov_trait_n
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

for(i in 1:length(relationship_name)){
VAR_STR[i,1]=paste0("$VAR_STR ",i," COR ASCII     ",relationship_path,"/",relationship_name[i])
}}


#$SOLUTION
#if(dmu_module!="rjmc"){
SOLUTION[1,1]="$RESIDUALS ASCII"
SOLUTION[2,1]="$SOLUTION"
#}
#$PRIOR

if(include_social_effect==FALSE){  #social effect 不知道该如何构建PRIOR文件，

PRIOR[1,1]="$PRIOR"
if(!is.null(provided_prior_file_name)&!is.null(provided_prior_file_path)){
given_prior=fread(paste0(provided_prior_file_path,"/",provided_prior_file_name),data.table=F)
given_prior=as.matrix(given_prior)
PRIOR[2:(nrow(given_prior)+1),1:ncol(given_prior)]=given_prior
}else{

iteration_pos=1
for(i in 1:(length(random_lev_number)+1)){  #随机效应的数量 

	for(j in 1:sum(random_total_lev_number%in%i)){
	
		for(k in 1:j){
			iteration_pos=iteration_pos+1
			PRIOR[iteration_pos,1]=i
			PRIOR[iteration_pos,2]=j
			PRIOR[iteration_pos,3]=k
			if(j!=k){
			PRIOR[iteration_pos,4]=0.5
			}else {PRIOR[iteration_pos,4]=1}
	}}}
}	
}

#残差限定

for(i in 1:length(residual_cov_trait)){

residual_cov_combination=residual_cov_trait[[i]]
residual_number=residual_number

if(!is.null(residual_cov_combination)){

residual_cov_combination_pos=match(residual_cov_combination,target_trait_name)

trait1_pos=residual_cov_combination_pos[1]
trait2_pos=residual_cov_combination_pos[2]
PRIOR=PRIOR[!(PRIOR[,1]%in%residual_number&PRIOR[,2]%in%trait1_pos&PRIOR[,3]%in%trait2_pos),]
PRIOR=PRIOR[!(PRIOR[,1]%in%residual_number&PRIOR[,2]%in%trait2_pos&PRIOR[,3]%in%trait1_pos),]		
}}


#不显示 PRIOR结果， 因为针对随机回归模型，PRIOR结果存在些许问题，使用DMU默认值即可
#如果不是人为提供PRIOR文件的话，PRIOR设置为NULL
if(!is.null(provided_prior_file_name)&!is.null(provided_prior_file_path)){

}else{
PRIOR=NULL
}


# module_parameter

if(dmu_module=="dmuai"){

PARAMETER[1,1]="$DMUAI"
PARAMETER[2,1]=EM_iter
PARAMETER[3,1]=iteration_criteria
PARAMETER[4,1]="1.0d-6"
PARAMETER[5,1]="1"
PARAMETER[6,1]="0"

}else if(dmu_module=="dmu5"){
PARAMETER[1,1]="$DMU5"
PARAMETER[2,2]=paste0(dmu5_iter," 1e-9")
PARAMETER[3,3]="512"

}else if(dmu_module=="rjmc"){
PARAMETER[1,1]="$RJMC"
PARAMETER[2,1]="3 0"
PARAMETER[3,1]=paste(gibbs_seeds,collapse=" ")
PARAMETER[4,1]=paste(gibbs_n_burnin,gibbs_n_rounds,gibbs_n_gaps,collapse=" ")

if(!is.null(PRIOR)){
n_var=length(unique(PRIOR[,1]))
}else{
n_var=length(unique(random_lev))+1
}
PARAMETER[5:(5+n_var-1),1]=paste(1:n_var,6*length(target_trait_name))
}






#输出结果
#DIR=rbind(COMMENT,ANALYSE,DATA,VARIABLE,VAR_STR,MODEL,SOLUTION,PRIOR,PARAMETER)
#DIR=DIR[!DIR[,1]%in%"",]

COMMENT=COMMENT[!COMMENT[,1]%in%"",]
ANALYSE=ANALYSE[!ANALYSE[,1]%in%"",]
DATA=DATA[!DATA[,1]%in%"",]
VARIABLE=VARIABLE[!VARIABLE[,1]%in%"",]
VAR_STR=VAR_STR[!VAR_STR[,1]%in%"",]
MODEL=MODEL[!MODEL[,1]%in%"",]
SOLUTION=SOLUTION[!SOLUTION[,1]%in%"",]
PRIOR=PRIOR[!PRIOR[,1]%in%"",]
PARAMETER=PARAMETER[!PARAMETER[,1]%in%"",]

DIR=rbind(COMMENT,"",ANALYSE,"",DATA,"",VARIABLE,"",VAR_STR,"",MODEL,"",SOLUTION,"",PRIOR,"",PARAMETER)

setwd(output_DIR_path)
if(is.null(output_DIR_name)){output_DIR_name="Trait.DIR"}
write.table(DIR,output_DIR_name,quote=F,row.names=F,col.names=F,sep=" ")
return(list(DIR=union_dir(DIR),random_effect_name=random_effect_name))
}


union_dir<-function(data){
union_data=matrix(NA,nrow=nrow(data),ncol=1)
for(i in 1:nrow(data)){
union_data[i,1]=paste(data[i,],collapse=" ")
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
Variable_n=length(phe_col_names)
real_n=Variable_n-integer_n
Integer_col_names=phe_col_names[1:integer_n]
Real_col_names=phe_col_names[(integer_n+1):Variable_n]

res_Integer_col_names=setdiff(Integer_col_names,c(genetic_effect_name,social_effect_group_name)) #除遗传和group之外的整型变量名称
res_Real_col_names=Real_col_names #剩余的实型变量名称


#统计最大的group_size
group_size=as.data.frame(table(phe[,social_effect_group_name]),stringsAsFactors=F)
group_size[,1]=as.numeric(group_size[,1])
max_size=max(group_size[,2])-1
colnames(group_size)=c("Group","Size")

if(max_size==0){stop("The number of indivuals in each group is equal to 1, please check your group record!!!")}


#构建group_ind_set和group的对应列表
group_ind_set <- vector(mode = "list", length = nrow(group_size))
group_i_set=as.numeric(group_size$Group)
for(i in 1:nrow(group_size)){
group_ind_set[[i]]=phe[phe[,social_effect_group_name]%in%group_i_set[i],genetic_effect_name]
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
phe=phe[,c(genetic_effect_name,res_Integer_col_names,social_effect_group_name,integer_group_names,
            res_Real_col_names,real_group_names)]

return(list(phe=phe,integer_n=integer_n+length(integer_group_names),
			integer_group_names=integer_group_names,real_group_names=real_group_names))			  
}
