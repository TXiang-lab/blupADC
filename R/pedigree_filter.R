#
trace_pedigree<-function(input_pedigree=NULL,  #第一列为个体号,第二列为父亲,第三列为母亲
						 input_pedigree_path=NULL,
						 input_pedigree_name=NULL,
						 multi_col=FALSE,#不对系谱数据进行重新整理
						 missing_value=NULL,
						 dup_error_check=TRUE,
						 sex_error_check=TRUE,
						 breed_error_check=FALSE,
						 birth_date_error_check=FALSE,
						 parent_error_check=TRUE,
						 trace_id=NULL,   #默认追溯所有系谱中的个体
						 trace_sibs=FALSE,#是否追溯同胞
						 trace_fullsibs=FALSE, #追溯全同胞
						 max_fullsibs=NULL,  #追溯全同胞时，全同胞最大数量
						 trace_direction="backward", #backward 为反向
						 #trace_offspring=FALSE, #追溯个体的子代
						 trace_generation=NULL, #追溯的代数	
						 trace_birth_date=NULL, #仅追溯出生日期在给定日期之后的个体
						 output_rename_pedigree=TRUE, #输出rename后的系谱
						 output_pedigree_path=NULL,
						 output_pedigree_name=NULL,
						 output_pedigree_tree=FALSE,
						 pedigree_tree_depth=3,
						 display_message=TRUE,
						 priority_rename_id=NULL,
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
pedigree=data.table::fread(paste0(input_pedigree_path,"/",input_pedigree_name),data.table=F)
}else{
pedigree=input_pedigree
}


if(ncol(pedigree)==3){
if(display_message==TRUE){	
cat("Peidgree provided has three columns,please make sure the format of pedigree_data has three columns: Offspring Sire Dam  \n")
}
}else if(ncol(pedigree)==4){
if(display_message==TRUE){
cat("Peidgree provided has four columns,please make sure the format of pedigree_data has four columns: Offspring Sire Dam  Birth_Date \n")
}
}else if(multi_col==TRUE){
if(display_message==TRUE){
cat("Peidgree provided has multiple columns,please make sure the format of pedigree_data similar to: Offspring Sire Dam  SireSire SireDam ......  \n")
}
pedigree=pedigree_format_convertion(pedigree)
}else{
stop("Error:peidgree_provided is not standard format!")
}

#将缺失值变为NA
if(!is.null(missing_value)){
pedigree[pedigree==missing_value]=NA
}else{
pedigree[pedigree=="-9999"]=NA
pedigree[pedigree=="0"]=NA
pedigree[pedigree==""]=NA
}

##将系谱变成长系谱，确保第一列包含所有个体
#id_F=na.omit(unique(pedigree[!pedigree[,2]%in%pedigree[,1],2]))
#id_M=na.omit(unique(pedigree[!pedigree[,3]%in%pedigree[,1],3]))
#if(length(id_F)==0){id_F=NA}
#if(length(id_M)==0){id_M=NA}
#
#pedigree=rbind(as.matrix(pedigree[,1:3]),cbind(id_F,NA,NA),cbind(id_M,NA,NA))
#
#if(ncol(input_pedigree)==4){
#pedigree=cbind(pedigree,input_pedigree[match(pedigree[,1],input_pedigree[,1]),4])
#colnames(pedigree)[4]="Birth_Date"}

ped=as.matrix(pedigree)
if(!mode(ped)=="character"){ped=apply(ped,2,as.character)}
colnames(ped)[1:3]=c("Offspring","Sire","Dam")
ped=ped[!is.na(ped[,1]),]
#检查系谱错误
if(display_message==TRUE){
cat("Start checking pedigree......\n")
}
ped_record_check_data=ped_record_check(ped=ped, dup_error_check=dup_error_check, sex_error_check=sex_error_check,
						                        breed_error_check=breed_error_check, birth_date_error_check=birth_date_error_check,
												parent_error_check=parent_error_check,display_message=display_message)
if(display_message==TRUE){cat("Finish checking pedigree......\n")}
ped=ped_record_check_data$ped

error_duplicated_id=ped_record_check_data$error_duplicated_id
error_sex_id=ped_record_check_data$error_sex_id
error_breed_id=ped_record_check_data$error_breed_id
error_birth_date_id=ped_record_check_data$error_birth_date_id
error_parent_id=ped_record_check_data$error_parent_id

rename_ped=NULL;rename_phenotype=NULL
#系谱排序, 数字化
if(!is.null(trace_birth_date)){
if(display_message==TRUE){cat("Please make sure_pedigree data has four columns \n")}
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
#trace_id_offspring=ped[ped[,2]%in%trace_id|ped[,3]%in%trace_id,1]
#trace_id=unique(c(trace_id,trace_id_offspring))
#}
#
#trace_id=unique(na.omit(c(trace_id,trace_id_parent))) #追溯给定个体的父母
#ped=ped[ped[,1]%in%trace_id,]
#}
ped[ped==0]=NA
#指定个体追溯系谱
if(!is.null(trace_id)){                     #根据给定的个体id, trace所有和这些个体id有联系的个体的系谱(默认trace这些个体以及它们的父母)
if(display_message==TRUE){cat("Tracing porvided id......\n")}
i=0
trace_type="Match"
trace_id_set=trace_id

trace_id_sibs=NULL
#找寻 trace_id 这一代的同胞
if(trace_fullsibs==TRUE){
trace_id_sibs=match_fullsibs(trace_id,ped,max_fullsibs)
}


while(trace_type=="Match"){
i=i+1
if(!is.null(trace_generation)){
if(i>=trace_generation)break
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
if(trace_direction=="backward"){ped[is.na(ped)]="0"}
rename_pedigree=ped

base_result=single_pedigree_cpp(ped)

IND_base=base_result[[1]]
IND_middle=base_result[[2]]
IND_offspring=base_result[[3]]






if(trace_direction=="backward"){

#########################################################################################
#2022.3.17                                                                              #
#在IND_offspring中，有些个体是另一些个体父母的同胞，具体问题复现可以根据 example_ped2来观察  #
#解决办法：找到IND_middle个体的同胞，将其从IND_offspring中去除                          #
IND_middle_sibs=match_all_sibs(IND_middle,ped)                                                    #
offsprint_sibs=IND_offspring[IND_offspring%in%IND_middle_sibs]                                      #
if(length(offsprint_sibs)>0){                                                                               #
IND_middle=c(IND_middle,offsprint_sibs)                                                                   #
IND_offspring=setdiff(IND_offspring,offsprint_sibs)                                                    #
}                                                                                                           #
#########################################################################################

final_ped=data.frame(Offspring=c(IND_base,IND_middle,IND_offspring),Generation="NA",stringsAsFactors=F)
final_ped[final_ped[,1]%in%IND_base,2]=paste0("A",sprintf("%5s",0))
final_ped[final_ped[,1]%in%IND_offspring,2]=paste0("B",sprintf("%5s",1001))

#接下来的目的是将 IND_middle 中的个体划分出  父母与子代
tmp_ped=ped
i=0
IND_base_set=IND_base
IND_middle_set=IND_middle
IND_offspring_set=IND_offspring
#test_tmp_base_result=tmp_base_result #当循环之间的结果是一样的时候，通过这个来跳出循环
while(length(IND_middle)>0){
if(i>1000){stop("Found generations of provided pedigree larger than 1000, please check your pedigree carefully!")}
i=i+1
#tmp_ped=tmp_ped[tmp_ped[,2]%in%IND_middle|tmp_ped[,3]%in%IND_middle,]
tmp_ped=tmp_ped[tmp_ped[,1]%in%IND_middle,]
tmp_base_result=single_pedigree(tmp_ped)
#if(identical(test_tmp_base_result,tmp_base_result)){break}#当循环之间的结果是一样的时候，通过这个来跳出循环
#test_tmp_base_result=tmp_base_result
tmp_IND_base=setdiff(tmp_base_result[[1]],IND_base_set);IND_base_set=c(IND_base_set,tmp_IND_base)
tmp_IND_middle=tmp_base_result[[2]];
tmp_IND_offspring=tmp_base_result[[3]];

final_ped[final_ped[,1]%in%tmp_IND_base,2]=paste0("A",sprintf("%5s",i))
final_ped[final_ped[,1]%in%tmp_IND_offspring,2]=paste0("B",sprintf("%5s",1000-i))

IND_middle=tmp_IND_middle
}

final_ped$Generation=as.numeric(as.factor(final_ped$Generation))

final_ped$Generation=final_ped$Generation-1

}else if(trace_direction=="forward"){

final_ped=data.frame(Offspring=c(IND_base,IND_middle,IND_offspring),stringsAsFactors=F)
final_ped$Sire=ped[match(final_ped[,1],ped[,1]),2]
final_ped$Dam=ped[match(final_ped[,1],ped[,1]),3]
IND_base_set=IND_base
pos=get_offspring_generation_cpp(final_ped,IND_base_set)
final_ped$Generation=pos
final_ped=final_ped[,c(1,4)]


}else{
stop("Non-standard input of trace_direction parameter")
}
final_ped=final_ped[final_ped[,1]!=0,]
generation=max(final_ped$Generation)-min(final_ped$Generation)+1
final_ped$Generation=final_ped$Generation-(min(final_ped$Generation)-1)

if(display_message==TRUE){cat(paste0("Found ",generation," generations in provided pedigree \n"))}

final_ped=final_ped[!duplicated(final_ped[,1]),]
final_ped=final_ped[!final_ped[,1]%in%NA,]
options (warn =-1)	
if(NA%in%as.numeric(final_ped[,1])){
final_ped=final_ped[order(final_ped[,2],final_ped[,1]),]
}else{
final_ped=final_ped[order(final_ped[,2],as.numeric(final_ped[,1])),]
}
options (warn =1)	

if(!is.null(trace_generation)&is.null(trace_id)){
trace_generation=as.numeric(trace_generation)
max_generation=max(as.numeric(final_ped[,2]))
min_generation=max_generation-trace_generation+1

if(min_generation< 0){
if(display_message==TRUE){cat(paste0("Attention:the max generation of this_pedigree is ",max_generation," \n"))}
min_generation=0
}

final_ped=final_ped[final_ped[,2]%in%c(min_generation:max_generation),]
final_ped[,2]=final_ped[,2]-min_generation
}

#set these animals have smallest renamed_id
if(!is.null(priority_rename_id)){
final_ped=rbind(final_ped[final_ped[,1]%in%priority_rename_id,],final_ped[!final_ped[,1]%in%priority_rename_id,])
}

final_ped=cbind(final_ped,1:nrow(final_ped))


colnames(final_ped)[3]="Offspring_Id"


for(i in 1:ncol(rename_pedigree)){

rename_pedigree[,i]=final_ped[match(rename_pedigree[,i],final_ped[,"Offspring"]),"Offspring_Id"]

}



final_ped=cbind(final_ped,match_parents(as.character(as.numeric(final_ped[,"Offspring_Id"])),rename_pedigree)[[1]])
final_ped=cbind(final_ped,1:nrow(final_ped));
colnames(final_ped)[4:6]=c("Sire_Id","Dam_Id","Order")
final_ped[,4]=as.numeric(as.character(final_ped[,4]));
final_ped[,5]=as.numeric(as.character(final_ped[,5]))


Unrename_ped=final_ped[,c(1,4,5,2)];Unrename_ped[,2]=final_ped[match(Unrename_ped[,2],final_ped[,3]),1]
Unrename_ped[,3]=final_ped[match(Unrename_ped[,3],final_ped[,3]),1]
Unrename_ped[is.na(Unrename_ped)]=0
colnames(Unrename_ped)=c("Offspring","Sire","Dam","Generation")

if(!is.null(output_pedigree_path)&!is.null(output_pedigree_name)){
setwd(output_pedigree_path)


if(output_rename_pedigree==TRUE){
final_ped[is.na(final_ped)]=0
data.table::fwrite(data.table(final_ped[,-c(1:2)]),paste0(output_pedigree_name,".txt"),col.names=F,quote=F,sep=" ",row.names=F)
tmp=final_ped[,c(1,3)]
colnames(tmp)=c("Original_Id","Renamed_Id")
data.table::fwrite(data.table(tmp),paste0(output_pedigree_name,"_renamed_key.txt"),col.names=T,quote=F,sep=" ",row.names=F)
}else{
if(display_message==TRUE){cat("Output ordered but non-renamed pedigree! \n")}
tmp=Unrename_ped
tmp[is.na(tmp)]=0
tmp[,4]=1:nrow(tmp)
data.table::fwrite(data.table(tmp),paste0(output_pedigree_name,".txt"),col.names=F,quote=F,sep=" ",row.names=F)
}
}
final_ped[is.na(final_ped)]=0
if(sum((final_ped$Sire_Id-final_ped$Offspring_Id)>0)>=1 |sum((final_ped$Dam_Id-final_ped$Offspring_Id)>0)>=1){
stop("Rename error: sire_id or dam_id larger than offspring_id \n")
}

if(min(final_ped[,2]==0))final_ped[,2]=final_ped[,2]+1
#追溯全系谱
#将系谱 rename排序后，
#依次从上到下，找个体的祖先，
#利用迭代的思想，因为系谱在上面的个体是系谱在下面个体的祖先，如果上面个体的祖先已知，那么下面个体的祖先就也是已知的。
full_generation_matrix=NULL
if(output_pedigree_tree==TRUE){
if(display_message==TRUE){cat("Constructing pedigree_tree......\n")}

#original_ped=final_ped
#original_ped[,4]=original_ped[match(original_ped[,4],original_ped[,3]),1]
#original_ped[,5]=original_ped[match(original_ped[,5],original_ped[,3]),1]
#original_ped[is.na(original_ped)]="0"
#original_ped=as.matrix(original_ped[,3:5])
#original_ped=apply(original_ped,2,as.character)

generation_names=get_generation_colnames(generation=pedigree_tree_depth)

Unrename_ped[is.na(Unrename_ped)]="0"
Unrename_ped=as.matrix(Unrename_ped)
Unrename_ped=apply(Unrename_ped,2,as.character)


full_generation_matrix=full_generation_conversion(generation_names=generation_names,ped=Unrename_ped[,1:3])

colnames(full_generation_matrix)=c("Offspring",generation_names)

}

#统计同胞数目
group_full_sibs=NULL
family_size_full_sibs=NULL
if(summary_sibs==TRUE){

#统计全同胞数目
full_sibs_pedigree=data.frame(Unrename_ped,stringsAsFactors=F)
full_sibs_pedigree$MF=paste0(full_sibs_pedigree[,2],"_",full_sibs_pedigree[,3])
full_sibs_pedigree=full_sibs_pedigree[!(full_sibs_pedigree[,2]%in%0|full_sibs_pedigree[,3]%in%0),]
full_sibs_com=aggregate(full_sibs_pedigree,by=list(full_sibs_pedigree$MF),length)
group_full_sibs=nrow(full_sibs_com)
family_size_full_sibs=mean(full_sibs_com[,2])
}



return(list(ped=Unrename_ped,rename_ped=final_ped,pedigree_tree=full_generation_matrix,
				   error_id_set=list(error_duplicated_id=error_duplicated_id,error_sex_id=error_sex_id,
			         error_breed_id=error_breed_id,error_birth_date_id=error_birth_date_id,error_parent_id=error_parent_id),
				sibs=list(group_full_sibs=group_full_sibs,family_size_full_sibs=family_size_full_sibs)))

}



#检查系谱错误
ped_record_check<-function(ped,  #matrix format
                           dup_error_check=TRUE, 
						   sex_error_check=TRUE,
						   breed_error_check=FALSE,
						   parent_error_check=TRUE,
						   birth_date_error_check=FALSE,
						   display_message=TRUE)
{
#将缺失值变为NA
#ped[ped=="-9999"]=NA
#ped[ped=="0"]=NA
#ped[ped==""]=NA
ped=as.matrix(ped)
error_duplicated_id=NULL;error_sex_id=NULL;error_breed_id=NULL;error_birth_date_id=NULL;error_parent_id=NULL;

#错误0：个体既出在子代又出现在父母列



#错误1：个体名相同，但是父母至少有一个不相同
if(dup_error_check==TRUE){
duplicated_id=unique(ped[duplicated(ped[,1]),1])

error_duplicated_id=NULL

if(length(duplicated_id)>0){
#需要改写成c++
for(id in duplicated_id){
duplicated_set=ped[ped[,1]%in%id,]
if(length(unique(duplicated_set[,2]))>=2| length(unique(duplicated_set[,3]))>=2){
error_duplicated_id=c(error_duplicated_id,id)
}
}

if(length(error_duplicated_id)>0){
if(display_message==TRUE){cat(paste0("Found ",length(error_duplicated_id)," duplicated id error records: offsprings with  same id but have different sire or dam records, records of sire and dam would be treated as missing value \n"))}
ped[ped[,1]%in%error_duplicated_id,2:3]=NA
}
ped=ped[!duplicated(ped[,1]),]  #选择非重复的行
}
}


#错误2  性别有误，个体既出现在父亲这一列，又出现在母亲这一列
if(sex_error_check==TRUE){
error_sex_id=unique(na.omit(c(ped[,2][ped[,2]%in%ped[,3]],ped[ped[,1]==ped[,2]|ped[,1]==ped[,3],1])))
if(length(error_sex_id)>0){

if(display_message==TRUE){
cat(paste0("Found ",length(error_sex_id)," sex error records: ids in the sire column also appear in the dam column or offspring column, these ids would be treated as missing value \n"))
}

ped[ped[,2]%in%error_sex_id,2:3]=NA
ped[ped[,3]%in%error_sex_id,2:3]=NA
}}

#错误3  品种错误，对于某些系谱记录来说，个体ID 和 父母ID 记录有品种信息，因此可根据品种记录找出父母和子女品种记录不一致的个体，仅针对纯种
#默认品种信息记录在Id的开始两位，eg. YY201201
if(breed_error_check==TRUE){

breed_individual=substr(ped[,1],1,2)
breed_sire=substr(ped[,2],1,2)
breed_dam=substr(ped[,3],1,2)

result_breed=(breed_individual==breed_sire) & (breed_individual==breed_dam) #逻辑值：父母的品种记录均和后代一致

if(sum(result_breed==FALSE)>=1){
if(display_message==TRUE){
cat("Found breed error: the breed of offspring, sire and dam are not equal  \n")
}
error_breed_id=ped[!result_breed,1]
ped[!result_breed,2:3]=NA
}
}

#错误4 根据出生日期的记录，判断后代的出生日期是否早于父母的出生日期，日期需要是数字格式，如：20180204
if(birth_date_error_check==TRUE){  

if(ncol(ped)<4){
stop("Error: peidgree provided doesn't have birth date records, please turn off birth_date_check function")
}

check_ped=ped[!(is.na(ped[,2])&is.na(ped[,3])),]  #去除父母为缺失的个体

ordered_offspring=check_ped[,4][order(as.numeric(check_ped[,4]))]

offspring_order=(1:length(ordered_offspring))[match(check_ped[,4],ordered_offspring)]

sire_order=offspring_order[match(check_ped[,2],check_ped[,1])]

dam_order=offspring_order[match(check_ped[,3],check_ped[,1])]

birth_date_result=(offspring_order>sire_order)|(offspring_order>dam_order) #逻辑值

if(sum(na.omit(birth_date_result==FALSE))>=1){
error_birth_date_id=check_ped[!birth_date_result,1]
error_birth_date_id=na.omit(error_birth_date_id)
ped[ped[,1]%in%error_birth_date_id,2:3]=NA
if(display_message==TRUE){
cat(paste0("Found ",length(error_birth_date_id)," birth date error records: the birth date of offspring is before than the birth date of sire and dam \n"))
}
}
}

#错误5： 父母
if(parent_error_check==TRUE){

offspring_id=ped[,1]
parent_id=match_parents(offspring_id,ped)[[1]]
sire_id=parent_id[,1]
dam_id=parent_id[,2]

sire_id_parents=match_parents(sire_id,ped)[[1]]
sire_id_sire=as.character(sire_id_parents[,1])
sire_id_dam=sire_id_parents[,2]

dam_id_parents=match_parents(dam_id,ped)[[1]]
dam_id_sire=dam_id_parents[,1]
dam_id_dam=dam_id_parents[,2]



pos1=sire_id_sire==offspring_id
pos2=sire_id_dam==offspring_id
pos3=dam_id_sire==offspring_id
pos4=dam_id_dam==offspring_id

if(TRUE%in%pos1){error_parent_id=c(error_parent_id,na.omit(offspring_id[pos1]))}
if(TRUE%in%pos2){error_parent_id=c(error_parent_id,na.omit(offspring_id[pos2]))}
if(TRUE%in%pos3){error_parent_id=c(error_parent_id,na.omit(offspring_id[pos3]))}
if(TRUE%in%pos4){error_parent_id=c(error_parent_id,na.omit(offspring_id[pos4]))}

if(length(error_parent_id)>0){
if(display_message==TRUE){
cat(paste0("Found ",length(error_parent_id)," parent id error records: the relationships between offsprings and parents have conflicts, records of sire and dam would be treated as missing value \n"))
}
ped[ped[,1]%in%error_parent_id,2:3]=NA
}
ped=ped[!duplicated(ped[,1]),]  #选择非重复的行
}


return(list(ped=ped,error_duplicated_id=error_duplicated_id,error_sex_id=error_sex_id,
			         error_breed_id=error_breed_id,error_birth_date_id=error_birth_date_id,error_parent_id=error_parent_id))
}


match_parents<-function(id,ped){ # ped 为三列向量（个体 父亲 母亲）， id 为一列向量 

parents=matrix(NA,ncol=2,nrow=length(id))

parents[,1]=ped[match(id,ped[,1]),2]
parents[,2]=ped[match(id,ped[,1]),3]
match_state=ifelse(sum(is.na(c(parents[,1],parents[,2])))==c(length(id)*2),"Complete_Unmatch","Match")
return(c(list(parents),list(match_state)))
}



#寻找个体的子代
match_offsprings_new<-function(id=NULL,ped=NULL){
offspring=NULL
if(!NA%in%id){
offspring=na.omit(unique(c(ped[ped[,2]%in%id,1],ped[ped[,3]%in%id,1])))
}
return(offspring)
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
get_generation_colnames<-function(generation=1){
generation_colnames=c("Sire","Dam")
total=generation_colnames

if(generation>=2){
for(i in 1:(generation-1)){
temp=NULL
for(i in 1:length(generation_colnames)){
temp=c(temp,paste0(generation_colnames[i],c("Sire","Dam")))
}
generation_colnames=temp
total=c(total,generation_colnames)
}}

return(total)
}


 
pedigree_format_convertion<-function(pedigree){

#check pedigree names
col=colnames(pedigree)
max_generation=round(nchar(col[1])/3)
max_generation_names=c("Offspring",get_generation_colnames(max_generation))


if(sum(col%in%max_generation_names==FALSE)>0){
error_colnames=col[!col%in%max_generation_names]
stop(paste0("Colnames of provided pedigree_data, including << ",error_colnames, ">> doesn't meet the requirement, please modify these columns names!"))
}

Total_result=list(pedigree[,1:3])
match_parents_column_name<-function(name,pedigree){
part=pedigree[,colnames(pedigree)%in%name]

if(paste0(name,"Sire")%in%colnames(pedigree)){
part=cbind(part,pedigree[,colnames(pedigree)%in%paste0(name,"Sire")])
}else{part=cbind(part,NA)}

if(paste0(name,"Dam")%in%colnames(pedigree)){
part=cbind(part,pedigree[,colnames(pedigree)%in%paste0(name,"Dam")])
}else{part=cbind(part,NA)}
colnames(part)=c("Offspring","Sire","Dam")
return(part)
}

for(i in colnames(pedigree)){
if(paste0(i,"Sire")%in%colnames(pedigree) | paste0(i,"Dam")%in%colnames(pedigree)){
result=list(match_parents_column_name(i,pedigree))
Total_result=c(Total_result,result)
}
}
pedigree=do.call(rbind,Total_result)
pedigree=pedigree[!duplicated(paste0(pedigree[,1],pedigree[,2],pedigree[,3])),]
if(display_message==TRUE){cat("Complete multi-coulums pedigree_data convert into standard 3 columns pedigree_data! \n")}
return(pedigree)
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
for(i in 1:max_fullsibs){
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



#找寻个体的全同胞合半同胞
match_all_sibs<-function(id=NULL,ped=NULL){
id_F=ped[ped[,1]%in%id,2]
id_M=ped[ped[,1]%in%id,3]
if(length(id_F)==0){id_F=NA}
if(length(id_M)==0){id_M=NA}
id_P=na.omit(c(id_F,id_M))
sibs=ped[ped[,2]%in%id_P|ped[,3]%in%id_P,1]

return(setdiff(sibs,id))
}