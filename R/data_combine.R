
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
input=get_input_data_type(input_data_type=input_data_type,input_data_hmp=input_data_hmp,input_data_plink_ped=input_data_plink_ped,
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
for(i in 1:length(input_data_name)){
data_map=fread(paste0(path,"/",input_data_name[i],".map"),data.table=F,header=F)
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
for(i in 1:length(input_data_name)){
cat("Reading plink data......\n")
data_ped_subset=fread(paste0(path,"/",input_data_name[i],".ped"),data.table=F,header=F,nThread=cpu_cores)
data_map_subset=fread(paste0(path,"/",input_data_name[i],".map"),data.table=F,header=F)
cat(paste0("Plink data file path",plink_iter,": ",path," ;\nPlink data file name",plink_iter,": ",input_data_name[i],"\n"))
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
for(i in 1:length(input_data_name)){
cat("Reading hapmap data......\n")
data_hmp=fread(paste0(path,"/",input_data_name[i],".hmp.txt"),data.table=F,header=T)
cat(paste0("Hapmap data file path",hmp_iter,": ",path," ;\nHapmap data file name",hmp_iter,": ",input_data_name[i],"\n"))
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


