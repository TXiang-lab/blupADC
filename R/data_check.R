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
input=get_input_data_type(input_data_type=input_data_type,input_data_hmp=input_data_hmp,input_data_plink_ped=input_data_plink_ped,
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

overlap_matrix=NULL
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
overlap_matrix=numeric_overlap_cpp(numeric=data_numeric_subset,dis_progress=TRUE,overlap_threshold=duplication_threshold,cpu_cores=cpu_cores)
overlap_matrix=overlap_matrix[overlap_matrix[,3]>duplication_threshold,]

#根据overlap_matrix的初步结果，利用所有SNP数据,再重新进行overlap计算
if(nrow(overlap_matrix)>500){stop("Plesase select appropriate numbers of selected snps in detecting overlap ")}
if(nrow(overlap_matrix)>=1){
for(i in 1:nrow(overlap_matrix)){
overlap_matrix[i,3]=sum(data_numeric[overlap_matrix[i,1],]==data_numeric[overlap_matrix[i,2],])/ncol(data_numeric)
}}

overlap_matrix=overlap_matrix[as.numeric(overlap_matrix[,3])>=duplication_threshold,]

if(nrow(overlap_matrix)>=1){
message(paste0("Found ",nrow(overlap_matrix)," duplicate genotypes in your chip data"))
IND1=IND_geno[overlap_matrix[,1]]
IND2=IND_geno[overlap_matrix[,2]]
overlap_matrix[,1]=IND1
overlap_matrix[,2]=IND2
if(!is.null(output_data_path)){
setwd(output_data_path)
write.table(overlap_matrix,"Duplicated_Genotype_Percentage.txt",quote=F,row.names=F,,col.names=F,sep="\t")}
overlap_IND=unique(c(overlap_matrix[,1],overlap_matrix[,2]))

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
return(list(sum_data=sum_data,duplicate_genotype=overlap_matrix,breed=pca_outlier))
}
}



