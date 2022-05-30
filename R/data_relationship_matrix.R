
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


if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
setwd(output_matrix_path)
}

#判断输入类型
if("Homozygous"%in%inbred_type){kinship_type=c(kinship_type,"G_A")}
if("G_Diag"%in%inbred_type){kinship_type=c(kinship_type,"G_A")}
if("Pedigree"%in%inbred_type){kinship_type=c(kinship_type,"P_A")}
if("H_diag"%in%inbred_type){kinship_type=c(kinship_type,"H_A")}


#基因型数据
if("G_A"%in%kinship_type | "G_Ainv"%in%kinship_type |"G_D"%in%kinship_type |"G_Dinv"%in%kinship_type| "H_A"%in%kinship_type| "H_Ainv"%in%kinship_type|"H_D"%in%kinship_type|"H_Dinv"%in%kinship_type){ 
#程序自动判断输入类型
input=get_input_data_type(input_data_type=input_data_type,input_data_hmp=input_data_hmp,input_data_plink_ped=input_data_plink_ped,
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
stop("Error:In terms of constring col_three matrix, provided pedigree id include character, please set IND_rename=TRUE !")
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

if(sum(IND_geno%in%IND_pedigree)==0){stop("Provided genotype id is not accordance with the pedigree id!")}
if(IND_rename==TRUE){
IND_geno=num_ped[match(IND_geno,IND_pedigree),1]
IND_pedigree=num_ped[,1]
}

if("col_three"%in%output_matrix_type){
if(NA %in% as.numeric(IND_geno)){
stop("Error:In terms of constring col_three matrix,, provided genotype id include character, please set the genotype id as integer or provide pedigree for recoding the genotype id!")
}
IND_geno=as.numeric(IND_geno)
IND_pedigree=as.numeric(IND_pedigree)
}

}else if(is.null(input_pedigree)&!is.null(input_data_numeric)){

IND_geno=rownames(input_data_numeric)
if("col_three"%in%output_matrix_type){
if(NA %in% as.numeric(IND_geno)){
stop("Error:In terms of constring col_three matrix, provided genotype id include character, please set the genotype id as integer or provide pedigree for recoding the genotype id!")
}
IND_geno=as.numeric(IND_geno)
}
}


G_A=NULL;G_Ainv=NULL;G_D=NULL;G_Dinv=NULL;P_A=NULL;P_Ainv=NULL;P_D=NULL;P_Dinv=NULL;H_A=NULL;H_Ainv=NULL;H_D=NULL;H_Dinv=NULL;
homo_inbred=NULL;diag_inbred=NULL;pedigree_inbred=NULL;H_inbred=NULL

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


if(!is.null(P_A)){

if("col_three" %in% output_matrix_type){
P_A_three=matrix_col3(P_A,IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_pedigree),"IND_pedigree.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(P_A_three),paste0(output_matrix_name,"P_A_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(P_A_three);gc();}
}

if("col_all" %in% output_matrix_type){
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_pedigree),"IND_pedigree.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(P_A),paste0(output_matrix_name,"P_A_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(P_A);gc();}
}

}


if(!is.null(P_Ainv)){

if("col_three" %in% output_matrix_type){
P_Ainv_three=matrix_col3(P_Ainv,IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_pedigree),"IND_pedigree.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(P_Ainv_three),paste0(output_matrix_name,"P_Ainv_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(P_Ainv_three);gc();}
}

if("col_all" %in% output_matrix_type){
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_pedigree),"IND_pedigree.txt",quote=F,row.names=F,col.names=F,sep="\t")
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
if(!is.null(P_D)){

if("col_three" %in% output_matrix_type){
P_D_three=matrix_col3(P_D,IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_pedigree),"IND_pedigree.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(P_D_three),paste0(output_matrix_name,"P_D_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(P_D_three);gc();}
}

if("col_all" %in% output_matrix_type){
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_pedigree),"IND_pedigree.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(P_D),paste0(output_matrix_name,"P_D_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(P_D);gc();}
}

}


if(!is.null(P_Dinv)){

if("col_three" %in% output_matrix_type){
P_Dinv_three=matrix_col3(P_Dinv,IND_geno=as.numeric(IND_pedigree),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_pedigree),"IND_pedigree.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(P_Dinv_three),paste0(output_matrix_name,"P_Dinv_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(P_Dinv_three);gc();}
}

if("col_all" %in% output_matrix_type){
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_pedigree),"IND_pedigree.txt",quote=F,row.names=F,col.names=F,sep="\t")
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

if(!is.null(G_A)){

if("col_three" %in% output_matrix_type){
G_A_three=matrix_col3(G_A,IND_geno=as.numeric(IND_geno),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_geno),"IND_geno.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(G_A_three),paste0(output_matrix_name,"G_A_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(G_A_three);gc();}
}

if("col_all" %in% output_matrix_type){
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_geno),"IND_geno.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(G_A),paste0(output_matrix_name,"G_A_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(G_A);gc();}
}

}


if(!is.null(G_Ainv)){

if("col_three" %in% output_matrix_type){
G_Ainv_three=matrix_col3(G_Ainv,IND_geno=as.numeric(IND_geno),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_geno),"IND_geno.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(G_Ainv_three),paste0(output_matrix_name,"G_Ainv_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(G_Ainv_three);gc();}
}

if("col_all" %in% output_matrix_type){
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_geno),"IND_geno.txt",quote=F,row.names=F,col.names=F,sep="\t")
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

if(!is.null(G_D)){

if("col_three" %in% output_matrix_type){
G_D_three=matrix_col3(G_D,IND_geno=as.numeric(IND_geno),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_geno),"IND_geno.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(G_D_three),paste0(output_matrix_name,"G_D_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(G_D_three);gc();}
}

if("col_all" %in% output_matrix_type){
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_geno),"IND_geno.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(G_D),paste0(output_matrix_name,"G_D_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(G_D);gc();}
}

}


if(!is.null(G_Dinv)){

if("col_three" %in% output_matrix_type){
G_Dinv_three=matrix_col3(G_Dinv,IND_geno=as.numeric(IND_geno),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_geno),"IND_geno.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(G_Dinv_three),paste0(output_matrix_name,"G_Dinv_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(G_Dinv_three);gc();}
}

if("col_all" %in% output_matrix_type){
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_geno),"IND_geno.txt",quote=F,row.names=F,col.names=F,sep="\t")
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
H_inbred=data.frame(Id=IND_Additive,Diag_inbred=H_A_inbred,stringsAsFactors=F)
H_A=H_Aresult$H
H_Ainv=NULL
rm(H_Aresult);gc();
}else if("H_A"%in%kinship_type&("H_Ainv"%in%kinship_type)){
H_Aresult=makeHA_cpp(num_ped,input_data_numeric,IND_geno,pos_A11,pos_A22,pos_geno,pos_A,pos_H22,APY_algorithm,APY_eigen_threshold,APY_n_core,direct=TRUE,inverse=TRUE,omega=SSBLUP_omega)
H_A_inbred=H_Aresult$inbred
H_inbred=data.frame(Id=IND_Additive,Diag_inbred=H_A_inbred,stringsAsFactors=F)
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
H_inbred=data.frame(Id=IND_Additive,Diag_inbred=H_A_inbred,stringsAsFactors=F)
H_A=H_Aresult$H
H_Ainv=NULL
rm(H_Aresult);gc();
}else if("H_A"%in%kinship_type&("H_Ainv"%in%kinship_type)){
H_Aresult=makeHA_metafounder_cpp(num_ped,input_data_numeric,pos_A11,pos_A22,pos_geno,pos_A,pos_H22,direct=TRUE,inverse=TRUE,omega=SSBLUP_omega)
H_A_inbred=H_Aresult$inbred
H_inbred=data.frame(Id=IND_Additive,Diag_inbred=H_A_inbred,stringsAsFactors=F)
H_A=H_Aresult$H
H_Ainv=H_Aresult$Hinv
rm(H_Aresult);gc();
}else if(!"H_A"%in%kinship_type&("H_Ainv"%in%kinship_type)){
H_Ainv=makeHA_metafounder_cpp(num_ped,input_data_numeric,pos_A11,pos_A22,pos_geno,pos_A,pos_H22,direct=FALSE,inverse=TRUE,omega=SSBLUP_omega)$Hinv
H_A=NULL
H_A_inbred=NULL
}

}

if(!is.null(H_A)){

if("col_three" %in% output_matrix_type){
H_A_three=matrix_col3(H_A,IND_geno=as.numeric(IND_Additive),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_Additive),"IND_SSBLUP.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(H_A_three),paste0(output_matrix_name,"H_A_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(H_A_three);gc();}
}

if("col_all" %in% output_matrix_type){
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_Additive),"IND_SSBLUP.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(H_A),paste0(output_matrix_name,"H_A_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(H_A);gc();}
}

}


if(!is.null(H_Ainv)){

if("col_three" %in% output_matrix_type){
H_Ainv_three=matrix_col3(H_Ainv,IND_geno=as.numeric(IND_Additive),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_Additive),"IND_SSBLUP.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(H_Ainv_three),paste0(output_matrix_name,"H_Ainv_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(H_Ainv_three);gc();}
}

if("col_all" %in% output_matrix_type){
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_Additive),"IND_SSBLUP.txt",quote=F,row.names=F,col.names=F,sep="\t")
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

if(!is.null(H_D)){

if("col_three" %in% output_matrix_type){
H_D_three=matrix_col3(H_D,IND_geno=as.numeric(IND_Additive),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_Additive),"IND_SSBLUP.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(H_D_three),paste0(output_matrix_name,"H_D_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(H_D_three);gc();}
}

if("col_all" %in% output_matrix_type){
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_Additive),"IND_SSBLUP.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(H_D),paste0(output_matrix_name,"H_D_col_all.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(H_D);gc();}
}

}


if(!is.null(H_Dinv)){

if("col_three" %in% output_matrix_type){
H_Dinv_three=matrix_col3(H_Dinv,IND_geno=as.numeric(IND_Additive),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col3_threshold) 
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_Additive),"IND_SSBLUP.txt",quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(H_Dinv_three),paste0(output_matrix_name,"H_Dinv_col_three.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}
if(return_result==FALSE){rm(H_Dinv_three);gc();}
}

if("col_all" %in% output_matrix_type){
if((!is.null(output_matrix_path)&!is.null(output_matrix_name))){
fwrite(data.table(IND_Additive),"IND_SSBLUP.txt",quote=F,row.names=F,col.names=F,sep="\t")
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
#H_inbred=data.frame(Id=c(IND_A11,IND_A22),H_inbreeding=H_inbred,stringsAsFactors=F)
if(!is.null(inbred_output_matrix_path)){
setwd(inbred_output_matrix_path)
cat("Saving inbreeding coefficients......\n")
if(Metafounder_algorithm==TRUE){
fwrite(data.table(H_inbred),"H_metafounder_inbreeding_coefficients.txt",quote=F,sep="\t",row.names=F,col.names=T)
}else {
fwrite(data.table(H_inbred),"H_inbreeding_coefficients.txt",quote=F,sep="\t",row.names=F,col.names=T)}
}
}else{
H_inbred=NULL
}

if(!exists("G_A")){G_A=NULL};if(!exists("G_Ainv")){G_Ainv=NULL};
if(!exists("G_D")){G_D=NULL};if(!exists("G_Dinv")){G_Dinv=NULL};
if(!exists("P_A")){P_A=NULL};if(!exists("P_Ainv")){P_Ainv=NULL};
if(!exists("P_D")){P_D=NULL};if(!exists("P_Dinv")){P_Dinv=NULL};
if(!exists("H_A")){H_A=NULL};if(!exists("H_Ainv")){H_Ainv=NULL};
if(!exists("H_D")){H_D=NULL};if(!exists("H_Dinv")){H_Dinv=NULL};

if(!exists("G_A_three")){G_A_three=NULL};if(!exists("G_Ainv_three")){G_Ainv_three=NULL};
if(!exists("G_D_three")){G_D_three=NULL};if(!exists("G_Dinv_three")){G_Dinv_three=NULL};
if(!exists("P_A_three")){P_A_three=NULL};if(!exists("P_Ainv_three")){P_Ainv_three=NULL};
if(!exists("P_D_three")){P_D_three=NULL};if(!exists("P_Dinv_three")){P_Dinv_three=NULL};
if(!exists("H_A_three")){H_A_three=NULL};if(!exists("H_Ainv_three")){H_Ainv_three=NULL};
if(!exists("H_D_three")){H_D_three=NULL};if(!exists("H_Dinv_three")){H_Dinv_three=NULL};


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
if(!is.null(H_D)){rownames(H_D)=colnames(H_D)=IND_Additive}
if(!is.null(H_Dinv)){rownames(H_Dinv)=colnames(H_Dinv)=IND_Additive}


if(return_result==TRUE){
return(list(G_A=list(A=G_A,Ainv=G_Ainv,A_col3=G_A_three,Ainv_col3=G_Ainv_three),
			G_D=list(D=G_D,Dinv=G_Dinv,D_col3=G_D_three,Dinv_col3=G_Dinv_three),
			P_A=list(A=P_A,Ainv=P_Ainv,A_col3=P_A_three,Ainv_col3=P_Ainv_three),
			P_D=list(D=P_D,Dinv=P_Dinv,D_col3=P_D_three,Dinv_col3=P_Dinv_three),
		    H_A=list(A=H_A,Ainv=H_Ainv,A_col3=H_A_three,Ainv_col3=H_Ainv_three),
			H_D=list(D=H_D,Dinv=H_Dinv,D_col3=H_D_three,Dinv_col3=H_Dinv_three),
		    Inbred=list(Homozygous=homo_inbred,G_diag=diag_inbred,Pedigree=pedigree_inbred,H_diag=H_inbred),
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
P_A_three=matrix_col3_memory(P_A@address,paste0(bigmemory_data_name,"_P_A"),bigmemory_data_path,
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
P_Ainv_three=matrix_col3_memory(pBigMat=P_Ainv@address,paste0(bigmemory_data_name,"_P_Ainv"),bigmemory_data_path,
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
P_D_three=matrix_col3_memory(pBigMat=P_D@address,paste0(bigmemory_data_name,"_P_D"),bigmemory_data_path,
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
P_Dinv_three=matrix_col3_memory(pBigMat=P_Dinv@address,paste0(bigmemory_data_name,"_P_Dinv"),bigmemory_data_path,
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
G_A_three=matrix_col3_memory(pBigMat=G_A@address,paste0(bigmemory_data_name,"_G_A"),bigmemory_data_path,
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
G_Ainv_three=matrix_col3_memory(pBigMat=G_Ainv@address,paste0(bigmemory_data_name,"_G_Ainv"),bigmemory_data_path,
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
G_D_three=matrix_col3_memory(pBigMat=G_D@address,paste0(bigmemory_data_name,"_G_D"),bigmemory_data_path,
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
G_Dinv_three=matrix_col3_memory(pBigMat=G_Dinv@address,paste0(bigmemory_data_name,"_G_Dinv"),bigmemory_data_path,
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
delete_bigmemory_file("temp_A22_inv",bigmemory_data_name,bigmemory_data_path,FALSE);
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
H_inbred=data.frame(Id=IND_Additive,Diag_inbred=H_A_inbred,stringsAsFactors=F)
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
H_inbred=data.frame(Id=IND_Additive,Diag_inbred=H_A_inbred,stringsAsFactors=F)
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
H_inbred=data.frame(Id=IND_Additive,Diag_inbred=H_A_inbred,stringsAsFactors=F)
H_A=H_Aresult$H
H_Ainv=NULL
rm(H_Aresult);gc();
}else if("H_A"%in%kinship_type&("H_Ainv"%in%kinship_type)){
H_Aresult=K_matrix_cal_memory(numeric_address=input_data_numeric@address,bigmemory_data_name,bigmemory_data_path,
					IND_geno=IND_geno,Pedigree=num_ped,pos_A11=pos_A11,pos_A22=pos_A22,
					pos_geno=pos_geno,pos_A=pos_A,pos_H22=pos_H22,
					H_A_direct=TRUE,omega=SSBLUP_omega,type=6,inverse=TRUE)
H_A_inbred=H_Aresult$inbred
H_inbred=data.frame(Id=IND_Additive,Diag_inbred=H_A_inbred,stringsAsFactors=F)
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
delete_bigmemory_file("temp_A22_inv",bigmemory_data_name,bigmemory_data_path,FALSE);
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
H_A_three=matrix_col3_memory(pBigMat=H_A@address,bigmemory_data_name=paste0(bigmemory_data_name,"_H_A"),bigmemory_data_path,
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
H_Ainv_three=matrix_col3_memory(pBigMat=H_Ainv@address,bigmemory_data_name=paste0(bigmemory_data_name,"_H_Ainv"),bigmemory_data_path,
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
#H_inbred=data.frame(Id=c(IND_A11,IND_A22),H_inbreeding=H_inbred,stringsAsFactors=F)
if(!is.null(inbred_output_matrix_path)){
setwd(inbred_output_matrix_path)
cat("Saving inbreeding coefficients......\n")
if(Metafounder_algorithm==TRUE){
fwrite(data.table(H_inbred),"H_metafounder_inbreeding_coefficients.txt",quote=F,sep="\t",row.names=F,col.names=T)
}else {
fwrite(data.table(H_inbred),"H_inbreeding_coefficients.txt",quote=F,sep="\t",row.names=F,col.names=T)}
}
}else{
H_inbred=NULL
}

if(!exists("G_A")){G_A=NULL};if(!exists("G_Ainv")){G_Ainv=NULL};
if(!exists("G_D")){G_D=NULL};if(!exists("G_Dinv")){G_Dinv=NULL};
if(!exists("P_A")){P_A=NULL};if(!exists("P_Ainv")){P_Ainv=NULL};
if(!exists("P_D")){P_D=NULL};if(!exists("P_Dinv")){P_Dinv=NULL};
if(!exists("H_A")){H_A=NULL};if(!exists("H_Ainv")){H_Ainv=NULL};
if(!exists("H_D")){H_D=NULL};if(!exists("H_Dinv")){H_Dinv=NULL};

if(!exists("G_A_three")){G_A_three=NULL};if(!exists("G_Ainv_three")){G_Ainv_three=NULL};
if(!exists("G_D_three")){G_D_three=NULL};if(!exists("G_Dinv_three")){G_Dinv_three=NULL};
if(!exists("P_A_three")){P_A_three=NULL};if(!exists("P_Ainv_three")){P_Ainv_three=NULL};
if(!exists("P_D_three")){P_D_three=NULL};if(!exists("P_Dinv_three")){P_Dinv_three=NULL};
if(!exists("H_A_three")){H_A_three=NULL};if(!exists("H_Ainv_three")){H_Ainv_three=NULL};
if(!exists("H_D_three")){H_D_three=NULL};if(!exists("H_Dinv_three")){H_Dinv_three=NULL};


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
return(list(G_A=list(A=G_A,Ainv=G_Ainv,A_col3=G_A_three,Ainv_col3=G_Ainv_three),
			G_D=list(D=G_D,Dinv=G_Dinv,D_col3=G_D_three,Dinv_col3=G_Dinv_three),
			P_A=list(A=P_A,Ainv=P_Ainv,A_col3=P_A_three,Ainv_col3=P_Ainv_three),
			P_D=list(D=P_D,Dinv=P_Dinv,D_col3=P_D_three,Dinv_col3=P_Dinv_three),
		    H_A=list(A=H_A,Ainv=H_Ainv,A_col3=H_A_three,Ainv_col3=H_Ainv_three),
			H_D=list(D=H_D,Dinv=H_Dinv,D_col3=H_D_three,Dinv_col3=H_Dinv_three),
		    Inbred=list(Homozygous=homo_inbred,G_diag=diag_inbred,Pedigree=pedigree_inbred,H_diag=H_inbred),
			phased_block=data_block))
}

}

}



#solve要比 ginv快，如果solve报错，再调用 MASS::ginv求广义逆
my.solve<-function(Matrix){

	Matrix_inv<- try(chol2inv(chol(Matrix)), silent = TRUE)
	if(class(Matrix_inv) == "try-error"){
	cat("Added small value to diagnoal to ensure matrix invertible \n")
	Matrix_inv <- try(chol2inv(chol(Matrix+diag(0.0001,nrow(Matrix)))),silent=TRUE)
	if(class(Matrix_inv) == "try-error"){
	cat("Using GINV To Calculate Inverse Matrix \n")
	Matrix_inv <- my.solve_old(Matrix)
	}
	}
	return(Matrix_inv)
}



get_input_data_type<-function(input_data_type=NULL,input_data_hmp=NULL,input_data_plink_ped=NULL,input_data_plink_map=NULL,
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
