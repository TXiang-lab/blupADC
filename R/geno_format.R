
geno_format<-function(             
		input_data_type=NULL, #"Plink" , "Hapmap" , "VCF" , "BLUPF90","Numeric"  "Haplotype"
		input_data_path=NULL,
		input_data_name=NULL,  #hapmap 为全名，vcf 和 plink 默认不包括后缀，为 name.vcf name.ped name.map 
		input_data_hmp=NULL,
		input_data_plink_ped=NULL,
		input_data_plink_map=NULL,
		input_data_blupf90=NULL,
		input_data_blupf90_map=NULL, #格式和 Haplotype-map格式一样
		input_data_numeric=NULL,	 
		input_data_numeric_map=NULL, #格式和 Haplotype-map格式一样
		input_data_vcf=NULL,
		input_data_haplotype_hap=NULL,
		input_data_haplotype_map=NULL,
		input_data_haplotype_sample=NULL,
		bigmemory_cal=FALSE,
		bigmemory_data_type="integer",
		bigmemory_data_path=getwd(),
		bigmemory_data_name="blupADC",
		phased_symbol="/",
		phased_genotype=FALSE,
		haplotype_window_nSNP=NULL,
		haplotype_window_kb=NULL,
		haplotype_window_block=NULL,	
		miss_base=NULL,		
		miss_base_num=0,
		output_data_name="blupADC_geno_convertion",
		output_data_type=NULL,  # "Plink" , "Hapmap" , "VCF"  "BLUPF90"  "Numeric"  "Haplotype"
		output_data_path=NULL,
		return_result=FALSE,
		cpu_cores=1  #调用的cpu数目，用于加速计算
		){	
options (warn =-1)	
library(data.table)
library(bigmemory)
options(bigmemory.allow.dimnames=TRUE)
options(bigmemory.typecast.warning=FALSE)
if(!is.null(output_data_path)){
if(!file.exists(output_data_path)){
dir.create(output_data_path,recursive=T)
cat("Output data path is not exist,software will constructs this path  automatically \n")
}
setwd(output_data_path)
}

output_data_type=unique(output_data_type)
if(bigmemory_cal==TRUE){return_result=TRUE}
#程序自动判断输入类型
	
input=get_input_data_type(input_data_type=input_data_type,input_data_hmp=input_data_hmp,input_data_plink_ped=input_data_plink_ped,
                                          input_data_plink_map=input_data_plink_map,input_data_blupf90=input_data_blupf90,
								    input_data_numeric=input_data_numeric,input_data_vcf=input_data_vcf,
									input_data_haplotype_hap=input_data_haplotype_hap,input_data_haplotype_map=input_data_haplotype_map,
									input_data_haplotype_sample=input_data_haplotype_sample,
								     miss_base=miss_base)

input_data_type=input$input_data_type
miss_base=input$miss_base

if(identical(input_data_type,output_data_type)){return(NULL)}
##############################################初始化

data_hmp=NULL;data_ped=NULL;data_map=NULL;data_blupf90=NULL;data_numeric=NULL;data_vcf=NULL
data_haplotype_hap=NULL;data_haplotype_map=NULL;data_haplotype_sample=NULL;output_block=NULL
haplotype_allele=NULL;

#判断是否存在本地文件
if(!is.null(input_data_path)&!is.null(input_data_name)){
if(input_data_type=="Plink"){
if(!file.exists(paste0(input_data_path,"/",input_data_name,".ped"))|!file.exists(paste0(input_data_path,"/",input_data_name,".map"))){
stop("Couldn't found provided input data!")}

}else if(input_data_type=="Haplotype"){
if(!file.exists(paste0(input_data_path,"/",input_data_name,".hap"))|!file.exists(paste0(input_data_path,"/",input_data_name,".sample"))|!file.exists(paste0(input_data_path,"/",input_data_name,".map"))){
stop("Couldn't found provided input data!")}
}else{
if(!file.exists(paste0(input_data_path,"/",input_data_name))){stop("Couldn't found provided input data!")}
}
}

#判断输入的bigmeory类型
if(!bigmemory_data_type%in%(c("double","float","integer","short","char"))){

stop("Unknow input bigmemory_data_type, bigmemory_data_type should be double or float or integer or short or char ")

}


####################################input data type is Hapmap#################################### 
if(input_data_type=="Hapmap"){
if(!is.null(input_data_path)&!is.null(input_data_name)){
cat("Start read the Hapmap format genotype data \n")
input_data_hmp=fread(paste0(input_data_path,"/",input_data_name),data.table=F,header=T)
cat("Complete read the Hapmap format genotype data  \n")}

if(!is.matrix(input_data_hmp)){input_data_hmp=as.matrix(input_data_hmp)}
input_data_hmp[,3]=as.numeric(input_data_hmp[,3])
input_data_hmp[,4]=as.numeric(input_data_hmp[,4])

IND_name=colnames(input_data_hmp)[-c(1:11)]
n_snp=nrow(input_data_hmp)
n_ind=ncol(input_data_hmp)-11

output_type_number=get_output_type_number(input_data_type=input_data_type,
										  output_data_type=output_data_type,
										  bigmemory_cal=bigmemory_cal,
										  phased_genotype=phased_genotype)

if(length(output_type_number)==0){stop("Please provided standard output data type!")}

for(type in output_type_number){
if(exists(".Random.seed")){fixed_openMP(.Random.seed);set.seed(19980204)}

if(type==1){
data_ped=hapmap_convertion(input_data_hmp,type,miss_base,phased_symbol,miss_base_num,cpu_cores)$ped;
data_map=input_data_hmp[,c(3,1,4,4)];data_map[,3]=as.numeric(data_map[,3])/1000000;
data_map[,3]=as.numeric(data_map[,3])
data_map[,4]=as.numeric(data_map[,4])
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_ped),paste0(output_data_name,".ped"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(data_map),paste0(output_data_name,".map"),quote=F,row.names=F,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_map,data_ped);gc();}
}

if(type==2){data_vcf=hapmap_convertion(input_data_hmp,type,miss_base,phased_symbol,miss_base_num,cpu_cores)$vcf;
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(get_vcf_header()),paste0(output_data_name,".vcf"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(data_vcf),paste0(output_data_name,".vcf"),append=TRUE,quote=F,row.names=F,col.names=T,sep="\t")}
if(return_result==FALSE){rm(data_vcf);gc();}
}

if(type==3){data_numeric=hapmap_convertion(input_data_hmp,type,miss_base,phased_symbol,miss_base_num,cpu_cores)$numeric
rownames(data_numeric)=IND_name
if(!is.null(output_data_path)&!is.null(output_data_name)){
fwrite(data.frame(data_numeric),paste0(output_data_name,".numeric"),quote=F,row.names=T,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_numeric);gc();}
}

if(type==4){data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=hapmap_convertion(input_data_hmp,type,miss_base,phased_symbol,miss_base_num,cpu_cores)$blupf90,stringsAsFactors=F)
if(!is.null(output_data_path)&!is.null(output_data_name)){
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=F,sep=" ")}
if(return_result==FALSE){rm(data_blupf90);gc();}
}

#构建bigmemory对象 
if(type==13){

if(exists("pBigMat_num_address")){rm(pBigMat_num_address)}
gc()
delete_bigmemory_file("numeric",bigmemory_data_name,bigmemory_data_path,FALSE);

pBigMat_num_address=make_bigmemory_object_address_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
																
hapmap_to_numeric_memory_cpp(input_data_hmp,pBigMat_num_address,miss_base,miss_base_num,cpu_cores,type=bigmemory_data_type)									 
data_numeric=make_bigmemory_object_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)	
options(bigmemory.allow.dimnames=TRUE)
rownames(data_numeric)=IND_name
colnames(data_numeric)=paste0("SNP",1:ncol(data_numeric))																
cat(paste0("bigmemory-Numeric data is save as ",bigmemory_data_path,"/",bigmemory_data_name,"_numeric \n"))
																
}


}
rm(input_data_hmp);gc();
}


####################################input data type is Plink#################################### 
if(input_data_type=="Plink"){
if(!is.null(input_data_path)&!is.null(input_data_name)){
cat("Start read the Plink format genotype data \n")
input_data_plink_ped=fread(paste0(input_data_path,"/",input_data_name,".ped"),data.table=F,header=F)
input_data_plink_map=fread(paste0(input_data_path,"/",input_data_name,".map"),data.table=F,header=F)
cat("Complete read the Plink format genotype data \n")}

IND_name=input_data_plink_ped[,2]
n_ind=length(IND_name)
n_snp=nrow(input_data_plink_map)

if(!is.matrix(input_data_plink_map)){input_data_plink_map=as.matrix(input_data_plink_map)}
if(!is.matrix(input_data_plink_ped)){input_data_plink_ped=as.matrix(input_data_plink_ped)}
input_data_plink_map[,1]=as.numeric(input_data_plink_map[,1])
input_data_plink_map[,3]=as.numeric(input_data_plink_map[,3])
input_data_plink_map[,4]=as.numeric(input_data_plink_map[,4])

output_type_number=get_output_type_number(input_data_type=input_data_type,
										  output_data_type=output_data_type,
										  bigmemory_cal=bigmemory_cal,
										  phased_genotype=phased_genotype)

if(length(output_type_number)==0){stop("Please provided standard output data type!")}

for(type in output_type_number){
if(exists(".Random.seed")){fixed_openMP(.Random.seed);set.seed(19980204)}
if(type==1){data_hmp=plink_convertion(input_data_plink_ped,input_data_plink_map,	type, miss_base,phased_symbol,miss_base_num,cpu_cores)$hmp
data_hmp[,1]=input_data_plink_map[,2]
data_hmp[,3]=input_data_plink_map[,1]
data_hmp[,4]=input_data_plink_map[,4]
data_hmp[,3]=as.numeric(data_hmp[,3])
data_hmp[,4]=as.numeric(data_hmp[,4])
data_hmp[,c(2,5:11)]="NA"
colnames(data_hmp)=c("rs#","alleles","chrom","pos","strand","assembly","center","protLSID","assayLSID","panelLSID","QCcode",IND_name)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.table(data_hmp),paste0(output_data_name,".hmp.txt"),quote=F,row.names=F,col.names=T,sep="\t")}
if(return_result==FALSE){rm(data_hmp);gc();}
}

if(type==2){data_vcf=plink_convertion(input_data_plink_ped,input_data_plink_map,	type, miss_base,phased_symbol,miss_base_num,cpu_cores)$vcf
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(get_vcf_header()),paste0(output_data_name,".vcf"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(data_vcf),paste0(output_data_name,".vcf"),append=TRUE,quote=F,row.names=F,col.names=T,sep="\t")}
if(return_result==FALSE){rm(data_vcf);gc();}
}

if(type==3){data_numeric=plink_convertion(input_data_plink_ped,input_data_plink_map,	type, miss_base,phased_symbol,miss_base_num,cpu_cores)$numeric
rownames(data_numeric)=IND_name
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_numeric),paste0(output_data_name,".numeric"),quote=F,row.names=T,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_numeric);gc();}
}

if(type==4){data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=plink_convertion(input_data_plink_ped,input_data_plink_map,	type, miss_base,phased_symbol,miss_base_num,cpu_cores)$blupf90,stringsAsFactors=F)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=F,sep=" ")}
if(return_result==FALSE){rm(data_blupf90);gc();}
}


#构建bigmemory对象 
if(type==13){

if(exists("pBigMat_num_address")){rm(pBigMat_num_address)}
gc()
delete_bigmemory_file("numeric",bigmemory_data_name,bigmemory_data_path,FALSE);

pBigMat_num_address=make_bigmemory_object_address_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
																
plink_to_numeric_memory_cpp(input_data_plink_ped,input_data_plink_map,pBigMat_num_address,cpu_cores,miss_base,miss_base_num,type=bigmemory_data_type)								 
data_numeric=make_bigmemory_object_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)	
options(bigmemory.allow.dimnames=TRUE)
rownames(data_numeric)=IND_name
colnames(data_numeric)=paste0("SNP",1:ncol(data_numeric))																
cat(paste0("bigmemory-Numeric data is save as ",bigmemory_data_path,"/",bigmemory_data_name,"_numeric \n"))
																
}


}
rm(input_data_plink_ped,input_data_plink_map);gc();
}


####################################input data type is VCF#################################### 
if(input_data_type=="VCF"){
if(!is.null(input_data_path)&!is.null(input_data_name)){
cat("Start read the VCF format genotype data \n")
input_data_vcf=fread(paste0(input_data_path,"/",input_data_name),data.table=F)
cat("Complete read the VCF format genotype data \n")}

if(!is.matrix(input_data_vcf)){input_data_vcf=as.matrix(input_data_vcf)}

input_data_vcf[,1]=as.numeric(input_data_vcf[,1])
input_data_vcf[,2]=as.numeric(input_data_vcf[,2])

IND_name=colnames(input_data_vcf)[-c(1:9)]
input_data_haplotype_map=input_data_vcf[,c(1,3,2,4,5)]


#构建block
if(!is.null(haplotype_window_nSNP)|!is.null(haplotype_window_kb)|!is.null(haplotype_window_block)){
block_result=get_haplotype_block(input_data_haplotype_map,haplotype_window_nSNP,
				                    haplotype_window_kb,haplotype_window_block)
block_hap=block_result$block_rcpp
block_start=block_hap[,1];block_end=block_hap[,2]
output_block=block_result$block
}else{
block_start=1:5
block_end=1:5
output_block=NULL
}

output_type_number=get_output_type_number(input_data_type=input_data_type,
										  output_data_type=output_data_type,
										  bigmemory_cal=bigmemory_cal,
										  phased_genotype=phased_genotype)
n_ind=ncol(input_data_vcf)-9
n_snp=nrow(input_data_vcf)
n_window=length(block_start)
data_haplotype_map=input_data_vcf[,c(1,3,2,4,5)];IND_name=colnames(input_data_vcf)[-c(1:9)]

if(length(output_type_number)==0){stop("Please provided standard output data type!")}
for(type in output_type_number){
if(exists(".Random.seed")){fixed_openMP(.Random.seed);set.seed(19980204)}
if(type==1){data_hmp=vcf_convertion(input_data_vcf,block_start,block_end,type,miss_base,miss_base_num,cpu_cores)$hmp
data_hmp[,1]=input_data_vcf[,3]
data_hmp[,3]=input_data_vcf[,1]
data_hmp[,4]=input_data_vcf[,2]
data_hmp[,c(2,5:11)]="NA"
data_hmp[,3]=as.numeric(data_hmp[,3])
data_hmp[,4]=as.numeric(data_hmp[,4])
colnames(data_hmp)=c("rs#","alleles","chrom","pos","strand","assembly","center","protLSID","assayLSID","panelLSID","QCcode",IND_name)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.table(data_hmp),paste0(output_data_name,".hmp.txt"),quote=F,row.names=F,col.names=T,sep="\t")}
if(return_result==FALSE){rm(data_hmp);gc();}
}

if(type==2){data_ped=vcf_convertion(input_data_vcf,block_start,block_end,type,miss_base,miss_base_num,cpu_cores)$ped
data_map=input_data_vcf[,c(1,3,2,2)];data_map[,3]=(as.numeric(data_map[,3])/1000000)
data_map[,3]=as.numeric(data_map[,3])
data_map[,4]=as.numeric(data_map[,4])
data_ped[,3:6]=0
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_ped),paste0(output_data_name,".ped"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(data_map),paste0(output_data_name,".map"),quote=F,row.names=F,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_map,data_ped);gc();}
}


if(type==8){
tmp=vcf_convertion(input_data_vcf,block_start,block_end,type,miss_base,miss_base_num,cpu_cores)
data_ped=tmp$ped
haplotype_allele=tmp$SNP
origin_haplotype_allele=haplotype_allele
rm(tmp);gc();

haplotype_allele_num=sapply(haplotype_allele,length)
block_pos_allele=output_block$block_pos
block_pos_allele$chr=input_data_haplotype_map[block_pos_allele[,2],1]
block_pos_allele$pos=paste0(input_data_haplotype_map[block_pos_allele[,2],3],"_",
						     input_data_haplotype_map[block_pos_allele[,3],3])

haplotype_allele_chr=rep(block_pos_allele$chr,times=haplotype_allele_num)
haplotype_allele_pos=rep(block_pos_allele$pos,times=haplotype_allele_num)

haplotype_allele=do.call(c,haplotype_allele)
haplotype_allele=data.frame(Chr=haplotype_allele_chr,
							SNP=paste0(1:length(haplotype_allele)),
							cM=0,
							Pos=haplotype_allele_pos,
							Allele=haplotype_allele)

data_map=haplotype_allele[,1:4]
haplotype_allele=origin_haplotype_allele
data_ped[,3:6]=0
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_ped),paste0(output_data_name,".ped"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(data_map),paste0(output_data_name,".map"),quote=F,row.names=F,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_map,data_ped);gc();}
}

if(type==3){data_haplotype_hap=vcf_convertion(input_data_vcf,block_start,block_end,type,miss_base,miss_base_num,cpu_cores)$hap
data_haplotype_map=input_data_vcf[,c(1,3,2,4,5)];IND_name=colnames(input_data_vcf)[-c(1:9)]
data_haplotype_sample=data.frame(IND_name)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_haplotype_hap),paste0(output_data_name,"_haplotype.hap"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(data_haplotype_map),paste0(output_data_name,"_haplotype.map"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(data_haplotype_sample),paste0(output_data_name,"_haplotype.sample"),quote=F,row.names=F,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_haplotype_hap,data_haplotype_map,data_haplotype_sample);gc();}
}

if(type==4){

tmp=vcf_convertion(input_data_vcf,block_start,block_end,type,miss_base,miss_base_num,cpu_cores)
haplotype_allele=tmp$SNP
data_numeric=tmp$numeric
rm(tmp);gc();
rownames(data_numeric)=IND_name

if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_numeric),paste0(output_data_name,".numeric"),quote=F,row.names=T,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_numeric);gc();}
}

if(type==5){
tmp=vcf_convertion(input_data_vcf,block_start,block_end,type,miss_base,miss_base_num,cpu_cores)
haplotype_allele=tmp$SNP
data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=tmp$blupf90,stringsAsFactors=F)
rm(tmp);gc();
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=F,sep=" ")}
if(return_result==FALSE){rm(data_blupf90);gc();}
}

if(type==6){data_numeric=vcf_convertion(input_data_vcf,block_start,block_end,type,miss_base,miss_base_num,cpu_cores)$numeric
rownames(data_numeric)=IND_name
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_numeric),paste0(output_data_name,".numeric"),quote=F,row.names=T,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_numeric);gc();}
}

if(type==7){data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=vcf_convertion(input_data_vcf,block_start,block_end,type,miss_base,miss_base_num,cpu_cores)$blupf90,stringsAsFactors=F)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=F,sep=" ")}
if(return_result==FALSE){rm(data_blupf90);gc();}
}
#构建bigmemory对象 
if(type%in%c(13,14,15)){
delete_bigmemory_file("haplotype_hap",bigmemory_data_name,bigmemory_data_path,FALSE);
pBigMat_hap_address=make_bigmemory_object_address_cpp(nrow=n_snp,ncol=n_ind*2,file_name=paste0(bigmemory_data_name,"_haplotype_hap"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
}

if(type%in%c(14,15)){
delete_bigmemory_file("haplotype_list",bigmemory_data_name,bigmemory_data_path,FALSE);

if(type==14){
delete_bigmemory_file("numeric",bigmemory_data_name,bigmemory_data_path,TRUE)
pBigMat_list_address=make_bigmemory_object_address_cpp(nrow=n_window,ncol=n_ind*2,file_name=paste0(bigmemory_data_name,"_haplotype_list"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
}
}

if(type%in%c(16)){
delete_bigmemory_file("numeric",bigmemory_data_name,bigmemory_data_path,TRUE);

pBigMat_num_address=make_bigmemory_object_address_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
}


if(type==13){
phased_vcf_to_haplotype_memory_cpp(input_data_vcf, pBigMat_hap_address, cpu_cores)
rm(pBigMat_hap_address);gc();
}


if(type==14){
haplotype_allele=phased_vcf_to_numeric_memory_cpp(pBigMat_hap_address,pBigMat_list_address,block_start,block_end,
									 numeric_file_name=paste0(bigmemory_data_name,"_numeric"),numeric_file_path=bigmemory_data_path,
									 input_data_vcf,cpu_cores)
rm(pBigMat_hap_address,pBigMat_list_address);gc();									 
}

if(type==15){	
tmp=phased_vcf_blupf90_memory_cpp(pBigMat_hap_address,pBigMat_list_address,block_start,block_end,input_data_vcf,cpu_cores)
haplotype_allele=tmp$haplotype_allele
data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=tmp$allele_string_vector,stringsAsFactors=F)									 
rm(tmp,pBigMat_hap_address,pBigMat_list_address);gc();	
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=T,sep=" ")}
}

if(type==16){unphased_vcf_to_numeric_memory_cpp(pBigMat_num_address,input_data_vcf,cpu_cores,miss_base,miss_base_num)								 
rm(pBigMat_num_address);gc();	
}
}
rm(input_data_vcf);gc();
#删除bigmemory中间文件

delete_bigmemory_file("haplotype_list",bigmemory_data_name,bigmemory_data_path,FALSE);

if(!13%in%output_type_number){

delete_bigmemory_file("haplotype_hap",bigmemory_data_name,bigmemory_data_path,TRUE);

}


if(13%in%output_type_number){
data_haplotype_hap=make_bigmemory_object_cpp(nrow=n_snp,ncol=n_ind*2,file_name=paste0(bigmemory_data_name,"_haplotype_hap"),
                                                               file_path=bigmemory_data_path,type=bigmemory_data_type)
data_haplotype_sample=data.frame(IND_name)
cat(paste0("bigmemory-Haplotype-hap data is save as ",bigmemory_data_path,"/",bigmemory_data_name,"_haplotype_hap \n"))
}

if(14%in%output_type_number|16%in%output_type_number){									 
data_numeric=make_bigmemory_object_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)	
options(bigmemory.allow.dimnames=TRUE)
rownames(data_numeric)=IND_name
colnames(data_numeric)=paste0("SNP",1:ncol(data_numeric))																
cat(paste0("bigmemory-Numeric data is save as ",bigmemory_data_path,"/",bigmemory_data_name,"_numeric \n"))
}

}
####################################input data type is Haplotype#################################### 
if(input_data_type=="Haplotype"){
if(!is.null(input_data_path)&!is.null(input_data_name)){
cat("Start read the Haplotype format genotype data \n")
input_data_haplotype_hap=fread(paste0(input_data_path,"/",input_data_name,".hap"),header=F,data.table=F)
input_data_haplotype_sample=fread(paste0(input_data_path,"/",input_data_name,".sample"),header=F,data.table=F)
input_data_haplotype_map=fread(paste0(input_data_path,"/",input_data_name,".map"),header=F,data.table=F)
input_data_haplotype_map[,1]=as.numeric(input_data_haplotype_map[,1])
input_data_haplotype_map[,3]=as.numeric(input_data_haplotype_map[,3])
cat("Complete read the Haplotype format genotype data \n")}

if(is.data.frame(input_data_haplotype_hap)){input_data_haplotype_hap=DataFrame_to_arma(input_data_haplotype_hap)}
if(typeof(input_data_haplotype_hap)!="integer"){input_data_haplotype_hap=NumericMatrix_to_arma(input_data_haplotype_hap)}

IND_name=as.character(input_data_haplotype_sample[,1])
#构建block
if((phased_genotype==TRUE)&all(is.null(haplotype_window_nSNP),is.null(haplotype_window_kb),is.null(haplotype_window_block))){
stop("Once set phased_genotype==TRUE,  haplotype_window_nSNP or haplotype_window_kb or haplotype_window_block must be provided!")
}

if(!is.null(haplotype_window_nSNP)|!is.null(haplotype_window_kb)|!is.null(haplotype_window_block)){
block_result=get_haplotype_block(input_data_haplotype_map,haplotype_window_nSNP,
				                    haplotype_window_kb,haplotype_window_block)
block_hap=block_result$block_rcpp
block_start=block_hap[,1];block_end=block_hap[,2]
output_block=block_result$block
}else{
block_start=1:5
block_end=1:5
output_block=NULL
}

n_ind=ncol(input_data_haplotype_hap)/2
n_snp=nrow(input_data_haplotype_hap)
n_window=length(block_start)


output_type_number=get_output_type_number(input_data_type=input_data_type,
										  output_data_type=output_data_type,
										  bigmemory_cal=bigmemory_cal,
										  phased_genotype=phased_genotype)

if(length(output_type_number)==0){stop("Please provided standard output data type!")}

for(type in output_type_number){
if(exists(".Random.seed")){fixed_openMP(.Random.seed);set.seed(19980204)}
if(type==1){
tmp=haplotype_convertion(input_data_haplotype_hap,block_start,block_end,IND_name,type,cpu_cores)
haplotype_allele=tmp$SNP
data_numeric=tmp$numeric
rm(tmp);gc();
rownames(data_numeric)=IND_name
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_numeric),paste0(output_data_name,".numeric"),quote=F,row.names=T,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_numeric);gc();}
}

if(type==2){
tmp=haplotype_convertion(input_data_haplotype_hap,block_start,block_end,IND_name,type,cpu_cores)
haplotype_allele=tmp$SNP
data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=tmp$blupf90,stringsAsFactors=F)
rm(tmp);gc();
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=F,sep=" ")}
if(return_result==FALSE){rm(data_blupf90);gc();}
}

if(type==3){

tmp=haplotype_convertion(input_data_haplotype_hap,block_start,block_end,IND_name,type,cpu_cores)
haplotype_allele=tmp$SNP
data_numeric=tmp$numeric
rm(tmp);gc();
rownames(data_numeric)=IND_name
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_numeric),paste0(output_data_name,".numeric"),quote=F,row.names=T,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_numeric);gc();}
}

if(type==4){
tmp=haplotype_convertion(input_data_haplotype_hap,block_start,block_end,IND_name,type,cpu_cores)
haplotype_allele=tmp$SNP
data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=tmp$blupf90,stringsAsFactors=F)
rm(tmp);gc();
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=F,sep=" ")}
if(return_result==FALSE){rm(data_blupf90);gc();}
}

#构建bigmemory对象 
if(type%in%c(11,12,13,14)){

if(exists("pBigMat_hap")){rm(pBigMat_hap)}
gc()
delete_bigmemory_file("haplotype_hap",bigmemory_data_name,bigmemory_data_path,FALSE);
pBigMat_hap=make_bigmemory_object_cpp(nrow=n_snp,ncol=n_ind*2,file_name=paste0(bigmemory_data_name,"_haplotype_hap"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
pBigMat_hap[]=input_data_haplotype_hap
rm(input_data_haplotype_hap,pBigMat_hap);gc();
pBigMat_hap_address=make_bigmemory_object_address_cpp(nrow=n_snp,ncol=n_ind*2,file_name=paste0(bigmemory_data_name,"_haplotype_hap"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)																
}

if(type%in%c(11,12)){

if(exists("pBigMat_list_address")){rm(pBigMat_list_address)}
gc()
delete_bigmemory_file("haplotype_list",bigmemory_data_name,bigmemory_data_path,FALSE);
delete_bigmemory_file("numeric",bigmemory_data_name,bigmemory_data_path,FALSE);

pBigMat_list_address=make_bigmemory_object_address_cpp(nrow=n_window,ncol=n_ind*2,file_name=paste0(bigmemory_data_name,"_haplotype_list"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
}

if(type%in%c(13)){

if(exists("pBigMat_num_address")){rm(pBigMat_num_address)}
gc()
delete_bigmemory_file("numeric",bigmemory_data_name,bigmemory_data_path,FALSE);

pBigMat_num_address=make_bigmemory_object_address_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
}


if(type==11){
haplotype_allele=phased_haplotype_to_numeric_memory_cpp(pBigMat_hap_address,pBigMat_list_address,
									   numeric_file_name=paste0(bigmemory_data_name,"_numeric"),
									   numeric_file_path=bigmemory_data_path,
                                             block_start,block_end,cpu_cores)

data_numeric=make_bigmemory_object_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
options(bigmemory.allow.dimnames=TRUE)
rownames(data_numeric)=IND_name
colnames(data_numeric)=paste0("SNP",1:ncol(data_numeric))																		
cat(paste0("bigmemory-Numeric data is save as ",bigmemory_data_path,"/",bigmemory_data_name,"_numeric \n"))
}


if(type==12){

tmp=phased_haplotype_to_blupf90_memory_cpp(block_start,block_end,pBigMat_hap_address,pBigMat_list_address,IND_name,cpu_cores)
haplotype_allele=tmp$SNP

data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=tmp$allele_string_vector,stringsAsFactors=F)								 
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=T,sep=" ")}
}

if(type==13){
unphased_haplotype_to_numeric_memory_cpp(pBigMat_hap_address,pBigMat_num_address,cpu_cores,miss_base,miss_base_num)

data_numeric=make_bigmemory_object_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
options(bigmemory.allow.dimnames=TRUE)
rownames(data_numeric)=IND_name
colnames(data_numeric)=paste0("SNP",1:ncol(data_numeric))																		
cat(paste0("bigmemory-Numeric data is save as ",bigmemory_data_path,"/",bigmemory_data_name,"_numeric \n"))
}

if(type==14){
data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=unphased_haplotype_to_blupf90_memory_cpp(pBigMat_hap_address,IND_name,cpu_cores),stringsAsFactors=F)								 
if(!is.null(output_data_path)&!is.null(output_data_name)){
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=T,sep=" ")}
}

}
}



####################################input data type is Numeric#################################### 
if(input_data_type=="Numeric"){
if(!is.null(input_data_path)&!is.null(input_data_name)){
cat("Start read the Numeric format genotype data \n")
input_data_numeric=fread(paste0(input_data_path,"/",input_data_name),data.table=F,header=F)
IND_name=as.character(input_data_numeric[,1])
input_data_numeric=as.matrix(input_data_numeric[,-1])
if(file.exists(paste0(input_data_path,"/",input_data_name,".map"))){
cat("Start read the Numeric_map data \n")
input_data_numeric_map=fread(paste0(input_data_path,"/",input_data_name,".map"),data.table=F,header=F)
}
cat("Complete read the Numeric format genotype data \n")}

output_type_number=get_output_type_number(input_data_type=input_data_type,
										  output_data_type=output_data_type,
										  bigmemory_cal=bigmemory_cal,
										  phased_genotype=phased_genotype)

if(length(output_type_number)==0){stop("Please provided standard output data type!")}


if(!is.null(input_data_numeric_map)){
if(ncol(input_data_numeric_map)==3){
cat("Provided input_data_numeric_map constains 3 columns, please make sure the format is: SNP_name, Chromosome, Position \n")
input_data_numeric_map=input_data_numeric_map[,c(2,1,3)]
input_data_numeric_map=cbind(input_data_numeric_map,"A","T")
}else if(ncol(input_data_numeric_map)==5){
cat("Provided input_data_numeric_map constains 5 columns, please make sure the format is: Chromosome,SNP_name, Position, Ref, Alt\n")
}else{
stop("Provided input_data_numeric_map didn't meet the requirement format, it should be 3 columns or 5 columns!")
}

input_data_numeric_map=as.matrix(input_data_numeric_map)
input_data_numeric_map[,1]=as.numeric(input_data_numeric_map[,1])
input_data_numeric_map[,3]=as.numeric(input_data_numeric_map[,3])
}


for(type in output_type_number){

if(type==1){data_blupf90=data.frame(V1=sprintf(paste0("%-",max(nchar(IND_name))+3,"s"),IND_name),V2=numeric_to_blupf90_cpp(input_data_numeric,IND_name,cpu_cores),stringsAsFactors=F)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data_blupf90,paste0(output_data_name,".blupf90"),quote=F,row.names=F,col.names=F,sep=" ")}
if(return_result==FALSE){rm(data_blupf90);gc();}
}

if(type==2){
if(is.null(input_data_numeric_map)){stop("Please provide Numeric map information!")}
data_ped=numeric_to_ped_cpp(IND_name,input_data_numeric_map,input_data_numeric,cpu_cores,miss_base)
data_map=input_data_numeric_map[,c(1,2,3,3)];data_map[,3]=as.numeric(data_map[,3])/1000000;
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_ped),paste0(output_data_name,".ped"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(data_map),paste0(output_data_name,".map"),quote=F,row.names=F,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_map,data_ped);gc();}

}


if(type==3){
if(is.null(input_data_numeric_map)){stop("Please provide Numeric map information!")}
input_data_numeric_map=as.matrix(input_data_numeric_map)
data_hmp=numeric_to_hapmap_cpp(IND_name,input_data_numeric_map,input_data_numeric,cpu_cores,miss_base)
data_hmp[,c(2,5:11)]="NA"
colnames(data_hmp)=c("rs#","alleles","chrom","pos","strand","assembly","center","protLSID","assayLSID","panelLSID","QCcode",IND_name)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.table(data_hmp),paste0(output_data_name,".hmp.txt"),quote=F,row.names=F,col.names=T,sep="\t")}
if(return_result==FALSE){rm(data_hmp);gc();}
}


if(type==4){
if(is.null(input_data_numeric_map)){stop("Please provide Numeric map information!")}
data_vcf=numeric_to_vcf_cpp(IND_name,input_data_numeric_map,input_data_numeric,phased_symbol,cpu_cores)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(get_vcf_header()),paste0(output_data_name,".vcf"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(data_vcf),paste0(output_data_name,".vcf"),append=TRUE,quote=F,row.names=F,col.names=T,sep="\t")}
if(return_result==FALSE){rm(data_vcf);gc();}
}

}
rm(input_data_numeric);gc();
}


####################################input data type is BLUPF90#################################### 
if(input_data_type=="BLUPF90"){
if(!is.null(input_data_path)&!is.null(input_data_name)){
cat("Start read the BLUPF90 format genotype data \n")
input_data_blupf90=fread(paste0(input_data_path,"/",input_data_name),data.table=F,header=F)
if(file.exists(paste0(input_data_path,"/",input_data_name,".map"))){
cat("Start read the BLUPF90_map data \n")
input_data_blupf90_map=fread(paste0(input_data_path,"/",input_data_name,".map"),data.table=F,header=F)
}
cat("Complete read the BLUPF90 format genotype data \n")}
IND_name=as.character(input_data_blupf90[,1]);input_data_blupf90=input_data_blupf90[,2];gc();
n_ind=length(IND_name)
n_snp=nchar(input_data_blupf90[1])
output_type_number=get_output_type_number(input_data_type=input_data_type,
										  output_data_type=output_data_type,
										  bigmemory_cal=bigmemory_cal,
										  phased_genotype=phased_genotype)

if(length(output_type_number)==0){stop("Please provided standard output data type!")}

if(!is.null(input_data_blupf90_map)){
if(ncol(input_data_blupf90_map)==3){
cat("Provided input_data_blupf90_map constains 3 columns(identical to BLUPF90 map format), please make sure the format is: SNP_name, Chromosome, Position \n")
input_data_blupf90_map=input_data_blupf90_map[,c(2,1,3)]
input_data_blupf90_map=cbind(input_data_blupf90_map,"A","T")
}else if(ncol(input_data_blupf90_map)==5){
cat("Provided input_data_blupf90_map constains 5 columns, please make sure the format is: Chromosome,SNP_name, Position, Ref, Alt\n")
}else{
stop("Provided input_data_blupf90_map didn't meet the requirement format, it should be 3 columns or 5 columns!")
}

input_data_blupf90_map=as.matrix(input_data_blupf90_map)
input_data_blupf90_map[,1]=as.numeric(input_data_blupf90_map[,1])
input_data_blupf90_map[,3]=as.numeric(input_data_blupf90_map[,3])
}


for(type in output_type_number){
if(exists(".Random.seed")){fixed_openMP(.Random.seed);set.seed(19980204)}

if(type==1){data_numeric=blupf90_to_numeric_cpp(input_data_blupf90,cpu_cores)
rownames(data_numeric)=IND_name
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_numeric),paste0(output_data_name,".numeric"),quote=F,row.names=T,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_numeric);gc();}
}


if(type==2){
if(is.null(input_data_blupf90_map)){stop("Please provide BLUPF90 map information!")}
data_ped=blupf90_to_ped_cpp(IND_name,input_data_blupf90_map,input_data_blupf90,cpu_cores,miss_base)
data_map=input_data_blupf90_map[,c(1,2,3,3)];data_map[,3]=as.numeric(data_map[,3])/1000000;
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(data_ped),paste0(output_data_name,".ped"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.frame(data_map),paste0(output_data_name,".map"),quote=F,row.names=F,col.names=F,sep="\t")}
if(return_result==FALSE){rm(data_map,data_ped);gc();}

}


if(type==3){
if(is.null(input_data_blupf90_map)){stop("Please provide BLUPF90 map information!")}
data_hmp=blupf90_to_hapmap_cpp(IND_name,input_data_blupf90_map,input_data_blupf90,cpu_cores,miss_base)
data_hmp[,c(2,5:11)]="NA"
colnames(data_hmp)=c("rs#","alleles","chrom","pos","strand","assembly","center","protLSID","assayLSID","panelLSID","QCcode",IND_name)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.table(data_hmp),paste0(output_data_name,".hmp.txt"),quote=F,row.names=F,col.names=T,sep="\t")}
if(return_result==FALSE){rm(data_hmp);gc();}
}


if(type==4){
if(is.null(input_data_blupf90_map)){stop("Please provide BLUPF90 map information!")}
data_vcf=blupf90_to_vcf_cpp(IND_name,input_data_blupf90_map,input_data_blupf90,phased_symbol,cpu_cores)
if(!is.null(output_data_path)&!is.null(output_data_name)){
cat("Start writing output data...... \n")
fwrite(data.frame(get_vcf_header()),paste0(output_data_name,".vcf"),quote=F,row.names=F,col.names=F,sep="\t")
fwrite(data.table(data_vcf),paste0(output_data_name,".vcf"),append=TRUE,quote=F,row.names=F,col.names=T,sep="\t")}
if(return_result==FALSE){rm(data_vcf);gc();}
}


#构建bigmemory对象 
if(type==11){

if(exists("pBigMat_num_address")){rm(pBigMat_num_address)}
gc()
delete_bigmemory_file("numeric",bigmemory_data_name,bigmemory_data_path,FALSE);

pBigMat_num_address=make_bigmemory_object_address_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)
																
blupf90_to_numeric_memory_cpp(input_data_blupf90,pBigMat_num_address,cpu_cores,type=bigmemory_data_type)								 
data_numeric=make_bigmemory_object_cpp(nrow=n_ind,ncol=n_snp,file_name=paste0(bigmemory_data_name,"_numeric"),
                                                                file_path=bigmemory_data_path,type=bigmemory_data_type)	
options(bigmemory.allow.dimnames=TRUE)
rownames(data_numeric)=IND_name
colnames(data_numeric)=paste0("SNP",1:ncol(data_numeric))																
cat(paste0("bigmemory-Numeric data is save as ",bigmemory_data_path,"/",bigmemory_data_name,"_numeric \n"))
																
}


}
rm(input_data_blupf90);gc();
}

options (warn =1)
##############################################

if(!exists("data_hmp")){data_hmp=NULL}
if(!exists("data_ped")){data_ped=NULL}
if(!exists("data_map")){data_map=NULL}
if(!exists("data_blupf90")){data_blupf90=NULL}
if(!exists("data_numeric")){data_numeric=NULL}
if(!exists("data_vcf")){data_vcf=NULL}
if(!exists("data_haplotype_hap")){data_haplotype_hap=NULL}
if(!exists("data_haplotype_map")){data_haplotype_map=NULL}
if(!exists("data_haplotype_sample")){data_haplotype_sample=NULL}
if(!exists("output_block")){output_block=NULL}


if(!is.null(haplotype_allele)&!is.null(output_block)){

haplotype_allele_num=sapply(haplotype_allele,length)
block_pos_allele=output_block$block_pos
block_pos_allele$chr=input_data_haplotype_map[block_pos_allele[,2],1]
block_pos_allele$pos=paste0(input_data_haplotype_map[block_pos_allele[,2],3],"_",
						     input_data_haplotype_map[block_pos_allele[,3],3])

haplotype_allele_chr=rep(block_pos_allele$chr,times=haplotype_allele_num)
haplotype_allele_pos=rep(block_pos_allele$pos,times=haplotype_allele_num)

haplotype_allele=do.call(c,haplotype_allele)
haplotype_allele=data.frame(Chr=haplotype_allele_chr,
							SNP=paste0(1:length(haplotype_allele)),
							cM=0,
							Pos=haplotype_allele_pos,
							Allele=haplotype_allele)
output_block=c(output_block,list(allele=haplotype_allele))							
}

if(return_result==TRUE){return(list(hmp=data_hmp,
                                            ped=data_ped,
								     map=data_map,
									 blupf90=data_blupf90,
									 numeric=data_numeric,
									 vcf=data_vcf,
									 haplotype_hap=data_haplotype_hap,
									 haplotype_map=data_haplotype_map,
									 haplotype_sample=data_haplotype_sample,
									 phased_block=output_block))}
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


get_vcf_header<-function(date=Sys.time()){
H1="##fileformat=VCFv4.2"
H2=paste0("##filedate=",date)
H3="##source=R-package:blupADC"
H4="##INFO=<ID=AF,Number=A,Type=Float"
H5="##INFO=<ID=DR2,Number=A,Type=Float"
H6="##INFO=<ID=IMP,Number=0,Type=Flag"
H7="##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype"
H8="##FORMAT=<ID=DS,Number=A,Type=Float"
return(rbind(H1,H2,H3,H4,H5,H6,H7,H8))
}


get_output_type_number<-function(input_data_type=NULL,
								 output_data_type=NULL,
								 bigmemory_cal=FALSE,
								 phased_genotype=FALSE)
								 {
data_type_set=c("Hapmap","Plink","VCF","Haplotype","Numeric","BLUPF90")
output_data_type=unique(output_data_type)
input_data_type=unique(input_data_type)

if(length(unique(input_data_type))!=1){stop("Provided input_data_type  more than two or equal zero!")}
if(sum(!input_data_type%in%data_type_set)>=1){stop("Provided input_data_type contain non-standard type!")}
if(sum(!output_data_type%in%data_type_set)>=1){stop("Provided output_data_type contain non-standard type!")}

output_data_type=output_data_type[!output_data_type%in%input_data_type]

type_number=NULL

if(input_data_type=="Plink"){

   if(bigmemory_cal==FALSE){
		if("Hapmap"%in%output_data_type){type_number=c(type_number,1)}
		if("VCF"%in%output_data_type){type_number=c(type_number,2)}	
		if("Numeric"%in%output_data_type){type_number=c(type_number,3)}
		if("BLUPF90"%in%output_data_type){type_number=c(type_number,4)}   
   }else{
		
		if("Hapmap"%in%output_data_type){type_number=c(type_number,1)}
		if("VCF"%in%output_data_type){type_number=c(type_number,2)}	
		if("Numeric"%in%output_data_type){type_number=c(type_number,13)}
		if("BLUPF90"%in%output_data_type){type_number=c(type_number,4)}    		
   }

}else if(input_data_type=="Hapmap"){

   if(bigmemory_cal==FALSE){
		if("Plink"%in%output_data_type){type_number=c(type_number,1)}
		if("VCF"%in%output_data_type){type_number=c(type_number,2)}	
		if("Numeric"%in%output_data_type){type_number=c(type_number,3)}
		if("BLUPF90"%in%output_data_type){type_number=c(type_number,4)}   
   }else{
		if("Plink"%in%output_data_type){type_number=c(type_number,1)}
		if("VCF"%in%output_data_type){type_number=c(type_number,2)}	
		if("Numeric"%in%output_data_type){type_number=c(type_number,13)}
		if("BLUPF90"%in%output_data_type){type_number=c(type_number,4)}   		
   }

}else if(input_data_type=="VCF"){

   if(bigmemory_cal==FALSE){
		if("Hapmap"%in%output_data_type){type_number=c(type_number,1)}
		if("Haplotype"%in%output_data_type){type_number=c(type_number,3)}	
		if("Numeric"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,4)}
		if("BLUPF90"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,5)}
		if("Numeric"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,6)}
		if("BLUPF90"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,7)} 
		if("Plink"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,8)} 
		if("Plink"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,2)} 
   }else{		
		if("Hapmap"%in%output_data_type){type_number=c(type_number,1)}
		if("Plink"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,8)} 
		if("Plink"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,2)} 
		if("Haplotype"%in%output_data_type){type_number=c(type_number,13)}	
		if("Numeric"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,14)}
		if("BLUPF90"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,15)}
		if("Numeric"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,16)}
		if("BLUPF90"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,7)}  
   }

}else if(input_data_type=="Haplotype"){

   if(bigmemory_cal==FALSE){
		if("Numeric"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,1)}  
		if("BLUPF90"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,2)} 			
		if("Numeric"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,3)}		
		if("BLUPF90"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,4)} 
	
   }else{ 
            
		if("Numeric"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,11)}  
		if("BLUPF90"%in%output_data_type&phased_genotype==TRUE){type_number=c(type_number,12)} 			
		if("Numeric"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,13)}		
		if("BLUPF90"%in%output_data_type&phased_genotype==FALSE){type_number=c(type_number,14)} 
   
   }
   
}else if(input_data_type=="Numeric"){


		if("Plink"%in%output_data_type|"Hapmap"%in%output_data_type|"VCF"%in%output_data_type){
	
		cat("Please make sure provide input_data_numeric_map, which should be 3 columns(SNP_name, Chromosome, Position) or 5 columns(Chromosome,SNP_name, Position, Ref, Alt) \n")

		}

   if(bigmemory_cal==FALSE){
		if("BLUPF90"%in%output_data_type){type_number=c(type_number,1)}
		if("Plink"%in%output_data_type){type_number=c(type_number,2)}
		if("Hapmap"%in%output_data_type){type_number=c(type_number,3)}
		if("VCF"%in%output_data_type){type_number=c(type_number,4)}		
	
   }else{ 
         if("BLUPF90"%in%output_data_type){type_number=c(type_number,1)}
		if("Plink"%in%output_data_type){type_number=c(type_number,2)}
		if("Hapmap"%in%output_data_type){type_number=c(type_number,3)}
		if("VCF"%in%output_data_type){type_number=c(type_number,4)}		  
   
   }
   
}else if(input_data_type=="BLUPF90"){

		if("Plink"%in%output_data_type|"Hapmap"%in%output_data_type|"VCF"%in%output_data_type){
	
		cat("Please make sure  provide input_data_blupf90_map, which should be 3 columns(SNP_name, Chromosome, Position) or 5 columns(Chromosome,SNP_name, Position, Ref, Alt) \n")

		}

   if(bigmemory_cal==FALSE){
		if("Numeric"%in%output_data_type){type_number=c(type_number,1)}
		if("Plink"%in%output_data_type){type_number=c(type_number,2)}
		if("Hapmap"%in%output_data_type){type_number=c(type_number,3)}
		if("VCF"%in%output_data_type){type_number=c(type_number,4)}
	
   }else{ 
         if("Numeric"%in%output_data_type){type_number=c(type_number,11)} 
		if("Plink"%in%output_data_type){type_number=c(type_number,2)}
		if("Hapmap"%in%output_data_type){type_number=c(type_number,3)}
		if("VCF"%in%output_data_type){type_number=c(type_number,4)}		 
   
   }
}
return(type_number)
}








#获取单倍型 block
get_haplotype_block<-function(data_map=NULL,
							  haplotype_window_nSNP=NULL,
							  haplotype_window_kb=NULL,
							  haplotype_window_block=NULL){

#同时提供三种参数中至少两种会报错
if(sum(c(!is.null(haplotype_window_nSNP),!is.null(haplotype_window_kb),!is.null(haplotype_window_block)))>=2){stop("Provided too much parateters: haplotype_window_nSNP, haplotype_window_kb, haplotype_window_block!")}
if(sum(c(!is.null(haplotype_window_nSNP),!is.null(haplotype_window_kb),!is.null(haplotype_window_block)))==0){stop("Provided no parateters: haplotype_window_nSNP, haplotype_window_kb, haplotype_window_block!")}

#获取SNP位置-考虑过了染色体
data_map=data.frame(data_map,stringsAsFactors=F)
data_map[,1]=as.numeric(data_map[,1])
data_map[,3]=as.numeric(data_map[,3])
chr_set=unique(data_map[,1])
chr_set_snp_n=0
for(i in 1:length(chr_set)){
tmp_data_map=data_map[data_map[,1]==chr_set[i],]
chr_set_snp_n=c(chr_set_snp_n,nrow(tmp_data_map))
}
chr_set_snp_n=cumsum(chr_set_snp_n)
chr_set_snp_n=chr_set_snp_n[-length(chr_set_snp_n)]

if(!is.null(haplotype_window_block)){

block_start=haplotype_window_block[,1]-1  #window_block为两列数据，均为位置信息
block_end=haplotype_window_block[,2]-1

}else if(!is.null(haplotype_window_nSNP)){ #根据SNP数目划分block

block_start=NULL
block_end=NULL
for(i in 1:length(chr_set)){
tmp_data_map=data_map[data_map[,1]==chr_set[i],]
block_start_tmp=seq(1,nrow(tmp_data_map),haplotype_window_nSNP)-1
block_end_tmp=seq(1,nrow(tmp_data_map),haplotype_window_nSNP)-2
block_end_tmp=c(block_end_tmp,nrow(tmp_data_map)-1)
block_end_tmp=block_end_tmp[-1]

#添加位置信息
block_start_tmp=block_start_tmp+chr_set_snp_n[i]
block_end_tmp=block_end_tmp+chr_set_snp_n[i]
block_start=c(block_start,block_start_tmp)
block_end=c(block_end,block_end_tmp)
}

}else if(!is.null(haplotype_window_kb)){  #根据物理距离划分block

block_start=NULL
block_end=NULL
for(i in 1:length(chr_set)){
tmp_data_map=data_map[data_map[,1]==chr_set[i],]
block=seq(min(tmp_data_map[,3]),max(tmp_data_map[,3]),haplotype_window_kb*1000)

block_1=block
block_2=block-1;
block_2=block_2[-1]; #去除最后一列
block_2=c(block_2,max(tmp_data_map[,3]))

#R-function too slow
#block_start_tmp=NULL
#block_end_tmp=NULL
#for(j in 1:length(block)){
#block_start_tmp=c(block_start_tmp,min(which(tmp_data_map[,3]>=block_1[j])))
#block_end_tmp=c(block_end_tmp,max(which(tmp_data_map[,3]<=block_2[j])))
#}

#Rcpp function 
block_result=define_block_window_kb_cpp(block_1,block_2,tmp_data_map[,3])
block_start_tmp=block_result[[1]]
block_end_tmp=block_result[[2]]

pos_status=block_start_tmp<=block_end_tmp
block_start_tmp=block_start_tmp[pos_status]
block_end_tmp=block_end_tmp[pos_status]
block_start_tmp=block_start_tmp-1
block_end_tmp=block_end_tmp-1
#添加位置信息
block_start_tmp=block_start_tmp+chr_set_snp_n[i]
block_end_tmp=block_end_tmp+chr_set_snp_n[i]
block_start=c(block_start,block_start_tmp)
block_end=c(block_end,block_end_tmp)
}
}
#进行单倍型构建分析
#写出block的信息
snp_map=as.character(data_map[,2])
block_pos=data.frame(block_num=paste0("block",1:length(block_start)),window_start_pos=block_start+1,window_end_pos=block_end+1)
block_snp=data.frame(block_num=paste0("block",1:length(block_start)),window_start_snp=snp_map[block_start+1],window_end_snp=snp_map[block_end+1])
return(list(block=list(block_pos=block_pos,block_snp=block_snp),
			block_rcpp=cbind(block_start,block_end) #c++程序的输入文件
			))
}

fixed_openMP<-function(x){
  Sx <- deparse(substitute(x))
  rm(list=Sx,envir=sys.frame(-1))
}

get_bigmemory_address<-function(data){

	return(data@address)
}

delete_bigmemory_file<-function(matrix_type,
								    bigmemory_data_name,
						             bigmemory_data_path,
							         message=FALSE){
								   
		
		file1=paste0(bigmemory_data_path,"/",bigmemory_data_name,"_",matrix_type,".desc")
		file2=paste0(bigmemory_data_path,"/",bigmemory_data_name,"_",matrix_type,".bk")

	    if(file.exists(file1)&file.exists(file2)){
			if(message==TRUE){
			cat(paste0("Found ",matrix_type,".desc & ",matrix_type,".bk files in path:", bigmemory_data_path,"\n"))
			cat("Software will delete these files automatically!\n")
			}
			unlink(file2)
			unlink(file1)
			
		}
	    if(file.exists(file1)){unlink(file1)}
		 if(file.exists(file2)){unlink(file2)}

		
}
