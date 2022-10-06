geno_qc_impute<-function(
      input_data_hmp=NULL,       # for a given hapmap not read from  R
      input_data_name=NULL,      # get data from the work directory which don't need to read by R, but data must put in the input_data_path
      input_data_type=NULL,  # there are three type of data, include  hapmap 、plink 、 vcf 、Blupf90 the default is hapmap
	 input_data_plink_map=NULL,  # for a given map not read from  R
	 input_data_plink_ped=NULL,  # for a given ped not read from  R
	 input_data_blupf90=NULL,
	 input_data_path=NULL,     	# the directory of genotype data 
	 input_data_vcf=NULL,
	 input_data_numeric=NULL,
	 input_data_haplotype_hap=NULL,
	 input_data_haplotype_map=NULL,
	 input_data_haplotype_sample=NULL,
	 input_data_numeric_map=NULL,
	 input_data_blupf90_map=NULL,	 
	 bigmemory_cal=FALSE,
	 bigmemory_data_type="integer",
	 bigmemory_data_path=getwd(),
	 bigmemory_data_name="blupADC",
	 phased_symbol="|",
	 phased_genotype=FALSE,
	 haplotype_window_nSNP=NULL,
	 haplotype_window_kb=NULL,
	 haplotype_window_block=NULL,	 
	 output_data_path=NULL,
	 output_data_name=NULL,
	 output_data_type="Plink",
	 data_analysis_method="QC_Imputation",   # there are three data_analysis_method , Plink 、 Beagle 、 QC_Imputation	 
	 cpu_cores=1, 
      miss_base="0",          # the default missing value of genotype 
	 qc_snp_rate=0.1,
	 qc_ind_rate=0.1,
	 qc_maf=0.05,
	 qc_hwe=0.0000001,
	 extra_parameter=NULL,
	 plink_software_path=ifelse(as.character(Sys.info()["sysname"])=="Windows",system.file("extdata/bin_windows", package = "blupADC"),system.file("extdata/bin_linux", package = "blupADC")),
	 plink_software_name="plink", 
	 chr_set=NULL, # 识别多个染色体-plink, 最多可识别多个95个 	 
	 keep_inds_set=NULL,   #字符串向量
	 keep_snps_set=NULL,   #字符串向量
	 keep_chroms_set=NULL, #数值向量
	 beagle_software_path=ifelse(as.character(Sys.info()["sysname"])=="Windows",system.file("extdata/bin_windows", package = "blupADC"),system.file("extdata/bin_linux", package = "blupADC")),
	 beagle_software_name="beagle.5.2.jar",
	 beagle_ped_path=NULL,  
	 beagle_ped_name=NULL,
	 beagle_ref_data_path=NULL, #Using Reference For Inputaion
	 beagle_ref_data_name=NULL,
	 Java_Space=NULL,	 
	 sigle_duo_trio=NULL  # becuase if you want to use pedigree, trio and duos would make your computation more slowly , the standard_format is  sigle_duo_trio=c(" singlescale=20  duoscale=20  trioscale=20 ")	
){  # determin whether output the hmp.file 
	#author:  Quanshun Mei 
	#date: 2019.09.09
	#purpose: conduct QC and Imputation 
	#data: the standard genotype file, include  .vcf 、.ped .map 、 .hmp ，
	#missing: the symbols for the missing genotype
	#SNPcall, INDcall, maf, hwe: the conditions that will be used during using Plink
	#path默认的格式为  ：  eg： beaglepath_4.1="/home/qsmei/beagle_4.1"
	# the version of beagle include beagle.4.0.jar 、beagle.4.1.jar 、beagle.5.1.jar，   only beagle.4.0.jar can add pedigree into imputation 

if(is.null(output_data_name)){stop("Please specify your output data name!!!")}
if(is.null(output_data_path)){stop("Please specify your output data path!!!")}
if(!is.null(chr_set)){chr_set=paste0(" --chr-set ",chr_set," ")}

library(data.table)

input=get_input_data_type(input_data_type=input_data_type,input_data_hmp=input_data_hmp,input_data_plink_ped=input_data_plink_ped,
                                          input_data_plink_map=input_data_plink_map,input_data_blupf90=input_data_blupf90,
								    input_data_numeric=input_data_numeric,input_data_vcf=input_data_vcf,
									input_data_haplotype_hap=input_data_haplotype_hap,input_data_haplotype_map=input_data_haplotype_map,
									input_data_haplotype_sample=input_data_haplotype_sample,
								     miss_base=miss_base)

input_data_type=input$input_data_type
miss_base=input$miss_base

#当提供输入数据为vcf格式，且不需要质控的时候，直接开始填充
if(!(input_data_type=="VCF"&data_analysis_method=="Imputation")){


#如果输入的基因型文件为本地Plink格式的基因型文件，直接拷贝基因型文件即可
if(input_data_type=="Plink"&!is.null(input_data_path)){
cat("Start Read The Plink Format Genotype Data \n")
system(paste0("cp -r ",input_data_path,"/",input_data_name,".ped   ",output_data_path,"/temp_plink.ped"))
system(paste0("cp -r ",input_data_path,"/",input_data_name,".map   ",output_data_path,"/temp_plink.map"))

input_data_ped=fread(paste0(output_data_path,"/temp_plink.ped"),select=1,data.table=F)
input_data_map=fread(paste0(output_data_path,"/temp_plink.map"),data.table=F)
cat("Complete Read The Plink Format Genotype Data \n")

}else if(input_data_type=="Plink"&!is.null(input_data_plink_ped)){
input_data_ped=input_data_plink_ped
input_data_map=input_data_plink_map
fwrite(data.table(input_data_plink_ped),paste0(output_data_path,"/temp_plink.ped"),quote=F,row.names=F,col.names=F,sep=" ")
fwrite(data.table(input_data_plink_map),paste0(output_data_path,"/temp_plink.map"),quote=F,row.names=F,col.names=F,sep=" ")
rm(input_data_plink_ped,input_data_plink_map);
gc();

}else{

#将输入数据转换为 Plink格式
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
	input_data_numeric_map=input_data_numeric_map,
	input_data_blupf90_map=input_data_blupf90_map,		
		miss_base=miss_base,
		cpu_cores=cpu_cores,
		output_data_type="Plink",
		output_data_name="temp_plink",
		output_data_path=output_data_path,	
		return_result=TRUE)	

input_data_ped=sum_data$ped
input_data_map=sum_data$map 
rm(sum_data);gc();

}

IND_n=nrow(input_data_ped)
SNP_n=nrow(input_data_map)
rm(input_data_ped,input_data_map);gc();


setwd(output_data_path)
output_data_path=getwd()

os_type=as.character(Sys.info()["sysname"])

if(os_type=="Linux"){

#针对软件路径含有空格的情况，如果路径存在空格，将软件copy 到输出路径下
if(length(grep(" ",plink_software_path))>0){
#cat("Attention: Phe path include space!")
file.copy(paste0(plink_software_path,"/",plink_software_name),getwd())
plink_software_path=getwd()
}

}else{

if(length(grep(" ",plink_software_path))>0){
#cat("Attention: Phe path include space!")
file.copy(paste0(plink_software_path,"/",plink_software_name,".exe"),getwd())
plink_software_path=getwd()
}

}


if(length(grep(" ",beagle_software_path))>0){
#cat("Attention: Phe path include space!")
file.copy(paste0(beagle_software_path,"/",beagle_software_name),getwd())
beagle_software_path=getwd()
}




#Part1:调用Plink,对Plink格式的数据进行质控
if(data_analysis_method%in%c("QC","QC_Imputation")){

output_data_name=paste0(output_data_name,"_QC")
if(!is.null(plink_software_path)){plink_software_path=paste0(plink_software_path,"/")}

temp=system(paste0(plink_software_path,plink_software_name),ignore.stdout=T)
if(temp==127){
stop("Software couldn't find plink in your computer!")
}

system(paste0(plink_software_path,plink_software_name," --noweb --file ",output_data_path,"/","temp_plink",
	        " --missing-genotype ",miss_base," --geno ",qc_snp_rate," --mind ",
			qc_ind_rate," --maf ",qc_maf," --hwe ",qc_hwe," --make-just-fam  ",chr_set," ",extra_parameter," ",
			" --allow-extra-chr --write-snplist --out ",	
			output_data_path,"/",output_data_name," --recode vcf-iid"),ignore.stdout=T)

if(file.exists(paste0(output_data_name,".fam"))){
IND_QC=fread(paste0(output_data_name,".fam"),header=F,data.table=F)
IND_QC_n=nrow(IND_QC)
}else{stop("All individuals have been filtered, please check your data carefully!")}

if(file.exists(paste0(output_data_name,".snplist"))){
SNP_QC=fread(paste0(output_data_name,".snplist"),header=F,data.table=F)
SNP_QC_n=nrow(SNP_QC)
}else{stop("All snps have been filtered, please check your data carefully!")}

cat("Quality control is done \n")
cat(paste0("Before QC ,there are ", SNP_n," SNPs; After QC ,there are ",SNP_QC_n," SNPs remained \n"))
cat(paste0("Before QC ,there are ", IND_n," indviduals; After QC ,there are ",IND_QC_n," indviduals remained \n"))



}else{
output_data_name=output_data_name
if(!is.null(plink_software_path)){plink_software_path=paste0(plink_software_path,"/")}
system(paste0(plink_software_path,plink_software_name," --noweb --file ",output_data_path,"/","temp_plink",
	        " --missing-genotype ",miss_base,chr_set," ",extra_parameter," ",
			" --allow-extra-chr  --out ",
			output_data_path,"/",output_data_name," --recode vcf-iid"),ignore.stdout=T)
}

#删除temp_plink data
system(paste0("rm -rf temp_plink.map")) 
system(paste0("rm -rf temp_plink.ped")) 


#删除copy过来的软件
if(file.exists(paste0(plink_software_name))){system(paste0("rm -rf ",plink_software_name)) }
if(file.exists(paste0(plink_software_name,".exe"))){system(paste0("rm -rf ",paste0(plink_software_name,".exe")))}



if(file.exists(paste0(output_data_name,".nosex"))){system(paste0("rm -rf ",output_data_name,".nosex")) }
#if(file.exists(paste0(output_data_name,".log"))){system(paste0("rm -rf ",output_data_name,".log")) }
file.rename(paste0(output_data_name,".log"),paste0(output_data_name,"_plink.log"))
if(file.exists(paste0(output_data_name,".snplist"))){system(paste0("rm -rf ",output_data_name,".snplist")) }

}else{

if(input_data_type=="VCF"&data_analysis_method=="Imputation"){
file.copy(paste0(input_data_path,"/",input_data_name),paste0(output_data_path,"/",output_data_name,".vcf"))
#system(paste0("cp -r ",input_data_path,"/",input_data_name,"   ",output_data_path,"/",output_data_name,".vcf"))
}

setwd(output_data_path)
output_data_path=getwd()

if(length(grep(" ",beagle_software_path))>0){
#cat("Attention: Phe path include space!")
file.copy(paste0(beagle_software_path,"/",beagle_software_name),getwd())
beagle_software_path=getwd()
}

}

if(file.exists(paste0(output_data_name,".fam"))){unlink(paste0(output_data_name,".fam"))}

#Part2:调用 Beagle 对数据进行填充
if(data_analysis_method%in%c("QC_Imputation","Imputation")){ 
beagle_software=paste0(beagle_software_path,"/",beagle_software_name)  
input_imputation_data=paste0(output_data_path,"/",output_data_name,".vcf")
output_imputation_data=paste0(output_data_path,"/",output_data_name,"_Imputation")
output_data_name=paste0(output_data_name,"_Imputation")



if(!is.null(beagle_ped_path)){
cat("Please make sure the version of Beagle in this analysi is 4.0 for which included pedigree_ module \n")
beagle_pedigree=paste0(" ped=",beagle_ped_path,"/",beagle_ped_name)}else{beagle_pedigree=NULL} 
				 
if(!is.null(beagle_ref_data_path)){
beagle_reference=paste0(" ref=",beagle_ref_data_path,"/",beagle_ref_data_name)}else{beagle_reference=NULL}	
			 
cat("Using Beagle to perform  imputation ......\n")

#调用 beagle

temp1=system("java",ignore.stdout = T,ignore.stderr = T)
if(temp1==127){
stop("Java doesn't install on you computer, please intall Java before!")
}

temp=system(paste0("java ",Java_Space," -jar ", beagle_software),ignore.stdout=T)
if(temp==1){
stop("Software couldn't find beagle in your computer!")
}

system(paste0("java ",Java_Space," -jar ", beagle_software,beagle_reference,beagle_pedigree,sigle_duo_trio,
			" gt=",input_imputation_data," nthreads=",cpu_cores," out=",output_imputation_data))

file.rename(paste0(output_imputation_data,".log"),paste0(output_imputation_data,"_beagle_output.log"))
#删除填充的input data
system(paste0("rm -rf ",input_imputation_data))
#解压 .gz 格式			
system(paste0("gzip  -d ",output_imputation_data,".vcf.gz")) #将.gz文件解压为 .vcf  ,便于进行 vcftools转换
}
 

#Part3:将质控或质控或填充好的数据进行输出, 默认输出的数据为： output_imputation_data .vcf  
 

 sum_data=geno_format(      
          input_data_name=paste0(output_data_name,".vcf"),  
	     input_data_path=output_data_path,		  
          input_data_type="VCF",  	
		miss_base=miss_base,
		cpu_cores=cpu_cores,
		bigmemory_cal=bigmemory_cal,
		bigmemory_data_type=bigmemory_data_type,
		bigmemory_data_path=bigmemory_data_path,
		bigmemory_data_name=bigmemory_data_name,
		phased_symbol=phased_symbol,
		phased_genotype=phased_genotype,
		haplotype_window_nSNP=haplotype_window_nSNP,
		haplotype_window_kb=haplotype_window_kb,
		haplotype_window_block=haplotype_window_block,		
		output_data_type=output_data_type,
		output_data_name=output_data_name,
		output_data_path=output_data_path,			
		return_result=FALSE)
		
if(file.exists(paste0(beagle_software_name))){system(paste0("rm -rf ",beagle_software_name))}
if(file.exists(paste0(output_data_name,".nosex"))){system(paste0("rm -rf ",output_data_name,".nosex")) }
if(file.exists(paste0(output_data_name,".fam"))){system(paste0("rm -rf ",output_data_name,".fam")) }
if(!"VCF"%in%output_data_type){system(paste0("rm -rf ",output_data_path,"/",output_data_name,".vcf"))}
}