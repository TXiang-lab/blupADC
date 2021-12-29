library(GHap)
library(blupADC)
library(profmem)

target_path=getwd()


for(i in c(1)){
rm(phased_result);gc();


mem_blupADC_bigmemory <- profmem({
#test-blupADC 
library(blupADC)
phased_result=geno_format(
         input_data_path="/root/JY_data/software/temp/GHap",      # input data path 
         input_data_name="blupADC_human",  # input data name,for vcf data 
         input_data_type="Haplotype",          # input data type
         phased_genotype=TRUE,           # whether the vcf data has been phased
		 bigmemory_cal=TRUE,             # format conversion via bigmemory object
         bigmemory_data_path=target_path,    # path of bigmemory data 
         bigmemory_data_name="test2", #name of bigmemory data 
         haplotype_window_nSNP=5,        # according to nSNP define block,
         output_data_type=c("Numeric"),# output data format
         return_result=TRUE,             #save result as a R environment variable
         cpu_cores=i                     # number of cpu 
                  )
rm(phased_result);gc();				  
})


mem_blupADC_normal <- profmem({
#test-blupADC 
library(blupADC)
phased_result=geno_format(
         input_data_path="/root/JY_data/software/temp/GHap",      # input data path 
         input_data_name="blupADC_human",  # input data name,for vcf data 
         input_data_type="Haplotype",          # input data type
         phased_genotype=TRUE,           # whether the vcf data has been phased
         haplotype_window_nSNP=5,        # according to nSNP define block,
         output_data_type=c("Numeric"),# output data format
         return_result=FALSE,             #save result as a R environment variable
		 output_data_path=target_path,
		 output_data_name="test1",
         cpu_cores=i                     # number of cpu 
                  )
rm(phased_result);gc();				  
				  
})

mem_hap <- profmem({
ghap.compress(input.file = "human", out.file = "human")
phase <- ghap.loadphase("human")
blocks.mkr <- ghap.blockgen(phase, windowsize = 5, slide = 5, unit = "marker")
ghap.haplotyping(phase, blocks.mkr, outfile = "human", binary = FALSE, ncores = i)
rm(phase,blocks.mkr);gc();
file.remove("human.hapalleles")
file.remove("human.hapgenotypes")
file.remove("human.hapsamples")
file.remove("human.phaseb")
})

setwd(target_path)
saveRDS(mem_hap,paste0("cpu",i,"_mem_hap.rds"))
saveRDS(mem_blupADC_normal,paste0("cpu",i,"_mem_blupADC_normal.rds"))
saveRDS(mem_blupADC_bigmemory,paste0("cpu",i,"_mem_blupADC_bigmemory.rds"))


mem_hap=data.frame( mem_hap="mem_hap",
					cpu_cores=i,
								 total_mem_MB=sum(na.omit(mem_hap$bytes))/(1024^2),
								 max_mem_MB=max(na.omit(mem_hap$bytes))/(1024^2))

mem_blupADC_bigmemory=data.frame( mem_blupADC_bigmemory="mem_blupADC_bigmemory",cpu_cores=i,
								 total_mem_MB=sum(na.omit(mem_blupADC_bigmemory$bytes))/(1024^2),
								 max_mem_MB=max(na.omit(mem_blupADC_bigmemory$bytes))/(1024^2))

mem_blupADC_normal=data.frame( mem_blupADC_normal="mem_blupADC_normal",cpu_cores=i,
								 total_mem_MB=sum(na.omit(mem_blupADC_normal$bytes))/(1024^2),
								 max_mem_MB=max(na.omit(mem_blupADC_normal$bytes))/(1024^2))


print(paste0("Cpu cores = ",i))
print(mem_hap)
print(mem_blupADC_normal)
print(mem_blupADC_bigmemory)

}			