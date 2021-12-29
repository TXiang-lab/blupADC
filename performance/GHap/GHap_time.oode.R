library(GHap)
library(blupADC)
library(profmem)
library(microbenchmark)
target_path=getwd()

for(i in c(1:20)){

GHap_time_stats <- microbenchmark(

blupADC_bigmemory = {
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
},
blupADC_normal = {
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
}, 

GHap={
ghap.compress(input.file = "human", out.file = "human")
phase <- ghap.loadphase("human")
blocks.mkr <- ghap.blockgen(phase, windowsize = 5, slide = 5, unit = "marker")
ghap.haplotyping(phase, blocks.mkr, outfile = "human", binary = FALSE, ncores = i)
rm(phase,blocks.mkr);gc();
file.remove("human.hapalleles")
file.remove("human.hapgenotypes")
file.remove("human.hapsamples")
file.remove("human.phaseb")
} 
,
times =1, 
unit = "s"
)
print(paste0("Cpu-cores= ",i))
print(GHap_time_stats)


saveRDS(GHap_time_stats,paste0("cpu_",i,"_GHap_time_stats.rds"))


}
