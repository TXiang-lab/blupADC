library(blupADC)
library(profmem)
library(nadiv)
library(microbenchmark)
library(bigmemory)
options(bigmemory.allow.dimnames=TRUE)

ped=fread("ped.txt",data.table=F)
for(i in c(1)){

mem_blupADC_bigmemory =profmem({
file.remove("blupADC_P_D.bk")
file.remove("blupADC_P_D.desc")
B=cal_kinship(input_pedigree=ped,          #provided hapmap data object
				bigmemory_cal=TRUE,
				gene_dropping=TRUE,
				gene_dropping_iteration=1,
				cpu_cores=i,
                kinship_type=c("P_D"),      #type of  kinship matrix
                return_result=TRUE)
file.remove("blupADC_P_D.bk")
file.remove("blupADC_P_D.desc")				
rm(B);gc();					
})

mem_blupADC_normal =profmem({
A=cal_kinship(input_pedigree=ped,          #provided hapmap data object
				gene_dropping=TRUE,
				gene_dropping_iteration=1,
				cpu_cores=i,
                kinship_type=c("P_D"),      #type of  kinship matrix
                return_result=TRUE)
rm(A);gc();				
})
mem_nadiv=profmem({
D=makeDsim(
ped,
N=1,
parallel = TRUE,
ncores = i,
invertD = FALSE,
calcSE = FALSE,
returnA = FALSE
)
rm(D);gc();	
}) 

saveRDS(mem_nadiv,paste0("cpu_",i,"_mem_nadiv.rds"))
saveRDS(mem_blupADC_normal,paste0("cpu_",i,"_mem_blupADC_normal.rds"))
saveRDS(mem_blupADC_bigmemory,paste0("cpu_",i,"_mem_blupADC_bigmemory.rds"))


mem_nadiv=data.frame( mem_nadiv="mem_nadiv",
					cpu_cores=i,
								 total_mem_MB=sum(na.omit(mem_nadiv$bytes))/(1024^2),
								 max_mem_MB=max(na.omit(mem_nadiv$bytes))/(1024^2))

mem_blupADC_bigmemory=data.frame( mem_blupADC_bigmemory="mem_blupADC_bigmemory",cpu_cores=i,
								 total_mem_MB=sum(na.omit(mem_blupADC_bigmemory$bytes))/(1024^2),
								 max_mem_MB=max(na.omit(mem_blupADC_bigmemory$bytes))/(1024^2))

mem_blupADC_normal=data.frame( mem_blupADC_normal="mem_blupADC_normal",cpu_cores=i,
								 total_mem_MB=sum(na.omit(mem_blupADC_normal$bytes))/(1024^2),
								 max_mem_MB=max(na.omit(mem_blupADC_normal$bytes))/(1024^2))

print(paste0("Cpu cores = ",i))
print(mem_nadiv)
print(mem_blupADC_normal)
print(mem_blupADC_bigmemory)
}