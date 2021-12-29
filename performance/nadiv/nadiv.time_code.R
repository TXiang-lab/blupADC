library(blupADC)
library(profmem)
library(nadiv)
library(microbenchmark)
library(bigmemory)
options(bigmemory.allow.dimnames=TRUE)
ped=fread("ped.txt",data.table=F)

for(i in c(1:20)){

makeDsim_time_stats <- microbenchmark(

blupADC_bigmemory = {
file.remove("blupADC_P_D.bk")
file.remove("blupADC_P_D.desc")
B=cal_kinship(input_pedigree=ped,          #provided hapmap data object
				bigmemory_cal=TRUE,
				gene_dropping=TRUE,
				gene_dropping_iteration=1000,
				cpu_cores=i,
                kinship_type=c("P_D"),      #type of  kinship matrix
                return_result=TRUE)
file.remove("blupADC_P_D.bk")
file.remove("blupADC_P_D.desc")				
rm(B);gc();					
},
blupADC_normal = {
A=cal_kinship(input_pedigree=ped,          #provided hapmap data object
				gene_dropping=TRUE,
				gene_dropping_iteration=1000,
				cpu_cores=i,
                kinship_type=c("P_D"),      #type of  kinship matrix
                return_result=TRUE)
rm(A);gc();				
}, 

nadiv={
D=makeDsim(
ped,
N=1000,
parallel = TRUE,
ncores = i,
invertD = FALSE,
calcSE = FALSE,
returnA = FALSE
)
rm(D);gc();	
} 
,
times =1, 
unit = "s"
)
print(paste0("Cpu-cores= ",i))
print(makeDsim_time_stats)
saveRDS(makeDsim_time_stats,paste0("cpu_",i,"_makeDsim_time_stats.rds"))
}