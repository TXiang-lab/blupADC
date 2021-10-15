# blupADC <img src="https://img.shields.io/badge/Issues-%2B-brightgreen.svg" /><img src="https://img.shields.io/badge/license-GPL3.0-blue.svg" /> <img src="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210617165506.png" alt="logo-blupADC"  height="250" align="right"/>   
### R package for animal and plant breeding
## Contents

-   [OVERVIEW](#overview)

-   [GETTING STARTED](#getting-started)

    -   [Installation](#installation)
    -   [Features](#features)

-   [USAGE](#usage)

------------------------------------------------------------------------
**Documents support two-language([English](https://qsmei.netlify.app/post/feature-0-overview/overview/) and [Chinese](https://qsmei.netlify.app/zh/post/feature-0-overview/overview/)).** 

### OVERVIEW

`blupADC` is an useful and powerful tool for handling genomic data and pedigree data in animal and plant breeding(**traditional blup and genomic selection**).  In the design of this package, most of data analysis problems in breeding have been considered, and  the speed of calculation is also the key point. In terms of the speed,  the core functions of this package are coded by c++ (`Rcpp` and `RcppArmadillo `) , and it also supports  parallel calculation (by applying `openMP` programming) and big data calculation(by importing `bigmemory ` package). 

`blupADC` provides many useful functions for the whole steps for animal and plant breeding, including pedigree analysis(**trace pedigree, rename pedigree, and correct pedigree errors**), genotype data format conversion(supports **Hapmap, Plink, Blupf90, Numeric,  VCF and Haplotype** format), genotype data quality control and imputation, construction of kinship matrix(**pedigree, genomic  and single-step**),and genetic evaluation( by interfacing with two famous breeding softwares, **DMU** and **BLUPF90**  in an easy way). 

Finally, we kindly provides an easier way of applying `blupADC`, which is a free  website([shinyapp](http://47.95.251.15:443/blupADC/)).  Several functions are still under development.  But the pitfall of this website is that it can't handle big data. 

😊 Good Luck Charlie ! 
If you have suggestion or question, please contact: hzau_qsmei@163.com ! 
### 👨‍💻 Citation 

Quanshun Mei, Chuanke Fu, Jieling Li, Shuhong Zhao, and Tao Xiang. "blupADC: An R package and shiny toolkit for comprehensive genetic data analysis in animal and plant breeding." *bioRxiv* (2021), **doi:** https://doi.org/10.1101/2021.09.09.459557
## New features 

### 1.0.3

- Incorporate  maternal effect, permanent effect, random regression effect, and social genetic  effect models in  the genetic evaluation by DMU (2021.8.24)

### 1.0.4

- Incorporate  haplotype format conversion ,haplotype-based numeric matrix construction and haplotype-based additive relationship matrix construction (2021.10.8)
- Import  bigmemory object in matrix save and calculation  for handling big data(2021.10.8)

## GETTING STARTED

### 🙊Installation

`blupADC` links to R packages `Rcpp`, `RcppArmadillo` , `data.table` and  `bigmemory` .  These dependencies should be installed before installing `blupADC`.  

```R
install.packages(c("Rcpp", "RcppArmadillo","data.table","bigmemory"))
```

**👉 Note: In the analysis of DMU  and BLUPF90 , we need to download software DMU  ([DMU download website](https://dmu.ghpc.au.dk/dmu/))  and BLUPF90 previously ([BLUPF90 download website](http://nce.ads.uga.edu/html/projects/programs/)). For convenience, we have encapsulated  the basic module of DMU and BLUPF90 in package `blupADC`.**  

 **For commercial use of DMU and BLUPF90,  user must contact the author of DMU and BLUPF90 !!!** 

#### Install blupADC on Linux 

```R
packageurl <- "https://github.com/TXiang-lab/blupADC/releases/download/V1.0.4/blupADC_1.0.4_R_x86_64-pc-linux-gnu.tar.gz"
install.packages(packageurl,repos=NULL,method="libcurl")
```

#### Install blupADC on Windows

```R
packageurl <- "https://github.com/TXiang-lab/blupADC/releases/download/V1.0.4/blupADC_1.0.4.zip"
install.packages(packageurl,repos=NULL)
```

👉 **Note:If the connection with github is not good(such as in China), user can download as below:**  

#### Install blupADC on Linux 

```R
packageurl <- "https://gitee.com/qsmei/blupADC/attach_files/851170/download/blupADC_1.0.4_R_x86_64-pc-linux-gnu.tar.gz"
install.packages(packageurl,repos=NULL,method="libcurl")
```

#### Install blupADC on Windows

```R
packageurl<-"https://gitee.com/qsmei/blupADC/attach_files/851169/download/blupADC_1.0.4.zip"
install.packages(packageurl,repos=NULL)
```


After installed successfully, the `blupADC` package can be loaded by typing

``` {.r}
library(blupADC)
```

**Note**: In terms of the relationship matrix construction, we highly recommend Microsoft R Open(faster than traditional R many times)

### 🙊Features

-   Feature 1. Genomic data format conversion
-   Feature 2. Genomic data quality control and genotype imputation
-   Feature 3. Breed composition analysis and duplication detection of genomic data
-   Feature 4. Pedigree tracing, and analysis
-   Feature 5. Pedigree visualization
-   Feature 6. Relationship matrix construction(A,G, and H) 
-   Feature 7. Genetic evaluation with DMU
-   Feature 8. Genetic  evaluation with BLUPF90

## Usage

`blupADC` provides several datasets objects, including `data_hmp`, `origin_pedigree`.

In addition, `blupADC` provides several files which are saved in `~/blupADC/extdata`. We can get the path of these files by typing

``` {.r}
system.file("extdata", package = "blupADC") # path of provided files
```

#### Feature 1. Genomic data format conversion ([see more details](https://qsmei.netlify.app/post/blupadc/))

``` R
library(blupADC)
format_result=geno_format(
    	input_data_hmp=example_data_hmp,  # provided data variable
        output_data_type=c("Plink","BLUPF90","Numeric"),# output data format
    	output_data_path=getwd(),   #output data path      
    	output_data_name="blupADC", #output data name    
        return_result = TRUE,       #save result in R environment
        cpu_cores=1                 # number of cpu 
                  )

#convert phased VCF data to haplotype format and  haplotype-based numeric format
library(blupADC)
data_path=system.file("extdata", package = "blupADC")  #  path of example files 
phased=geno_format(
         input_data_path=data_path,      # input data path 
         input_data_name="example.vcf",  # input data name,for vcf data
         input_data_type="VCF",          # input data type
         phased_genotype=TRUE,           # whether the vcf data has been phased
         haplotype_window_nSNP=5,        # according to nSNP define haplotype-block,
    	 bigmemory_cal=TRUE,             # format conversion via bigmemory object
    	 bigmemory_data_path=getwd(),    # path of bigmemory data 
    	 bigmemory_data_name="test_blupADC", #name of bigmemory data 
         output_data_type=c("Haplotype","Numeric"),# output data format
         return_result=TRUE,             #save result in R environment
         cpu_cores=1                     # number of cpu 
                  )
```

#### Feature 2. Genomic data quality control and genotype imputation ([see more details](https://qsmei.netlify.app/post/feature-2-qc_imputation/qc_imputation/))

``` R
library(blupADC)
geno_qc_impute(
            input_data_hmp=example_data_hmp,        #provided data variable
            data_analysis_method="QC_Imputation",   #analysis method type,QC + imputatoin
            output_data_path=getwd(),               #output data path
            output_data_name="YY_data",             #output data name
            output_data_type="VCF"                  #output data format 
            )                       
```

#### Feature 3. Breed composition analysis and duplication detection of genomic data ([see more details](https://qsmei.netlify.app/post/feature-3-overlap_pca/blupadc/))

``` R
library(blupADC)
check_result=geno_check(
                  input_data_hmp=example_PCA_data_hmp,   #provided hapmap data object
                  duplication_check=FALSE,       #whether check the duplication of genotype
                  breed_check=TRUE,               # whether check the record of breed
                  breed_record=example_PCA_Breed, # provided breed record
                  output_data_path=getwd(),       #output path
                  return_result=TRUE              #save result as a R environment variable
                  )
```

#### Feature 4. Pedigree tracing, analysis ([see more details](https://qsmei.netlify.app/post/feature-4-trace_pedigree/pedigree/))

``` R
library(blupADC)
pedigree_result=trace_pedigree(
                input_pedigree=example_ped1,   #provided pedigree data variable
                trace_generation=3,            # trace generation
                output_pedigree_tree=T         # output pedigree tree
                )  
```

#### Feature 5. Pedigree visualization ([see more details](https://qsmei.netlify.app/post/feature-5-visualize_pedigree/pedigree/))

``` R
library(blupADC)
plot=ggped(
       input_pedigree=example_ped2,
       trace_id=c("121"),
       trace_sibs=TRUE   #whether plot the sibs of subset-id  
        ) 
```

#### Feature 6. Relationship matrix construction(A,G, and H)  ([see more details](https://qsmei.netlify.app/post/feature-6-kinship_matrix/relationship_matrix/))

``` R
library(blupADC)
data_path=system.file("extdata", package = "blupADC")  #  path of example files 
kinship_result=cal_kinship(
        		input_data_path=data_path,      # input data path 
        		input_data_name="example.vcf",  # input data name,for vcf data
         		input_data_type="VCF",          # input data type
    			kinship_type=c("G_A","G_D"),      #type of  kinship matrix
    			dominance_type=c("genotypic"),    #type of dominance effect
    			inbred_type=c("Homozygous"),      #type of inbreeding coefficients
    			return_result=TRUE)               #save result as a R environment variable         
```

#### Feature 7. Genetic evaluation with DMU ([see more details](https://qsmei.netlify.app/post/feature-7-run_dmu/run_dmu/))

``` R
library(blupADC)
data_path=system.file("extdata", package = "blupADC")  #  path of example files 
  
run_DMU(
        phe_col_names=c("Id","Mean","Sex","Herd_Year_Season","Litter","Trait1","Trait2","Age"), # colnames of phenotype 
        target_trait_name=list(c("Trait1")),                     #trait name 
        fixed_effect_name=list(c("Sex","Herd_Year_Season")),     #fixed effect name
        random_effect_name=list(c("Id","Litter")),               #random effect name
        covariate_effect_name=NULL,                              #covariate effect name
        phe_path=data_path,                          #path of phenotype file
        phe_name="phenotype.txt",                    #name of phenotype file
        integer_n=5,                                 #number of integer variable 
        analysis_model="PBLUP_A",                    #model of genetic evaluation
        dmu_module="dmuai",                          #modeule of estimating variance components 
        relationship_path=data_path,                 #path of relationship file 
        relationship_name="pedigree.txt",            #name of relationship file 
        output_result_path=getwd()                   # output path 
        )
```

#### Feature 8. Genetic evaluation with BLUPF90 ([see more details](https://qsmei.netlify.app/post/feature-8-run_blupf90/blupf90/))

``` R
library(blupADC)
data_path=system.file("extdata", package = "blupADC")  #  path of example files 
  
run_BLUPF90(
        phe_col_names=c("Id","Mean","Sex","Herd_Year_Season","Litter","Trait1","Trait2","Age"), # colnames of phenotype 
        target_trait_name=list(c("Trait1")),                     #trait name 
        fixed_effect_name=list(c("Sex","Herd_Year_Season")),     #fixed effect name
        random_effect_name=list(c("Id","Litter")),               #random effect name
        covariate_effect_name=NULL,                              #covariate effect name
        phe_path=data_path,                          #path of phenotype file
        phe_name="phenotype.txt",                    #name of phenotype file
        analysis_model="PBLUP_A",                    #model of genetic evaluation
        relationship_path=data_path,                 #path of relationship file 
        relationship_name="pedigree.txt",            #name of relationship file 
        output_result_path=getwd()                   # output path 
        )   
```
