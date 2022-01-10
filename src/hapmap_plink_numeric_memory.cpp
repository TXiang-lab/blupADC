#include "shared_function.h"

// Numgeno  in Rcpp
template <typename T>
void hapmap_to_numeric_memory_cpp(CharacterMatrix & data_hmp, 
						   XPtr<BigMatrix> pMat_num,
						   std::string miss_base="N", 
						   int miss_base_num=0, 
						   int cpu_cores=1,
						   std::string type="integer"){  //输入数据，行为SNP, 列为个体

	omp_set_num_threads(cpu_cores);
	arma::Mat<T> data_numeric((T*) pMat_num->matrix(), pMat_num -> nrow(), pMat_num -> ncol(), false);
	std::string tmp_allele,tmp_allele1,tmp_allele2;	
	int n_snp=data_hmp.nrow(),n_ind=data_hmp.ncol()-11,i,j;
	CharacterVector snp_vec;
	std::string factor_result,Ref,Alt,allele,allele1,allele2;
	CharacterMatrix Ref_type(n_snp,2);	
	CharacterVector snp1,snp2;
	Rcout<<"bigmemory-Start Hapmap to Numeric data format conversion......"<<endl;
		
	for(i=0;i<n_snp;i++){	
		snp_vec=data_hmp.row(i);
		snp_vec.erase(0,11);
		factor_result=pair_base_factor_cpp(snp_vec,miss_base=miss_base);
		Ref=factor_result[0];
		Alt=factor_result[1];
		Ref_type(i,0)=Ref;
		Ref_type(i,1)=Alt;
	}

    bool status_progress=true,allele1_status,allele2_status;
	Progress p(n_snp*n_ind,status_progress);	
	#pragma omp parallel for private(i,j,Ref,allele,allele1,allele2,allele1_status,allele2_status)	
	for(i=0;i<n_snp;i++){		
	
		Ref=Ref_type(i,0);
		
		for(j=0;j<n_ind;j++){
			p.increment();
			allele=data_hmp(i,j+11);
			allele1=allele[0];
			allele2=allele[1];
			
			allele1_status=(allele1==Ref);
			allele2_status=(allele2==Ref);
			if(allele1_status&&allele2_status){					
				data_numeric(j,i)=0;		
			}else if(allele1==miss_base|allele2==miss_base){				
				data_numeric(j,i)=miss_base_num;  //将缺失值视为隐性纯合，设为0						
			}else if(allele1_status||allele2_status){	
				data_numeric(j,i)=1;					
			}else {					
				data_numeric(j,i)=2;
			}			
		}
	}
	Rcout<<" "<<endl;	
	Rcout<<"bigmemory-Complete Hapmap to Numeric data format conversion!"<<endl;
}




// [[Rcpp::export]]
void hapmap_to_numeric_memory_cpp(CharacterMatrix & data_hmp, 
						   SEXP pBigMat_num,
						   std::string miss_base="N", 
						   int miss_base_num=0, 
						   int cpu_cores=1,
						   std::string type="short"){
							   
    Rcpp::XPtr<BigMatrix> pMat_num(pBigMat_num);
    int type_number;
	if(type=="char"){
		type_number=1;
	}else if(type=="short"){
		type_number=2;
	}else if(type=="integer"){
		type_number=4;
	}else if(type=="float"){	
		type_number=6;	
	}else if(type=="double"){
		type_number=8;
	}
	
	
    switch(type_number) {
    case 1:
        return hapmap_to_numeric_memory_cpp<char>(data_hmp,pMat_num,miss_base,miss_base_num,cpu_cores,type);
    case 2:
        return hapmap_to_numeric_memory_cpp<short>(data_hmp,pMat_num,miss_base,miss_base_num,cpu_cores,type);
    case 4:
        return hapmap_to_numeric_memory_cpp<int>(data_hmp,pMat_num,miss_base,miss_base_num,cpu_cores,type);
    case 6:
        return hapmap_to_numeric_memory_cpp<float>(data_hmp,pMat_num,miss_base,miss_base_num,cpu_cores,type);
    case 8:
        return hapmap_to_numeric_memory_cpp<double>(data_hmp,pMat_num,miss_base,miss_base_num,cpu_cores,type);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}



template <typename T>
void plink_to_numeric_memory_cpp(CharacterMatrix & data_ped,
									 CharacterMatrix data_map,
									 XPtr<BigMatrix> pMat_num,
									 int cpu_cores=5, 
									 std::string miss_base="0", 
									 int miss_base_num=0,
						             std::string type="short"){  //输入数据，行为SNP, 列为个体
	
			omp_set_num_threads(cpu_cores);	
			arma::Mat<T> data_numeric((T*) pMat_num->matrix(), pMat_num -> nrow(), pMat_num -> ncol(), false);
			//转换数据框为字符串矩阵
            int i,j,k,n_snp=(data_ped.ncol()-6)/2,n_ind=data_ped.nrow();
			CharacterVector temp_vec,temp_vec1,temp_vec2,SNP_name=data_map.column(0);
			CharacterMatrix Ref_type(n_snp,2);

			std::string base_SNP1,base_SNP2,temp_char,temp_base_SNP1,temp_base_SNP2,Ref,Alt;

	        Rcout<<"bigmemory-Start Plink to numeric(0,1,2) format conversion...... "<<endl;			
	        for(k=0;k<n_snp;k++){	 
	        	temp_vec1=data_ped.column(2*k+6);
				temp_vec2=data_ped.column(2*k+7);
				temp_base_SNP1=single_base_factor_cpp(temp_vec1,miss_base=miss_base);
				temp_base_SNP2=single_base_factor_cpp(temp_vec2,miss_base=miss_base);

		         if(temp_base_SNP1<temp_base_SNP2){
		         	Ref=temp_base_SNP1[0];
		         	Alt=temp_base_SNP2[1];		
		         }else{
		         	Ref=temp_base_SNP2[0];
		         	Alt=temp_base_SNP1[1];								
		         }				
				Ref_type(k,0)=Ref;
				Ref_type(k,1)=Alt;
	        }				
			bool allele1_status,allele2_status;			
			std::string allele1,allele2;
             Progress p(n_snp*n_ind,true);	
#pragma omp parallel for private(i,j,base_SNP1,allele1_status,allele2_status)				 
			for( j=0; j<n_snp;j++){	
				base_SNP1=Ref_type(j,0);
				for( i=0; i<n_ind;i++){					
					p.increment();
					allele1=data_ped(i,2*j+6);
					allele2=data_ped(i,2*j+7);
					
					allele1_status=(allele1==base_SNP1);
					allele2_status=(allele2==base_SNP1);
					if(allele1_status&&allele2_status){					
						data_numeric(i,j)=0;		
					}else if(allele1==miss_base|allele2==miss_base){				
						data_numeric(i,j)=miss_base_num;  //将缺失值视为隐性纯合，设为0						
					}else if(allele1_status||allele2_status){	
						data_numeric(i,j)=1;					
					}else {					
						data_numeric(i,j)=2;
					}
			}}
			Rcout<<" "<<endl;
	        Rcout<<"bigmemory-Complete Plink to numeric(0,1,2) format conversion...... "<<endl;			
}


// [[Rcpp::export]]
void plink_to_numeric_memory_cpp(CharacterMatrix & data_ped,
						  CharacterMatrix data_map,
						  SEXP pBigMat_num,
						  int cpu_cores=5, 
						  std::string miss_base="0", 
						  int miss_base_num=0,
						    std::string type="short"){
							   
    Rcpp::XPtr<BigMatrix> pMat_num(pBigMat_num);
    int type_number;
	if(type=="char"){
		type_number=1;
	}else if(type=="short"){
		type_number=2;
	}else if(type=="integer"){
		type_number=4;
	}else if(type=="float"){	
		type_number=6;	
	}else if(type=="double"){
		type_number=8;
	}
	
	
	
    switch(type_number) {
    case 1:
        return plink_to_numeric_memory_cpp<char>(data_ped,data_map,pMat_num,cpu_cores,miss_base,miss_base_num,type);
    case 2:
        return plink_to_numeric_memory_cpp<short>(data_ped,data_map,pMat_num,cpu_cores,miss_base,miss_base_num,type);
    case 4:
        return plink_to_numeric_memory_cpp<int>(data_ped,data_map,pMat_num,cpu_cores,miss_base,miss_base_num,type);
    case 6:
        return plink_to_numeric_memory_cpp<float>(data_ped,data_map,pMat_num,cpu_cores,miss_base,miss_base_num,type);
    case 8:
        return plink_to_numeric_memory_cpp<double>(data_ped,data_map,pMat_num,cpu_cores,miss_base,miss_base_num,type);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}



template <typename T>
void blupf90_to_numeric_memory_cpp(std::vector<std::string> & data_blupf90,
							XPtr<BigMatrix> pMat_num,
							int cpu_cores=1,
							std::string type="short"){  //输入数据，行为SNP, 列为个体		
			
			omp_set_num_threads(cpu_cores);
			arma::Mat<T> data_numeric((T*) pMat_num->matrix(), pMat_num -> nrow(), pMat_num -> ncol(), false);
			Rcout<<"Start BLUPF90 to Numeric format conversion......"<<endl;	
			int i,j,k;	
			std::string tmp_ind_string,allele_string;
			int n_ind=data_blupf90.size(),n_snp=data_blupf90[0].size();		
			
			#pragma omp parallel for private(i,j,tmp_ind_string,allele_string)
			for(i=0;i<n_ind;i++){
				
				tmp_ind_string=data_blupf90[i];
				
				for(j=0;j<n_snp;j++){
				
					allele_string=tmp_ind_string[j];
					
					data_numeric(i,j)=std::stoi(allele_string);
				}
			}		
			
			Rcout<<" "<<endl;			
			Rcout<<"Complete BLUPF90 to Numeric format conversion!"<<endl;	
}




// [[Rcpp::export]]
void blupf90_to_numeric_memory_cpp(std::vector<std::string> & data_blupf90,
							SEXP pBigMat_num,
							int cpu_cores=1,
							std::string type="short"){
							   
    Rcpp::XPtr<BigMatrix> pMat_num(pBigMat_num);
    int type_number;
	if(type=="char"){
		type_number=1;
	}else if(type=="short"){
		type_number=2;
	}else if(type=="integer"){
		type_number=4;
	}else if(type=="float"){	
		type_number=6;	
	}else if(type=="double"){
		type_number=8;
	}
	
	
	
    switch(type_number) {
    case 1:
        return blupf90_to_numeric_memory_cpp<char>(data_blupf90,pMat_num,cpu_cores,type);
    case 2:
        return blupf90_to_numeric_memory_cpp<short>(data_blupf90,pMat_num,cpu_cores,type);
    case 4:
        return blupf90_to_numeric_memory_cpp<int>(data_blupf90,pMat_num,cpu_cores,type);
    case 6:
        return blupf90_to_numeric_memory_cpp<float>(data_blupf90,pMat_num,cpu_cores,type);
    case 8:
        return blupf90_to_numeric_memory_cpp<double>(data_blupf90,pMat_num,cpu_cores,type);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}