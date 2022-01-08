#include "shared_function.h"

template <typename T>
std::string paste_vec_memory(arma::Col<T> data){
	
	int data_size=data.size();
	std::string tmp;
	int tmp1;
	for (int i=0;i<data_size;i++){	
		tmp1=data[i];
		tmp +=std::to_string(tmp1);
	}
    return tmp;
}


template <typename T>
std::string	get_blupf90_allele_string_phased_vcf_memory(arma::Col<T> allele_string_col1, arma::Col<T> allele_string_col2,
												 int max_length,std::vector<int> cumsum_haplo_type_num, int n_snp){
	
	std::string allele,allele1,allele2,allele_string(n_snp+max_length,'0'),ind_pos1_string,ind_pos2_string;
	int ind_pos1,ind_pos2,tmp,ind_pos1_num,ind_pos2_num;

	for(int i=0;i<allele_string_col1.size();i++){	
		tmp=cumsum_haplo_type_num[i];
		
		ind_pos1=allele_string_col1[i]+tmp+max_length;
	    ind_pos2=allele_string_col2[i]+tmp+max_length;	
		
		if(ind_pos1==ind_pos2){
		   allele_string[ind_pos1]=2+'0';
		}else{
		   allele_string[ind_pos1]=1+'0';
		   allele_string[ind_pos2]=1+'0';		
		}
	}
	for(int j=0;j<max_length;j++){	
		allele_string[j]=' ';	
	}
	return allele_string;
}


template <typename T>
std::vector<int> get_haplotype_set_memory(arma::Mat<T> & data_haplotype){
		
		int col_n=data_haplotype.n_cols;	
		std::string chr;	
		arma::Col<T> tmp;
		std::vector<std::string> strVec(col_n);
		
		for (int i=0;i<col_n;i++){	
			tmp=data_haplotype.col(i);	
			chr=paste_vec_memory(tmp);
			strVec[i]=chr;
			}
		std::vector<std::string> strVec1=strVec;	
		sort(strVec.begin(), strVec.end());	 //排序为了后续的unique			
		auto last=std::unique(strVec.begin(),strVec.end());	
		strVec.erase(last, strVec.end()); 
		
		int size_set=strVec1.size(),index;
		std::string tmp3;
		std::vector<int>  index_set(size_set); 
		
		for(int j=0;j<size_set;j++){
			vector<string>::iterator itr1 = std::find(strVec.begin(), strVec.end(), strVec1[j]);
			index = std::distance(std::begin(strVec),itr1);
			index_set[j]=index;
		}
		return index_set;
	}


template <typename T>
void allele_get_haplotype_set_short_memory(arma::Mat<T> & data_haplotype, Rcpp::List & haplotype_allele, int i_pos){
		
		int col_n=data_haplotype.n_cols;	
		std::string chr;	
		arma::Col<T> tmp;
		std::vector<std::string> strVec(col_n);
		
		for (int i=0;i<col_n;i++){	
			tmp=data_haplotype.col(i);	
			chr=paste_vec_memory(tmp);
			strVec[i]=chr;
			}
		std::vector<std::string> strVec1=strVec;	
		sort(strVec.begin(), strVec.end());	 //排序为了后续的unique			
		auto last=std::unique(strVec.begin(),strVec.end());	
		strVec.erase(last, strVec.end()); 
		
		int size_set=strVec1.size(),index;
		std::string tmp3;
		std::vector<int>  index_set(size_set); 
		
		for(int j=0;j<size_set;j++){
			vector<string>::iterator itr1 = std::find(strVec.begin(), strVec.end(), strVec1[j]);
			index = std::distance(std::begin(strVec),itr1);
			index_set[j]=index;
		}
				
		haplotype_allele[i_pos]=strVec;
	
	}

//phased_vcf_to_haplotype_memory_cpp
template <typename T>
void phased_vcf_to_haplotype_memory_cpp(const CharacterMatrix & data_vcf,XPtr<BigMatrix> pMat,int cpu_cores=1){

	Rcout<<"Start VCF to Haplotype data format conversion......"<<endl;		
	omp_set_num_threads(cpu_cores);
	arma::Mat<T> data_hap((T*) pMat->matrix(), pMat -> nrow(), pMat -> ncol(), false);
	int i,j,n_snp=data_vcf.nrow(),n_ind=data_vcf.ncol()-9;
	std::string allele1,allele2;
	std::string allele;	
	bool status_progress=true;
	Progress p(n_snp*n_ind,status_progress);		
	#pragma omp parallel for private(i,j,allele,allele1,allele2)
	for(i=0;i<n_snp;i++){
		for(j=0;j<n_ind;j++){		
			allele=data_vcf(i,j+9);
			allele1=allele[0];
			allele2=allele[2];
			data_hap(i,2*j)=std::stoi(allele1);
			data_hap(i,2*j+1)=std::stoi(allele2);
			p.increment();
		}
	}
	Rcout<<""<<endl;
	Rcout<<"Complete VCF to Haplotype data format conversion!"<<endl;	
}


// [[Rcpp::export]]
void phased_vcf_to_haplotype_memory_cpp(const CharacterMatrix & data_vcf,SEXP pBigMat_hap,int cpu_cores=1) {
    Rcpp::XPtr<BigMatrix> pMat(pBigMat_hap);
    
    switch(pMat->matrix_type()) {
    case 1:
        return phased_vcf_to_haplotype_memory_cpp<char>(data_vcf, pMat, cpu_cores);
    case 2:
        return phased_vcf_to_haplotype_memory_cpp<short>(data_vcf, pMat, cpu_cores);
    case 4:
        return phased_vcf_to_haplotype_memory_cpp<int>(data_vcf, pMat, cpu_cores);
    case 6:
        return phased_vcf_to_haplotype_memory_cpp<float>(data_vcf, pMat, cpu_cores);		
    case 8:
        return phased_vcf_to_haplotype_memory_cpp<double>(data_vcf, pMat, cpu_cores);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}


// phased_vcf to numeric 
template <typename T>
Rcpp::List phased_vcf_to_numeric_memory_cpp(
									  XPtr<BigMatrix> pMat_hap,
									  XPtr<BigMatrix> pMat_list,
									  std::string numeric_file_name,
									  std::string numeric_file_path,
                                      std::vector<int> block_start,
									  std::vector<int> block_end,
									  const CharacterMatrix & data_vcf,
									  int cpu_cores=1){


	arma::Mat<T> data_hap((T*) pMat_hap->matrix(), pMat_hap -> nrow(), pMat_hap -> ncol(), false);
	arma::Mat<T> haplo_pos_list3((T*) pMat_list->matrix(), pMat_list -> nrow(), pMat_list -> ncol(), false);



	omp_set_num_threads(cpu_cores);
	
	int i,j,n_snp=data_vcf.nrow(),n_ind=data_vcf.ncol()-9,m,n,tmp;
	std::string allele1,allele2;
	std::string allele;	
	Progress p0(n_snp*n_ind,true);
	Rcout<<"bigmemory-Start Phased haplotype data format conversion......"<<endl;
	#pragma omp parallel for private(i,j,allele,allele1,allele2)
	for(i=0;i<n_snp;i++){
		for(j=0;j<n_ind;j++){	
			p0.increment();
			allele=data_vcf(i,j+9);
			allele1=allele[0];
			allele2=allele[2];
			data_hap(i,2*j)=std::stoi(allele1);
			data_hap(i,2*j+1)=std::stoi(allele2);
		}
	}
	Rcout<<""<<endl;
	Rcout<<"bigmemory-Complete Phased haplotype data format conversion!"<<endl;
//hapmap to numeric 
	Rcout<<"bigmemory-Start construct haplotype list...... "<<endl;	
	int window_n=block_start.size();
	Function generate_list("vector");
	Rcpp::List haplotype_allele=generate_list("list", window_n);
	std::vector<int> haplo_pos,haplo_type_num(window_n),cumsum_haplo_type_num;
	arma::Mat<T> data_haplotype_subset;
	// get the overview set of haplo_data
	for(i=0; i < window_n; i++){
		data_haplotype_subset=data_hap.rows(block_start[i],block_end[i]);
		allele_get_haplotype_set_short_memory(data_haplotype_subset,haplotype_allele,i);
    }
	Progress p(window_n*n_ind,true);					
	#pragma omp parallel for private(i,j,data_haplotype_subset,haplo_pos)
	for(i=0; i < window_n; i++){
		data_haplotype_subset=data_hap.rows(block_start[i],block_end[i]);
		haplo_pos=get_haplotype_set_memory(data_haplotype_subset);
		haplo_type_num[i]=*max_element(haplo_pos.begin(),haplo_pos.end())+1;
		
		for(j=0;j<n_ind;j++){
			p.increment();
			haplo_pos_list3(i,2*j)=haplo_pos[2*j];
			haplo_pos_list3(i,2*j+1)=haplo_pos[2*j+1];
			
		}
	}
	Rcout<<" "<<endl;
	Rcout<<"bigmemory-Complete construct haplotype list!"<<endl;
	int sum_haplo_type_num=std::accumulate(haplo_type_num.begin(),haplo_type_num.end(),0);
	int ind_pos1,ind_pos2;  // all types of haplot
	
	cumulativeSum(haplo_type_num,cumsum_haplo_type_num);
	
	int type_num=pMat_hap->matrix_type();
	std::string type;
	if(type_num==1){
		type="char";
	}else if(type_num==2){
		type="short";
	}else if(type_num==4){
		type="integer";
	}else if(type_num==8){	
		type="double";
	}else{
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
	
	SEXP pBigMat3=make_bigmemory_object_address_cpp(n_ind,sum_haplo_type_num,numeric_file_name,numeric_file_path,type);
	
	Rcpp::XPtr<BigMatrix> pMat_num(pBigMat3);
	
	arma::Mat<T> data_numeric((T*) pMat_num->matrix(), pMat_num -> nrow(), pMat_num -> ncol(), false);
	Rcout<<"bigmemory-Start phased-haplotype-numeric(0,1,2) format conversion...... "<<endl;	
	Progress p2(window_n*n_ind,true);
		
	#pragma omp parallel for private(m,n,tmp,ind_pos1,ind_pos2)
	for(m=0;m<window_n;m++){
		
		tmp=cumsum_haplo_type_num[m];
		
		for(n=0;n<n_ind;n++){
			p2.increment();
			ind_pos1=haplo_pos_list3(m,2*n)+tmp;
		    ind_pos2=haplo_pos_list3(m,2*n+1)+tmp;
		    data_numeric(n,ind_pos1)=data_numeric(n,ind_pos1)+1;
		    data_numeric(n,ind_pos2)=data_numeric(n,ind_pos2)+1;
		}
	}
    Rcout<<""<<endl;			
	Rcout<<"bigmemory-Complete phased-haplotype-numeric(0,1,2) format conversion!"<<endl;
		return haplotype_allele;	
}


// [[Rcpp::export]]
Rcpp::List phased_vcf_to_numeric_memory_cpp(SEXP pBigMat_hap,
								  SEXP pBigMat_list,
								  std::string numeric_file_name,
								  std::string numeric_file_path,											  
								  std::vector<int> block_start,
								  std::vector<int> block_end,
								  const CharacterMatrix & data_vcf,
								  int cpu_cores=1){
    Rcpp::XPtr<BigMatrix> pMat_hap(pBigMat_hap);
	Rcpp::XPtr<BigMatrix> pMat_list(pBigMat_list);
    
    switch(pMat_hap->matrix_type()) {
    case 1:
        return phased_vcf_to_numeric_memory_cpp<char>(pMat_hap,pMat_list,numeric_file_name,numeric_file_path,block_start,block_end,data_vcf,cpu_cores);
    case 2:
        return phased_vcf_to_numeric_memory_cpp<short>(pMat_hap,pMat_list,numeric_file_name,numeric_file_path,block_start,block_end,data_vcf,cpu_cores);
    case 4:
        return phased_vcf_to_numeric_memory_cpp<int>(pMat_hap,pMat_list,numeric_file_name,numeric_file_path,block_start,block_end,data_vcf,cpu_cores);
    case 6:
        return phased_vcf_to_numeric_memory_cpp<float>(pMat_hap,pMat_list,numeric_file_name,numeric_file_path,block_start,block_end,data_vcf,cpu_cores);
    case 8:
        return phased_vcf_to_numeric_memory_cpp<double>(pMat_hap,pMat_list,numeric_file_name,numeric_file_path,block_start,block_end,data_vcf,cpu_cores);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}



//unphased-vcf to numeric
template <typename T>
void unphased_vcf_to_numeric_memory_cpp(XPtr<BigMatrix> pMat_num,
												const CharacterMatrix & data_vcf,
												int cpu_cores=1,
												std::string miss_base=".",
												int miss_base_num=0){
	
	Rcout<<"Start VCF to Numeric data format conversion......"<<endl;		
	arma::Mat<T> data_numeric((T*) pMat_num->matrix(), pMat_num -> nrow(), pMat_num -> ncol(), false);
	
	omp_set_num_threads(cpu_cores);	
	int i,j,n_snp=data_vcf.nrow(),n_ind=data_vcf.ncol()-9,allele_num;
	std::string allele1,allele2;
	std::string allele;
	Progress p(n_snp*n_ind,true);		
	#pragma omp parallel for private(i,j,allele,allele1,allele2,allele_num)
	for(i=0;i<n_snp;i++){
		for(j=0;j<n_ind;j++){
			p.increment();
			allele=data_vcf(i,j+9);
			allele1=allele[0];
			allele2=allele[2];
			
			if(allele1==miss_base|allele2==miss_base){
			data_numeric(j,i)=miss_base_num;
			}else{
			data_numeric(j,i)=std::stoi(allele1)+std::stoi(allele2);		
			}
		}
	}
	Rcout<<""<<endl;	
	Rcout<<"Complete VCF to Numeric data format conversion!"<<endl;	

}

// [[Rcpp::export]]
void unphased_vcf_to_numeric_memory_cpp(SEXP pBigMat_num,
									const CharacterMatrix & data_vcf,
									int cpu_cores=1,
									std::string miss_base=".",
									int miss_base_num=0){
										
    Rcpp::XPtr<BigMatrix> pMat_num(pBigMat_num);

    
    switch(pMat_num->matrix_type()) {
    case 1:
        return unphased_vcf_to_numeric_memory_cpp<char>(pMat_num,data_vcf,cpu_cores,miss_base,miss_base_num);
    case 2:
        return unphased_vcf_to_numeric_memory_cpp<short>(pMat_num,data_vcf,cpu_cores,miss_base,miss_base_num);
    case 4:
        return unphased_vcf_to_numeric_memory_cpp<int>(pMat_num,data_vcf,cpu_cores,miss_base,miss_base_num);
    case 6:
        return unphased_vcf_to_numeric_memory_cpp<float>(pMat_num,data_vcf,cpu_cores,miss_base,miss_base_num);
    case 8:
        return unphased_vcf_to_numeric_memory_cpp<double>(pMat_num,data_vcf,cpu_cores,miss_base,miss_base_num);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}

// phased_vcf to numeric 
template <typename T>
Rcpp::List phased_vcf_blupf90_memory_cpp(
									  XPtr<BigMatrix> pMat_hap,
									  XPtr<BigMatrix> pMat_list,
                                      std::vector<int> block_start,
									  std::vector<int> block_end,
									  const CharacterMatrix & data_vcf,
									  int cpu_cores=1){


	arma::Mat<T> data_hap((T*) pMat_hap->matrix(), pMat_hap -> nrow(), pMat_hap -> ncol(), false);
	arma::Mat<T> haplo_pos_list3((T*) pMat_list->matrix(), pMat_list -> nrow(), pMat_list -> ncol(), false);



	omp_set_num_threads(cpu_cores);
	
	int i,j,n_snp=data_vcf.nrow(),n_ind=data_vcf.ncol()-9,m,n,tmp;
	std::string allele1,allele2;
	std::string allele;	
	CharacterVector IND_name=colnames(data_vcf);
	IND_name.erase(0,9);
	Rcout<<"bigmemory-Start Phased haplotype data format conversion......"<<endl;
	#pragma omp parallel for private(i,j,allele,allele1,allele2)
	for(i=0;i<n_snp;i++){
		for(j=0;j<n_ind;j++){		
			allele=data_vcf(i,j+9);
			allele1=allele[0];
			allele2=allele[2];
			data_hap(i,2*j)=std::stoi(allele1);
			data_hap(i,2*j+1)=std::stoi(allele2);
		}
	}
	Rcout<<""<<endl;
	Rcout<<"bigmemory-Complete Phased haplotype data format conversion!"<<endl;
//hapmap to numeric 
	Rcout<<"bigmemory-Start phased-haplotype-BLUPF90 format conversion...... "<<endl;	
	int window_n=block_start.size();
	Function generate_list("vector");
	Rcpp::List haplotype_allele=generate_list("list", window_n);	
	std::vector<int> haplo_pos,haplo_type_num(window_n),cumsum_haplo_type_num;
	arma::Mat<T> data_haplotype_subset;
	// get the overview set of haplo_data
	for(i=0; i < window_n; i++){
		data_haplotype_subset=data_hap.rows(block_start[i],block_end[i]);
		allele_get_haplotype_set_short_memory(data_haplotype_subset,haplotype_allele,i);
    }
	Progress p(window_n*n_ind,true);					
	#pragma omp parallel for private(i,j,data_haplotype_subset,haplo_pos)
	for(i=0; i < window_n; i++){
		data_haplotype_subset=data_hap.rows(block_start[i],block_end[i]);
		haplo_pos=get_haplotype_set_memory(data_haplotype_subset);
		haplo_type_num[i]=*max_element(haplo_pos.begin(),haplo_pos.end())+1;
		
		for(j=0;j<n_ind;j++){
			p.increment();
			haplo_pos_list3(i,2*j)=haplo_pos[2*j];
			haplo_pos_list3(i,2*j+1)=haplo_pos[2*j+1];
			
		}
	}
	int sum_haplo_type_num=std::accumulate(haplo_type_num.begin(),haplo_type_num.end(),0);
	int ind_pos1,ind_pos2;  // all types of haplot
	
	cumulativeSum(haplo_type_num,cumsum_haplo_type_num);
	
	int n_snp_hap=sum_haplo_type_num;
	std::string ind_name;
	std::vector<int> ind_name_size(IND_name.size());
	//确定个体名称中最大字符串长度
	for(i=0;i<IND_name.size();i++){
		ind_name=IND_name[i];
		ind_name_size[i]=ind_name.length();
	}
	int max_length=*max_element(ind_name_size.begin(),ind_name_size.end())+3;	
	max_length=0;
	std::string allele_string(n_snp_hap+max_length,' ');
	std::vector<std::string> allele_string_vector(n_ind);
	arma::Col<T> allele_string_col1,allele_string_col2;
	#pragma omp parallel for private(j,allele_string_col1,allele_string_col2,allele_string)
	for(j=0;j<n_ind;j++){
	
		allele_string_col1=haplo_pos_list3.col(2*j);
		allele_string_col2=haplo_pos_list3.col(2*j+1);
		
		allele_string=get_blupf90_allele_string_phased_vcf_memory(allele_string_col1,allele_string_col2,max_length,cumsum_haplo_type_num,n_snp_hap);
	
		allele_string_vector[j]=allele_string;	
	}
	
	//添加个体名称
	//for(j=0;j<n_ind;j++){		
	//	ind_name=IND_name[j];	
	//	for(int k=0;k<ind_name.size();k++){
	//		allele_string_vector[j][k]=ind_name[k];
	//	}
	//}	
	
    Rcout<<""<<endl;		
	Rcout<<"bigmemory-Complete phased-haplotype-BLUPF90 format conversion!"<<endl;
		return List::create(Named("allele_string_vector") = allele_string_vector,
							Named("haplotype_allele")=haplotype_allele);	
		
}


// [[Rcpp::export]]
Rcpp::List phased_vcf_blupf90_memory_cpp(SEXP pBigMat_hap,
								  SEXP pBigMat_list,										  
								  std::vector<int> block_start,
								  std::vector<int> block_end,
								  const CharacterMatrix & data_vcf,
								  int cpu_cores=1){
    Rcpp::XPtr<BigMatrix> pMat_hap(pBigMat_hap);
	Rcpp::XPtr<BigMatrix> pMat_list(pBigMat_list);
    
    switch(pMat_hap->matrix_type()) {
    case 1:
        return phased_vcf_blupf90_memory_cpp<char>(pMat_hap,pMat_list,block_start,block_end,data_vcf,cpu_cores);
    case 2:
        return phased_vcf_blupf90_memory_cpp<short>(pMat_hap,pMat_list,block_start,block_end,data_vcf,cpu_cores);
    case 4:
        return phased_vcf_blupf90_memory_cpp<int>(pMat_hap,pMat_list,block_start,block_end,data_vcf,cpu_cores);
    case 6:
        return phased_vcf_blupf90_memory_cpp<float>(pMat_hap,pMat_list,block_start,block_end,data_vcf,cpu_cores);		
    case 8:
        return phased_vcf_blupf90_memory_cpp<double>(pMat_hap,pMat_list,block_start,block_end,data_vcf,cpu_cores);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}