#include "shared_function.h"

// [[Rcpp::export]]
void  cumulativeSum(const std::vector<int> input, std::vector<int>& result) {
	result.push_back(0);
	result.push_back(input[0]);
	for (int i = 2; i < input.size(); i++) {
		result.push_back(result[i - 1] + input[i-1]);
	}
	
}


// [[Rcpp::export]]
std::string paste_vec_short( arma::Col<int> &data){
	
	int data_size=data.size();
	std::string tmp;
	int tmp1;
	for (int i=0;i<data_size;i++){	
		tmp1=data[i];
		tmp +=std::to_string(tmp1);
	}
    return tmp;
}

// [[Rcpp::export]]
std::vector<int> get_haplotype_set_short( arma::Mat<int> & data_haplotype){
		
		int col_n=data_haplotype.n_cols;	
		std::string chr;	
		arma::Col<int> tmp;
		std::vector<std::string> strVec(col_n);
		
		for (int i=0;i<col_n;i++){	
			tmp=data_haplotype.col(i);	
			chr=paste_vec_short(tmp);
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


//获取haplotype to numeric的基因型
// [[Rcpp::export]]
void allele_get_haplotype_set_short( arma::Mat<int> & data_haplotype,Rcpp::List & haplotype_allele, int i_pos){
		
		int col_n=data_haplotype.n_cols;	
		std::string chr;	
		arma::Col<int> tmp;
		std::vector<std::string> strVec(col_n);
		
		for (int i=0;i<col_n;i++){	
			tmp=data_haplotype.col(i);	
			chr=paste_vec_short(tmp);
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

// [[Rcpp::export]]
List define_block_window_kb_cpp(std::vector<int> block1,std::vector<int> block2,std::vector<int> tmp_data_map_3){
	
	std::vector<int> block_start_tmp(block1.size());
	std::vector<int> block_end_tmp(block1.size());
	int tmp1,tmp2;
	
	for(int i=0;i<block1.size();i++){	
		
		tmp1=block1[i];
		tmp2=block2[i];

		for(int j=0;j<tmp_data_map_3.size();j++){			
			if(tmp_data_map_3[j]>=tmp1){
			   block_start_tmp[i]=j;
			   break;
			}		
		}
		
		for(int k=0;k<tmp_data_map_3.size();k++){			
			if(tmp_data_map_3[k]>tmp2){
			   block_end_tmp[i]=k-1;
			   break;
			}	
		}				
	}
	if(tmp_data_map_3[tmp_data_map_3.size()-1]==block2[block1.size()-1]){block_end_tmp[block_end_tmp.size()-1]=tmp_data_map_3.size()-1;}
	
	for(int i=0;i<=block_start_tmp.size()-1;++i){block_start_tmp[i]=block_start_tmp[i]+1;}
	for(int i=0;i<=block_end_tmp.size()-1;++i){block_end_tmp[i]=block_end_tmp[i]+1;}
	return List::create(Named("block_start_tmp") = block_start_tmp,_["block_end_tmp"] = block_end_tmp);
	
}

// [[Rcpp::export]]
std::string single_base_factor_cpp(CharacterVector data, std::string miss_base){  
	int N_size=data.size();
	std::string SNP,a,A,target;
	for(int j=0;j<N_size;j++){
	SNP=data[j];
	if(SNP!=miss_base){		
	a=SNP[0];
	A=SNP[0];
	break;
	}
	if(j==N_size-1){
		return miss_base+miss_base; // All bases are missing 
	}
	}
	
	for(int i=0;i<N_size;i++){		
		target=data[i];		
		if(target!=a && target!=miss_base){
			
			if(a<target){
				A=target;
				break;
			}else if(a>target){
				a=target;
				break;
			}
			}
		}

	return a+A 	; // 第一个字符为ref,第二个字符为 alt
	 }




// [[Rcpp::export]]
std::string pair_base_factor_cpp(CharacterVector data, std::string miss_base){  //  主要用于处理hapmap的行数据，获取Reference allele 和 Alternative allele

	int N_size=data.size();
	std::string miss_base_set=miss_base+miss_base;   // 由缺失的碱基组成的set
	std::string SNP,a,A,aa,target,target1,target2;
	for(int j=0;j<N_size;j++){
	SNP=data[j];
	if(SNP!=miss_base_set){		
	a=SNP[0];
	A=SNP[0];
	aa=a+a;
	break;
	}
	}
	for(int i=0;i<N_size;i++){		
		target=data[i];		
		if(target!=aa && target!=miss_base_set){
			target1=target[0];
			target2=target[1];
			
			if(a<target1){
				A=target1;
				break;
			}else if(a<target2){
				A=target2;
				break;
			}else if(a>target1){
				a=target1;
				break;
			}else if(a>target2){	
				a=target2;	
				break;
			}
		}
	}
	CharacterVector base_pair(2);
	base_pair[0]=a;
	base_pair[1]=A;
	return a+A 	; // 第一个字符为ref,第二个字符为 alt
	 }
	 

// [[Rcpp::export]]
std::string	get_blupf90_allele_string_phased_vcf(arma::Col<int> allele_string_col1, arma::Col<int> allele_string_col2,
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



// [[Rcpp::export]]
std::string	get_blupf90_allele_string_numeric(arma::Row<int> allele_string_row,int max_length){

	std::string allele_string(allele_string_row.size()+max_length,' ');
	int allele_num;
	for(int i=0;i<allele_string_row.size();i++){		
	allele_num=allele_string_row[i];
	allele_string[i+max_length]=allele_num+'0';		
	}
	
	return allele_string;
}

// [[Rcpp::export]]
SEXP make_bigmemory_object_cpp(int nrow, int ncol,std::string file_name,std::string file_path, std::string type){

  Environment stats("package:bigmemory");
  Function filebacked_big_matrix = stats["filebacked.big.matrix"];
  Function attach_big_matrix = stats["attach.big.matrix"];

  std::string backingfile=file_name+".bk",descriptorfile=file_name+".desc",backingpath=file_path;
  std::string file1=file_path+"/"+file_name+".bk";
  std::string file2=file_path+"/"+file_name+".desc";

	fstream fs1,fs2;
	fs1.open(file1, ios::in);
	fs2.open(file2, ios::in);
	if(!fs1&!fs2){
	
	filebacked_big_matrix(Named("nrow", nrow),
						  Named("ncol", ncol),
						  Named("type", type),
						  Named("backingfile", backingfile),
						  Named("backingpath", backingpath),
						  Named("descriptorfile", descriptorfile));	
		

	}
	return (attach_big_matrix(file2));
}


// [[Rcpp::export]]
SEXP make_bigmemory_object_address_cpp(int nrow, int ncol,std::string file_name,std::string file_path, std::string type){

  Environment stats("package:bigmemory");
  Function filebacked_big_matrix = stats["filebacked.big.matrix"];
  Function attach_big_matrix = stats["attach.big.matrix"];
  Function attr_cpp("attr");

  std::string backingfile=file_name+".bk",descriptorfile=file_name+".desc",backingpath=file_path;
  std::string file1=file_path+"/"+file_name+".bk";
  std::string file2=file_path+"/"+file_name+".desc";

	fstream fs1,fs2;
	fs1.open(file1, ios::in);
	fs2.open(file2, ios::in);
	if(!fs1&!fs2){
	
	filebacked_big_matrix(Named("nrow", nrow),
						  Named("ncol", ncol),
						  Named("type", type),
						  Named("backingfile", backingfile),
						  Named("backingpath", backingpath),
						  Named("descriptorfile", descriptorfile));	
		

	}
	return attr_cpp(attach_big_matrix(file2),"address");
}	

// [[Rcpp::export]]
arma::Mat<int> NumericMatrix_to_arma(NumericMatrix & data_numeric){ 


	return Rcpp::as<arma::Mat<int>>(data_numeric);

}

// [[Rcpp::export]]
arma::Mat<int> DataFrame_to_arma(DataFrame & data_numeric){ 


	return Rcpp::as<arma::Mat<int>>(internal::convert_using_rfunction(data_numeric, "as.matrix"));

}


// [[Rcpp::export]]
IntegerVector get_offspring_generation_cpp(DataFrame ped, CharacterVector IND_base){
	
CharacterVector IND_base_set=IND_base,Offspring=ped[0],Sire=ped[1],Dam=ped[2];
IntegerVector pos,pos1,pos2,generation(ped.nrow()),IND_base_set_pos;
LogicalVector status1,status2,status;

	int i=0;
while(IND_base_set.size()>0){

	i=i+1;
	status1=in(Sire,IND_base_set);

	status2=in(Dam,IND_base_set);

	status=status1|status2;	

	IND_base_set=Offspring[status];

	IND_base_set=na_omit(IND_base_set);
	IND_base_set_pos=match(IND_base_set,Offspring)-1;
	generation[IND_base_set_pos]=i;

}

return generation;
}

// [[Rcpp::export]]
List single_pedigree_cpp(CharacterMatrix ped){
	
CharacterVector Offspring=ped.column(0),Sire=ped.column(1),Dam=ped.column(2);
LogicalVector status1,status2,status,s1;

status1=is_na(Sire);
status2=is_na(Dam);
status=status1 & status2;
CharacterVector f0=Offspring[status];

Offspring=na_omit(Offspring);
Sire=na_omit(Sire);
Dam=na_omit(Dam);

CharacterVector n=unique(Offspring),P=union_(Sire,Dam);
P=unique(P);

CharacterVector f1=setdiff(P,n);
f1=union_(f0,f1);
f1=unique(f1);
CharacterVector f2=setdiff(P,f1);
CharacterVector fn=setdiff(n,f1);
fn=setdiff(fn,f2);

return List::create(Named("a") = sort_unique(f1),
					Named("b") = sort_unique(f2),
					Named("c") = sort_unique(fn));
}





// [[Rcpp::export]]
std::string	get_blupf90_allele_string_unphased_haplotype(arma::Col<int> allele_string_col1, arma::Col<int> allele_string_col2,
												 int max_length){

	std::string allele_string(allele_string_col1.size()+max_length,' ');
	int allele_num,allele1_num,allele2_num;
	for(int i=0;i<allele_string_col1.size();i++){		
	allele1_num=allele_string_col1[i];
	allele2_num=allele_string_col2[i];	
	
	allele_num=allele1_num+allele2_num;
	allele_string[i+max_length]=allele_num+'0';		
	}
	
	return allele_string;
}



// [[Rcpp::export]]
CharacterVector union_cpp(CharacterVector X, CharacterVector Y){ //合并字符串向量
	CharacterVector Z(X.size()+Y.size());
		
	for(int i=0;i<X.size();i++){
		Z[i]=X[i];
	}
	
	for(int i=0;i<Y.size();i++){
		Z[i+X.size()]=Y[i];
	}
	return Z;
}


// [[Rcpp::export]]
void delete_bigmemory_file_cpp(std::string  matrix_type,
								std::string bigmemory_data_name,
						       std::string bigmemory_data_path,
							   bool message){
								   

		std::string file1,file2;
		fstream fs1,fs2;
		
		file1=bigmemory_data_path+"/"+bigmemory_data_name+"_"+matrix_type+".desc";
		file2=bigmemory_data_path+"/"+bigmemory_data_name+"_"+matrix_type+".bk";

		fs1.open(file1, ios::in);
	    fs2.open(file2, ios::in);

	    if(fs1&&fs2){
			if(message==true){
			Rcout<<"Found "+matrix_type+".desc & "+matrix_type+".bk files in path:"+ bigmemory_data_path<<endl;
			Rcout<<"Software will delete these files automatically!"<<endl;
			}
			remove(file1.c_str());
			remove(file2.c_str());	
		}
								   
								   
}


// [[Rcpp::export]]
arma::Mat<int> get_allele(arma::Mat<int> Pedigree){
    int n_ind=Pedigree.n_rows;	
	arma::Col<int> Sire=Pedigree.col(1);
	arma::Col<int> Dam=Pedigree.col(2);
	arma::Mat<int> allele(n_ind,3);
	int Sire_i_pos,Dam_i_pos;
	arma::uvec founder_s=find(Sire==0);  //父亲为0的个体
	arma::uvec founder_d=find(Dam==0);   //母亲为0的个体
   	
	allele.col(0)=Pedigree.col(0);
	allele(founder_s, arma::uvec {1})=conv_to<arma::Mat<int>>::from(founder_s)+1; //为这些个体的父亲分配allele
	allele(founder_d, arma::uvec {2})=conv_to<arma::Mat<int>>::from(founder_d)+1+founder_s.size(); //为这些个体的母亲分配allele
                                                                                  //这些个体相当于基础群的个体
																				  //所有子代都是从他们的基因型里面	
//为所有的个体分配allle
	for(int i=0; i<n_ind;i++){
    Sire_i_pos=Sire[i]-1;		
	Dam_i_pos=Dam[i]-1;		
	
	if(Sire[i]>0){
	int k = rand()%2+1; //从个体父亲的allele中随机抽取一个传给下一代(1表示抽父亲的父亲，2表示抽父亲的母亲)
	allele(i,1)=allele(Sire_i_pos,k);	
	}

	if(Dam[i]>0){
	int k = rand()%2+1;  //从个体母亲的allele中随机抽取一个传给下一代(1表示抽母亲的父亲，2表示抽母亲的母亲)
	allele(i,2)=allele(Dam_i_pos,k);	
	}
	
	}
	
	
	return allele;
}



// convert matrix into three column format
// IND_geno must be integer format for accelerating
// [[Rcpp::export]]
arma::Mat<double> matrix_col3_old(arma::Mat<double> & G,arma::Col<int> IND_geno, bool det, 
                              int cpu_cores, double threshold) { 

	omp_set_num_threads(cpu_cores);
	
	int n_ind=IND_geno.size(), pos,i,j;
	double det_value ;
	int max_row=n_ind*(n_ind+1)/2+1;
	arma::Mat<double> col3_matrix(n_ind*(n_ind+1)/2+1,3);
	col3_matrix.fill(0);
	#pragma omp parallel for private(pos)
	//不需要指定 private
	for( i=0;i<n_ind;i++){	
		for( j=i;j<n_ind;j++){

		pos=max_row-(((n_ind-i)*(n_ind-i-1)/2+n_ind-j-1)+1);		
		col3_matrix(pos,0)=IND_geno[i];	
		col3_matrix(pos,1)=IND_geno[j];			
		col3_matrix(pos,2)=G(i,j);						
		}
	}	
	if(det==true){
	// get real part from complex 
	det_value=log_det(G).real(); 
	col3_matrix(0,2)=det_value;
	}
    
	arma::Col<double> col3=col3_matrix.col(2);
	arma::uvec row_condition=find(abs(col3)>threshold);  //矩阵元素小于0.0001,，可以视作无贡献，可删除
	arma::uvec col={0,1,2};
	col3_matrix=col3_matrix.submat(row_condition,col);
	return col3_matrix;
}


// convert matrix into three column format
// IND_geno must be integer format for accelerating
// [[Rcpp::export]]
SEXP matrix_col3_memory_old(SEXP pBigMat,std::string bigmemory_data_name,std::string bigmemory_data_path,
						arma::Col<int> IND_geno, bool det, int cpu_cores, double threshold) { 

//Rcout<<"Save bigmemory-col_3 matrix as "<<bigmemory_data_path+"/"+bigmemory_data_name+"_col3"+" ......"<<endl;
	Rcpp::XPtr<BigMatrix> pMat(pBigMat); 
	arma::Mat<double> G((double*) pMat->matrix(), pMat -> nrow(), pMat -> ncol(), false); 

	omp_set_num_threads(cpu_cores);

	long n_ind=IND_geno.size(), pos=0,i,j;
	double det_value ;
	
	delete_bigmemory_file_cpp("col3",bigmemory_data_name,bigmemory_data_path,true);	
		
	
	SEXP pBigMat_col3;

	if(det==TRUE){
	
	pBigMat_col3=make_bigmemory_object_address_cpp(n_ind*(n_ind+1)/2+1,3,bigmemory_data_name+"_col3",bigmemory_data_path,"double");			
	Rcpp::XPtr<BigMatrix> pMat_col3(pBigMat_col3);
	arma::Mat<double> col3_matrix((double*) pMat_col3->matrix(), pMat_col3 -> nrow(), pMat_col3 -> ncol(), false);	
	col3_matrix.fill(0);
	int max_row=n_ind*(n_ind+1)/2+1;
	#pragma omp parallel for private(pos)
	//不需要指定 private
	for( i=0;i<n_ind;i++){	
		for( j=i;j<n_ind;j++){
	
		pos=max_row-(((n_ind-i)*(n_ind-i-1)/2+n_ind-j-1)+1);
		col3_matrix(pos,0)=IND_geno[i];	
		col3_matrix(pos,1)=IND_geno[j];			
		col3_matrix(pos,2)=G(i,j);						
		}
	}	
	// get real part from complex 
	det_value=log_det(G).real(); 
	col3_matrix(0,2)=det_value;
	}else{
		
	pBigMat_col3=make_bigmemory_object_address_cpp(n_ind*(n_ind+1)/2,3,bigmemory_data_name+"_col3",bigmemory_data_path,"double");			
	Rcpp::XPtr<BigMatrix> pMat_col3(pBigMat_col3);
	arma::Mat<double> col3_matrix((double*) pMat_col3->matrix(), pMat_col3 -> nrow(), pMat_col3 -> ncol(), false);	
	col3_matrix.fill(0);
	#pragma omp parallel for private(pos)
	//不需要指定 private
	for( i=0;i<n_ind;i++){	
		for( j=i;j<n_ind;j++){

		pos=((n_ind-i)*(n_ind-i-1)/2+n_ind-j-1);
		col3_matrix(pos,0)=IND_geno[i];	
		col3_matrix(pos,1)=IND_geno[j];			
		col3_matrix(pos,2)=G(i,j);						
		}
	}		
		
	}
	return pBigMat_col3;
}



// [[Rcpp::export]]
SEXP matrix_col3_memory_alt(arma::Mat<double> G,std::string bigmemory_data_name,std::string bigmemory_data_path,
						arma::Col<int> IND_geno, bool det, int cpu_cores, double threshold) { 
	omp_set_num_threads(cpu_cores);

	long n_ind=IND_geno.size(), pos=0,i,j;
	double det_value ;
	
	delete_bigmemory_file_cpp("col3",bigmemory_data_name,bigmemory_data_path,true);	
		
	
	SEXP pBigMat_col3;

	if(det==TRUE){
	
	pBigMat_col3=make_bigmemory_object_address_cpp(n_ind*(n_ind+1)/2+1,3,bigmemory_data_name+"_col3",bigmemory_data_path,"double");			
	Rcpp::XPtr<BigMatrix> pMat_col3(pBigMat_col3);
	arma::Mat<double> col3_matrix((double*) pMat_col3->matrix(), pMat_col3 -> nrow(), pMat_col3 -> ncol(), false);	
	col3_matrix.fill(0);
	int max_row=n_ind*(n_ind+1)/2+1;
	#pragma omp parallel for private(pos)
	//不需要指定 private
	for( i=0;i<n_ind;i++){	
		for( j=i;j<n_ind;j++){
	
		pos=max_row-(((n_ind-i)*(n_ind-i-1)/2+n_ind-j-1)+1);
		col3_matrix(pos,0)=IND_geno[i];	
		col3_matrix(pos,1)=IND_geno[j];			
		col3_matrix(pos,2)=G(i,j);						
		}
	}	
	// get real part from complex 
	det_value=log_det(G).real(); 
	col3_matrix(0,2)=det_value;
	}else{
		
	pBigMat_col3=make_bigmemory_object_address_cpp(n_ind*(n_ind+1)/2,3,bigmemory_data_name+"_col3",bigmemory_data_path,"double");			
	Rcpp::XPtr<BigMatrix> pMat_col3(pBigMat_col3);
	arma::Mat<double> col3_matrix((double*) pMat_col3->matrix(), pMat_col3 -> nrow(), pMat_col3 -> ncol(), false);	
	col3_matrix.fill(0);
	#pragma omp parallel for private(pos)
	//不需要指定 private
	for( i=0;i<n_ind;i++){	
		for( j=i;j<n_ind;j++){

		pos=((n_ind-i)*(n_ind-i-1)/2+n_ind-j-1);
		col3_matrix(pos,0)=IND_geno[i];	
		col3_matrix(pos,1)=IND_geno[j];			
		col3_matrix(pos,2)=G(i,j);						
		}
	}		
		
	}
	return pBigMat_col3;
}



// [[Rcpp::export]]
arma::Mat<double> matrix_col3(arma::Mat<double> & G,arma::Col<int> IND_geno,bool det=false,int cpu_cores=1,double threshold=0) { 

	Rcout<<"Output col3  matrix type......"<<endl;
	omp_set_num_threads(cpu_cores);	
	uvec upper_indices = trimatl_ind(size(G));
	vec upper_part = G(upper_indices);
	long ind_n=G.n_rows;
	long max_row=ind_n*(ind_n+1)/2,ele_size=upper_part.size();

	vec ind_row(ele_size),ind_col(ele_size);
	arma::uvec row_condition=find(abs(upper_part)>threshold);  //矩阵元素小于0.0001,，可以视作无贡献，可删除
	arma::Mat<double> col3_matrix;

	col3_matrix=arma::mat(row_condition.size(),3,fill::zeros);
	long i,j,pos;
	#pragma omp parallel for private(i,j,pos)
	////不需要指定 private
	for(i=0;i<ind_n;i++){	
		for( j=i;j<ind_n;j++){		
			pos=max_row-(((ind_n-i)*(ind_n-i-1)/2+ind_n-j-1)+1);
			ind_row[pos]=IND_geno[i];
			ind_col[pos]=IND_geno[j];
		}
	}	
	col3_matrix.col(0)=ind_row(row_condition);
	col3_matrix.col(1)=ind_col(row_condition);
	col3_matrix.col(2)=upper_part(row_condition);

	if(det==true){
	double det_value=log_det(G).real(); 
	arma::mat det_mat(1,3,fill::zeros);
	det_mat(0,2)=det_value;
	col3_matrix=join_cols(det_mat,col3_matrix);
	}
	return col3_matrix;
}



// convert matrix into three column format
// IND_geno must be integer format for accelerating
// [[Rcpp::export]]
SEXP matrix_col3_memory(SEXP pBigMat,std::string bigmemory_data_name,std::string bigmemory_data_path,
						arma::Col<int> IND_geno, bool det, int cpu_cores, double threshold) { 

	omp_set_num_threads(cpu_cores);	
	Rcpp::XPtr<BigMatrix> pMat(pBigMat); 
	arma::Mat<double> G((double*) pMat->matrix(), pMat -> nrow(), pMat -> ncol(), false); 

	uvec upper_indices = trimatl_ind(size(G));
	vec upper_part = G(upper_indices);
	long ind_n=G.n_rows;
	long max_row=ind_n*(ind_n+1)/2,ele_size=upper_part.size();

	vec ind_row(ele_size),ind_col(ele_size);
	arma::uvec row_condition=find(abs(upper_part)>threshold);  //矩阵元素小于0.0001,，可以视作无贡献，可删除
	arma::Mat<double> col3_matrix;

	col3_matrix=arma::mat(row_condition.size(),3,fill::zeros);
	long i,j,pos;
	#pragma omp parallel for private(i,j,pos)
	////不需要指定 private
	for(i=0;i<ind_n;i++){	
		for( j=i;j<ind_n;j++){		
			pos=max_row-(((ind_n-i)*(ind_n-i-1)/2+ind_n-j-1)+1);
			ind_row[pos]=IND_geno[i];
			ind_col[pos]=IND_geno[j];
		}
	}	
	col3_matrix.col(0)=ind_row(row_condition);
	col3_matrix.col(1)=ind_col(row_condition);
	col3_matrix.col(2)=upper_part(row_condition);

	if(det==true){
	double det_value=log_det(G).real(); 
	arma::mat det_mat(1,3,fill::zeros);
	det_mat(0,2)=det_value;
	col3_matrix=join_cols(det_mat,col3_matrix);
	}

	
	delete_bigmemory_file_cpp("col3",bigmemory_data_name,bigmemory_data_path,true);	
			
	SEXP pBigMat_col3;	
	pBigMat_col3=make_bigmemory_object_address_cpp(col3_matrix.n_rows,3,bigmemory_data_name+"_col3",bigmemory_data_path,"double");			
	Rcpp::XPtr<BigMatrix> pMat_col3(pBigMat_col3);
	arma::Mat<double> col3_matrix_memory((double*) pMat_col3->matrix(), pMat_col3 -> nrow(), pMat_col3 -> ncol(), false);
	col3_matrix_memory=col3_matrix;
	
	
	SEXP pBigMat_col3_bigmemory_object=make_bigmemory_object_cpp(col3_matrix.n_rows,3,bigmemory_data_name+"_col3",bigmemory_data_path,"double");	
	
	return pBigMat_col3_bigmemory_object;
}
