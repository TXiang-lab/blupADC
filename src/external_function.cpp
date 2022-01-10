#include "shared_function.h"

// [[Rcpp::export]]
LogicalVector judge_character_string(CharacterVector left, String b){
int n_ind=left.size();
LogicalVector result(n_ind);
for(int i=0;i<n_ind;i++){
result[i]=(left[i]==b);
}
return result;
}


// get full_generation pedigree
// [[Rcpp::export]]
CharacterMatrix full_generation_conversion(CharacterVector generation_names, CharacterMatrix ped){ //ped need renamed

int n_ind=ped.nrow(),Col_n=generation_names.size();
CharacterVector IND_pedigree=ped.column(0);
CharacterVector id,sire,dam;

CharacterMatrix full_generation(n_ind,Col_n+1);
full_generation.fill("0");

full_generation.column(0)=ped.column(0);
full_generation.column(1)=ped.column(1);
full_generation.column(2)=ped.column(2);

Function paste0_cpp("paste0");

for(int i=0;i<n_ind;i++){

id=full_generation(i,0),sire=full_generation(i,1),dam=full_generation(i,2);

IntegerVector pos_sire_vector=match(sire,IND_pedigree);
IntegerVector pos_dam_vector=match(dam,IND_pedigree);


// Father 
if(!is_true(all(is_na(pos_sire_vector)))){ //父亲存在于第一列

int pos_sire=pos_sire_vector[0]-1;

CharacterVector record=full_generation.row(pos_sire);
record.erase(0);

LogicalVector logic=judge_character_string(record,"0");

if(!(is_true(all(logic)))){  //父亲的祖先不都是缺失的

CharacterVector record_names=paste0_cpp("Sire",generation_names[!logic]);
CharacterVector record_non_NA=record[!logic];
IntegerVector pos=match(record_names,generation_names);
IntegerVector pos1=pos[!is_na(pos)];
CharacterVector record_non_NA_final=record_non_NA[!is_na(pos)];

for(int j=0;j<record_non_NA_final.size();j++){
	int k=pos1[j];
	full_generation(i,k)=record_non_NA_final[j];
}

}

}

// Father 
if(!is_true(all(is_na(pos_dam_vector)))){ //父亲存在于第一列

int pos_dam=pos_dam_vector[0]-1;

CharacterVector record=full_generation.row(pos_dam);
record.erase(0);

LogicalVector logic=judge_character_string(record,"0");

if(!(is_true(all(logic)))){  //父亲的祖先不都是缺失的

CharacterVector record_names=paste0_cpp("Dam",generation_names[!logic]);
CharacterVector record_non_NA=record[!logic];
IntegerVector pos=match(record_names,generation_names);
IntegerVector pos1=pos[!is_na(pos)];
CharacterVector record_non_NA_final=record_non_NA[!is_na(pos)];

for(int j=0;j<record_non_NA_final.size();j++){
	int k=pos1[j];
	full_generation(i,k)=record_non_NA_final[j];
}

}

}

}

return full_generation;
}


// [[Rcpp::export]]
arma::mat  numeric_overlap_cpp(arma::mat numeric,double overlap_threshold=0.9,int cpu_cores=1, bool dis_progress=true){
//对两个for循环使用 openMP报错，就使用一个for循环，问题也能解决，虽然可能会慢	
			//omp_set_num_threads(cpu_cores);
	int n_ind=numeric.n_rows, SNP_n=numeric.n_cols, iter=-1,i,j,j_sorce;
	double zero_vec_number,standard_error=SNP_n*(1-overlap_threshold);
	arma::mat all_one,result_matrix(n_ind*(n_ind+1)/2,3);
	result_matrix.fill(0);
	arma::uvec zero_vec;	
	all_one.ones(SNP_n,1);
	arma::mat score=numeric*all_one;
	rowvec IND1,IND2,result;	
    Progress p(n_ind*(n_ind+1)/2,dis_progress);	
	for( i=0;i<(n_ind-1);i++){	
	IND1=numeric.row(i);	
	int i_sorce=score(i,0);	
			
//#pragma omp parallel for private(j,j_sorce,IND2,result,zero_vec,zero_vec_number,iter)			
		for( j=i+1;j<n_ind;j++){	
		if ( ! Progress::check_abort()) {
		p.increment(); 	
		}
		 j_sorce=score(j,0);	
        if(j_sorce>i_sorce-standard_error && j_sorce<i_sorce+standard_error){	
		iter=iter+1;			
		IND2=numeric.row(j);
		result=IND1-IND2;
		zero_vec = find(result >= 0);
		zero_vec_number=zero_vec.size();	
        if(zero_vec_number/SNP_n>overlap_threshold){	
		result_matrix(iter,0)=i+1;
		result_matrix(iter,1)=j+1;
		result_matrix(iter,2)=zero_vec_number/SNP_n;		
		} 
		}   
	 }
	 }	 
	
	return result_matrix;
}



//generate social genetic effect phenotype
// [[Rcpp::export]]
List get_rest_id(List ref, NumericVector ref_index, NumericMatrix phe, int column_n, int group_column , int id_column) {	
	
	NumericVector group_data=phe.column(group_column);  //group的类别
	NumericVector id_data=phe.column(id_column);  //group的类别
	
	NumericVector target_id(1);
	NumericVector group_i(1);

	int ind_size=phe.nrow();   // 个体数

	//创建整型列的group矩阵
	NumericMatrix group_integer(ind_size,column_n);
	std::fill(group_integer.begin(), group_integer.end(),999); //初始化相同的值

	//创建实型列的group矩阵
	IntegerMatrix group_real(ind_size,column_n);	
	
	for(int i=0; i < ind_size; i++){
		
	group_i[0]=group_data[i];   //当前行的group信息
	target_id[0]=id_data[i];    //当前行的个体信息
	
	IntegerVector pos=match(group_i,ref_index)-1;;  //group等于该列的数据
	NumericVector group_id_set=ref[pos[0]];    //该group下所有的个体
	NumericVector rest_id_set=setdiff(group_id_set,target_id); //剩余的个体

	int rest_n=rest_id_set.size();
	
	IntegerVector rest_id_temp=rep(1,rest_n);
	NumericVector rest_id_indicate=as<NumericVector>(rest_id_temp);
	
	for(int j=0; j < rest_n; j++){
	
	group_integer(i,j)=rest_id_set[j];
	group_real(i,j)=1;	
	}
	}
	return List::create(Named("integer_group") = group_integer,_["real_group"] = group_real);
	
}


