#include "shared_function.h"

// [[Rcpp::export]]
arma::Mat<int> phased_haplotype_to_numeric_cpp(std::vector<int> block_start,
											   std::vector<int> block_end, 
											   arma::Mat<int> &data_hap,
											   Rcpp::List &haplotype_allele,
											   int cpu_cores=1){  //输入数据，行为SNP, 列为个体		
			omp_set_num_threads(cpu_cores);
			Rcout<<"Please make sure haplotype data has been phased!"<<endl;	
			Rcout<<"Start phased-haplotype-numeric(0,1,2) format conversion...... "<<endl;	
			int window_n=block_start.size();
			int i,j,m,n,tmp;	
			int n_ind=(data_hap.n_cols)/2;
			std::vector<int> haplo_pos,haplo_type_num(window_n),cumsum_haplo_type_num;
			arma::Mat<int> data_haplotype_subset;
			arma::Mat<int> haplo_pos_list3(window_n,n_ind*2);

			// get the overview set of haplo_data
			for(i=0; i < window_n; i++){
				data_haplotype_subset=data_hap.rows(block_start[i],block_end[i]);
				allele_get_haplotype_set_short(data_haplotype_subset,haplotype_allele,i);
            }
			bool status_progress=true;
			Progress p(window_n*n_ind,status_progress);					
			#pragma omp parallel for private(i,j,data_haplotype_subset,haplo_pos)
			for(i=0; i < window_n; i++){
				data_haplotype_subset=data_hap.rows(block_start[i],block_end[i]);
				haplo_pos=get_haplotype_set_short(data_haplotype_subset);
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

			arma::Mat<int> data_numeric(n_ind,sum_haplo_type_num);	
			data_numeric.fill(0);
			#pragma omp parallel for private(m,n,tmp,ind_pos1,ind_pos2)
			for(m=0;m<window_n;m++){
				
				tmp=cumsum_haplo_type_num[m];
				
				for(n=0;n<n_ind;n++){
					
					ind_pos1=haplo_pos_list3(m,2*n)+tmp;
				    ind_pos2=haplo_pos_list3(m,2*n+1)+tmp;
				    data_numeric(n,ind_pos1)=data_numeric(n,ind_pos1)+1;
				    data_numeric(n,ind_pos2)=data_numeric(n,ind_pos2)+1;
				}
			}
    Rcout<<" "<<endl;			
	Rcout<<"Complete phased-haplotype-numeric(0,1,2) format conversion!"<<endl;	
			return data_numeric;
			}

// [[Rcpp::export]]
std::vector<std::string> phased_haplotype_to_blupf90_cpp(std::vector<int> block_start,std::vector<int> block_end, arma::Mat<int> &data_hap,
														 CharacterVector IND_name,
														 Rcpp::List &haplotype_allele,
														 int cpu_cores=1){  //输入数据，行为SNP, 列为个体		
			omp_set_num_threads(cpu_cores);
			Rcout<<"Please make sure haplotype data has been phased!"<<endl;	
			Rcout<<"Start phased-haplotype-BLUPF90 format conversion...... "<<endl;	
			int window_n=block_start.size();
			int i,j,k,m,n;	
			int n_ind=(data_hap.n_cols)/2;
			std::vector<int> haplo_pos,haplo_type_num(window_n),cumsum_haplo_type_num;
			arma::Mat<int> data_haplotype_subset;
			arma::Mat<int> haplo_pos_list3(window_n,n_ind*2);
			// get the overview set of haplo_data
			for(i=0; i < window_n; i++){
				data_haplotype_subset=data_hap.rows(block_start[i],block_end[i]);
				allele_get_haplotype_set_short(data_haplotype_subset,haplotype_allele,i);
            }
			bool status_progress=true;
			Progress p(window_n*n_ind,status_progress);					
			#pragma omp parallel for private(i,j,data_haplotype_subset,haplo_pos)
			for(i=0; i < window_n; i++){
				data_haplotype_subset=data_hap.rows(block_start[i],block_end[i]);
				haplo_pos=get_haplotype_set_short(data_haplotype_subset);
				haplo_type_num[i]=*max_element(haplo_pos.begin(),haplo_pos.end())+1;
				
				for(j=0;j<n_ind;j++){
					p.increment();
					haplo_pos_list3(i,2*j)=haplo_pos[2*j];
					haplo_pos_list3(i,2*j+1)=haplo_pos[2*j+1];
					
				}
			}
	int sum_haplo_type_num=std::accumulate(haplo_type_num.begin(),haplo_type_num.end(),0); //number of total snps(columns of data_numeric)
	int ind_pos1,ind_pos2;  // all types of haplot
	cumulativeSum(haplo_type_num,cumsum_haplo_type_num);
	
	int n_snp=sum_haplo_type_num,tmp;
	std::string ind_name;
	std::vector<int> ind_name_size(IND_name.size());
	//确定个体名称中最大字符串长度
	for(i=0;i<IND_name.size();i++){
		ind_name=IND_name[i];
		ind_name_size[i]=ind_name.length();
	}
	int max_length=*max_element(ind_name_size.begin(),ind_name_size.end())+3;	
	max_length=0;
	std::string allele1,allele2;	
	std::string allele,allele_string(n_snp+max_length,' ');
	std::vector<std::string> allele_string_vector(n_ind);
	arma::Col<int> allele_string_col1,allele_string_col2;
	#pragma omp parallel for private(j,allele_string_col1,allele_string_col2,allele_string)
	for(j=0;j<n_ind;j++){
	
		allele_string_col1=haplo_pos_list3.col(2*j);
		allele_string_col2=haplo_pos_list3.col(2*j+1);
		
		allele_string=get_blupf90_allele_string_phased_vcf(allele_string_col1,allele_string_col2,max_length,cumsum_haplo_type_num,n_snp);
	
		allele_string_vector[j]=allele_string;	
	}
	
	//添加个体名称
	//for(j=0;j<n_ind;j++){		
	//	ind_name=IND_name[j];	
	//	for(k=0;k<ind_name.size();k++){
	//		allele_string_vector[j][k]=ind_name[k];
	//	}
	//}
	
    Rcout<<" "<<endl;			
	Rcout<<"Complete phased-haplotype-BLUPF90 format conversion!"<<endl;	
			return allele_string_vector;
}



// [[Rcpp::export]]
arma::Mat<int> unphased_haplotype_to_numeric_cpp(arma::Mat<int> &data_hap,int cpu_cores=1){  //输入数据，行为SNP, 列为个体		
			omp_set_num_threads(cpu_cores);
			Rcout<<"Please make sure haplotype data has been phased!"<<endl;	
			Rcout<<"Start unphased-haplotype-numeric(0,1,2) format conversion...... "<<endl;	
			
			int i,j,allele1,allele2;	
			int n_ind=data_hap.n_cols/2,n_snp=data_hap.n_rows;
			arma::Mat<int> data_numeric(n_ind,n_snp);	
			data_numeric.fill(0);
			#pragma omp parallel for private(i,j,allele1,allele2)
			for(i=0;i<n_snp;i++){
				for(j=0;j<n_ind;j++){					
					allele1=data_hap(i,2*j);
					allele2=data_hap(i,2*j+1);
				    data_numeric(j,i)=allele1+allele2;
				}
			}
    Rcout<<" "<<endl;			
	Rcout<<"Complete phased-haplotype-numeric(0,1,2) format conversion!"<<endl;	
			return data_numeric;
			}
			



// [[Rcpp::export]]
std::vector<std::string> unphased_haplotype_to_blupf90_cpp(arma::Mat<int> & data_hap,CharacterVector IND_name, int cpu_cores=1){  //输入数据，行为SNP, 列为个体		
			omp_set_num_threads(cpu_cores);
			Rcout<<"Start unphased-haplotype-BLUPF90 format conversion......"<<endl;	
			
			int i,j,k;	
			int n_ind=data_hap.n_cols/2,n_snp=data_hap.n_rows;
	        std::string ind_name;
	        std::vector<int> ind_name_size(IND_name.size());
	        //确定个体名称中最大字符串长度
	        for(i=0;i<IND_name.size();i++){
	        	ind_name=IND_name[i];
	        	ind_name_size[i]=ind_name.length();
	        }
	        int max_length=*max_element(ind_name_size.begin(),ind_name_size.end())+3;	
	        max_length=0;
	        std::string allele,allele_string(n_snp+max_length,' ');
	        std::vector<std::string> allele_string_vector(n_ind);
	        arma::Col<int> allele_string_col1,allele_string_col2;
	        #pragma omp parallel for private(i,allele_string_col1,allele_string_col2,allele_string)
			for(i=0;i<n_ind;i++){	
			
				allele_string_col1=data_hap.col(2*i);
				allele_string_col2=data_hap.col(2*i+1);
				
	        	allele_string=get_blupf90_allele_string_unphased_haplotype(allele_string_col1,allele_string_col2,max_length);
	        
	        	allele_string_vector[i]=allele_string;								
					
				}
	        
	        //添加个体名称
	        //for(j=0;j<n_ind;j++){		
	        //	ind_name=IND_name[j];	
	        //	for(k=0;k<ind_name.size();k++){
	        //		allele_string_vector[j][k]=ind_name[k];
	        //	}
	        //}			
			
    Rcout<<" "<<endl;			
	Rcout<<"Complete unphased-haplotype-BLUPF90 format conversion!"<<endl;	
			return allele_string_vector;
			}			
			

// type 1: haplotype; type 2: numeric ;
// [[Rcpp::export]]
List  haplotype_convertion(
					arma::Mat<int> & data_hap, 	
					std::vector<int> block_start,
					std::vector<int> block_end,
					CharacterVector IND_name,
					int type=1, 	
					int cpu_cores=1){
		
		arma::Mat<int> data_numeric;
		std::vector<std::string> data_blupf90;
		int window_n=block_start.size();
		Function generate_list("vector");
		Rcpp::List haplotype_allele=generate_list("list", window_n);
		
		switch(type){

		case 1:{        //type:  phased haplotype to numeric 
		
		data_numeric=phased_haplotype_to_numeric_cpp(block_start,block_end,data_hap,haplotype_allele,cpu_cores);
			
			break;}	
			
		case 2:{        //type:  phased haplotype to blupf90 
		
		data_blupf90=phased_haplotype_to_blupf90_cpp(block_start,block_end,data_hap,IND_name,haplotype_allele,cpu_cores);
		
			break;}	
		
		case 3:{       //type:  unphased haplotype to numeric 
		
		data_numeric=unphased_haplotype_to_numeric_cpp(data_hap,cpu_cores);
		
			break;}
			
		case 4:{        //type:unphased haplotype to blupf90 

		data_blupf90=unphased_haplotype_to_blupf90_cpp(data_hap,IND_name,cpu_cores);
			
			break;}

		default:
		
        throw Rcpp::exception("unknown input data type!");			
			
		}
		
		return List::create(Named("numeric") = data_numeric,
							Named("blupf90") = data_numeric,
							Named("SNP")=haplotype_allele);
}			
			