#include "shared_function.h"


// [[Rcpp::export]]
std::vector<std::string> numeric_to_blupf90_cpp(arma::Mat<int> &data_numeric,CharacterVector IND_name, int cpu_cores=1){  //输入数据，行为SNP, 列为个体		
			omp_set_num_threads(cpu_cores);
			Rcout<<"Start Numeric to BLUPF90 format conversion......"<<endl;	
			int i,j,k;	
			int n_ind=data_numeric.n_rows,n_snp=data_numeric.n_cols;
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
	        arma::Row<int> allele_string_row;
	        #pragma omp parallel for private(i,allele_string_row,allele_string)
			for(i=0;i<n_ind;i++){	
			
				allele_string_row=data_numeric.row(i);
				
	        	allele_string=get_blupf90_allele_string_numeric(allele_string_row,max_length);
	        
	        	allele_string_vector[i]=allele_string;								
					
				}
		
			
    Rcout<<" "<<endl;			
	Rcout<<"Complete Numeric to BLUPF90 format conversion!"<<endl;	
			return allele_string_vector;
			}	


// [[Rcpp::export]]
CharacterMatrix numeric_to_hapmap_cpp(CharacterVector IND_name, 
									  CharacterMatrix &data_numeric_map,
									   arma::Mat<int> &data_numeric,
										int cpu_cores=1,
										std::string miss_base="N"){  //输入数据，
			omp_set_num_threads(cpu_cores);
			Rcout<<"Start Numeric to Hapmap format conversion......"<<endl;	
			int i,j,k;	
			int n_ind=data_numeric.n_rows,n_snp=data_numeric.n_cols;
			CharacterVector tmp1=data_numeric_map.column(3),tmp2=data_numeric_map.column(4);
			std::vector<std::string> Ref=Rcpp::as<std::vector<std::string> >(tmp1);
			std::vector<std::string> Alt=Rcpp::as<std::vector<std::string> >(tmp2);	
			std::string allele_aa,allele_Aa,allele_AA,allele_NN=miss_base+miss_base;
			CharacterMatrix data_hmp(n_snp,n_ind+11);
			int allele_numeric;			
	        #pragma omp parallel for private(i,j,allele_numeric,allele_aa,allele_Aa,allele_AA)
			for(i=0;i<n_ind;i++){	
	
				
				for(j=0;j<n_snp;j++){
					
					allele_numeric=data_numeric(i,j);
					
					if(allele_numeric==0){
						allele_aa=Ref[j]+Ref[j];
						data_hmp(j,i+11)=allele_aa;
					}else if(allele_numeric==1){
						allele_Aa=Ref[j]+Alt[j];
						data_hmp(j,i+11)=allele_Aa;
					}else if(allele_numeric==2){
						allele_AA=Alt[j]+Alt[j];	
						data_hmp(j,i+11)=allele_AA;	
					}else{	
						data_hmp(j,i+11)=allele_NN;							
					}		
					
				}
			}	

			
    Rcout<<" "<<endl;			
	Rcout<<"Complete Numeric to Hapmap format conversion!"<<endl;	
	data_hmp.column(0)=data_numeric_map.column(1);
	data_hmp.column(2)=data_numeric_map.column(0);
	data_hmp.column(3)=data_numeric_map.column(2);
	CharacterVector title={"rs#","alleles","chrom","pos","strand","center","protLSID","assembly","assayLSID","panelLSID","QCcode"};
	Function C_cpp("c");
	CharacterVector title_IND_name=C_cpp(title,IND_name);
	colnames(data_hmp)=title_IND_name;
	return data_hmp;
}


// [[Rcpp::export]]
CharacterMatrix numeric_to_ped_cpp(CharacterVector IND_name,
									CharacterMatrix &data_numeric_map,
									arma::Mat<int> &data_numeric,
									int cpu_cores=1,
									std::string miss_base="N"){  //输入数据，
			omp_set_num_threads(cpu_cores);
			Rcout<<"Start Numeric to Plink format conversion......"<<endl;	
			int i,j,k;	
			int n_ind=data_numeric.n_rows,n_snp=data_numeric.n_cols;
	        std::string ind_name;
	        std::vector<int> ind_name_size(IND_name.size());
			CharacterVector tmp1=data_numeric_map.column(3),tmp2=data_numeric_map.column(4);
			std::vector<std::string> Ref=Rcpp::as<std::vector<std::string> >(tmp1);
			std::vector<std::string> Alt=Rcpp::as<std::vector<std::string> >(tmp2);	
			CharacterMatrix data_ped(n_ind,n_snp*2+6);
			int allele_numeric;
			#pragma omp parallel for private(i,j,allele_numeric)
			for(i=0;i<n_ind;i++){
				
				for(j=0;j<n_snp;j++){
				
					allele_numeric=data_numeric(i,j);
					
					if(allele_numeric==0){
						data_ped(i,6+2*j)=Ref[j];
						data_ped(i,6+2*j+1)=Ref[j];
					}else if(allele_numeric==1){
						data_ped(i,6+2*j)=Ref[j];
						data_ped(i,6+2*j+1)=Alt[j];
					}else if(allele_numeric==2){
						data_ped(i,6+2*j)=Alt[j];
						data_ped(i,6+2*j+1)=Alt[j];	
					}else{	
						data_ped(i,6+2*j)=miss_base;
						data_ped(i,6+2*j+1)=miss_base;						
					}
					
				}
			}		
			
	Rcout<<" "<<endl;			
	Rcout<<"Complete BLUPF90 to Plink data format conversion!"<<endl;			
	data_ped.column(0)=IND_name;
	data_ped.column(1)=IND_name;

CharacterVector zero1={"0"};
CharacterVector zero=rep(zero1,n_snp);

	data_ped.column(2)=zero;
	data_ped.column(3)=zero;
	data_ped.column(4)=zero;
	data_ped.column(5)=zero;	
	
	return  data_ped ;	
}


// [[Rcpp::export]]
CharacterMatrix numeric_to_vcf_cpp(CharacterVector IND_name,
									CharacterMatrix &data_numeric_map,
									arma::Mat<int> &data_numeric,
									std::string phased_symbol="/",
									int cpu_cores=1){  //输入数据，
			omp_set_num_threads(cpu_cores);
			Rcout<<"Start Numeric to VCF format conversion......"<<endl;	
			int i,j,k;	
			int n_ind=data_numeric.n_rows,n_snp=data_numeric.n_cols;
			CharacterMatrix data_vcf(n_snp,n_ind+9);
			CharacterVector tmp1=data_numeric_map.column(3),tmp2=data_numeric_map.column(4);
			data_vcf.column(3)=tmp1;
			data_vcf.column(4)=tmp2;
			std::vector<std::string> Ref=Rcpp::as<std::vector<std::string> >(tmp1);
			std::vector<std::string> Alt=Rcpp::as<std::vector<std::string> >(tmp2);	
			std::string allele_aa,allele_Aa,allele_AA;

			int allele_numeric;	
			std::string allele_aa_num="0"+phased_symbol+"0";
			std::string allele_Aa_num="1"+phased_symbol+"0";
			std::string allele_AA_num="1"+phased_symbol+"1";
			std::string allele_NN_num="."+phased_symbol+".";
			
			#pragma omp parallel for private(i,j,allele_numeric)
			for(i=0;i<n_ind;i++){
				
				for(j=0;j<n_snp;j++){
				
					allele_numeric=data_numeric(i,j);
					
					if(allele_numeric==0){
						data_vcf(j,i+9)=allele_aa_num;
					}else if(allele_numeric==1){
						data_vcf(j,i+9)=allele_Aa_num;
					}else if(allele_numeric==2){
						data_vcf(j,i+9)=allele_AA_num;	
					}else{
						data_vcf(j,i+9)=allele_NN_num;	
					}
					
				}
			}		
			
Rcout<<" "<<endl;	
Rcout<<"Complete Numeric to VCF data format conversion!"<<endl;
CharacterVector chr=data_numeric_map.column(0),snp=data_numeric_map.column(1),pos=data_numeric_map.column(2);
//CharacterVector qual=rep("1",n_snp),filter=rep("PASS",n_snp),info=rep("2",n_snp),format=rep("GT",n_snp);
CharacterVector qual1={"."},filter1={"PASS"},info1={"."},format1={"GT"};
CharacterVector qual=rep(qual1,n_snp),filter=rep(filter1,n_snp),info=rep(info1,n_snp),format=rep(format1,n_snp);
data_vcf.column(0)=chr;
data_vcf.column(1)=pos;
data_vcf.column(2)=snp;
data_vcf.column(5)=qual;
data_vcf.column(6)=filter;
data_vcf.column(7)=info;
data_vcf.column(8)=format;
	
IND_name.push_front("FORMAT");
IND_name.push_front("INFO");
IND_name.push_front("FILTER");
IND_name.push_front("QUAL");
IND_name.push_front("ALT");
IND_name.push_front("REF");
IND_name.push_front("ID");
IND_name.push_front("POS");
IND_name.push_front("#CHROM");
colnames(data_vcf)=IND_name;
    return data_vcf;
}

// [[Rcpp::export]]
arma::Mat<int> blupf90_to_numeric_cpp(std::vector<std::string> & data_blupf90,int cpu_cores=1){  //输入数据，行为SNP, 列为个体		
			
			omp_set_num_threads(cpu_cores);
			Rcout<<"Start BLUPF90 to Numeric format conversion......"<<endl;	
			int i,j,k;	
			std::string tmp_ind_string,allele_string;
			int n_ind=data_blupf90.size(),n_snp=data_blupf90[0].size();
			arma::Mat<int> data_numeric(n_ind,n_snp);			
			
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
			return data_numeric;
}



// [[Rcpp::export]]
CharacterMatrix blupf90_to_hapmap_cpp(CharacterVector IND_name,
									  CharacterMatrix &data_blupf90_map,
                                      std::vector<std::string> & data_blupf90,
									  int cpu_cores=1,
									  std::string miss_base="N"){  //输入数据，	
			
			omp_set_num_threads(cpu_cores);
			Rcout<<"Start BLUPF90 to Hapmap format conversion......"<<endl;	
			int i,j,k;	
			std::string tmp_ind_string,allele_string;
			int n_ind=data_blupf90.size(),n_snp=data_blupf90[0].size();
			arma::Mat<int> data_numeric(n_ind,n_snp);	
			CharacterVector tmp1=data_blupf90_map.column(3),tmp2=data_blupf90_map.column(4);
			std::vector<std::string> Ref=Rcpp::as<std::vector<std::string> >(tmp1);
			std::vector<std::string> Alt=Rcpp::as<std::vector<std::string> >(tmp2);
			std::string allele_aa,allele_Aa,allele_AA,allele_NN=miss_base+miss_base;
			CharacterMatrix data_hmp(n_snp,n_ind+11);
			int allele_numeric;
			#pragma omp parallel for private(i,j,tmp_ind_string,allele_string,allele_numeric,allele_aa,allele_Aa,allele_AA)
			for(i=0;i<n_ind;i++){
				
				tmp_ind_string=data_blupf90[i];
				
				for(j=0;j<n_snp;j++){
				
					allele_string=tmp_ind_string[j];
					allele_numeric=std::stoi(allele_string);				
					if(allele_numeric==0){
						allele_aa=Ref[j]+Ref[j];
						data_hmp(j,i+11)=allele_aa;
					}else if(allele_numeric==1){
						allele_Aa=Ref[j]+Alt[j];
						data_hmp(j,i+11)=allele_Aa;
					}else if(allele_numeric==2){
						allele_AA=Alt[j]+Alt[j];	
						data_hmp(j,i+11)=allele_AA;	
					}else{
						data_hmp(j,i+11)=allele_NN;
					}
					
				}
			}		
			
	Rcout<<" "<<endl;			
	Rcout<<"Complete BLUPF90 to Hapmap format conversion!"<<endl;
			
	data_hmp.column(0)=data_blupf90_map.column(1);
	data_hmp.column(2)=data_blupf90_map.column(0);
	data_hmp.column(3)=data_blupf90_map.column(2);
	CharacterVector title={"rs#","alleles","chrom","pos","strand","center","protLSID","assembly","assayLSID","panelLSID","QCcode"};
	Function C_cpp("c");
	CharacterVector title_IND_name=C_cpp(title,IND_name);
	colnames(data_hmp)=title_IND_name;
	return data_hmp;			
}


// [[Rcpp::export]]
CharacterMatrix blupf90_to_ped_cpp(CharacterVector IND_name,
								   CharacterMatrix &data_blupf90_map,
                                   std::vector<std::string> & data_blupf90,
								   int cpu_cores=1,
								   std::string miss_base="N"){  
			omp_set_num_threads(cpu_cores);
			Rcout<<"Start BLUPF90 to Plink format conversion......"<<endl;	
			int i,j,k;	
			std::string tmp_ind_string,allele_string;
			int n_ind=data_blupf90.size(),n_snp=data_blupf90[0].size();
			arma::Mat<int> data_numeric(n_ind,n_snp);			
			CharacterVector tmp1=data_blupf90_map.column(3),tmp2=data_blupf90_map.column(4);
			std::vector<std::string> Ref=Rcpp::as<std::vector<std::string> >(tmp1);
			std::vector<std::string> Alt=Rcpp::as<std::vector<std::string> >(tmp2);
			CharacterMatrix data_ped(n_ind,n_snp*2+6);
			int allele_numeric;
			#pragma omp parallel for private(i,j,tmp_ind_string,allele_string,allele_numeric)
			for(i=0;i<n_ind;i++){
				
				tmp_ind_string=data_blupf90[i];
				
				for(j=0;j<n_snp;j++){
				
					allele_string=tmp_ind_string[j];
					allele_numeric=std::stoi(allele_string);				
					if(allele_numeric==0){
						data_ped(i,6+2*j)=Ref[j];
						data_ped(i,6+2*j+1)=Ref[j];
					}else if(allele_numeric==1){
						data_ped(i,6+2*j)=Ref[j];
						data_ped(i,6+2*j+1)=Alt[j];
					}else if(allele_numeric==2){
						data_ped(i,6+2*j)=Alt[j];
						data_ped(i,6+2*j+1)=Alt[j];	
					}else{
						data_ped(i,6+2*j)=miss_base;
						data_ped(i,6+2*j+1)=miss_base;		
					}
					
				}
			}		
			
	Rcout<<" "<<endl;			
	Rcout<<"Complete BLUPF90 to Plink data format conversion!"<<endl;			
	data_ped.column(0)=IND_name;
	data_ped.column(1)=IND_name;

CharacterVector zero1={"0"};
CharacterVector zero=rep(zero1,n_snp);

	data_ped.column(2)=zero;
	data_ped.column(3)=zero;
	data_ped.column(4)=zero;
	data_ped.column(5)=zero;	
	
	return  data_ped ;			
}



// [[Rcpp::export]]
CharacterMatrix blupf90_to_vcf_cpp(CharacterVector IND_name,
								   CharacterMatrix &data_blupf90_map,
                                   std::vector<std::string> & data_blupf90,
								   std::string phased_symbol="/",
								   int cpu_cores=1){  //输入数据
			omp_set_num_threads(cpu_cores);
			Rcout<<"Start BLUPF90 to VCF format conversion......"<<endl;	
			int i,j,k;	
			std::string tmp_ind_string,allele_string;
			int n_ind=data_blupf90.size(),n_snp=data_blupf90[0].size();	
			CharacterMatrix data_vcf(n_snp,n_ind+9);			
			CharacterVector tmp1=data_blupf90_map.column(3),tmp2=data_blupf90_map.column(4);
			std::vector<std::string> Ref=Rcpp::as<std::vector<std::string> >(tmp1);
			std::vector<std::string> Alt=Rcpp::as<std::vector<std::string> >(tmp2);
			data_vcf.column(3)=tmp1;
			data_vcf.column(4)=tmp2;

			int allele_numeric;	
			std::string allele_aa_num="0"+phased_symbol+"0";
			std::string allele_Aa_num="1"+phased_symbol+"0";
			std::string allele_AA_num="1"+phased_symbol+"1";
			std::string allele_NN_num="."+phased_symbol+".";
			
			#pragma omp parallel for private(i,j,tmp_ind_string,allele_string,allele_numeric)
			for(i=0;i<n_ind;i++){
				
				tmp_ind_string=data_blupf90[i];
				
				for(j=0;j<n_snp;j++){
				
					allele_string=tmp_ind_string[j];
					allele_numeric=std::stoi(allele_string);	
					
					if(allele_numeric==0){
						data_vcf(j,i+9)=allele_aa_num;
					}else if(allele_numeric==1){
						data_vcf(j,i+9)=allele_Aa_num;
					}else if(allele_numeric==2){
						data_vcf(j,i+9)=allele_AA_num;	
					}else{
						data_vcf(j,i+9)=allele_NN_num;	
					}
					
				}
			}		
			
Rcout<<" "<<endl;	
Rcout<<"Complete BLUPF90 to VCF data format conversion!"<<endl;
CharacterVector chr=data_blupf90_map.column(0),snp=data_blupf90_map.column(1),pos=data_blupf90_map.column(2);
//CharacterVector qual=rep("1",n_snp),filter=rep("PASS",n_snp),info=rep("2",n_snp),format=rep("GT",n_snp);
CharacterVector qual1={"."},filter1={"PASS"},info1={"."},format1={"GT"};
CharacterVector qual=rep(qual1,n_snp),filter=rep(filter1,n_snp),info=rep(info1,n_snp),format=rep(format1,n_snp);
data_vcf.column(0)=chr;
data_vcf.column(1)=pos;
data_vcf.column(2)=snp;
data_vcf.column(5)=qual;
data_vcf.column(6)=filter;
data_vcf.column(7)=info;
data_vcf.column(8)=format;
	
IND_name.push_front("FORMAT");
IND_name.push_front("INFO");
IND_name.push_front("FILTER");
IND_name.push_front("QUAL");
IND_name.push_front("ALT");
IND_name.push_front("REF");
IND_name.push_front("ID");
IND_name.push_front("POS");
IND_name.push_front("#CHROM");
colnames(data_vcf)=IND_name;
    return data_vcf;
}
