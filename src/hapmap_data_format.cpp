#include "shared_function.h"

//Convert Plink_to_Hapmap 	
// [[Rcpp::export]]
CharacterMatrix hapmap_to_ped_cpp(CharacterMatrix & data_hmp,int cpu_cores=5){
	
	omp_set_num_threads(cpu_cores);
    int n_ind=data_hmp.ncol()-11,n_snp=data_hmp.nrow(),i,j;
	CharacterVector temp_vec1,IND_name=colnames(data_hmp);
	IND_name.erase(0,11);
	std::string char1,char2,char3;
	CharacterMatrix data_ped(n_ind,n_snp*2+6);
	Rcout<<"Start Hapmap to Plink data format conversion......"<<endl;	
	Progress p(n_snp*n_ind,true);	
	#pragma omp parallel for private(i,j,char1,char2,char3)
	for( i=0; i < n_ind; i++){
	for( j=0;j<n_snp;j++){
		p.increment();
	    char1=data_hmp(j,i+11);
		char2=char1[0]; //疑问：如果不定义char2,直接令 data_ped(i,2*j)=char1[0]， openMP会报错，具体原因待查
		char3=char1[1];
		data_ped(i,6+2*j)=char2;
		data_ped(i,6+2*j+1)=char3;
		}
		}
	Rcout<<"Complete Hapmap to Plink data format conversion!"<<endl;			
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


// Numgeno  in Rcpp
// [[Rcpp::export]]
arma::Mat<int> hapmap_to_numeric_cpp(CharacterMatrix & data_hmp, std::string miss_base="N", int miss_base_num=0, int cpu_cores=1){  //输入数据，行为SNP, 列为个体

	omp_set_num_threads(cpu_cores);
	std::string tmp_allele,tmp_allele1,tmp_allele2;	
	int n_snp=data_hmp.nrow(),n_ind=data_hmp.ncol()-11,i,j;
	arma::Mat<int> data_numeric(n_ind,n_snp);
	data_numeric.fill(0);
	CharacterVector snp_vec;
	std::string factor_result,Ref,Alt,allele,allele1,allele2;
	CharacterMatrix Ref_type(n_snp,2);	
	CharacterVector snp1,snp2;
	Rcout<<"Start Hapmap to Numeric data format conversion......"<<endl;
		
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
	Progress p(n_snp*n_ind,true);	
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
	Rcout<<"Complete Hapmap to Numeric data format conversion!"<<endl;
    return data_numeric;
}



// [[Rcpp::export]]
std::string	get_blupf90_allele_string_hapmap(CharacterVector allele_string_col, CharacterMatrix Ref_type,
                                            std::string miss_base, int miss_base_num, int max_length){
	
				bool allele1_status,allele2_status;	
				int allele_num,n_snp=allele_string_col.size();
				std::string Ref,allele,allele1,allele2,allele_string(n_snp+max_length,' ');
				
				for(int i=0;i<n_snp;i++){
					Ref=Ref_type(i,0);
					allele=allele_string_col[i];
					allele1=allele[0];
					allele2=allele[1];
				
					allele1_status=(allele1==Ref);
					allele2_status=(allele2==Ref);
					if(allele1_status&&allele2_status){					
						allele_num=0;		
					}else if(allele1==miss_base|allele2==miss_base){				
						allele_num=miss_base_num;  //将缺失值视为隐性纯合，设为0						
					}else if(allele1_status||allele2_status){	
						allele_num=1;					
					}else {					
						allele_num=2;
					}
				allele_string[i+max_length]=allele_num+'0';						
				}
				return allele_string;
}


// Numgeno  in Rcpp
// [[Rcpp::export]]
std::vector<std::string> hapmap_to_blupf90_cpp(CharacterMatrix &data_hmp, std::string miss_base="N", int miss_base_num=0, int cpu_cores=1){  //输入数据，行为SNP, 列为个体

	omp_set_num_threads(cpu_cores);
	std::string tmp_allele,tmp_allele1,tmp_allele2;	
	int n_snp=data_hmp.nrow(),n_ind=data_hmp.ncol()-11,i,j,k;
	CharacterVector snp_vec;
	std::string factor_result,Ref,Alt,allele,allele1,allele2;
	CharacterMatrix Ref_type(n_snp,2);	
	CharacterVector snp1,snp2;
	Rcout<<"Start Hapmap to BLUPF90 data format conversion......"<<endl;
		
	for(i=0;i<n_snp;i++){	
		snp_vec=data_hmp.row(i);
		snp_vec.erase(0,11);
		factor_result=pair_base_factor_cpp(snp_vec,miss_base=miss_base);
		Ref=factor_result[0];
		Alt=factor_result[1];
		Ref_type(i,0)=Ref;
		Ref_type(i,1)=Alt;
	}


	        CharacterVector IND_name=colnames(data_hmp);
			IND_name.erase(0,11);
	        std::string ind_name;
	        std::vector<int> ind_name_size(IND_name.size());
	        //确定个体名称中最大字符串长度
	        for(i=0;i<IND_name.size();i++){
	        	ind_name=IND_name[i];
	        	ind_name_size[i]=ind_name.length();
	        }
	        int max_length=*max_element(ind_name_size.begin(),ind_name_size.end())+3;	
			max_length=0;
	        std::string allele_string(n_snp+max_length,' ');
	        std::vector<std::string> allele_string_vector(n_ind);
	        CharacterVector	allele_string_col;	
            Progress p(n_ind,true);	
#pragma omp parallel for private(i,allele_string_col,allele_string)				 
			for( i=0; i<n_ind;i++){	
			p.increment();
			allele_string_col=data_hmp.column(i+11);	
			allele_string=get_blupf90_allele_string_hapmap(allele_string_col,Ref_type,miss_base,miss_base_num,max_length);			
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
	Rcout<<"Complete Hapmap to BLUPF90 data format conversion!"<<endl;
    return allele_string_vector;
}



//hapmap to vcf
// [[Rcpp::export]]
CharacterMatrix hapmap_to_VCF_cpp(CharacterMatrix &data_hmp, int cpu_cores=1,std::string miss_base="N",std::string phased_symbol="/"){

	omp_set_num_threads(cpu_cores);
	
	int i,j,n_snp=data_hmp.nrow(),n_ind=data_hmp.ncol()-11;
	
	std::string allele1,allele2,allele,allele_num,allele1_num,allele2_num;
	std::string factor_result,tmp_allele,tmp_allele1,tmp_allele2,Ref,Alt;
	CharacterVector snp_vec;
	CharacterMatrix Ref_type(n_snp,2);
	CharacterMatrix data_vcf(n_snp,n_ind+9);	
	CharacterVector snp1,snp2;
	for(i=0;i<n_snp;i++){	
		snp_vec=data_hmp.row(i);
		snp_vec.erase(0,11);
		factor_result=pair_base_factor_cpp(snp_vec,miss_base=miss_base);
		Ref=factor_result[0];
		Alt=factor_result[1];
		Ref_type(i,0)=Ref;
		Ref_type(i,1)=Alt;
	}
	Rcout<<"Start Hapmap to VCF data format conversion......"<<endl;	
	Progress p(n_snp*n_ind,true);
	#pragma omp parallel for private(i,j,Ref,Alt,allele_num,allele1_num,allele2_num,allele,allele1,allele2)	
	for(i=0;i<n_snp;i++){	

		Ref=Ref_type(i,0);
		Alt=Ref_type(i,1);		
		data_vcf(i,3)=Ref;
		data_vcf(i,4)=Alt;
		
		for(j=0;j<n_ind;j++){
			p.increment();
			allele=data_hmp(i,j+11);
			allele1=allele[0];
			allele2=allele[1];
			
			if(allele1==Ref){
				allele1_num="0";
			}else if(allele1==Alt){
				allele1_num="1";
			}else{
				allele1_num=miss_base;
			}
			
			if(allele2==Ref){
				allele2_num="0";
			}else if(allele2==Alt){
				allele2_num="1";
			}else{
				allele2_num=miss_base;
			}
			
			allele_num=allele1_num+phased_symbol+allele2_num;
			
			data_vcf(i,j+9)=allele_num;
		}
	}
	Rcout<<" "<<endl;	
	Rcout<<"Complete Hapmap to VCF data format conversion!"<<endl;
CharacterVector chr=data_hmp.column(2),snp=data_hmp.column(0),pos=data_hmp.column(3);
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
	
CharacterVector IND_name=colnames(data_hmp);
IND_name.erase(0,11);
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


//hapmap to vcf
// [[Rcpp::export]]
CharacterMatrix hapmap_to_VCF_cpp1(CharacterMatrix & data_hmp, int cpu_cores=1,std::string miss_base=".",std::string phased_symbol="/"){

	omp_set_num_threads(cpu_cores);
	
	int i,j,n_snp=data_hmp.nrow(),n_ind=data_hmp.ncol()-11;
	
	std::string allele1,allele2,allele,allele_num,allele1_num,allele2_num;
	std::string factor_result,tmp_allele,tmp_allele1,tmp_allele2,Ref,Alt;
	CharacterVector snp_vec;
	CharacterMatrix Ref_type(n_snp,2);
	CharacterMatrix data_vcf(n_snp,n_ind+9);	
	CharacterVector snp1,snp2;
	Rcout<<"Start Hapmap to VCF data format conversion......"<<endl;
	Progress p(n_snp*n_ind,true);	

	for(i=0;i<n_snp;i++){	
		snp_vec=data_hmp.row(i);
		snp_vec.erase(0,11);
		factor_result=pair_base_factor_cpp(snp_vec,miss_base=miss_base);
		Ref=factor_result[0];
		Alt=factor_result[1];
		Ref_type(i,0)=Ref;
		Ref_type(i,1)=Alt;
	}

	std::string ref_ref="0"+phased_symbol+"0";
	std::string ref_alt="0"+phased_symbol+"1";
	std::string alt_ref="1"+phased_symbol+"0";
	std::string alt_alt="1"+phased_symbol+"1";
	std::string miss_allle=miss_base+phased_symbol+miss_base;
	bool allele1_status,allele2_status;
	#pragma omp parallel for private(i,j,Ref,Alt,allele_num,allele1_status,allele2_status,allele,allele1,allele2)	
	for(i=0;i<n_snp;i++){	

		Ref=Ref_type(i,0);
		Alt=Ref_type(i,1);		
		data_vcf(i,3)=Ref;
		data_vcf(i,4)=Alt;
		
		for(j=0;j<n_ind;j++){
			p.increment();
			allele=data_hmp(i,j+11);
			allele1=allele[0];
			allele2=allele[1];

			allele1_status=(allele1==Ref);
			allele2_status=(allele2==Ref);
			if(allele1_status&&allele2_status){					
				allele_num=ref_ref;		
			}else if(allele1==miss_base|allele2==miss_base){				
				allele_num=miss_allle;
			}else if(!allele1_status&&!allele2_status){	
				allele_num=alt_alt;						
			}else if(allele1_status){	
				allele_num=ref_alt;		
			}else{	
				allele_num=alt_ref;					
			}
			
			data_vcf(i,j+9)=allele_num;
		}
	}
	Rcout<<" "<<endl;	
	Rcout<<"Complete Hapmap to VCF data format conversion!"<<endl;
CharacterVector chr=data_hmp.column(2),snp=data_hmp.column(0),pos=data_hmp.column(3);
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
	
CharacterVector IND_name=colnames(data_hmp);
IND_name.erase(0,11);
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

// type 1: haplotype; type 2: numeric ;
// [[Rcpp::export]]
List  hapmap_convertion(
				    CharacterMatrix & data_hmp, 	
					int type=1, 
					std::string miss_base="N",
					std::string phased_symbol="/",
					int miss_base_num=5,
					int cpu_cores=1){
		
		arma::Mat<int> data_numeric;
		CharacterMatrix data_ped,data_vcf;
		std::vector<std::string> data_blupf90;
		
		switch(type){
			
		case 1:{        //type:  hapmap_to_ped_cpp
		
		data_ped=hapmap_to_ped_cpp(data_hmp,cpu_cores);
		
			break;}


		case 2:{        //type:  hapmap_to_VCF_cpp
		
		data_vcf=hapmap_to_VCF_cpp(data_hmp,cpu_cores,miss_base,phased_symbol);
			
			break;}	
			
			
		case 3:{        //type: hapmap_to_numeric_cpp
		
				
		data_numeric=hapmap_to_numeric_cpp(data_hmp,miss_base,miss_base_num,cpu_cores);
			
			break;}
			
		case 4:{        //type:  hapmap_to_blupf90_cpp
		
		data_blupf90=hapmap_to_blupf90_cpp(data_hmp,miss_base,miss_base_num,cpu_cores);
		
			break;}	
			
		default:
		
        throw Rcpp::exception("unknown input data type!");			
			
		}
		
		return List::create(Named("numeric") = data_numeric,
							Named("ped") = data_ped,
							Named("vcf")=data_vcf,
							Named("blupf90") = data_blupf90);
}