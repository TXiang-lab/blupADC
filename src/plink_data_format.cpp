#include "shared_function.h"



//Convert Plink_to_Hapmap 	
// [[Rcpp::export]]
CharacterMatrix ped_to_hapmap_cpp(CharacterMatrix & data_ped,CharacterMatrix data_map,int cpu_cores=5){
	
	omp_set_num_threads(cpu_cores);
	int n_ind=data_ped.nrow(); // 个体数
	int n_snp=(data_ped.ncol()-6)/2; //SNP数量	
	std::string char1,char2,char3;
	CharacterMatrix data_hmp(n_snp,n_ind+11);
	Rcout<<"Start Plink to Hapmap data format convertion......"<<endl;
	int i,j;
	#pragma omp parallel for private(i,j,char1,char2,char3)
	for(j=0;j<n_snp;j++){

	for(i=0; i < n_ind; i++){
	    char1=data_ped(i,2*j+6);
		char2=data_ped(i,2*j+7);
		char3=char1+char2;
		data_hmp(j,i+11)=char3;
		}	
		} 
	Rcout<<"Complete Plink to Hapmap data format convertion!"<<endl;
	
	data_hmp.column(0)=data_map.column(1);
	data_hmp.column(2)=data_map.column(0);
	data_hmp.column(3)=data_map.column(3);
	CharacterVector title={"rs#","alleles","chrom","pos","strand","center","protLSID","assembly","assayLSID","panelLSID","QCcode"};
	CharacterVector IND_name=data_ped.column(1);
	Function C_cpp("c");
	CharacterVector title_IND_name=C_cpp(title,IND_name);
	colnames(data_hmp)=title_IND_name;
	return data_hmp;	
}


//vcf to numeric
// [[Rcpp::export]]
CharacterMatrix plink_to_vcf_cpp(CharacterMatrix & data_ped,CharacterMatrix data_map, int cpu_cores=1,std::string miss_base="0",std::string phased_symbol="/"){

	omp_set_num_threads(cpu_cores);
	
	int i,j,n_snp=data_map.nrow(),n_ind=data_ped.nrow();
	
	std::string allele1,allele2,allele;
	std::string temp_SNP1,temp_SNP2,temp_base_SNP1,temp_base_SNP2,Ref,Alt;
	CharacterMatrix data_vcf(n_snp,n_ind+9),data_tmp_base_SNP(n_snp,2);	
	CharacterVector snp1,snp2;
	Rcout<<"Start Plink to VCF data format conversion......"<<endl;	
	for(i=0;i<n_snp;i++){	
		snp1=data_ped.column(2*i+6);
		snp2=data_ped.column(2*i+7);	
		temp_SNP1=single_base_factor_cpp(snp1,miss_base=miss_base);
		temp_SNP2=single_base_factor_cpp(snp2,miss_base=miss_base);
		data_tmp_base_SNP(i,0)=temp_SNP1;
		data_tmp_base_SNP(i,1)=temp_SNP2;
	}
	Progress p(n_snp*n_ind,true);	
	#pragma omp parallel for private(i,j,temp_base_SNP1,temp_base_SNP2,Ref,Alt,allele,allele1,allele2)	
	for(i=0;i<n_snp;i++){		
	
		temp_base_SNP1=data_tmp_base_SNP(i,0);
		temp_base_SNP2=data_tmp_base_SNP(i,1);
		
		if(temp_base_SNP1<temp_base_SNP2){
			Ref=temp_base_SNP1[0];
			Alt=temp_base_SNP2[1];		
		}else{
			Ref=temp_base_SNP2[0];
			Alt=temp_base_SNP1[1];								
		}
		
		data_vcf(i,3)=Ref;
		data_vcf(i,4)=Alt;
		
		for(j=0;j<n_ind;j++){
			p.increment();
			if(data_ped(j,2*i+6)==Ref){
				allele1="0";
			}else if(data_ped(j,2*i+6)==Alt){
				allele1="1";
			}else{
				allele1=miss_base;
			}

			if(data_ped(j,2*i+7)==Ref){
				allele2="0";
			}else if(data_ped(j,2*i+7)==Alt){
				allele2="1";
			}else{
				allele2=miss_base;
			}
			
			allele=allele1+phased_symbol+allele2;
			
			data_vcf(i,j+9)=allele;
		}
	}
	Rcout<<" "<<endl;	
	Rcout<<"Complete Plink to VCF data format conversion!"<<endl;
CharacterVector chr=data_map.column(0),snp=data_map.column(1),pos=data_map.column(3);
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
	
CharacterVector IND_name=data_ped.column(0);
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
arma::Mat<int> plink_to_numeric_cpp(CharacterMatrix & data_ped,CharacterMatrix data_map,int cpu_cores=5, std::string miss_base="0", int miss_base_num=0){  //输入数据，行为SNP, 列为个体
	
			omp_set_num_threads(cpu_cores);	
			//转换数据框为字符串矩阵
            int i,j,k,n_snp=(data_ped.ncol()-6)/2,n_ind=data_ped.nrow();
			CharacterVector temp_vec,temp_vec1,temp_vec2,SNP_name=data_map.column(0);
			CharacterMatrix Ref_type(n_snp,2);

			std::string base_SNP1,base_SNP2,temp_char,temp_base_SNP1,temp_base_SNP2,Ref,Alt;

	        Rcout<<"Start Plink to numeric(0,1,2) format conversion...... "<<endl;			
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
		    arma::Mat<int> data_numeric(n_ind,n_snp);  //行为个体,列为SNP
			data_numeric.fill(0);
			std::string allele1,allele2;
             Progress p(n_snp*n_ind,true);	
#pragma omp parallel for private(i,j,base_SNP1,allele1,allele2,allele1_status,allele2_status)				 
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
	        Rcout<<"Complete Plink to numeric(0,1,2) format conversion...... "<<endl;			
			return data_numeric;
}



// [[Rcpp::export]]
std::string	get_blupf90_allele_string_plink(CharacterVector & allele_string_row, CharacterMatrix & Ref_type,
                                            std::string miss_base, int miss_base_num, int max_length){
	
				bool allele1_status,allele2_status;	
				int allele_num,n_snp=(allele_string_row.size()-6)/2;
				std::string base_SNP1,allele1,allele2,allele_string(n_snp+max_length,' ');
				
				for(int i=0;i<n_snp;i++){
					base_SNP1=Ref_type(i,0);
					allele1=allele_string_row[2*i+6];
					allele2=allele_string_row[2*i+7];
				
					allele1_status=(allele1==base_SNP1);
					allele2_status=(allele2==base_SNP1);
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

// [[Rcpp::export]]
std::vector<std::string> plink_to_blupf90_cpp(CharacterMatrix & data_ped,CharacterMatrix data_map,int cpu_cores=5, std::string miss_base="0", int miss_base_num=0){  //输入数据，行为SNP, 列为个体
	
			omp_set_num_threads(cpu_cores);	
			//转换数据框为字符串矩阵
            int i,j,k,n_snp=(data_ped.ncol()-6)/2,n_ind=data_ped.nrow();
			CharacterVector temp_vec,temp_vec1,temp_vec2,SNP_name=data_map.column(0);
			CharacterMatrix Ref_type(n_snp,2);

			std::string base_SNP1,base_SNP2,temp_char,temp_base_SNP1,temp_base_SNP2,Ref,Alt;

	        Rcout<<"Start Plink to numeric(0,1,2) format conversion...... "<<endl;		
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
	        CharacterVector IND_name=data_ped.column(0);
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
	        CharacterVector	allele_string_row;	
            Progress p(n_snp*n_ind,false);	
#pragma omp parallel for private(i,allele_string_row,allele_string)				 
			for( i=0; i<n_ind;i++){	
			p.increment();
			allele_string_row=data_ped.row(i);	
			allele_string=get_blupf90_allele_string_plink(allele_string_row,Ref_type,miss_base,miss_base_num,max_length);			
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
	        Rcout<<"Complete Plink to numeric(0,1,2) format conversion...... "<<endl;			
			return allele_string_vector;
}


// type 1: haplotype; type 2: numeric ;
// [[Rcpp::export]]
List  plink_convertion(
					CharacterMatrix & data_ped,
					CharacterMatrix & data_map,					
					int type=1, 
					std::string miss_base="N",
					std::string phased_symbol="/",
					int miss_base_num=5,
					int cpu_cores=1){
		
		arma::Mat<int> data_numeric;
		CharacterMatrix data_hmp,data_vcf;
		std::vector<std::string> data_blupf90;
		
		switch(type){
			
		case 1:{        //type:  plink to Hapmap 
		
		data_hmp=ped_to_hapmap_cpp(data_ped,data_map,cpu_cores);
		
			break;}


		case 2:{        //type:  plink to vcf
		
		data_vcf=plink_to_vcf_cpp(data_ped,data_map,cpu_cores,miss_base,phased_symbol);
			
			break;}	
			
			
		case 3:{        //type: plink to numeric
		
				
		data_numeric=plink_to_numeric_cpp(data_ped,data_map,cpu_cores,miss_base,miss_base_num);
			
			break;}
			
		case 4:{        //type:  plink to blupf90
		
		data_blupf90=plink_to_blupf90_cpp(data_ped,data_map,cpu_cores, miss_base,miss_base_num);
		
			break;}	
			
		default:
		
        throw Rcpp::exception("unknown input data type!");			
			
		}
		
		return List::create(Named("numeric") = data_numeric,
							Named("hmp") = data_hmp,
							Named("vcf")=data_vcf,
							Named("blupf90") = data_blupf90);
}