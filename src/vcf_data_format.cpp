#include "shared_function.h"

//phased vcf to haplotype 
// [[Rcpp::export]]
arma::Mat<int> phased_vcf_to_haplotype_cpp(CharacterMatrix & data_vcf,int cpu_cores=1){
    Rcout<<"Start Phased-VCF to Haplotype data format conversion......"<<endl;	
	omp_set_num_threads(cpu_cores);
	
	int i,j,n_snp=data_vcf.nrow(),n_ind=data_vcf.ncol()-9;
	std::string allele1,allele2;
	std::string allele;
	arma::Mat<int> data_hap(n_snp,n_ind*2,fill::zeros);	
	Progress p(n_snp*n_ind,true);		
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
	Rcout<<" "<<endl;
    Rcout<<"Complete Phased-VCF to  Haplotype data format conversion!"<<endl;		
    return data_hap;
}

	
// phased_vcf to numeric 
// [[Rcpp::export]]
arma::Mat<int> phased_vcf_to_numeric_cpp(std::vector<int> block_start,
									  std::vector<int> block_end,
									  CharacterMatrix & data_vcf,
									  Rcpp::List &haplotype_allele,
									  int cpu_cores=1){

	omp_set_num_threads(cpu_cores);
	
	int i,j,n_snp=data_vcf.nrow(),n_ind=data_vcf.ncol()-9,m,n,tmp;
	std::string allele1,allele2;
	std::string allele;
	arma::Mat<int> data_hap(n_snp,n_ind*2,fill::zeros);	
	Rcout<<"Start Phased VCF data format conversion......"<<endl;
	Progress p0(n_snp*n_ind,true);
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
	Rcout<<" "<<endl;
	Rcout<<"Complete Phased VCF data format conversion!"<<endl;
//hapmap to numeric 
	Rcout<<"Start construct haplotype list...... "<<endl;	
	int window_n=block_start.size();
	std::vector<int> haplo_pos,haplo_type_num(window_n),cumsum_haplo_type_num;
	arma::Mat<int> data_haplotype_subset;
	arma::Mat<int> haplo_pos_list3(window_n,n_ind*2,fill::zeros);
	// get the overview set of haplo_data
	for(i=0; i < window_n; i++){
		data_haplotype_subset=data_hap.rows(block_start[i],block_end[i]);
		allele_get_haplotype_set_short(data_haplotype_subset,haplotype_allele,i);

	}
	Progress p1(window_n*n_ind,true);
					
	#pragma omp parallel for private(i,j,data_haplotype_subset,haplo_pos)
	for(i=0; i < window_n; i++){
		data_haplotype_subset=data_hap.rows(block_start[i],block_end[i]);
		haplo_pos=get_haplotype_set_short(data_haplotype_subset);
		haplo_type_num[i]=*max_element(haplo_pos.begin(),haplo_pos.end())+1;
		
		for(j=0;j<n_ind;j++){
			p1.increment();
			haplo_pos_list3(i,2*j)=haplo_pos[2*j];
			haplo_pos_list3(i,2*j+1)=haplo_pos[2*j+1];
			
		}
	}
	Rcout<<" "<<endl;
	Rcout<<"Complete construct haplotype list!"<<endl;		
	int sum_haplo_type_num=std::accumulate(haplo_type_num.begin(),haplo_type_num.end(),0);
	int ind_pos1,ind_pos2;  // all types of haplot
	
	cumulativeSum(haplo_type_num,cumsum_haplo_type_num);
	
	arma::Mat<int> data_numeric(n_ind,sum_haplo_type_num,fill::zeros);
	Rcout<<"Start haplotype-numeric(0,1,2) format conversion...... "<<endl;		
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
    Rcout<<" "<<endl;	
	Rcout<<"Complete haplotype-numeric(0,1,2) format conversion!"<<endl;	
			
    return data_numeric;
}



// phased_vcf to numeric 
// [[Rcpp::export]]
std::vector<std::string>  phased_vcf_to_blupf90_cpp(std::vector<int> block_start,
									       std::vector<int> block_end,
									       CharacterMatrix & data_vcf,
										   Rcpp::List &haplotype_allele,
									       int cpu_cores=1){

	omp_set_num_threads(cpu_cores);
	
	int i,j,n_snp=data_vcf.nrow(),n_ind=data_vcf.ncol()-9,m,n,tmp;
	CharacterVector IND_name=colnames(data_vcf);
	IND_name.erase(0,9);
	std::string allele1,allele2;
	std::string allele;
	arma::Mat<int> data_hap(n_snp,n_ind*2,fill::zeros);	
	Rcout<<"Start Phased VCF data format conversion......"<<endl;
	Progress p1(n_snp*n_ind,true);
	#pragma omp parallel for private(i,j,allele,allele1,allele2)
	for(i=0;i<n_snp;i++){
		for(j=0;j<n_ind;j++){		
			p1.increment();
			allele=data_vcf(i,j+9);
			allele1=allele[0];
			allele2=allele[2];
			data_hap(i,2*j)=std::stoi(allele1);
			data_hap(i,2*j+1)=std::stoi(allele2);
		}
	}
	Rcout<<""<<endl;
	Rcout<<"Complete Phased VCF data format conversion!"<<endl;
//hapmap to numeric 
	Rcout<<"Start construct haplotype list...... "<<endl;		
	int window_n=block_start.size();
	std::vector<int> haplo_pos,haplo_type_num(window_n),cumsum_haplo_type_num;
	arma::Mat<int> data_haplotype_subset;
	arma::Mat<int> haplo_pos_list3(window_n,n_ind*2,fill::zeros);
	// get the overview set of haplo_data
	for(i=0; i < window_n; i++){
		data_haplotype_subset=data_hap.rows(block_start[i],block_end[i]);
		allele_get_haplotype_set_short(data_haplotype_subset,haplotype_allele,i);

	}	
	Progress p(window_n*n_ind,true);					
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
	Rcout<<" "<<endl;
	Rcout<<"Complete construct haplotype list!"<<endl;	
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
	arma::Col<int> allele_string_col1,allele_string_col2;
	Progress p2(n_ind,true);
	Rcout<<"Start phased-haplotype-BLUPF90 format conversion......"<<endl;
	#pragma omp parallel for private(j,allele_string_col1,allele_string_col2,allele_string)
	for(j=0;j<n_ind;j++){
		p2.increment();
		allele_string_col1=haplo_pos_list3.col(2*j);
		allele_string_col2=haplo_pos_list3.col(2*j+1);
		
		allele_string=get_blupf90_allele_string_phased_vcf(allele_string_col1,allele_string_col2,max_length,cumsum_haplo_type_num,n_snp_hap);
	
		allele_string_vector[j]=allele_string;	
	}
	
	//添加个体名称
	//for(j=0;j<n_ind;j++){		
	//	ind_name=IND_name[j];	
	//	for(int k=0;k<ind_name.size();k++){
	//		allele_string_vector[j][k]=ind_name[k];
	//	}
	//}	
	
    Rcout<<" "<<endl;	
	Rcout<<"Complete phased-haplotype-BLUPF90 format conversion!"<<endl;	
			return allele_string_vector;
}


// phased_vcf to Plink 
// [[Rcpp::export]]
CharacterMatrix phased_vcf_to_Plink_cpp(std::vector<int> block_start,
									  std::vector<int> block_end,
									  CharacterMatrix & data_vcf,
									  Rcpp::List &haplotype_allele,
									  int cpu_cores=1){

	omp_set_num_threads(cpu_cores);
	
	int i,j,n_snp=data_vcf.nrow(),n_ind=data_vcf.ncol()-9,m,n,tmp;
	std::string allele1,allele2;
	std::string allele;
	arma::Mat<int> data_hap(n_snp,n_ind*2,fill::zeros);	
	Rcout<<"Start Phased VCF data format conversion......"<<endl;
	Progress p0(n_snp*n_ind,true);
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
	Rcout<<" "<<endl;
	Rcout<<"Complete Phased VCF data format conversion!"<<endl;
//hapmap to numeric 
	Rcout<<"Start construct haplotype list...... "<<endl;	
	int window_n=block_start.size();
	std::vector<int> haplo_pos,haplo_type_num(window_n),cumsum_haplo_type_num;
	arma::Mat<int> data_haplotype_subset;
	arma::Mat<int> haplo_pos_list3(window_n,n_ind*2,fill::zeros);
	// get the overview set of haplo_data
	for(i=0; i < window_n; i++){
		data_haplotype_subset=data_hap.rows(block_start[i],block_end[i]);
		allele_get_haplotype_set_short(data_haplotype_subset,haplotype_allele,i);

	}	
	Progress p1(window_n*n_ind,true);					
	#pragma omp parallel for private(i,j,data_haplotype_subset,haplo_pos)
	for(i=0; i < window_n; i++){
		data_haplotype_subset=data_hap.rows(block_start[i],block_end[i]);
		haplo_pos=get_haplotype_set_short(data_haplotype_subset);
		haplo_type_num[i]=*max_element(haplo_pos.begin(),haplo_pos.end())+1;
		
		for(j=0;j<n_ind;j++){
			p1.increment();
			haplo_pos_list3(i,2*j)=haplo_pos[2*j];
			haplo_pos_list3(i,2*j+1)=haplo_pos[2*j+1];
			
		}
	}
	Rcout<<" "<<endl;
	Rcout<<"Complete construct haplotype list!"<<endl;		
	int sum_haplo_type_num=std::accumulate(haplo_type_num.begin(),haplo_type_num.end(),0);
	int ind_pos1,ind_pos2;  // all types of haplot
	
	cumulativeSum(haplo_type_num,cumsum_haplo_type_num);
	
	arma::Mat<int> data_numeric(n_ind,sum_haplo_type_num,fill::zeros);
	//Rcout<<"Start haplotype-numeric(0,1,2) format conversion...... "<<endl;		
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
    //Rcout<<" "<<endl;	
	//Rcout<<"Complete haplotype-numeric(0,1,2) format conversion!"<<endl;	


CharacterMatrix data_ped(n_ind,sum_haplo_type_num*2+6);	

Rcout<<"Start Phased-VCF-to-Plink format conversion...... "<<endl;	

std::string  Ref="A",Alt="T";
int num_allele;
	Progress p(sum_haplo_type_num*n_ind,true);	
	#pragma omp parallel for private(i,j,num_allele)
	for(i=0;i<sum_haplo_type_num;i++){
		for(j=0;j<n_ind;j++){
			p.increment();
			num_allele=data_numeric(j,i);
			
			if(num_allele==0){
				data_ped(j,2*i+6)=Ref;				
				data_ped(j,2*i+6+1)=Ref;			
			}else if(num_allele==1){		
				data_ped(j,2*i+6)=Ref;				
				data_ped(j,2*i+6+1)=Alt;					
			}else{				
				data_ped(j,2*i+6)=Alt;				
				data_ped(j,2*i+6+1)=Alt;									
			}	
			

		}
	}


Rcout<<"Complete Phased-VCF-to-Plink format conversion!"<<endl;	

	CharacterVector SNP_name=data_vcf.column(2),IND_name=colnames(data_vcf);
	IND_name.erase(0,9);
	data_ped.column(0)=IND_name;
	data_ped.column(1)=IND_name;
	
    return data_ped;
}


//unphased-vcf to Plink
// [[Rcpp::export]]
CharacterMatrix vcf_to_plink_cpp(CharacterMatrix & data_vcf,int cpu_cores=1){
	Rcout<<"Start VCF to Plink data format conversion......"<<endl;	
	omp_set_num_threads(cpu_cores);	
	int i,j,n_snp=data_vcf.nrow(),n_ind=data_vcf.ncol()-9;
	std::string allele1,allele2;
	std::string allele,Ref,Alt;
	CharacterMatrix data_ped(n_ind,n_snp*2+6);	
	Progress p(n_snp*n_ind,true);	
	#pragma omp parallel for private(i,j,Ref,Alt,allele,allele1,allele2)
	for(i=0;i<n_snp;i++){
		Ref=data_vcf(i,3);
		Alt=data_vcf(i,4);
		for(j=0;j<n_ind;j++){
			p.increment();
			allele=data_vcf(i,j+9);
			allele1=allele[0];
			allele2=allele[2];
			
			if(allele1=="0"){
				data_ped(j,2*i+6)=Ref;				
			}else{
				data_ped(j,2*i+6)=Alt;		
			}
			
			if(allele2=="0"){
				data_ped(j,2*i+6+1)=Ref;				
			}else{
				data_ped(j,2*i+6+1)=Alt;
			}

		}
	}
	Rcout<<" "<<endl;
	Rcout<<"Complete VCF to Plink data format conversion!"<<endl;		
	CharacterVector SNP_name=data_vcf.column(2),IND_name=colnames(data_vcf);
	IND_name.erase(0,9);
	data_ped.column(0)=IND_name;
	data_ped.column(1)=IND_name;
    return data_ped;
}



//unphased-vcf to hapmap
// [[Rcpp::export]]
CharacterMatrix vcf_to_hapmap_cpp(CharacterMatrix & data_vcf,int cpu_cores=1,std::string miss_base="0"){
	
	omp_set_num_threads(cpu_cores);

	int i,j,n_snp=data_vcf.nrow(),n_ind=data_vcf.ncol()-9;
	std::string allele1,allele2;
	std::string allele,Ref,Alt;
	std::string tmp_allele,tmp1,tmp2,tmp3,tmp4,tmp5;
	CharacterMatrix data_hmp(n_snp,n_ind+11);	
	Rcout<<"Start VCF to Hapmap data format conversion......"<<endl;	
	bool ref1_status,ref2_status;	
	std::vector<std::string> allele_set1(n_snp),allele_set2(n_snp),allele_set3(n_snp),allele_set4(n_snp),allele_set5(n_snp);
	//get hapmap allele 
	for(i=0;i<n_snp;i++){
		Ref=data_vcf(i,3);
		Alt=data_vcf(i,4);		
		allele_set1[i]=Ref+Ref;
		allele_set2[i]=Alt+Alt;
		allele_set3[i]=Alt+Ref;
		allele_set4[i]=Ref+Alt;	
		allele_set5[i]=miss_base+miss_base;
	}
	Progress p(n_snp*n_ind,true);	
	#pragma omp parallel for schedule(dynamic) firstprivate(tmp1,tmp2,tmp3,tmp4,tmp5)  private(i,j,allele,allele1,allele2,ref1_status,ref2_status)

	for(i=0;i<n_snp;i++){
		tmp1=allele_set1[i];
		tmp2=allele_set2[i];
		tmp3=allele_set3[i];
		tmp4=allele_set4[i];
		tmp5=allele_set5[i];
		for(j=0;j<n_ind;j++){
			p.increment();
			allele=data_vcf(i,j+9);
			allele1=allele[0];
			allele2=allele[2];
			
			ref1_status=(allele1=="0");
			ref2_status=(allele2=="0");	

			if(ref1_status&&ref2_status){
			   data_hmp(i,j+11)=tmp1;
			}else if((allele1==miss_base)||(allele2==miss_base)){
			   data_hmp(i,j+11)=tmp5;
			}else if((!ref1_status)&&(!ref2_status)){				
			   data_hmp(i,j+11)=tmp2;				
			}else if(!ref1_status){
			   data_hmp(i,j+11)=tmp3;
			}else{
			   data_hmp(i,j+11)=tmp4;
			}
		
		}
	}
	Rcout<<" "<<endl;
	Rcout<<"Complete VCF to Hapmap data format conversion!"<<endl;

    return data_hmp;
}


//unphased-vcf to numeric
// [[Rcpp::export]]
arma::Mat<int> unphased_vcf_to_numeric_cpp(CharacterMatrix & data_vcf,int cpu_cores=1,std::string miss_base=".",int miss_base_num=0){
	Rcout<<"Start VCF to Numeric data format conversion......"<<endl;		
	omp_set_num_threads(cpu_cores);	
	int i,j,n_snp=data_vcf.nrow(),n_ind=data_vcf.ncol()-9,allele_num;
	std::string allele1,allele2;
	std::string allele;
	arma::Mat<int> data_numeric(n_ind,n_snp,fill::zeros);	
	Progress p(n_snp*n_ind,true);		
	#pragma omp parallel for private(i,j,allele,allele1,allele2,allele_num)
	for(i=0;i<n_snp;i++){
		for(j=0;j<n_ind;j++){
			p.increment();
			allele=data_vcf(i,j+9);
			allele1=allele[0];
			allele2=allele[2];
			
			if(allele1==miss_base||allele2==miss_base){
			data_numeric(j,i)=miss_base_num;
			}else{
			data_numeric(j,i)=std::stoi(allele1)+std::stoi(allele2);		
			}
		}
	}
	Rcout<<" "<<endl;
	Rcout<<"Complete VCF to Numeric data format conversion!"<<endl;	

    return data_numeric;
}


//unphased-vcf to BLUPF90
// [[Rcpp::export]]
std::vector<std::string> unphased_vcf_to_blupf90_unopenMP(CharacterMatrix & data_vcf,int cpu_cores=1,std::string miss_base="5",int miss_base_num=5){
	Rcout<<"Start VCF to BLUPF90 data format conversion......"<<endl;		
	omp_set_num_threads(cpu_cores);	
	int i,j,k,n_snp=data_vcf.nrow(),n_ind=data_vcf.ncol()-9,allele_num;
	
	CharacterVector IND_name=colnames(data_vcf);
	IND_name.erase(0,9);
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
	std::vector<std::string> allele_string_vector(n_ind,allele_string);
	Progress p(n_snp*n_ind,true);	
	#pragma omp parallel for private(i,j,allele,allele1,allele2,allele_num)
	for(j=0;j<n_ind;j++){
		
		for(i=0;i<n_snp;i++){
			p.increment();
			allele=data_vcf(i,j+9);
			allele1=allele[0];
			allele2=allele[2];		
			if(allele1==miss_base|allele2==miss_base){
				allele_num=miss_base_num;
			}else{
			    allele_num=std::stoi(allele1)+std::stoi(allele2);		
			}
			allele_string_vector[j][i+max_length]=allele_num+'0';
		}
	}
	
	//添加个体名称
	//for(j=0;j<n_ind;j++){		
	//	ind_name=IND_name[j];	
	//	for(k=0;k<ind_name.size();k++){
	//		allele_string_vector[j][k]=ind_name[k];
	//	}
	//}
	
	Rcout<<" "<<endl;
	Rcout<<"Complete VCF to BLUPF90 data format conversion!"<<endl;	
    return allele_string_vector;
}


// [[Rcpp::export]]
std::string	get_blupf90_allele_string_unphased_vcf(CharacterVector allele_string_row, std::string miss_base, int miss_base_num, int max_length){
	
	std::string allele,allele1,allele2,allele_string(allele_string_row.size()+max_length,' ');
	int allele_num;
	for(int i=0;i<allele_string_row.size();i++){		
	allele=allele_string_row[i];	
	allele1=allele[0];
	allele2=allele[2];	
	if(allele1==miss_base|allele2==miss_base){
	   allele_num=miss_base_num;
	}else{
	   allele_num=std::stoi(allele1)+std::stoi(allele2);		
		}
	   allele_string[i+max_length]=allele_num+'0';		
	}
	return allele_string;
}


//unphased-vcf to BLUPF90
// [[Rcpp::export]]
std::vector<std::string> unphased_vcf_to_blupf90_cpp(CharacterMatrix & data_vcf,int cpu_cores=1,std::string miss_base="5",int miss_base_num=5){
	Rcout<<"Start VCF to BLUPF90 data format conversion......"<<endl;		
	omp_set_num_threads(cpu_cores);	
	int i,j,k,n_snp=data_vcf.nrow(),n_ind=data_vcf.ncol()-9,allele_num;
	
	CharacterVector IND_name=colnames(data_vcf);
	IND_name.erase(0,9);
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
	CharacterVector	allele_string_row;
	Progress p(n_ind,true);	
	#pragma omp parallel for private(j,allele_string,allele_string_row)
	for(j=0;j<n_ind;j++){
		p.increment();
		allele_string_row=data_vcf.column(j+9);
	
		allele_string=get_blupf90_allele_string_unphased_vcf(allele_string_row,miss_base, miss_base_num, max_length);

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
	Rcout<<"Complete VCF to BLUPF90 data format conversion!"<<endl;	
    return allele_string_vector;
}

// type 1: haplotype; type 2: numeric ;
// [[Rcpp::export]]
List  vcf_convertion(
					CharacterMatrix & data_vcf, 	
					std::vector<int> block_start,
					std::vector<int> block_end,
					int type=1, 
					std::string miss_base="N",
					int miss_base_num=5,
					int cpu_cores=1){
		
		arma::Mat<int> data_haplotype,data_numeric;
		CharacterMatrix data_ped,data_hmp;
		std::vector<std::string> data_blupf90;

		int window_n=block_start.size();
		Function generate_list("vector");
		Rcpp::List haplotype_allele=generate_list("list", window_n);
		
		switch(type){
	
		case 1:{       //type:  vcf to Hapmap 
		
		data_hmp=vcf_to_hapmap_cpp(data_vcf,cpu_cores,miss_base);
		
			break;}

		case 2:{        //type:  vcf to Plink 
		
		data_ped=vcf_to_plink_cpp(data_vcf,cpu_cores);

			break;}
			

		case 8:{        //type:  phased-vcf to Plink 
		
		data_ped=phased_vcf_to_Plink_cpp(block_start,block_end,data_vcf,haplotype_allele,cpu_cores);

			break;}
			
		case 3:{        //type: phased vcf to haplotype
				
		data_haplotype=phased_vcf_to_haplotype_cpp(data_vcf,cpu_cores);
			
			break;}
			
		case 4:{        //type: phased vcf to numeric 
		
		data_numeric=phased_vcf_to_numeric_cpp(block_start,block_end,data_vcf,haplotype_allele,cpu_cores);
		
			break;}	

		case 5:{	   //type:  phased vcf to BLUPF90 

		data_blupf90=phased_vcf_to_blupf90_cpp(block_start,block_end,data_vcf,haplotype_allele,cpu_cores);
		
			break;}					
			
		case 6:{        //type:  unphased vcf to numeric 
		
		data_numeric=unphased_vcf_to_numeric_cpp(data_vcf,cpu_cores,miss_base);
		
			break;}	

		case 7:{	   //type:  unphased vcf to BLUPF90 

		data_blupf90=unphased_vcf_to_blupf90_cpp(data_vcf,cpu_cores=1,miss_base,miss_base_num);
		
			break;}
			
		default:
		
        throw Rcpp::exception("unknown input data type!");			
			
		}
		
		return List::create(Named("numeric") = data_numeric,
							Named("hap") = data_haplotype,
							Named("ped") = data_ped,
							Named("hmp") = data_hmp,
							Named("blupf90") = data_blupf90,
							Named("SNP")=haplotype_allele);
}