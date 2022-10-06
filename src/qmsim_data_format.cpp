#include "shared_function.h"

// [[Rcpp::export]]
arma::Mat<int> qmsim_to_numeric_cpp(std::vector<std::string> & data_qmsim,int cpu_cores=1){  //输入数据，行为SNP, 列为个体		
			
			omp_set_num_threads(cpu_cores);
			Rcout<<"Start QMsim to Numeric format conversion......"<<endl;	
			int i,j,k;	
			std::string tmp_ind_string,allele_string;
			int n_ind=data_qmsim.size(),n_snp=data_qmsim[0].size();
			arma::Mat<int> data_numeric(n_ind,n_snp);			
			
			#pragma omp parallel for private(i,j,tmp_ind_string,allele_string)
			for(i=0;i<n_ind;i++){
				
				tmp_ind_string=data_qmsim[i];
				
				for(j=0;j<n_snp;j++){
				
					allele_string=tmp_ind_string[j];
					
					if(allele_string=="3"|allele_string=="4"|allele_string=="5"){
						allele_string="1"; // default set missing as heterzygot
					}
					
					data_numeric(i,j)=std::stoi(allele_string);
				}
			}		
			
			Rcout<<" "<<endl;			
			Rcout<<"Complete QMsim to Numeric format conversion!"<<endl;	
			return data_numeric;
}

// [[Rcpp::export]]
std::vector<std::string> qmsim_to_blupf90_cpp(std::vector<std::string> & data_qmsim,int cpu_cores=1){  //输入数据，行为SNP, 列为个体		
			
			omp_set_num_threads(cpu_cores);
			Rcout<<"Start QMsim to BLUPF90 format conversion......"<<endl;	
			int i,j,k;	
			std::string tmp_ind_string;
			int n_ind=data_qmsim.size(),n_snp=data_qmsim[0].size();
			std::vector<std::string> allele_string_vector(n_ind);		
			
	        #pragma omp parallel for private(i,tmp_ind_string)
			for(i=0;i<n_ind;i++){
				
				tmp_ind_string=data_qmsim[i];								
				std::replace(tmp_ind_string.begin(), tmp_ind_string.end(), char('4'), char('1'));
				std::replace(tmp_ind_string.begin(), tmp_ind_string.end(), char('3'), char('1'));
				allele_string_vector[i]=tmp_ind_string;				
				
			}		
			
			Rcout<<" "<<endl;			
			Rcout<<"Complete QMsim to BLUPF90 format conversion!"<<endl;	
			return allele_string_vector;
}


// [[Rcpp::export]]
CharacterMatrix qmsim_to_hapmap_cpp(CharacterVector IND_name,
									  CharacterMatrix &data_qmsim_map,
                                      std::vector<std::string> & data_qmsim,
									  int cpu_cores=1,
									  std::string miss_base="N"){  //输入数据，	
			
			omp_set_num_threads(cpu_cores);
			Rcout<<"Start QMsim to Hapmap format conversion......"<<endl;	
			int i,j,k;	
			std::string tmp_ind_string,allele_string;
			int n_ind=data_qmsim.size(),n_snp=data_qmsim[0].size();
			arma::Mat<int> data_numeric(n_ind,n_snp);	
			CharacterVector tmp1=data_qmsim_map.column(3),tmp2=data_qmsim_map.column(4);
			std::vector<std::string> Ref=Rcpp::as<std::vector<std::string> >(tmp1);
			std::vector<std::string> Alt=Rcpp::as<std::vector<std::string> >(tmp2);
			std::string allele_aa,allele_Aa,allele_AA,allele_NN=miss_base+miss_base;
			CharacterMatrix data_hmp(n_snp,n_ind+11);
			int allele_numeric;
			#pragma omp parallel for private(i,j,tmp_ind_string,allele_string,allele_numeric,allele_aa,allele_Aa,allele_AA)
			for(i=0;i<n_ind;i++){
				
				tmp_ind_string=data_qmsim[i];
				
				for(j=0;j<n_snp;j++){
				
					allele_string=tmp_ind_string[j];
					allele_numeric=std::stoi(allele_string);				
					if(allele_numeric==0){
						allele_aa=Ref[j]+Ref[j];
						data_hmp(j,i+11)=allele_aa;
					}else if(allele_numeric==3){
						allele_Aa=Ref[j]+Alt[j];
						data_hmp(j,i+11)=allele_Aa;
					}else if(allele_numeric==4){
						allele_Aa=Alt[j]+Ref[j];
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
	Rcout<<"Complete QMsim to Hapmap format conversion!"<<endl;
			
	data_hmp.column(0)=data_qmsim_map.column(0);
	data_hmp.column(2)=data_qmsim_map.column(1);
	data_hmp.column(3)=data_qmsim_map.column(2);
	CharacterVector title={"rs#","alleles","chrom","pos","strand","center","protLSID","assembly","assayLSID","panelLSID","QCcode"};
	Function C_cpp("c");
	CharacterVector title_IND_name=C_cpp(title,IND_name);
	colnames(data_hmp)=title_IND_name;
	return data_hmp;			
}


// [[Rcpp::export]]
CharacterMatrix qmsim_to_ped_cpp(CharacterVector IND_name,
								   CharacterMatrix &data_qmsim_map,
                                   std::vector<std::string> & data_qmsim,
								   int cpu_cores=1,
								   std::string miss_base="N"){  
			omp_set_num_threads(cpu_cores);
			Rcout<<"Start QMsim to Plink format conversion......"<<endl;	
			int i,j,k;	
			std::string tmp_ind_string,allele_string;
			int n_ind=data_qmsim.size(),n_snp=data_qmsim[0].size();
			arma::Mat<int> data_numeric(n_ind,n_snp);			
			CharacterVector tmp1=data_qmsim_map.column(3),tmp2=data_qmsim_map.column(4);
			std::vector<std::string> Ref=Rcpp::as<std::vector<std::string> >(tmp1);
			std::vector<std::string> Alt=Rcpp::as<std::vector<std::string> >(tmp2);
			CharacterMatrix data_ped(n_ind,n_snp*2+6);
			int allele_numeric;
			#pragma omp parallel for private(i,j,tmp_ind_string,allele_string,allele_numeric)
			for(i=0;i<n_ind;i++){
				
				tmp_ind_string=data_qmsim[i];
				
				for(j=0;j<n_snp;j++){
				
					allele_string=tmp_ind_string[j];
					allele_numeric=std::stoi(allele_string);				
					if(allele_numeric==0){
						data_ped(i,6+2*j)=Ref[j];
						data_ped(i,6+2*j+1)=Ref[j];
					}else if(allele_numeric==3){
						data_ped(i,6+2*j)=Ref[j];
						data_ped(i,6+2*j+1)=Alt[j];
					}else if(allele_numeric==4){
						data_ped(i,6+2*j)=Alt[j];
						data_ped(i,6+2*j+1)=Ref[j];
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
	Rcout<<"Complete QMsim to Plink data format conversion!"<<endl;			
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
CharacterMatrix qmsim_to_vcf_cpp(CharacterVector IND_name,
								   CharacterMatrix &data_qmsim_map,
                                   std::vector<std::string> & data_qmsim,
								   std::string phased_symbol="|",
								   int cpu_cores=1){  //输入数据
			omp_set_num_threads(cpu_cores);
			Rcout<<"Start QMsim to VCF format conversion......"<<endl;	
			int i,j,k;	
			std::string tmp_ind_string,allele_string;
			int n_ind=data_qmsim.size(),n_snp=data_qmsim[0].size();	
			CharacterMatrix data_vcf(n_snp,n_ind+9);			
			CharacterVector tmp1=data_qmsim_map.column(3),tmp2=data_qmsim_map.column(4);
			std::vector<std::string> Ref=Rcpp::as<std::vector<std::string> >(tmp1);
			std::vector<std::string> Alt=Rcpp::as<std::vector<std::string> >(tmp2);
			data_vcf.column(3)=tmp1;
			data_vcf.column(4)=tmp2;

			int allele_numeric;	
			std::string allele_aa_num="0"+phased_symbol+"0";
			std::string allele_Aa_num="1"+phased_symbol+"0";
			std::string allele_aA_num="0"+phased_symbol+"1";
			std::string allele_AA_num="1"+phased_symbol+"1";
			std::string allele_NN_num="."+phased_symbol+".";
			
			#pragma omp parallel for private(i,j,tmp_ind_string,allele_string,allele_numeric)
			for(i=0;i<n_ind;i++){
				
				tmp_ind_string=data_qmsim[i];
				
				for(j=0;j<n_snp;j++){
				
					allele_string=tmp_ind_string[j];
					allele_numeric=std::stoi(allele_string);	
					
					if(allele_numeric==0){
						data_vcf(j,i+9)=allele_aa_num;
					}else if(allele_numeric==3){
						data_vcf(j,i+9)=allele_aA_num;
					}else if(allele_numeric==4){
						data_vcf(j,i+9)=allele_Aa_num;
					}else if(allele_numeric==2){
						data_vcf(j,i+9)=allele_AA_num;	
					}else{
						data_vcf(j,i+9)=allele_NN_num;	
					}
					
				}
			}		
			
Rcout<<" "<<endl;	
Rcout<<"Complete QMsim to VCF data format conversion!"<<endl;
CharacterVector chr=data_qmsim_map.column(1),snp=data_qmsim_map.column(0),pos=data_qmsim_map.column(2);
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
arma::Mat<int> qmsim_to_boa_cpp(std::vector<std::string> & data_qmsim,int cpu_cores=1){  	
			
			omp_set_num_threads(cpu_cores);
			Rcout<<"Start QMsim to Boa format conversion......"<<endl;	
			int i,j,k;	
			std::string tmp_ind_string,allele_string;
			int n_ind=data_qmsim.size(),n_snp=data_qmsim[0].size();
			arma::Mat<int> data_boa(n_snp,n_ind*2);			
			
			#pragma omp parallel for private(i,j,tmp_ind_string,allele_string)
			for(i=0;i<n_ind;i++){				
				
				for(j=0;j<n_snp;j++){
				
						data_boa(j,2*i)=1;
						data_boa(j,2*i+1)=2;
					}
					
				
			}		
			
			Rcout<<" "<<endl;			
			Rcout<<"Complete QMsim to Boa format conversion!"<<endl;	
			return data_boa;
}

//user_define_alt_ref_phased_vcf
// [[Rcpp::export]]
void user_define_phased_vcf_to_cpp(CharacterMatrix & data_vcf,
									std::vector<int> pos_conflict_ref,
									 int cpu_cores=1){
    Rcout<<"Start Phased-VCF to User-defined(Ref,Alt) Phased-VCF data format conversion......"<<endl;	
	omp_set_num_threads(cpu_cores);
	
	int i,j,pos,n_confilct_snp=pos_conflict_ref.size(),n_ind=data_vcf.ncol()-9;
	std::string allele;
	//CharacterMatrix my_data_vcf=clone(data_vcf);	
	Progress p(n_confilct_snp*n_ind,true);		
	#pragma omp parallel for private(i,j,pos,allele)
	for(i=0;i<n_confilct_snp;i++){
		pos=pos_conflict_ref[i];
		//Ref=data_vcf(pos,3);
		//Alt=data_vcf(i,4);
		//new_Ref=data_map(pos,1);
		//new_Alt=data_map(i,2);
		//if(Ref!=new_Ref){ //Ref and new_Ref are not the same, then need to modify 
			for(j=0;j<n_ind;j++){
				p.increment();
				allele=data_vcf(pos,j+9);
	
				if(allele=="0|1"){
					data_vcf(pos,j+9)="1|0";
				}else if(allele=="1|0"){
					data_vcf(pos,j+9)="0|1";
				}else if(allele=="0|0"){
					data_vcf(pos,j+9)="1|1";
				}else if(allele=="1|1"){
					data_vcf(pos,j+9)="0|0";
				}
				

			}
		//}
	}
	Rcout<<" "<<endl;
    Rcout<<"Complete PPhased-VCF to User-defined(Ref,Alt) Phased-VCF data format conversion!"<<endl;		
    //return data_vcf;
}


//cleaning old-version vcf, e.g. omit non-sense data in genotype 
// [[Rcpp::export]]
void purify_vcf_cpp(CharacterMatrix & data_vcf,
					int cpu_cores=1){
    Rcout<<"Start purify VCF format data......"<<endl;	
	omp_set_num_threads(cpu_cores);
	
	int i,j,pos,n_snp=data_vcf.nrow(),n_ind=data_vcf.ncol()-9;
	std::string allele,phased_symbol,tmp;
	tmp=data_vcf(1,10);
	phased_symbol=tmp[1];
	//CharacterMatrix my_data_vcf=clone(data_vcf);	
	Progress p(n_snp*n_ind,true);		
	#pragma omp parallel for private(i,j,pos,allele)
	for(i=0;i<n_snp;i++){
			for(j=0;j<n_ind;j++){
				p.increment();
				allele=data_vcf(i,j+9);	
				data_vcf(i,j+9)=allele[0]+phased_symbol+allele[2];			
			}
	}
	Rcout<<" "<<endl;
    Rcout<<"Complete purify VCF format data!"<<endl;		
    //return data_vcf;
}