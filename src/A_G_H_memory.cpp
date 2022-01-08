#include "shared_function.h"


template <typename T>
SEXP bigmemory_double_type(SEXP pBigMat, std::string numeric_file_name, std::string numeric_file_path, int cpu_cores=1){
	
	Rcpp::XPtr<BigMatrix> pMat(pBigMat);
	int type_number=(pMat->matrix_type());
	arma::Mat<T> data_numeric((T*) pMat->matrix(), pMat -> nrow(), pMat -> ncol(), false);
	int n_ind=data_numeric.n_rows,snp_n=data_numeric.n_cols;
	if(type_number!=8){
	Rcout<<"bigmemory numeric data(double type) is save as "<<numeric_file_path+"/"+numeric_file_name+"_double_numeric"<<endl;
	
	SEXP pBigMat_double=make_bigmemory_object_address_cpp(pMat -> nrow(),pMat -> ncol(),numeric_file_name+"_double_numeric",numeric_file_path,"double");
	Rcpp::XPtr<BigMatrix> pMat_double(pBigMat_double);
	arma::Mat<double> M_double((double*) pMat_double->matrix(), pMat_double -> nrow(), pMat_double -> ncol(), false);	
		//M_double=conv_to<arma::Mat<double>>::from(data_numeric);
		int i,j;
		omp_set_num_threads(cpu_cores);
		#pragma omp parallel for private(i,j)		
		for(i=0;i<n_ind;i++){
			for(j=0;j<snp_n;j++){
				M_double(i,j)=data_numeric(i,j);
			}
		}
		//M_double(0,4)=1000;	
		return pBigMat_double;
	}else{
	Rcout<<"bigmemory numeric data(double type) is save as "<<numeric_file_path+"/"+numeric_file_name+"_double_numeric"<<endl;	
		return pBigMat;
	}
}


// [[Rcpp::export]]
SEXP bigmemory_double_type(SEXP pBigMat, std::string numeric_file_name, std::string numeric_file_path,int cpu_cores=1){
	
	Rcpp::XPtr<BigMatrix> pMat(pBigMat);
    switch(pMat->matrix_type()) {
    case 1:		
        Rcout<<"Input bigmemory object is char type, convert to double type...... "<<endl;		
		return bigmemory_double_type<char>(pBigMat,numeric_file_name,numeric_file_path,cpu_cores);
    case 2:
        Rcout<<"Input bigmemory object is short type, convert to double type...... "<<endl;
		return bigmemory_double_type<short>(pBigMat,numeric_file_name,numeric_file_path,cpu_cores);
    case 4:
        Rcout<<"Input bigmemory object is integer type, convert to double type...... "<<endl;
		return bigmemory_double_type<int>(pBigMat,numeric_file_name,numeric_file_path,cpu_cores);
    case 6:
        Rcout<<"Input bigmemory object is double type, convert to double type......"<<endl;
		return bigmemory_double_type<double>(pBigMat,numeric_file_name,numeric_file_path,cpu_cores);
    case 8:
        Rcout<<"Input bigmemory object is double type,no need to convet!"<<endl;
		return bigmemory_double_type<double>(pBigMat,numeric_file_name,numeric_file_path,cpu_cores);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }	
	
}







// [[Rcpp::export]]
SEXP bigmemory_object_convert(arma::Mat<double> & input_data_numeric,
							   std::string bigmemory_data_name,
						       std::string bigmemory_data_path){
								   
		int n_ind=input_data_numeric.n_rows,n_snp=input_data_numeric.n_cols;
		
		SEXP pBigMat_temp=make_bigmemory_object_address_cpp(n_ind,n_snp,
															bigmemory_data_name+"_convert_numeric",bigmemory_data_path,"double");
		Rcpp::XPtr<BigMatrix> pMat_temp(pBigMat_temp);
		arma::Mat<double> M1((double*) pMat_temp->matrix(), pMat_temp -> nrow(), pMat_temp -> ncol(), false);	
		M1.submat(span(0,n_ind-1),span(0,n_snp-1))=input_data_numeric;

		SEXP pBigMat_temp_object=make_bigmemory_object_cpp(n_ind,n_snp,
															bigmemory_data_name+"_convert_numeric",bigmemory_data_path,"double");		
		return pBigMat_temp_object;
								   
}

// [[Rcpp::export]]
arma::Mat<double> makeA_tmp_cpp(arma::Mat<int> Pedigree){   //稀疏矩阵运行速度会变慢


	//Rcout<<"Start constructing pedigree additive relationship matrix......"<<endl;  		
	arma::Mat<double> A_matrix(Pedigree.n_rows,Pedigree.n_rows);	
	A_matrix.fill(0);
	arma::Col<int> Animal=Pedigree.col(0);
	arma::Col<int> Sire=Pedigree.col(1);
	arma::Col<int> Dam=Pedigree.col(2);
	int IND=Pedigree.n_rows;
	bool status1,status2;
	
	for(int i=0; i < IND; i++){

	status1=(Sire[i]==0);
	status2=(Dam[i]==0);


	if(!status1&&!status2){  //父母均已知
		A_matrix(i,i)=1+0.5*(A_matrix(Sire[i]-1,Dam[i]-1));	
		for(int j=0;j<i;j++){
			if(i>j) A_matrix(j,i)=A_matrix(i,j)=0.5*(A_matrix(j,Sire[i]-1)+A_matrix(j,Dam[i]-1));
	    }
	}else if(status1&&status2){  //父母均未知
		A_matrix(i,i)=1;	
		for(int j=0;j<i;j++){
			if(i>j) A_matrix(j,i)=0;
	    }
	}else if(!status1&&status2){  //父已知，母未知	
		A_matrix(i,i)=1;
		for(int j=0;j<i;j++){
			if(i>j) A_matrix(j,i)=A_matrix(i,j)=0.5*(A_matrix(j,Sire[i]-1));
		}
	}else{  //父未知，母已知	
		A_matrix(i,i)=1;
		for(int j=0;j<i;j++){
			if(i>j) A_matrix(j,i)=A_matrix(i,j)=0.5*(A_matrix(j,Dam[i]-1));
	    }
	}
	}
	//Rcout<<"Complete constructing pedigree additive relationship matrix!"<<endl;  	
	return A_matrix;
	}

//计算加性效应矩阵	
// [[Rcpp::export]]
List G_matrix_memory_cpp(SEXP pBigMat,std::string bigmemory_data_name, std::string bigmemory_data_path,
							bool base=false,bool trace=false, bool metafounder=false,
							bool inverse=false){ //M矩阵
							
	Rcpp::XPtr<BigMatrix> pMat(pBigMat); 
	arma::Mat<double> M_temp((double*) pMat->matrix(), pMat -> nrow(), pMat -> ncol(), false); 

delete_bigmemory_file_cpp("temp_numeric",bigmemory_data_name,bigmemory_data_path,false);	
	
	
	//生成临时的M1矩阵，因为M1矩阵经过标准化后，内容会发生改变
	SEXP pBigMat_temp=make_bigmemory_object_address_cpp(pMat -> nrow(),pMat -> ncol(),bigmemory_data_name+"_temp_numeric",bigmemory_data_path,"double");
	Rcpp::XPtr<BigMatrix> pMat_temp(pBigMat_temp);
	arma::Mat<double> M1((double*) pMat_temp->matrix(), pMat_temp -> nrow(), pMat_temp -> ncol(), false);	
	int n_ind=M1.n_rows; // 个体数
	int n_snp=M1.n_cols; //SNP数量	
	M1.submat(span(0,n_ind-1),span(0,n_snp-1))=M_temp;

delete_bigmemory_file_cpp("G_A",bigmemory_data_name,bigmemory_data_path,true);	
	

	SEXP pBigMat_G=make_bigmemory_object_address_cpp(n_ind,n_ind,bigmemory_data_name+"_G_A",bigmemory_data_path,"double");
	Rcpp::XPtr<BigMatrix> pMat_G(pBigMat_G);
	arma::Mat<double> G((double*) pMat_G->matrix(), pMat_G -> nrow(), pMat_G -> ncol(), false);

	Rcout<<"bigmemory-Start constructing genomic additive relationship matrix......"<<endl;  
	rowvec p_frequency(n_snp),q_frequency(n_snp);

	arma::Mat<double> b;
	rowvec a;
	if(base==TRUE){
	p_frequency.fill(0.5);
	q_frequency.fill(0.5);
	for(int i=0; i < n_snp; i++){
	M1.col(i)=M1.col(i)-1;
	}	
	}else{
	for(int i=0; i < n_snp; i++){

	p_frequency[i]=mean(M1.col(i));
	M1.col(i)=M1.col(i)-mean(M1.col(i));
	}
	p_frequency=0.5*p_frequency;
	q_frequency=1-p_frequency;
	}

	if(trace==true){
	b=M1*M1.t();
	G=(2*n_ind*b)/(arma::trace(b));				
	}else if(metafounder==true){
	G=M1*M1.t()/(n_snp/2);			
	}else{
	a=2*p_frequency%q_frequency;
	G=M1*M1.t()/sum(a);
	}

Rcout<<"bigmemory-Complete constructing genomic additive relationship matrix!"<<endl; 

	SEXP G_bigmemory_object=make_bigmemory_object_cpp(n_ind,n_ind,bigmemory_data_name+"_G_A",bigmemory_data_path,"double");	
	SEXP G_inv_bigmemory_object=R_NilValue;
	
	if(inverse==true){

delete_bigmemory_file_cpp("G_Ainv",bigmemory_data_name,bigmemory_data_path,true);		
		
	SEXP pBigMat_Ginv=make_bigmemory_object_address_cpp(n_ind,n_ind,bigmemory_data_name+"_G_Ainv",bigmemory_data_path,"double");
	Rcpp::XPtr<BigMatrix> pMat_Ginv(pBigMat_Ginv);

	arma::Mat<double> Ginv((double*) pMat_Ginv->matrix(), pMat_Ginv -> nrow(), pMat_Ginv -> ncol(), false);	
	G.diag()=G.diag()+0.0001;	
	Ginv.submat(span(0,n_ind-1),span(0,n_ind-1))=inv(G);   //如果直接令 Ginv=inv(G), Ginv的值并没有改变，反而会等于G,不晓得原因
	G.diag()=G.diag()-0.0001;	
	
	G_inv_bigmemory_object=make_bigmemory_object_cpp(n_ind,n_ind,bigmemory_data_name+"_G_Ainv",bigmemory_data_path,"double");
Rcout<<"bigmemory-Complete constructing inverse-genomic additive relationship matrix!"<<endl; 	
	}
//删除临时的 temp_numeric.desc 和  .bk文件
delete_bigmemory_file_cpp("temp_numeric",bigmemory_data_name,bigmemory_data_path,false);
	
	return List::create(Named("G") = G_bigmemory_object,_["Ginv"] = G_inv_bigmemory_object);
	
	
	
	
}


//计算显性离差矩阵
// [[Rcpp::export]]
List Deviation_matrix_memory_cpp(SEXP pBigMat,std::string bigmemory_data_name, std::string bigmemory_data_path,
							     bool base=false,bool trace=false, bool inverse=false){


	Rcpp::XPtr<BigMatrix> pMat(pBigMat); 
	arma::Mat<double> M_temp((double*) pMat->matrix(), pMat -> nrow(), pMat -> ncol(), false); 

delete_bigmemory_file_cpp("temp_numeric",bigmemory_data_name,bigmemory_data_path,false);	
	
	
	SEXP pBigMat_temp=make_bigmemory_object_address_cpp(pMat -> nrow(),pMat -> ncol(),bigmemory_data_name+"_temp_numeric",bigmemory_data_path,"double");
	Rcpp::XPtr<BigMatrix> pMat_temp(pBigMat_temp);
	arma::Mat<double> M1((double*) pMat_temp->matrix(), pMat_temp -> nrow(), pMat_temp -> ncol(), false);	
	int n_ind=M1.n_rows; // 个体数
	int n_snp=M1.n_cols; //SNP数量	
	M1.submat(span(0,n_ind-1),span(0,n_snp-1))=M_temp;

delete_bigmemory_file_cpp("G_D",bigmemory_data_name,bigmemory_data_path,true);	
	

	SEXP pBigMat_D=make_bigmemory_object_address_cpp(n_ind,n_ind,bigmemory_data_name+"_G_D",bigmemory_data_path,"double");
	Rcpp::XPtr<BigMatrix> pMat_D(pBigMat_D);
	arma::Mat<double> Devi((double*) pMat_D->matrix(), pMat_D -> nrow(), pMat_D -> ncol(), false);
	
	Rcout<<"bigmemory-Start constructing genomic dominance(classical) relationship matrix......"<<endl;  

	double p_freq,q_freq;
	arma::Row<double> p_frequency(n_snp), q_frequency(n_snp);
	if(base==TRUE){
	p_frequency.fill(0.5);
	q_frequency.fill(0.5);
	for( int i=0; i < n_snp; i++){
	for( int j=0; j < n_ind; j++){
	if(M1(j,i)==2){M1(j,i)=-2*0.25;
	}else if(M1(j,i)==1){M1(j,i)=2*0.25;
	}else{M1(j,i)=-2*0.25;}	
	}}
    }else{	
	for( int i=0; i < n_snp; i++){
	p_freq=0.5*mean(M1.col(i));
	q_freq=1-p_freq;
	for( int j=0; j < n_ind; j++){
	if(M1(j,i)==2){M1(j,i)=-2*q_freq*q_freq;
	}else if(M1(j,i)==1){M1(j,i)=2*p_freq*q_freq;
	}else{M1(j,i)=-2*p_freq*p_freq;}
	p_frequency[i]=p_freq;
	}
	}
	 q_frequency=1-p_frequency;
	} 
	 
	arma::Row<double> a=4*pow(p_frequency,2)%pow(q_frequency,2); //行向量*行向量
	Devi=M1*M1.t()/sum(a);

	SEXP D_bigmemory_object=make_bigmemory_object_cpp(n_ind,n_ind,bigmemory_data_name+"_G_D",bigmemory_data_path,"double");
	SEXP D_inv_bigmemory_object=R_NilValue;

Rcout<<"bigmemory-Complete constructing genomic dominance(classical) relationship matrix!"<<endl; 

	if(inverse==TRUE){
		
delete_bigmemory_file_cpp("G_Dinv",bigmemory_data_name,bigmemory_data_path,true);	
			
	SEXP pBigMat_Devi_inv=make_bigmemory_object_address_cpp(n_ind,n_ind,bigmemory_data_name+"_G_Dinv",bigmemory_data_path,"double");
	Rcpp::XPtr<BigMatrix> pMat_Devi_inv(pBigMat_Devi_inv);
	arma::Mat<double> Devi_inv((double*) pMat_Devi_inv->matrix(), pMat_Devi_inv -> nrow(), pMat_Devi_inv -> ncol(), false);		
	Devi.diag()=Devi.diag()+0.0001;	
	Devi_inv.submat(span(0,n_ind-1),span(0,n_ind-1)) =arma::inv(Devi);	
	Devi.diag()=Devi.diag()-0.0001;	
	D_inv_bigmemory_object=make_bigmemory_object_cpp(n_ind,n_ind,bigmemory_data_name+"_G_Dinv",bigmemory_data_path,"double");
Rcout<<"bigmemory-Complete constructing inverse-genomic dominance(classical) relationship matrix!"<<endl; 	
	}
	
	//删除临时的 temp_numeric.desc 和  .bk文件
	std::string temp_name1=bigmemory_data_path+"/"+bigmemory_data_name+"_temp_numeric.desc";
	std::string temp_name2=bigmemory_data_path+"/"+bigmemory_data_name+"_temp_numeric.bk";
	remove(temp_name1.c_str());
	remove(temp_name2.c_str());	
	
	return List::create(Named("D") = D_bigmemory_object,_["Dinv"] = D_inv_bigmemory_object);	
}


//计算显性效应矩阵
// [[Rcpp::export]]
List  Dominance_matrix_memory_cpp(SEXP pBigMat,std::string bigmemory_data_name, std::string bigmemory_data_path,
							bool base=false,bool trace=false, bool inverse=false){


	Rcpp::XPtr<BigMatrix> pMat(pBigMat); 
	arma::Mat<double> M_temp((double*) pMat->matrix(), pMat -> nrow(), pMat -> ncol(), false); 
	int n_ind=M_temp.n_rows; // 个体数
	int n_snp=M_temp.n_cols; //SNP数量

delete_bigmemory_file_cpp("temp_numeric",bigmemory_data_name,bigmemory_data_path,false);
	
	SEXP pBigMat_temp=make_bigmemory_object_address_cpp(pMat -> nrow(),pMat -> ncol(),bigmemory_data_name+"_temp_numeric",bigmemory_data_path,"double");
	Rcpp::XPtr<BigMatrix> pMat_temp(pBigMat_temp);
	arma::Mat<double> M1((double*) pMat_temp->matrix(), pMat_temp -> nrow(), pMat_temp -> ncol(), false);		
	M1.submat(span(0,n_ind-1),span(0,n_snp-1))=M_temp;

delete_bigmemory_file_cpp("G_D",bigmemory_data_name,bigmemory_data_path,true);
	
	SEXP pBigMat_D=make_bigmemory_object_address_cpp(n_ind,n_ind,bigmemory_data_name+"_G_D",bigmemory_data_path,"double");
	Rcpp::XPtr<BigMatrix> pMat_D(pBigMat_D);
	arma::Mat<double> Domi((double*) pMat_D->matrix(), pMat_D -> nrow(), pMat_D -> ncol(), false);
	
	Rcout<<"bigmemory-Start constructing genomic dominance(genotypic) relationship matrix......"<<endl;  

	double p_freq,q_freq;
	arma::Row<double> p_frequency(n_snp), q_frequency(n_snp);

	if(base==TRUE){
	p_frequency.fill(0.5);
	q_frequency.fill(0.5);
	for( int i=0; i < n_snp; i++){
	for( int j=0; j < n_ind; j++){
	if(M1(j,i)==2){M1(j,i)=-2*0.25;
	}else if(M1(j,i)==1){M1(j,i)=1-2*0.25;
	}else{M1(j,i)=-2*0.25;}	
	}}
    }else{
	
	for( int i=0; i < n_snp; i++){
	p_freq=0.5*mean(M1.col(i));
	q_freq=1-p_freq;
	for( int j=0; j < n_ind; j++){
	if(M1(j,i)==2){M1(j,i)=-2*p_freq*q_freq;
	}else if(M1(j,i)==1){M1(j,i)=1-2*p_freq*q_freq;
	}else{M1(j,i)=-2*p_freq*q_freq;}
	p_frequency[i]=p_freq;
	}
	}
	 q_frequency=1-p_frequency;
	}


	arma::Row<double> a=2*p_frequency%q_frequency;
	Domi=M1*M1.t()/sum(a%(1-a));
Rcout<<"bigmemory-Complete constructing genomic dominance(genotypic) relationship matrix!"<<endl; 

	SEXP D_bigmemory_object=make_bigmemory_object_cpp(n_ind,n_ind,bigmemory_data_name+"_G_D",bigmemory_data_path,"double");
	SEXP D_inv_bigmemory_object=R_NilValue;


	if(inverse==TRUE){
		
delete_bigmemory_file_cpp("G_Dinv",bigmemory_data_name,bigmemory_data_path,true);	
	
	SEXP pBigMat_Domi_inv=make_bigmemory_object_address_cpp(n_ind,n_ind,bigmemory_data_name+"_G_Dinv",bigmemory_data_path,"double");
	Rcpp::XPtr<BigMatrix> pMat_Domi_inv(pBigMat_Domi_inv);
	arma::Mat<double> Domi_inv((double*) pMat_Domi_inv->matrix(), pMat_Domi_inv -> nrow(), pMat_Domi_inv -> ncol(), false);		
	Domi.diag()=Domi.diag()+0.0001;	
	Domi_inv.submat(span(0,n_ind-1),span(0,n_ind-1)) =arma::inv(Domi);	
	Domi.diag()=Domi.diag()-0.0001;	
	D_inv_bigmemory_object=make_bigmemory_object_cpp(n_ind,n_ind,bigmemory_data_name+"_G_Dinv",bigmemory_data_path,"double");
Rcout<<"bigmemory-Complete constructing inverse-genomic dominance(genotypic) relationship matrix!"<<endl; 	
	}
	
	//删除临时的 temp_numeric.desc 和  .bk文件
	std::string temp_name1=bigmemory_data_path+"/"+bigmemory_data_name+"_temp_numeric.desc";
	std::string temp_name2=bigmemory_data_path+"/"+bigmemory_data_name+"_temp_numeric.bk";
	remove(temp_name1.c_str());
	remove(temp_name2.c_str());		
	
	return List::create(Named("D") = D_bigmemory_object,_["Dinv"] = D_inv_bigmemory_object);
}	  



// [[Rcpp::export]]
void APY_memory_cpp(SEXP pBigMat_G,SEXP pBigMat_Ginv,CharacterVector IND_Proven, CharacterVector IND_G){ //G为所有个体的亲缘关系矩阵
																				   //APY算法	

	Rcpp::XPtr<BigMatrix> pMat_G(pBigMat_G);
	arma::Mat<double> G((double*) pMat_G->matrix(), pMat_G -> nrow(), pMat_G -> ncol(), false);

	Rcpp::XPtr<BigMatrix> pMat_Ginv(pBigMat_Ginv);
	arma::Mat<double> G_inv((double*) pMat_Ginv->matrix(), pMat_Ginv -> nrow(), pMat_Ginv -> ncol(), false);

																				   
	CharacterVector IND_Young=setdiff(IND_G,IND_Proven);
	IntegerVector IND_Young_pos1=match(IND_Young,IND_G)-1;
	IntegerVector IND_Proven_pos1=match(IND_Proven,IND_G)-1;
	uvec IND_Young_pos=as<arma::uvec>(IND_Young_pos1);
	uvec IND_Proven_pos=as<arma::uvec>(IND_Proven_pos1);
	int IND_Young_n=IND_Young_pos.n_rows;
	int IND_Proven_n=IND_Proven_pos.n_rows;
	int n_ind=G.n_rows;
	arma::Mat<double> Gpp_inv=inv(G(IND_Proven_pos,IND_Proven_pos));
	arma::Mat<double> P_young=G.submat(IND_Young_pos,IND_Proven_pos)*Gpp_inv;
	arma::Mat<double> G_ip=G.submat(IND_Young_pos,IND_Proven_pos);
	arma::Col<double> tmp(IND_Young_n);
	for ( int i=0;i<IND_Young_n; ++i){
		arma::Mat<double> tmp1=G_ip.row(i);
	    arma::Mat<double> tmp2=tmp1*Gpp_inv*tmp1.t();
		tmp[i]=tmp2(0,0);	
	}
	arma::Col<double> gii=G.diag();
	gii=gii(IND_Young_pos);
	arma::Col<double> M_inv=1/(gii-tmp);
	arma::Mat<double> Py_t_M_inv=P_young.t()*diagmat(M_inv);
	//for(int i=0;i<IND_Young_n;i++){
	//	Py_t_M_inv(i,i)=Py_t_M_inv(i,i)*M_inv(i); // Py_t_M_inv=t(P_young)*M_inv	
	//}
	
	//mat G_inv(n_ind,n_ind,fill::zeros); 由于不能主动释放G，因此为了节省空间，用G存储G_inv
	G_inv.submat(span(0,IND_Proven_n-1),span(0,IND_Proven_n-1))=Py_t_M_inv*P_young+Gpp_inv;
	G_inv.submat(span(0,IND_Proven_n-1),span(IND_Proven_n,n_ind-1))=-1*Py_t_M_inv;
	G_inv.submat(span(IND_Proven_n,n_ind-1),span(0,IND_Proven_n-1))=-1*Py_t_M_inv.t();
	G_inv.submat(span(IND_Proven_n,n_ind-1),span(IND_Proven_n,n_ind-1))=diagmat(M_inv);
	//reorder the position, 
	CharacterVector IND_G_final=union_cpp(IND_Proven,IND_Young);
	IntegerVector reorder_pos1=match(IND_G,IND_G_final)-1;
	uvec reorder_pos=as<arma::uvec>(reorder_pos1);	
	G_inv.submat(span(0,n_ind-1),span(0,n_ind-1))=G_inv.submat(reorder_pos,reorder_pos);
	
}

// [[Rcpp::export]]
List APY_inverse_memory_cpp(SEXP pBigMat,std::string bigmemory_data_name, std::string bigmemory_data_path,
					        CharacterVector IND_geno, 
					        double APY_eigen_threshold=0.95,
					        int APY_n_core=0,
					        bool re_inverse=false){

Rcout<<"bigmemory-Start constructing inverse-genomic additive relationship matrix(APY method)......"<<endl;
List tmp=G_matrix_memory_cpp(pBigMat,bigmemory_data_name,bigmemory_data_path,false,false,false,false);

SEXP pBigMat_G=make_bigmemory_object_address_cpp(10,10,bigmemory_data_name+"_G_A",bigmemory_data_path,"double");
Rcpp::XPtr<BigMatrix> pMat_G(pBigMat_G);
arma::Mat<double> G((double*) pMat_G->matrix(), pMat_G -> nrow(), pMat_G -> ncol(), false);
int n_ind=G.n_rows;

SEXP pBigMat_Ginv=make_bigmemory_object_address_cpp(n_ind,n_ind,bigmemory_data_name+"_G_Ainv",bigmemory_data_path,"double");

int N_selected,i;
if(APY_n_core!=0){
N_selected=APY_n_core;
}else{
arma::Col<double>  eigval = eig_sym( G ),tmp_eigval;
double sum_score=sum(eigval);
for(i=0;i<eigval.size();i++){
	tmp_eigval=eigval.subvec(0,i);
	if((sum(tmp_eigval)/sum_score)>=APY_eigen_threshold){
		break;
	}	
}
N_selected=i;
}
Rcout<<N_selected<<endl;
CharacterVector IND_Proven=sample(IND_geno,N_selected,false);

APY_memory_cpp(pBigMat_G,pBigMat_Ginv,IND_Proven,IND_geno);
SEXP G_inv_bigmemory_object=make_bigmemory_object_cpp(n_ind,n_ind,bigmemory_data_name+"_G_Ainv",bigmemory_data_path,"double");
SEXP G_bigmemory_object=R_NilValue;

Rcout<<"bigmemory-Complete constructing inverse-genomic additive relationship matrix(APY method)"<<endl;

if(re_inverse==true){
Rcpp::XPtr<BigMatrix> pMat_Ginv(pBigMat_Ginv);
arma::Mat<double> G_inv((double*) pMat_Ginv->matrix(), pMat_Ginv -> nrow(), pMat_Ginv -> ncol(), false);	
	
G_inv.diag()=G_inv.diag()+0.0001;	
G.submat(span(0,n_ind-1),span(0,n_ind-1))=arma::inv(G_inv);	
G_inv.diag()=G_inv.diag()-0.0001;	
G_bigmemory_object=make_bigmemory_object_cpp(n_ind,n_ind,bigmemory_data_name+"_G_A",bigmemory_data_path,"double");
}else{
	
//删除临时的 G_A.desc 和  .bk文件
std::string temp_name1=bigmemory_data_path+"/"+bigmemory_data_name+"_G_A.desc";
std::string temp_name2=bigmemory_data_path+"/"+bigmemory_data_name+"_G_A.bk";
remove(temp_name1.c_str());
remove(temp_name2.c_str());	

Rcout<<"bigmemory-Complete constructing inverse-(inverse-genomic additive relationship matrix(APY method))!"<<endl;
	
}	

return List::create(Named("G") = G_bigmemory_object,_["Ginv"] = G_inv_bigmemory_object);
 
}
	


// [[Rcpp::export]]
List makeA_memory_cpp(arma::Mat<int> Pedigree,
					  std::string bigmemory_data_name,
					  std::string bigmemory_data_path){   //稀疏矩阵运行速度会变慢

delete_bigmemory_file_cpp("P_A",bigmemory_data_name,bigmemory_data_path,true);

	SEXP pBigMat_A=make_bigmemory_object_address_cpp(Pedigree.n_rows,Pedigree.n_rows,
	                                                  bigmemory_data_name+"_P_A",
													  bigmemory_data_path,"double");
	Rcpp::XPtr<BigMatrix> pMat_A(pBigMat_A);
	arma::Mat<double> A_matrix((double*) pMat_A->matrix(), pMat_A -> nrow(), pMat_A -> ncol(), false);
	A_matrix.fill(0);
	arma::Col<int> Animal=Pedigree.col(0);
	arma::Col<int> Sire=Pedigree.col(1);
	arma::Col<int> Dam=Pedigree.col(2);
	int n_ind=Pedigree.n_rows;
	Rcout<<"bigmemory-Start constructing pedigree additive relationship matrix......"<<endl;  
	
	
	bool status1,status2;
	
	for(int i=0; i < n_ind; i++){

	status1=(Sire[i]==0);
	status2=(Dam[i]==0);


	if(!status1&&!status2){  //父母均已知
		A_matrix(i,i)=1+0.5*(A_matrix(Sire[i]-1,Dam[i]-1));	
		for(int j=0;j<i;j++){
			if(i>j) A_matrix(j,i)=A_matrix(i,j)=0.5*(A_matrix(j,Sire[i]-1)+A_matrix(j,Dam[i]-1));
	    }
	}else if(status1&&status2){  //父母均未知
		A_matrix(i,i)=1;	
		for(int j=0;j<i;j++){
			if(i>j) A_matrix(j,i)=0;
	    }
	}else if(!status1&&status2){  //父已知，母未知	
		A_matrix(i,i)=1;
		for(int j=0;j<i;j++){
			if(i>j) A_matrix(j,i)=A_matrix(i,j)=0.5*(A_matrix(j,Sire[i]-1));
		}
	}else{  //父未知，母已知	
		A_matrix(i,i)=1;
		for(int j=0;j<i;j++){
			if(i>j) A_matrix(j,i)=A_matrix(i,j)=0.5*(A_matrix(j,Dam[i]-1));
	    }
	}
	}
Rcout<<"bigmemory-Complete constructing pedigree additive relationship matrix!"<<endl; 
SEXP P_A_bigmemory_object=make_bigmemory_object_cpp(n_ind,n_ind,bigmemory_data_name+"_P_A",bigmemory_data_path,"double");		
	return List::create(Named("A") = P_A_bigmemory_object);
	}


// [[Rcpp::export]]
List makeAinv_memory_cpp(arma::Mat<int> Pedigree,
						std::string bigmemory_data_name,
					    std::string bigmemory_data_path,bool inbreeding=false){


	int n_ind=Pedigree.n_rows;		
	arma::Col<int> Animal=Pedigree.col(0);
	arma::Col<int> Sire=Pedigree.col(1);
	arma::Col<int> Dam=Pedigree.col(2);

delete_bigmemory_file_cpp("P_Ainv",bigmemory_data_name,bigmemory_data_path,true);
	
	SEXP pBigMat_Ainv=make_bigmemory_object_address_cpp(Pedigree.n_rows,Pedigree.n_rows,
	                                                  bigmemory_data_name+"_P_Ainv",
													  bigmemory_data_path,"double");
	Rcpp::XPtr<BigMatrix> pMat_Ainv(pBigMat_Ainv);
	arma::Mat<double> A_inv((double*) pMat_Ainv->matrix(), pMat_Ainv -> nrow(), pMat_Ainv -> ncol(), false);
	A_inv.fill(0);
		
	double tmp,d;
	Rcout<<"bigmemory-Start constructing inverse-pedigree additive relationship matrix!"<<endl;  

	
if(inbreeding==true){	
	arma::Mat<double> L(n_ind,n_ind);
	arma::Mat<double> A(n_ind,n_ind);
	for(int i=0; i < n_ind; i++){
	
        if(Sire[i] != 0 && Dam[i] != 0) // Both parents known 
        {
             tmp = 0.0;

			for(int j=0; j < i+1; j++){

			if(j!=i){L(i, j) = 0.5 * (L(Sire[i]-1,j) + L(Dam[i]-1,j));			
			}else{L(i, j) = sqrt(1.0-0.25*(A(Sire[i]-1,Sire[i]-1)+A(Dam[i]-1,Dam[i]-1)));}
			tmp=tmp+pow(L(i, j),2);	
		
		   }	
		A(i,i)=tmp;
		d=1.0/(pow(L(i, i),2)+0.0);
		A_inv(i,i)=A_inv(i,i)+d;
		A_inv(i,Dam[i]-1)=A_inv(Dam[i]-1,i)=A_inv(i,Dam[i]-1)-0.5*d;
		A_inv(i,Sire[i]-1)=A_inv(Sire[i]-1,i)=A_inv(i,Sire[i]-1)-0.5*d;
		A_inv(Dam[i]-1,Dam[i]-1)=A_inv(Dam[i]-1,Dam[i]-1)+0.25*d;
		A_inv(Sire[i]-1,Sire[i]-1)=A_inv(Sire[i]-1,Sire[i]-1)+0.25*d;
		A_inv(Sire[i]-1,Dam[i]-1)=A_inv(Sire[i]-1,Dam[i]-1)+0.25*d;
		A_inv(Dam[i]-1,Sire[i]-1)=A_inv(Dam[i]-1,Sire[i]-1)+0.25*d;	
		
	   }
	else if(Sire[i] != 0 && Dam[i] == 0) // Sire know, Dam unknown 
        {
            tmp = 0.0;
			for(int j=0; j < i+1; j++){
			if(j!=i){L(i, j) = 0.5 *L(Sire[i]-1,j);		
			}else{L(i, j) = sqrt(1.0-0.25*A(Sire[i]-1,Sire[i]-1));}
			tmp=tmp+pow(L(i,j),2);
		   }	
 
		A(i,i)=tmp;
	
		d=1.0/(pow(L(i, i),2));

		A_inv(i,i)=A_inv(i,i)+d;
		
		A_inv(i,Sire[i]-1)=A_inv(Sire[i]-1,i)=A_inv(i,Sire[i]-1)-0.5*d;
		A_inv(Sire[i]-1,Sire[i]-1)=A_inv(Sire[i]-1,Sire[i]-1)+0.25*d;
			
	   }		   
	else if(Sire[i] == 0 && Dam[i] != 0) // Sire unknow, Dam known 
        {
             tmp = 0.0;

			for(int j=0; j < i+1; j++){

			if(j!=i){L(i, j) = 0.5 *L(Dam[i]-1,j);		
			}else{L(i, j) = sqrt(1.0-0.25*A(Dam[i]-1,Dam[i]-1));}
			tmp=tmp+pow(L(i, j),2);	
		   }	
		
	   	A(i,i)=tmp;
		d=1.0/(pow(L(i, i),2)+0.0);
		
		A_inv(i,i)=A_inv(i,i)+d;
		A_inv(i,Dam[i]-1)=A_inv(i,Dam[i]-1)-0.5*d;
		A_inv(Dam[i]-1,i)=A_inv(Dam[i]-1,i)-0.5*d;
		A_inv(Dam[i]-1,Dam[i])=A_inv(Dam[i]-1,Dam[i]-1)+0.25*d;
		
	   }	
	  else{ //(tmp_Sire == 0 && tmp_Dam == 0); // Both parents unknown, 不能使用else if, 否则 Sire=1，Dam=0也会进行此循环中
        
            tmp = 0.0;
			for(int j=0; j < i+1; j++){	
			if(j!=i){L(i, j) = 0.0;		
			}else{L(i, j) =1.0;}		
		    tmp=tmp+pow(L(i, j),2);
		    }	
	   	A(i,i)=tmp; // tmp=A_ii
		d=1.0/(pow(L(i,i),2));
		A_inv(i,i)=A_inv(i,i)+d;
	
		}

}
}else{

	for(int i=0; i < n_ind; i++){

        if(Sire[i] == 0 && Dam[i] == 0){  // Both parents unknown 
	    d=1.0;
		A_inv(i,i)=A_inv(i,i)+d;
	
		}	

        if(Sire[i] != 0 && Dam[i] == 0) { // Sire know, Dam unknown    

		 d=4.0/3.0;
		A_inv(i,i)=A_inv(i,i)+d;	
		A_inv(i,Sire[i]-1)=A_inv(Sire[i]-1,i)=A_inv(i,Sire[i]-1)-0.5*d;
		A_inv(Sire[i]-1,Sire[i]-1)=A_inv(Sire[i]-1,Sire[i]-1)+0.25*d;

		}		   
        if(Sire[i] == 0 && Dam[i] != 0){  // Sire unknow, Dam known     
		d=4.0/3.0;
		A_inv(i,i)=A_inv(i,i)+d;
		A_inv(i,Dam[i]-1)=A_inv(Dam[i]-1,i)=A_inv(i,Dam[i]-1)-0.5*d;
		A_inv(Dam[i]-1,Dam[i]-1)=A_inv(Dam[i]-1,Dam[i]-1)+0.25*d	;

	}	
	
        if(Sire[i] != 0 && Dam[i] != 0){ // Both parents known 
		d=2.0;
		A_inv(i,i)=A_inv(i,i)+d;		
		A_inv(i,Dam[i]-1)=A_inv(Dam[i]-1,i)=A_inv(i,Dam[i]-1)-0.5*d;
		A_inv(i,Sire[i]-1)=A_inv(Sire[i]-1,i)=A_inv(i,Sire[i]-1)-0.5*d;
		A_inv(Dam[i]-1,Dam[i]-1)=A_inv(Dam[i]-1,Dam[i]-1)+0.25*d;
		A_inv(Sire[i]-1,Sire[i]-1)=A_inv(Sire[i]-1,Sire[i]-1)+0.25*d;
		A_inv(Sire[i]-1,Dam[i]-1)=A_inv(Sire[i]-1,Dam[i]-1)+0.25*d;
		A_inv(Dam[i]-1,Sire[i]-1)=A_inv(Dam[i]-1,Sire[i]-1)+0.25*d;

		}
	   }
}

Rcout<<"bigmemory-Complete constructing inverse-pedigree additive relationship matrix!"<<endl; 
SEXP P_Ainv_bigmemory_object=make_bigmemory_object_cpp(n_ind,n_ind,bigmemory_data_name+"_P_Ainv",bigmemory_data_path,"double");		
	return List::create(Named("Ainv") = P_Ainv_bigmemory_object);
}	


// [[Rcpp::export]]
DataFrame cal_homo_inbred_memory_cpp(SEXP pBigMat, CharacterVector IND_geno){

	Rcpp::XPtr<BigMatrix> pMat(pBigMat); 
	arma::Mat<double> data_genumeric((double*) pMat->matrix(), pMat -> nrow(), pMat -> ncol(), false); 
	
	int n_ind=data_genumeric.n_rows,SNP_n=data_genumeric.n_cols;
	NumericVector hetezygous(n_ind);
	
	for(int i=0; i < n_ind; i++){
	
		for(int j=0; j < SNP_n; j++){
        
			if(data_genumeric(i,j)==1){			
				hetezygous[i]+= 1;
			}
		}
	}
	hetezygous=1-(hetezygous/SNP_n);
		     
	 return DataFrame::create( Named("Id") = IND_geno , _["Homo_inbred"] = hetezygous );
}


// [[Rcpp::export]]
arma::Mat<double> makeInbreeding_memory_cpp(arma::Mat<int> Pedigree){
	int n_ind=Pedigree.n_rows;	
	arma::Col<int> Animal=Pedigree.col(0);
	arma::Col<int> Sire=Pedigree.col(1);
	arma::Col<int> Dam=Pedigree.col(2);
	arma::Mat<double> L(n_ind,n_ind);
	arma::Mat<double> A(n_ind,n_ind);	

	double tmp;

	for(int i=0; i < n_ind; i++){
	
        if(Sire[i] != 0 && Dam[i]  != 0) // Both parents known 
        {
             tmp = 0.0;

			for(int j=0; j < i+1; j++){

			if(j!=i){L(i, j) = 0.5 * (L(Sire[i]-1,j) + L(Dam[i]-1,j));			
			}else{L(i, j) = sqrt(1.0-0.25*(A(Sire[i]-1,Sire[i]-1)+A(Dam[i]-1,Dam[i]-1)));}
			tmp=tmp+pow(L(i, j),2);	
		
		   }	
		A(i,i)=tmp;
				
	   }
	else if(Sire[i] != 0 && Dam[i] == 0) // Sire know, Dam unknown 
        {
            tmp = 0.0;
			for(int j=0; j < i+1; j++){
			if(j!=i){L(i, j) = 0.5 *L(Sire[i]-1,j);		
			}else{L(i, j) = sqrt(1.0-0.25*A(Sire[i]-1,Sire[i]-1));}
			tmp=tmp+pow(L(i,j),2);
		   }	
 
		A(i,i)=tmp;
				
	   }		   
	else if(Sire[i] == 0 && Dam[i]  != 0) // Sire unknow, Dam known 
        {
             tmp = 0.0;

			for(int j=0; j < i+1; j++){

			if(j!=i){L(i, j) = 0.5 *L(Dam[i]-1,j);		
			}else{L(i, j) = sqrt(1.0-0.25*A(Dam[i]-1,Dam[i]-1));}
			tmp=tmp+pow(L(i, j),2);	
		   }	
		
	   	A(i,i)=tmp;
			
	   }	
	  else{ //(tmp_Sire == 0 && tmp_Dam == 0); // Both parents unknown, 不能使用else if, 否则 Sire=1，Dam=0也会进行此循环中
        
            tmp = 0.0;
			for(int j=0; j < i+1; j++){	
			if(j!=i){L(i, j) = 0.0;		
			}else{L(i, j) =1.0;}		
		    tmp=tmp+pow(L(i, j),2);
		    }	
	   	A(i,i)=tmp; // tmp=A_ii
		}

}
	 return A; 
}



// [[Rcpp::export]]
List makeD_memory_cpp(arma::Mat<int> Pedigree,std::string bigmemory_data_name,
					 std::string bigmemory_data_path, bool inverse=false){   //稀疏矩阵运行速度会变慢

	int Sire_i_pos,Sire_j_pos,Dam_i_pos,Dam_j_pos,n_ind=Pedigree.n_rows;
	delete_bigmemory_file_cpp("P_D",bigmemory_data_name,bigmemory_data_path,true);
	delete_bigmemory_file_cpp("P_A_temp",bigmemory_data_name,bigmemory_data_path,true);
		
	SEXP pBigMat_D=make_bigmemory_object_address_cpp(n_ind,n_ind,
	                                                 bigmemory_data_name+"_P_D",
													 bigmemory_data_path,"double");
	Rcpp::XPtr<BigMatrix> pMat_D(pBigMat_D);
	arma::Mat<double> D_matrix((double*) pMat_D->matrix(), pMat_D -> nrow(), pMat_D -> ncol(), false);
	D_matrix.fill(0);


	SEXP pBigMat_A=make_bigmemory_object_address_cpp(n_ind,n_ind,
	                                                  bigmemory_data_name+"_P_A_temp",
													  bigmemory_data_path,"double");
	Rcpp::XPtr<BigMatrix> pMat_A(pBigMat_A);
	arma::Mat<double> A_matrix((double*) pMat_A->matrix(), pMat_A -> nrow(), pMat_A -> ncol(), false);
	A_matrix.fill(0);
	
	arma::Col<int> Animal=Pedigree.col(0);
	arma::Col<int> Sire=Pedigree.col(1);
	arma::Col<int> Dam=Pedigree.col(2);

	Rcout<<"bigmemory-Start constructing pedigree dominance relationship matrix......"<<endl;  	

	for(int i=0; i < n_ind; i++){

	if(Sire[i]==0&&Dam[i]==0){  //父母均未知
		A_matrix(i,i)=1;	
		D_matrix(i,i)=1;
		for(int j=0;j<i;j++){
			if(i>j) {
				A_matrix(j,i)=0;
	    }
	}
	}

	if(Sire[i]!=0&&Dam[i]!=0){  //个体i的父母均已知
		Sire_i_pos=Sire[i]-1;	
		Dam_i_pos=Dam[i]-1;		
		A_matrix(i,i)=1+0.5*(A_matrix(Sire[i]-1,Dam[i]-1));	
		//D_matrix(i,i)=1;
		D_matrix(i,i)=2-A_matrix(i,i);// nadiv对角线的计算方式
		for(int j=0;j<i;j++){
		Sire_j_pos=Sire[j]-1;
		Dam_j_pos=Dam[j]-1;
			if(i>j) {
				A_matrix(j,i)=A_matrix(i,j)=0.5*(A_matrix(j,Sire_i_pos)+A_matrix(j,Dam_i_pos));
			if(Sire[j]!=0&&Dam[j]!=0){
						D_matrix(j,i)=D_matrix(i,j)=0.25*(A_matrix(Sire_j_pos,Dam_i_pos)*A_matrix(Dam_j_pos,Sire_i_pos)+A_matrix(Sire_j_pos,Sire_i_pos)*A_matrix(Dam_j_pos,Dam_i_pos));  
				}
			}
	    }
	}

	if(Sire[i]!=0&&Dam[i]==0){  //父已知，母未知	
		Sire_i_pos=Sire[i]-1;	
		A_matrix(i,i)=1;
		D_matrix(i,i)=1;  			
		for(int j=0;j<i;j++){
			if(i>j) A_matrix(j,i)=A_matrix(i,j)=0.5*(A_matrix(j,Sire[i]-1));
		}
	}
		
	if(Sire[i]==0&&Dam[i]!=0){  //父未知，母已知	
		A_matrix(i,i)=1;
		D_matrix(i,i)=1;		
		for(int j=0;j<i;j++){
			if(i>j) A_matrix(j,i)=A_matrix(i,j)=0.5*(A_matrix(j,Dam[i]-1));
	    }
	}
	}

delete_bigmemory_file_cpp("P_A_temp",bigmemory_data_name,bigmemory_data_path,false);

	SEXP P_D_bigmemory_object=make_bigmemory_object_cpp(n_ind,n_ind,bigmemory_data_name+"_P_D",bigmemory_data_path,"double");	
	SEXP P_D_inv_bigmemory_object=R_NilValue;
Rcout<<"bigmemory-Complete constructing pedigree dominance relationship matrix!"<<endl;  	
if(inverse==true){
	
delete_bigmemory_file_cpp("P_Dinv",bigmemory_data_name,bigmemory_data_path,true);

	SEXP pBigMat_Dinv=make_bigmemory_object_address_cpp(n_ind,n_ind,bigmemory_data_name+"_P_Dinv",bigmemory_data_path,"double");
	Rcpp::XPtr<BigMatrix> pMat_Dinv(pBigMat_Dinv);
	arma::Mat<double> Dinv_matrix((double*) pMat_Dinv->matrix(), pMat_Dinv -> nrow(), pMat_Dinv -> ncol(), false);
	
	D_matrix.diag()=D_matrix.diag()+0.0001;	
	Dinv_matrix.submat(span(0,n_ind-1),span(0,n_ind-1)) =arma::inv(D_matrix);	
	D_matrix.diag()=D_matrix.diag()-0.0001;	
	P_D_inv_bigmemory_object=make_bigmemory_object_cpp(n_ind,n_ind,bigmemory_data_name+"_P_Dinv",bigmemory_data_path,"double");
	
Rcout<<"bigmemory-Complete constructing inverse-pedigree dominance relationship matrix!"<<endl; 				
	}
	return List::create(Named("D") = P_D_bigmemory_object,_["Dinv"] = P_D_inv_bigmemory_object);
	}




// openMP, not very accurate
// [[Rcpp::export]]
List gene_dropping_memory_D(arma::Mat<int> Pedigree,
							std::string bigmemory_data_name, std::string bigmemory_data_path,
                            bool inverse=false, 
							int iteration=1000, bool diagnoal=true, int cpu_cores=1){
	
     omp_set_num_threads(cpu_cores);
//父母只要有任何一方缺失，显性亲缘关系就为0

	
	int n_ind=Pedigree.n_rows;	
	arma::Mat<int> allele(n_ind,3);

	delete_bigmemory_file_cpp("P_D",bigmemory_data_name,bigmemory_data_path,true);
		
	SEXP pBigMat_D=make_bigmemory_object_address_cpp(n_ind,n_ind,bigmemory_data_name+"_P_D",bigmemory_data_path,"double");
	Rcpp::XPtr<BigMatrix> pMat_D(pBigMat_D);
	arma::Mat<double> D_sim((double*) pMat_D->matrix(), pMat_D -> nrow(), pMat_D -> ncol(), false);
	D_sim.fill(0);	

	arma::Col<int> Sire=Pedigree.col(1);
	arma::Col<int> Dam=Pedigree.col(2);
	arma::Col<int> Sire_allele,Dam_allele;
	int Sire_allele_j,Sire_allele_k,Dam_allele_j,Dam_allele_k,i,j,k;
	bool X,Y;
	arma::uvec founder_s=find(Sire==0);  //父亲为0的个体
	arma::uvec founder_d=find(Dam==0);   //母亲为0的个体
	arma::uvec founder=intersect(founder_s,founder_d); //父母同时为0的个体的行
   	int  pos=founder.size(); //父母同时为0的个体所占的时间
    //int test_iteration=0;
	Rcout<<"bigmemory-Start constructing pedigree dominance relationship matrix(gene dropping)......"<<endl;  
	Rcout<<"The max iteration of gene_dropping is set as: "<<iteration<<endl;	
#pragma omp parallel for private(allele,Sire_allele,Dam_allele,Sire_allele_j,Dam_allele_j,Sire_allele_k,Dam_allele_k,X,Y)	
	for( i=0;i<iteration;i++){
	allele=get_allele(Pedigree);
	 Sire_allele=allele.col(1);
	 Dam_allele=allele.col(2);	
		for( j=pos;j<n_ind;j++){
			 Sire_allele_j=Sire_allele[j];
			 Dam_allele_j=Dam_allele[j];
			for( k=pos;k<j;k++){
				
			 Sire_allele_k=Sire_allele[k];
			 Dam_allele_k=Dam_allele[k];		
			 X=(Sire_allele_j==Sire_allele_k)&(Dam_allele_j==Dam_allele_k);
			 Y=(Sire_allele_j==Dam_allele_k)&(Dam_allele_j==Sire_allele_k);
			if(X|Y){
				//test_iteration =test_iteration+1;	test_iteration 是唯一的，按顺序输出的

				D_sim(j,k) =D_sim(k,j)=D_sim(j,k)+1;
			}		
	}}}
	D_sim.cols(founder_s).fill(0);
	D_sim.rows(founder_s).fill(0);
	D_sim.cols(founder_d).fill(0);
	D_sim.rows(founder_d).fill(0);
	D_sim=D_sim/iteration;
	if(diagnoal==true){  //将D的对角线转为换 (2-A的对角线)	
	//Rcout<<"Calculate Diagnoal Of Dominance Relationship Matrix"<<endl;
	arma::Mat<double> A=makeA_tmp_cpp(Pedigree);
	D_sim.diag()=2-A.diag();
	}

	SEXP P_D_bigmemory_object=make_bigmemory_object_cpp(n_ind,n_ind,bigmemory_data_name+"_P_D",bigmemory_data_path,"double");	
	SEXP P_D_inv_bigmemory_object=R_NilValue;	
	Rcout<<"bigmemory-Complete constructing pedigree dominance relationship matrix(gene dropping)!"<<endl;  	

	if(inverse==true){

	delete_bigmemory_file_cpp("P_Dinv",bigmemory_data_name,bigmemory_data_path,true);

	SEXP pBigMat_Dinv=make_bigmemory_object_address_cpp(n_ind,n_ind,bigmemory_data_name+"_P_Dinv",bigmemory_data_path,"double");
	Rcpp::XPtr<BigMatrix> pMat_Dinv(pBigMat_Dinv);
	arma::Mat<double> Dinv_matrix((double*) pMat_Dinv->matrix(), pMat_Dinv -> nrow(), pMat_Dinv -> ncol(), false);		
		
		D_sim.diag()=D_sim.diag()+0.0001;
		Dinv_matrix.submat(span(0,n_ind-1),span(0,n_ind-1))=arma::inv(D_sim);
		D_sim.diag()=D_sim.diag()-0.0001;	
	
	P_D_inv_bigmemory_object=make_bigmemory_object_cpp(n_ind,n_ind,bigmemory_data_name+"_P_Dinv",bigmemory_data_path,"double");
	
	Rcout<<"bigmemory-Complete constructing inverse-pedigree dominance relationship matrix(gene dropping)!"<<endl;			
	}
	return List::create(Named("D") = P_D_bigmemory_object,_["Dinv"] = P_D_inv_bigmemory_object);
	
}
	



// [[Rcpp::export]]
List makeHA_memory_cpp(arma::Mat<int> & Pedigree, SEXP pBigMat, 
					   CharacterVector IND_geno,
				       std::string bigmemory_data_name, 
					   std::string bigmemory_data_path,
						   arma::uvec pos_A11,
						   arma::uvec pos_A22,
						   arma::uvec pos_geno,
						   arma::uvec pos_A,
						   arma::uvec pos_H22,
						  bool APY_algorithm=false,
						  double APY_eigen_threshold=0.95,
						  int APY_n_core=0,
						   bool direct=false,
						   bool inverse=true,
						   double omega=0.05){	
						   

//获取M
	Rcpp::XPtr<BigMatrix> pMat(pBigMat); 
	arma::Mat<double> M((double*) pMat->matrix(), pMat -> nrow(), pMat -> ncol(), false); 


	int n_ind1=pos_A11.size(),n_ind2=pos_A22.size(),n_ind=pos_A11.size()+pos_A22.size(),n_ind_geno=pos_geno.size();
	int n_snp=M.n_cols,n_pedigree=Pedigree.n_rows;

//创建临时的P_A	
	delete_bigmemory_file_cpp("temp_P_A",bigmemory_data_name,bigmemory_data_path,false);
SEXP tmp1=makeA_memory_cpp(Pedigree,bigmemory_data_name+"_temp",bigmemory_data_path);
	
	SEXP pBigMat_A=make_bigmemory_object_address_cpp(n_pedigree,n_pedigree,bigmemory_data_name+"_temp_P_A",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_A(pBigMat_A);
	arma::Mat<double> P_A((double*) pMat_A->matrix(), pMat_A -> nrow(), pMat_A -> ncol(), false);

//创建临时的 A22
	delete_bigmemory_file_cpp("temp_A22",bigmemory_data_name,bigmemory_data_path,false);
	SEXP pBigMat_A22=make_bigmemory_object_address_cpp(n_ind2,n_ind2,bigmemory_data_name+"_temp_A22",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_A22(pBigMat_A22);
	arma::Mat<double> A22((double*) pMat_A22->matrix(), pMat_A22 -> nrow(), pMat_A22 -> ncol(), false);

//创建临时的 A22_inv
	delete_bigmemory_file_cpp("temp_A22_inv",bigmemory_data_name,bigmemory_data_path,false);
	SEXP pBigMat_A22_inv=make_bigmemory_object_address_cpp(n_ind2,n_ind2,bigmemory_data_name+"_temp_A22_inv",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_A22_inv(pBigMat_A22_inv);
	arma::Mat<double> A22_inv((double*) pMat_A22_inv->matrix(), pMat_A22_inv -> nrow(), pMat_A22_inv -> ncol(), false);

//创建临时的 M_new
	delete_bigmemory_file_cpp("temp_M_new",bigmemory_data_name,bigmemory_data_path,false);
	SEXP pBigMat_M_new=make_bigmemory_object_address_cpp(n_ind_geno,n_snp,bigmemory_data_name+"_temp_M_new",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_M_new(pBigMat_M_new);
	arma::Mat<double> M_new((double*) pMat_M_new->matrix(), pMat_M_new -> nrow(), pMat_M_new -> ncol(), false);
	
	M_new.submat(span(0,n_ind_geno-1),span(0,n_snp-1))=M.rows(pos_geno);


	SEXP H_A_bigmemory_object=R_NilValue;
	SEXP H_A_inv_bigmemory_object=R_NilValue;
	arma::Col<double> inbred;	
	
	if(APY_algorithm==false){

//创建临时的 A11
	delete_bigmemory_file_cpp("temp_A11",bigmemory_data_name,bigmemory_data_path,false);
	SEXP pBigMat_A11=make_bigmemory_object_address_cpp(n_ind1,n_ind1,bigmemory_data_name+"_temp_A11",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_A11(pBigMat_A11);
	arma::Mat<double> A11((double*) pMat_A11->matrix(), pMat_A11 -> nrow(), pMat_A11 -> ncol(), false);


//创建临时的 A12
	delete_bigmemory_file_cpp("temp_A12",bigmemory_data_name,bigmemory_data_path,false);
	SEXP pBigMat_A12=make_bigmemory_object_address_cpp(n_ind1,n_ind2,bigmemory_data_name+"_temp_A12",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_A12(pBigMat_A12);
	arma::Mat<double> A12((double*) pMat_A12->matrix(), pMat_A12 -> nrow(), pMat_A12 -> ncol(), false);


//创建临时的G_A	
	delete_bigmemory_file_cpp("temp_G_A",bigmemory_data_name,bigmemory_data_path,false);
List tmp2=G_matrix_memory_cpp(pBigMat_M_new,bigmemory_data_name+"_temp",bigmemory_data_path,false,false,false,false);
	
	SEXP pBigMat_G_A=make_bigmemory_object_address_cpp(n_ind_geno,n_ind_geno,bigmemory_data_name+"_temp_G_A",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_G_A(pBigMat_G_A);
	arma::Mat<double> G_A((double*) pMat_G_A->matrix(), pMat_G_A -> nrow(), pMat_G_A -> ncol(), false);

				Rcout<<"bigmemory-Start constructing H_Additive relationship matrix......"<<endl; 						   
				A11.submat(span(0,n_ind1-1),span(0,n_ind1-1))=P_A.submat(pos_A11,pos_A11);
				A22.submat(span(0,n_ind2-1),span(0,n_ind2-1))=P_A.submat(pos_A22,pos_A22);
				A22_inv.submat(span(0,n_ind2-1),span(0,n_ind2-1))=arma::inv(A22);
				A12.submat(span(0,n_ind1-1),span(0,n_ind2-1))=P_A.submat(pos_A11,pos_A22);
				
                double diag_GA=mean(G_A.diag());
				double offdiag_GA=mean(G_A.elem(trimatl_ind(size(G_A),-1)));
                Rcout<<"The mean of diagnal and off-diagnal of G_A is : "<<diag_GA<<" and "<<offdiag_GA<<endl;
				
				double diag_A22=mean(A22.diag());
				double offdiag_A22=mean(A22.elem(trimatl_ind(size(A22),-1)));
				Rcout<<"The mean of diagnal and off-diagnal of P_A22 is : "<<diag_A22<<" and "<<offdiag_A22<<endl;

				double beta=(diag_A22-offdiag_A22)/(diag_GA-offdiag_GA);
				double alpha=(offdiag_A22-beta*offdiag_GA);
				
                Rcout<<"Blending G_A to P_A , the adjusting parameter is: alpha= "<<alpha<<" , beta= "<<beta<<endl;
                G_A=alpha+beta*G_A;		
				Rcout<<"Adjusting G_A , the adjusting parameter omega is: omega= "<<omega<<endl;
				G_A=G_A*(1-omega)+A22*omega;
				
				
				
				if(direct==true){
				delete_bigmemory_file_cpp("H_A",bigmemory_data_name,bigmemory_data_path,true);
				SEXP pBigMat_H_A=make_bigmemory_object_address_cpp(n_ind,n_ind,bigmemory_data_name+"_H_A",bigmemory_data_path,"double");	
				Rcpp::XPtr<BigMatrix> pMat_H_A(pBigMat_H_A);
				arma::Mat<double> H_A((double*) pMat_H_A->matrix(), pMat_H_A -> nrow(), pMat_H_A -> ncol(), false);
					
				arma::Mat<double> tmp=A12*A22_inv*G_A;	
			
				H_A.submat(span(0,n_ind-1),span(0,n_ind-1))=P_A.submat(pos_A,pos_A);
				H_A.submat(span(0,n_ind1-1),span(0,n_ind1-1))=A11+A12*A22_inv*(G_A-A22)*A22_inv*(A12.t());
				H_A.submat(span(0,n_ind1-1),span(n_ind1,n_ind-1))=tmp;
				H_A.submat(span(n_ind1,n_ind-1),span(0,n_ind1-1))=tmp.t();
				H_A.submat(span(n_ind1,n_ind-1),span(n_ind1,n_ind-1))=G_A;
				inbred=H_A.diag()-1;
				H_A_bigmemory_object=make_bigmemory_object_cpp(n_ind1,n_ind2,bigmemory_data_name+"_H_A",bigmemory_data_path,"double");	
				Rcout<<"bigmemory-Complete constructing H_Additive relationship matrix!"<<endl; 
				}
								

				if(inverse==true){

				delete_bigmemory_file_cpp("H_Ainv",bigmemory_data_name,bigmemory_data_path,true);
				SEXP pBigMat_H_Ainv=make_bigmemory_object_address_cpp(n_ind,n_ind,bigmemory_data_name+"_H_Ainv",bigmemory_data_path,"double");	
				Rcpp::XPtr<BigMatrix> pMat_H_Ainv(pBigMat_H_Ainv);
				arma::Mat<double> H_Ainv((double*) pMat_H_Ainv->matrix(), pMat_H_Ainv -> nrow(), pMat_H_Ainv -> ncol(), false);
				
//创建临时的 G_Ainv
	delete_bigmemory_file_cpp("temp_G_Ainv",bigmemory_data_name,bigmemory_data_path,false);
	SEXP pBigMat_G_Ainv=make_bigmemory_object_address_cpp(n_ind_geno,n_ind_geno,bigmemory_data_name+"_temp_G_Ainv",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_G_Ainv(pBigMat_G_Ainv);
	arma::Mat<double> G_Ainv((double*) pMat_G_Ainv->matrix(), pMat_G_Ainv -> nrow(), pMat_G_Ainv -> ncol(), false);


//创建临时的 P_Ainv

	delete_bigmemory_file_cpp("temp_P_Ainv",bigmemory_data_name,bigmemory_data_path,false);
SEXP tmp2=makeAinv_memory_cpp(Pedigree,bigmemory_data_name+"_temp",bigmemory_data_path);
	
	SEXP pBigMat_Ainv=make_bigmemory_object_address_cpp(n_pedigree,n_pedigree,bigmemory_data_name+"_temp_P_Ainv",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_Ainv(pBigMat_Ainv);
	arma::Mat<double> P_Ainv((double*) pMat_Ainv->matrix(), pMat_Ainv -> nrow(), pMat_Ainv -> ncol(), false);
					
				G_Ainv.submat(span(0,n_ind_geno-1),span(0,n_ind_geno-1))=arma::inv(G_A);
				
				H_Ainv.submat(span(0,n_ind-1),span(0,n_ind-1))=P_Ainv.submat(pos_A,pos_A);
				H_Ainv.submat(pos_H22,pos_H22)=H_Ainv.submat(pos_H22,pos_H22)+G_Ainv-A22_inv;	
Rcout<<"bigmemory-Complete constructing inverse-H_Additive relationship matrix!"<<endl; 	
				H_A_inv_bigmemory_object=make_bigmemory_object_cpp(n_ind1,n_ind2,bigmemory_data_name+"_H_Ainv",bigmemory_data_path,"double");	
			
				}			
	}else{

				delete_bigmemory_file_cpp("H_Ainv",bigmemory_data_name,bigmemory_data_path,true);
				SEXP pBigMat_H_Ainv=make_bigmemory_object_address_cpp(n_ind,n_ind,bigmemory_data_name+"_H_Ainv",bigmemory_data_path,"double");	
				Rcpp::XPtr<BigMatrix> pMat_H_Ainv(pBigMat_H_Ainv);
				arma::Mat<double> H_Ainv((double*) pMat_H_Ainv->matrix(), pMat_H_Ainv -> nrow(), pMat_H_Ainv -> ncol(), false);

//创建临时的 P_Ainv

	delete_bigmemory_file_cpp("temp_P_Ainv",bigmemory_data_name,bigmemory_data_path,false);
SEXP tmp2=makeAinv_memory_cpp(Pedigree,bigmemory_data_name+"_temp",bigmemory_data_path);
	
	SEXP pBigMat_Ainv=make_bigmemory_object_address_cpp(n_pedigree,n_pedigree,bigmemory_data_name+"_temp_P_Ainv",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_Ainv(pBigMat_Ainv);
	arma::Mat<double> P_Ainv((double*) pMat_Ainv->matrix(), pMat_Ainv -> nrow(), pMat_Ainv -> ncol(), false);
			

				delete_bigmemory_file_cpp("temp_G_Ainv",bigmemory_data_name,bigmemory_data_path,false);
				APY_inverse_memory_cpp(pBigMat_M_new,bigmemory_data_name+"_temp",bigmemory_data_path,IND_geno,
										APY_eigen_threshold,APY_n_core,false);
					   
				SEXP pBigMat_G_Ainv=make_bigmemory_object_address_cpp(n_ind_geno,n_ind_geno,bigmemory_data_name+"_temp_G_Ainv",
																		bigmemory_data_path,"double");	
				Rcpp::XPtr<BigMatrix> pMat_G_Ainv(pBigMat_G_Ainv);
				arma::Mat<double> G_Ainv((double*) pMat_G_Ainv->matrix(), pMat_G_Ainv -> nrow(), pMat_G_Ainv -> ncol(), false);

				H_Ainv.submat(span(0,n_ind-1),span(0,n_ind-1))=P_Ainv.submat(pos_A,pos_A);
				H_Ainv.submat(pos_H22,pos_H22)=H_Ainv.submat(pos_H22,pos_H22)+G_Ainv-A22_inv;	
Rcout<<"bigmemory-Complete constructing inverse-H_Additive relationship matrix!"<<endl; 	
				H_A_inv_bigmemory_object=make_bigmemory_object_cpp(n_ind1,n_ind2,bigmemory_data_name+"_H_Ainv",bigmemory_data_path,"double");	
			
		
		
	}
delete_bigmemory_file_cpp("temp_P_A",bigmemory_data_name,bigmemory_data_path,false);
delete_bigmemory_file_cpp("temp_A22",bigmemory_data_name,bigmemory_data_path,false);
delete_bigmemory_file_cpp("temp_A22_inv",bigmemory_data_name,bigmemory_data_path,false);
delete_bigmemory_file_cpp("temp_A11",bigmemory_data_name,bigmemory_data_path,false);
delete_bigmemory_file_cpp("temp_A12",bigmemory_data_name,bigmemory_data_path,false);
delete_bigmemory_file_cpp("temp_M_new",bigmemory_data_name,bigmemory_data_path,false);
delete_bigmemory_file_cpp("temp_G_A",bigmemory_data_name,bigmemory_data_path,false);
delete_bigmemory_file_cpp("temp_G_Ainv",bigmemory_data_name,bigmemory_data_path,false);
delete_bigmemory_file_cpp("temp_P_Ainv",bigmemory_data_name,bigmemory_data_path,false);	
	
				return List::create(Named("H") = H_A_bigmemory_object,_["Hinv"] = H_A_inv_bigmemory_object,Named("inbred") =inbred);
}						   
	





// [[Rcpp::export]]
List makeHA_metafounder_memory_cpp(arma::Mat<int> & Pedigree, SEXP pBigMat, 
					               CharacterVector IND_geno,
				                   std::string bigmemory_data_name, 
					               std::string bigmemory_data_path,
						           arma::uvec pos_A11,
						           arma::uvec pos_A22,
						           arma::uvec pos_geno,
						           arma::uvec pos_A,
						           arma::uvec pos_H22,
						           bool direct=false,
						           bool inverse=true,
						           double omega=0.05){	
						   

//获取M
	Rcpp::XPtr<BigMatrix> pMat(pBigMat); 
	arma::Mat<double> M((double*) pMat->matrix(), pMat -> nrow(), pMat -> ncol(), false); 


	int n_ind1=pos_A11.size(),n_ind2=pos_A22.size(),n_ind=pos_A11.size()+pos_A22.size(),n_ind_geno=pos_geno.size();
	int n_snp=M.n_cols,n_pedigree=Pedigree.n_rows;

//创建临时的P_A	
	delete_bigmemory_file_cpp("temp_P_A",bigmemory_data_name,bigmemory_data_path,false);
SEXP tmp1=makeA_memory_cpp(Pedigree,bigmemory_data_name+"_temp",bigmemory_data_path);
	
	SEXP pBigMat_A=make_bigmemory_object_address_cpp(n_pedigree,n_pedigree,bigmemory_data_name+"_temp_P_A",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_A(pBigMat_A);
	arma::Mat<double> P_A((double*) pMat_A->matrix(), pMat_A -> nrow(), pMat_A -> ncol(), false);

//创建临时的 A22
	delete_bigmemory_file_cpp("temp_A22",bigmemory_data_name,bigmemory_data_path,false);
	SEXP pBigMat_A22=make_bigmemory_object_address_cpp(n_ind2,n_ind2,bigmemory_data_name+"_temp_A22",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_A22(pBigMat_A22);
	arma::Mat<double> A22((double*) pMat_A22->matrix(), pMat_A22 -> nrow(), pMat_A22 -> ncol(), false);

//创建临时的 A22_inv
	delete_bigmemory_file_cpp("temp_A22_inv",bigmemory_data_name,bigmemory_data_path,false);
	SEXP pBigMat_A22_inv=make_bigmemory_object_address_cpp(n_ind2,n_ind2,bigmemory_data_name+"_temp_A22_inv",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_A22_inv(pBigMat_A22_inv);
	arma::Mat<double> A22_inv((double*) pMat_A22_inv->matrix(), pMat_A22_inv -> nrow(), pMat_A22_inv -> ncol(), false);

//创建临时的 M_new
	delete_bigmemory_file_cpp("temp_M_new",bigmemory_data_name,bigmemory_data_path,false);
	SEXP pBigMat_M_new=make_bigmemory_object_address_cpp(n_ind_geno,n_snp,bigmemory_data_name+"_temp_M_new",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_M_new(pBigMat_M_new);
	arma::Mat<double> M_new((double*) pMat_M_new->matrix(), pMat_M_new -> nrow(), pMat_M_new -> ncol(), false);
	
	M_new.submat(span(0,n_ind_geno-1),span(0,n_snp-1))=M.rows(pos_geno);


	SEXP H_A_bigmemory_object=R_NilValue;
	SEXP H_A_inv_bigmemory_object=R_NilValue;

//创建临时的 A11
	delete_bigmemory_file_cpp("temp_A11",bigmemory_data_name,bigmemory_data_path,false);
	SEXP pBigMat_A11=make_bigmemory_object_address_cpp(n_ind1,n_ind1,bigmemory_data_name+"_temp_A11",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_A11(pBigMat_A11);
	arma::Mat<double> A11((double*) pMat_A11->matrix(), pMat_A11 -> nrow(), pMat_A11 -> ncol(), false);


//创建临时的 A12
	delete_bigmemory_file_cpp("temp_A12",bigmemory_data_name,bigmemory_data_path,false);
	SEXP pBigMat_A12=make_bigmemory_object_address_cpp(n_ind1,n_ind2,bigmemory_data_name+"_temp_A12",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_A12(pBigMat_A12);
	arma::Mat<double> A12((double*) pMat_A12->matrix(), pMat_A12 -> nrow(), pMat_A12 -> ncol(), false);


//创建临时的G_A	
	delete_bigmemory_file_cpp("temp_G_A",bigmemory_data_name,bigmemory_data_path,false);
List tmp2=G_matrix_memory_cpp(pBigMat_M_new,bigmemory_data_name+"_temp",bigmemory_data_path,true,false,true,false); //metafounder的G矩阵构建
	
	SEXP pBigMat_G_A=make_bigmemory_object_address_cpp(n_ind_geno,n_ind_geno,bigmemory_data_name+"_temp_G_A",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_G_A(pBigMat_G_A);
	arma::Mat<double> G_A((double*) pMat_G_A->matrix(), pMat_G_A -> nrow(), pMat_G_A -> ncol(), false);

				Rcout<<"bigmemory-Start constructing metafounder-H_Additive relationship matrix......"<<endl; 						   


				A22.submat(span(0,n_ind2-1),span(0,n_ind2-1))=P_A.submat(pos_A22,pos_A22);
				A22_inv.submat(span(0,n_ind2-1),span(0,n_ind2-1))=arma::inv(A22);
				arma::Mat<double> identity1(n_ind2,1,fill::ones);
				arma::Mat<double> identity2(n_ind,n_ind,fill::ones);			
				arma::Row<double> miu=arma::inv(identity1.t()*A22_inv*identity1)*identity1.t()*A22_inv*M_new;
				arma::Row<double> p=miu/2;
				double gama=8*arma::var(p);
				Rcout<<"Estimated gama="<<gama<<endl;
				double k=1-gama/2;
				P_A=P_A*(1-gama/2)+gama*identity2;
				A22.submat(span(0,n_ind2-1),span(0,n_ind2-1))=P_A.submat(pos_A22,pos_A22);
				A22_inv.submat(span(0,n_ind2-1),span(0,n_ind2-1))=arma::inv(A22);
				
				A11.submat(span(0,n_ind1-1),span(0,n_ind1-1))=P_A.submat(pos_A11,pos_A11);
				A12.submat(span(0,n_ind1-1),span(0,n_ind2-1))=P_A.submat(pos_A11,pos_A22);
						
				Rcout<<"Adjusting G_A , the adjusting parameter omega is: omega= "<<omega<<endl;
				G_A=G_A*(1-omega)+A22*omega;
				arma::Col<double> inbred1=A11.diag(),inbred2=A22.diag(),inbred;
				inbred=join_cols(inbred1,inbred2);
				
				if(direct==true){
				delete_bigmemory_file_cpp("H_A",bigmemory_data_name,bigmemory_data_path,true);
				SEXP pBigMat_H_A=make_bigmemory_object_address_cpp(n_ind,n_ind,bigmemory_data_name+"_H_A",bigmemory_data_path,"double");	
				Rcpp::XPtr<BigMatrix> pMat_H_A(pBigMat_H_A);
				arma::Mat<double> H_A((double*) pMat_H_A->matrix(), pMat_H_A -> nrow(), pMat_H_A -> ncol(), false);
					
				arma::Mat<double> tmp=A12*A22_inv*G_A;	
			
				H_A.submat(span(0,n_ind-1),span(0,n_ind-1))=P_A.submat(pos_A,pos_A);
				H_A.submat(span(0,n_ind1-1),span(0,n_ind1-1))=A11+A12*A22_inv*(G_A-A22)*A22_inv*(A12.t());
				H_A.submat(span(0,n_ind1-1),span(n_ind1,n_ind-1))=tmp;
				H_A.submat(span(n_ind1,n_ind-1),span(0,n_ind1-1))=tmp.t();
				H_A.submat(span(n_ind1,n_ind-1),span(n_ind1,n_ind-1))=G_A;
				H_A=H_A/k;   //使得metafounder的方差组分是 compatable
				H_A_bigmemory_object=make_bigmemory_object_cpp(n_ind1,n_ind2,bigmemory_data_name+"_H_A",bigmemory_data_path,"double");	
				Rcout<<"bigmemory-Complete constructing metafounder-H_Additive relationship matrix!"<<endl; 				
				}
								

				if(inverse==true){

				delete_bigmemory_file_cpp("H_Ainv",bigmemory_data_name,bigmemory_data_path,true);
				SEXP pBigMat_H_Ainv=make_bigmemory_object_address_cpp(n_ind,n_ind,bigmemory_data_name+"_H_Ainv",bigmemory_data_path,"double");	
				Rcpp::XPtr<BigMatrix> pMat_H_Ainv(pBigMat_H_Ainv);
				arma::Mat<double> H_Ainv((double*) pMat_H_Ainv->matrix(), pMat_H_Ainv -> nrow(), pMat_H_Ainv -> ncol(), false);
				
//创建临时的 G_Ainv
	delete_bigmemory_file_cpp("temp_G_Ainv",bigmemory_data_name,bigmemory_data_path,false);
	SEXP pBigMat_G_Ainv=make_bigmemory_object_address_cpp(n_ind_geno,n_ind_geno,bigmemory_data_name+"_temp_G_Ainv",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_G_Ainv(pBigMat_G_Ainv);
	arma::Mat<double> G_Ainv((double*) pMat_G_Ainv->matrix(), pMat_G_Ainv -> nrow(), pMat_G_Ainv -> ncol(), false);


//创建临时的 P_Ainv

	delete_bigmemory_file_cpp("temp_P_Ainv",bigmemory_data_name,bigmemory_data_path,false);
	
	SEXP pBigMat_Ainv=make_bigmemory_object_address_cpp(n_pedigree,n_pedigree,bigmemory_data_name+"_temp_P_Ainv",bigmemory_data_path,"double");	
	Rcpp::XPtr<BigMatrix> pMat_Ainv(pBigMat_Ainv);
	arma::Mat<double> P_Ainv((double*) pMat_Ainv->matrix(), pMat_Ainv -> nrow(), pMat_Ainv -> ncol(), false);
				P_Ainv.submat(span(0,n_pedigree-1),span(0,n_pedigree-1))=arma::inv(P_A);
				
				G_Ainv.submat(span(0,n_ind_geno-1),span(0,n_ind_geno-1))=arma::inv(G_A);
				
				H_Ainv.submat(span(0,n_ind-1),span(0,n_ind-1))=P_Ainv.submat(pos_A,pos_A);
				H_Ainv.submat(pos_H22,pos_H22)=H_Ainv.submat(pos_H22,pos_H22)+G_Ainv-A22_inv;	
Rcout<<"bigmemory-Complete constructing inverse-H_Additive relationship matrix!"<<endl; 	
				H_A_inv_bigmemory_object=make_bigmemory_object_cpp(n_ind1,n_ind2,bigmemory_data_name+"_H_Ainv",bigmemory_data_path,"double");	
				Rcout<<"bigmemory-Complete constructing inverse-metafounder-H_Additive relationship matrix!"<<endl; 			
				}		
				
delete_bigmemory_file_cpp("temp_P_A",bigmemory_data_name,bigmemory_data_path,false);
delete_bigmemory_file_cpp("temp_A22",bigmemory_data_name,bigmemory_data_path,false);
delete_bigmemory_file_cpp("temp_A22_inv",bigmemory_data_name,bigmemory_data_path,false);
delete_bigmemory_file_cpp("temp_A11",bigmemory_data_name,bigmemory_data_path,false);
delete_bigmemory_file_cpp("temp_A12",bigmemory_data_name,bigmemory_data_path,false);
delete_bigmemory_file_cpp("temp_M_new",bigmemory_data_name,bigmemory_data_path,false);
delete_bigmemory_file_cpp("temp_G_A",bigmemory_data_name,bigmemory_data_path,false);
delete_bigmemory_file_cpp("temp_G_Ainv",bigmemory_data_name,bigmemory_data_path,false);
delete_bigmemory_file_cpp("temp_P_Ainv",bigmemory_data_name,bigmemory_data_path,false);	
	
				return List::create(Named("H") = H_A_bigmemory_object,_["Hinv"] = H_A_inv_bigmemory_object,Named("inbred") =inbred);
}




// type 1: haplotype; type 2: numeric ;
// [[Rcpp::export]]
List  K_matrix_cal_memory(SEXP numeric_address,
				          std::string bigmemory_data_name, 
				          std::string bigmemory_data_path,
						  CharacterVector IND_geno,
						  arma::Mat<int> Pedigree,
						  arma::uvec pos_A11,
						  arma::uvec pos_A22,
						  arma::uvec pos_geno,
						  arma::uvec pos_A,
						  arma::uvec pos_H22,
						  bool H_A_direct=false,
						  double omega=0.05,	  
				          bool base=false,
				          bool trace=false,
						  bool metafounder=false,
				          bool inverse=false,
						  int gene_dropping_iteration=1000,
				          int  type=1,
				          int cpu_cores=1,
						  bool APY_algorithm=false,
						  double APY_eigen_threshold=0.95,
						  int APY_n_core=0,
						  bool re_inverse=false){


		
		SEXP data_numeric_double_address;
		
		switch(type){
	
		case 1:{       //type1:  genomic additive 

		delete_bigmemory_file_cpp("double_numeric",bigmemory_data_name,bigmemory_data_path,true);
		
		data_numeric_double_address=bigmemory_double_type(numeric_address, bigmemory_data_name,bigmemory_data_path,cpu_cores);
		
		
		return G_matrix_memory_cpp(data_numeric_double_address,bigmemory_data_name, bigmemory_data_path,
							       base,trace,metafounder,inverse);
				break;}			

		case 2:{       //type2:  genomic dominance(classical) 

		delete_bigmemory_file_cpp("double_numeric",bigmemory_data_name,bigmemory_data_path,true);
		
		data_numeric_double_address=bigmemory_double_type(numeric_address, bigmemory_data_name,bigmemory_data_path,cpu_cores);

		return Deviation_matrix_memory_cpp(data_numeric_double_address,bigmemory_data_name, bigmemory_data_path,
							               base,trace,inverse);
			break;}

		case 3:{       //type3:  genomic dominance(genotypic) 

		delete_bigmemory_file_cpp("double_numeric",bigmemory_data_name,bigmemory_data_path,true);
		
		data_numeric_double_address=bigmemory_double_type(numeric_address, bigmemory_data_name,bigmemory_data_path,cpu_cores);

		return  Dominance_matrix_memory_cpp(data_numeric_double_address,bigmemory_data_name, bigmemory_data_path,
							        base,trace,inverse);
			break;}

		case 4:{       //type4:  APY

		delete_bigmemory_file_cpp("double_numeric",bigmemory_data_name,bigmemory_data_path,true);
		
		data_numeric_double_address=bigmemory_double_type(numeric_address, bigmemory_data_name,bigmemory_data_path,cpu_cores);

		return APY_inverse_memory_cpp(data_numeric_double_address,bigmemory_data_name,bigmemory_data_path,IND_geno,
				                       APY_eigen_threshold,APY_n_core,re_inverse);
			break;}

		case 5:{       //type5:  H_A

		delete_bigmemory_file_cpp("double_numeric",bigmemory_data_name,bigmemory_data_path,true);
		
		data_numeric_double_address=bigmemory_double_type(numeric_address, bigmemory_data_name,bigmemory_data_path,cpu_cores);

		return 	makeHA_memory_cpp(Pedigree,data_numeric_double_address,IND_geno,bigmemory_data_name,bigmemory_data_path,
						          pos_A11,pos_A22,pos_geno,
						          pos_A,pos_H22,APY_algorithm,APY_eigen_threshold,APY_n_core,H_A_direct,inverse,omega);
			break;}

		case 6:{       //type6:  H_A_metafounder

		delete_bigmemory_file_cpp("double_numeric",bigmemory_data_name,bigmemory_data_path,true);
		
		data_numeric_double_address=bigmemory_double_type(numeric_address, bigmemory_data_name,bigmemory_data_path,cpu_cores);

		return 	 makeHA_metafounder_memory_cpp(Pedigree,data_numeric_double_address,IND_geno,bigmemory_data_name,bigmemory_data_path,
						       pos_A11,pos_A22,pos_geno,
						       pos_A,pos_H22,H_A_direct,inverse,omega);
			break;}

		case 7:{       //type7:  P_A

		return 	makeA_memory_cpp(Pedigree,bigmemory_data_name,bigmemory_data_path);
			break;}

		case 8:{       //type8:  P_Ainv

		return 	makeAinv_memory_cpp(Pedigree,bigmemory_data_name,bigmemory_data_path);
			break;}

		case 9:{       //type9:  P_D

		return 	makeD_memory_cpp(Pedigree,bigmemory_data_name,bigmemory_data_path,inverse);
			break;}

		case 10:{       //type10:  P_D(gene_dropping)

		return 	gene_dropping_memory_D(Pedigree,bigmemory_data_name,bigmemory_data_path,
									   inverse,gene_dropping_iteration,true, cpu_cores);
			break;}
			
		default:
		
        throw Rcpp::exception("unknown input data type!");			
			
		}
}