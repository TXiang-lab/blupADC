#include "shared_function.h"
#define	matrix_threshold  0.00001


// [[Rcpp::export]]
arma::mat cal_breed_pro(arma::Mat<int> Pedigree, int n_meta){   //前三列为个体号，父亲，母亲，必须确保前面几行为Metafounder
	int n_ind=Pedigree.n_rows;
	int Sire_id,Dam_id;
	arma::Col<int> Sire=Pedigree.col(1);
	arma::Col<int> Dam=Pedigree.col(2);

    arma::mat Breed_matrix(n_ind,n_meta);	
	Breed_matrix.fill(0);
	//Breed_matrix.col(0)=Pedigree.col(0);
	//Breed_matrix.col(n_breed)=Breed;

	
	for(int i=0; i < n_ind; i++){
		Sire_id=Sire[i];
		Dam_id=Dam[i];		
		
	
		if(Sire_id==0&&Dam_id==0){  //父母均未知,
		
			Breed_matrix(i,i)=1; //父母未知的情况下,此时品种必须是已知的
				
		}else{                     //父母均已知	

			
			for(int j=0;j<n_meta;j++){
			Breed_matrix(i,j)=(Breed_matrix(Sire_id-1,j)+Breed_matrix(Dam_id-1,j))/2; //父母均已知,计算品种j的比例				
		   }
	    }
	}

	return Breed_matrix; 
}

// [[Rcpp::export]]
arma::Mat<double> makeA_cpp(arma::Mat<int> Pedigree){   //稀疏矩阵运行速度会变慢


	Rcout<<"Start constructing pedigree additive relationship matrix......"<<endl;  		
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
	Rcout<<"Complete constructing pedigree additive relationship matrix!"<<endl;  	
	return A_matrix;
	}


// [[Rcpp::export]]
arma::Mat<double> makeA_meta_cpp(arma::Mat<int> Pedigree, arma::Mat<double> Gamma){   //稀疏矩阵运行速度会变慢


	Rcout<<"Start constructing pedigree additive relationship matrix......"<<endl;  		
	arma::Mat<double> A_matrix(Pedigree.n_rows,Pedigree.n_rows);	
	A_matrix.fill(0);
	arma::Col<int> Animal=Pedigree.col(0);
	arma::Col<int> Sire=Pedigree.col(1);
	arma::Col<int> Dam=Pedigree.col(2);
	int IND=Pedigree.n_rows;
	bool status1,status2;
	int n_meta=Gamma.n_rows;


A_matrix.submat(span(0,n_meta-1),span(0,n_meta-1))=Gamma;


	
	for(int i=n_meta; i < IND; i++){

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
	Rcout<<"Complete constructing pedigree additive relationship matrix!"<<endl;  	
	return A_matrix;
	}


//计算加性效应矩阵	
// [[Rcpp::export]]
List G_matrix_cpp(arma::Mat<int> & M, bool base=false,bool trace=false, bool metafounder=false, bool inverse=false){ //M矩阵
	
	Rcout<<"Start constructing genomic additive relationship matrix......"<<endl;  	
	arma::Mat<double> M1=conv_to<arma::Mat<double>>::from(M);
	int ind_size=M1.n_rows; // 个体数
	int snp_size=M1.n_cols; //SNP数量

	arma::Row<double> p_frequency(snp_size),q_frequency(snp_size);
	arma::Mat<double> G,G_inv,b;
	arma::Row<double> a;
	if(base==TRUE){
	p_frequency.fill(0.5);
	q_frequency.fill(0.5);
	for(int i=0; i < snp_size; i++){
	M1.col(i)=M1.col(i)-1;
	}	
	}else{
	for(int i=0; i < snp_size; i++){

	p_frequency[i]=mean(M1.col(i));
	M1.col(i)=M1.col(i)-mean(M1.col(i));
	}
	p_frequency=0.5*p_frequency;
	q_frequency=1-p_frequency;
	}

	if(trace==true){
	b=M1*M1.t();
	G=(2*ind_size*b)/(arma::trace(b));				
	}else if(metafounder==true){
	G=M1*M1.t()/(snp_size/2);			
	}else{
	a=2*p_frequency%q_frequency;
	G=M1*M1.t()/sum(a);
	}
	Rcout<<"Complete constructing genomic additive relationship matrix!"<<endl;  	
	if(inverse==TRUE){
	G.diag()=G.diag()+0.0001;	
	G_inv = inv(G);
	G.diag()=G.diag()-0.0001;
Rcout<<"Complete constructing inverse-genomic additive relationship matrix!"<<endl; 		
	}	
    return List::create(Named("G") = G,_["Ginv"] = G_inv);
	
}



//计算显性离差矩阵
// [[Rcpp::export]]
List Deviation_matrix_cpp(arma::Mat<int> & M, bool inverse=false,bool base=false, bool trace=false){

	Rcout<<"Start constructing genomic dominance(classical) relationship matrix......"<<endl;  	
	arma::Mat<double> M1=conv_to<arma::Mat<double>>::from(M);
	int ind_size=M1.n_rows; // 个体数
	int snp_size=M1.n_cols; //SNP数量
	double p_freq,q_freq;
	arma::Row<double> p_frequency(snp_size), q_frequency(snp_size);
	arma::Mat<double> Devi_inv;
	if(base==TRUE){
	p_frequency.fill(0.5);
	q_frequency.fill(0.5);
	for( int i=0; i < snp_size; i++){
	for( int j=0; j < ind_size; j++){
	if(M1(j,i)==2){M1(j,i)=-2*0.25;
	}else if(M1(j,i)==1){M1(j,i)=2*0.25;
	}else{M1(j,i)=-2*0.25;}	
	}}
    }else{	
	for( int i=0; i < snp_size; i++){
	p_freq=0.5*mean(M1.col(i));
	q_freq=1-p_freq;
	for( int j=0; j < ind_size; j++){
	if(M1(j,i)==2){M1(j,i)=-2*q_freq*q_freq;
	}else if(M1(j,i)==1){M1(j,i)=2*p_freq*q_freq;
	}else{M1(j,i)=-2*p_freq*p_freq;}
	p_frequency[i]=p_freq;
	}
	}
	 q_frequency=1-p_frequency;
	} 
	 
	arma::Row<double> a=4*pow(p_frequency,2)%pow(q_frequency,2); //行向量*行向量
	arma::Mat<double> Devi=M1*M1.t()/sum(a);
Rcout<<"Complete constructing genomic dominance(classical) relationship matrix!"<<endl; 	
	if(inverse==TRUE){	
	Devi.diag()=Devi.diag()+0.0001;	
	Devi_inv = inv(Devi);
	Devi.diag()=Devi.diag()-0.0001;
	
Rcout<<"Complete constructing inverse-genomic dominance(classical) relationship matrix!"<<endl; 		
	}
	
List L = List::create(Named("D") = Devi,_["Dinv"] = Devi_inv); 	
    return L;		
	
}


//计算显性效应矩阵
// [[Rcpp::export]]
List Dominance_matrix_cpp(arma::Mat<int> & M,bool inverse=true,bool base=false, bool trace=false){

	Rcout<<"Start constructing genomic dominance(genotypic) relationship matrix......"<<endl;  
	arma::Mat<double> M1=conv_to<arma::Mat<double>>::from(M);
	int ind_size=M1.n_rows; // 个体数
	int snp_size=M1.n_cols; //SNP数量
	double p_freq,q_freq;
	arma::Row<double> p_frequency(snp_size),q_frequency(snp_size);
	arma::Mat<double> Domi_inv;

	if(base==TRUE){
	p_frequency.fill(0.5);
	q_frequency.fill(0.5);
	for( int i=0; i < snp_size; i++){
	for( int j=0; j < ind_size; j++){
	if(M1(j,i)==2){M1(j,i)=-2*0.25;
	}else if(M1(j,i)==1){M1(j,i)=1-2*0.25;
	}else{M1(j,i)=-2*0.25;}	
	}}
    }else{
	
	for( int i=0; i < snp_size; i++){
	p_freq=0.5*mean(M1.col(i));
	q_freq=1-p_freq;
	for( int j=0; j < ind_size; j++){
	if(M1(j,i)==2){M1(j,i)=-2*p_freq*q_freq;
	}else if(M1(j,i)==1){M1(j,i)=1-2*p_freq*q_freq;
	}else{M1(j,i)=-2*p_freq*q_freq;}
	p_frequency[i]=p_freq;
	}
	}
	 q_frequency=1-p_frequency;
	}
	
	arma::Row<double> a=2*p_frequency%q_frequency;
	arma::Mat<double> Domi=M1*M1.t()/sum(a%(1-a));
Rcout<<"Complete constructing genomic dominance(genotypic) relationship matrix!"<<endl; 
	
	if(inverse==TRUE){	
	Domi.diag()=Domi.diag()+0.0001;	
	Domi_inv = inv(Domi);
	Domi.diag()=Domi.diag()-0.0001;
	
Rcout<<"Complete constructing inverse-genomic dominance(genotypic) relationship matrix!"<<endl; 		
	}
	
List L = List::create(Named("D") = Domi,_["Dinv"] = Domi_inv); 	
    return L;	
}		  



// [[Rcpp::export]]
arma::Mat<double> APY_cpp(arma::Mat<double> G,CharacterVector IND_Proven, CharacterVector IND_G){ //G为所有个体的亲缘关系矩阵
																				   //APY算法																				   
	CharacterVector IND_Young=setdiff(IND_G,IND_Proven);
	IntegerVector IND_Young_pos1=match(IND_Young,IND_G)-1;
	IntegerVector IND_Proven_pos1=match(IND_Proven,IND_G)-1;
	uvec IND_Young_pos=as<arma::uvec>(IND_Young_pos1);
	uvec IND_Proven_pos=as<arma::uvec>(IND_Proven_pos1);
	int IND_Young_n=IND_Young_pos.n_rows;
	int IND_Proven_n=IND_Proven_pos.n_rows;
	int IND_Total=G.n_rows;
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
	
	//mat G_inv(IND_Total,IND_Total,fill::zeros); 由于不能主动释放G，因此为了节省空间，用G存储G_inv
	G.submat(span(0,IND_Proven_n-1),span(0,IND_Proven_n-1))=Py_t_M_inv*P_young+Gpp_inv;
	G.submat(span(0,IND_Proven_n-1),span(IND_Proven_n,IND_Total-1))=-1*Py_t_M_inv;
	G.submat(span(IND_Proven_n,IND_Total-1),span(0,IND_Proven_n-1))=-1*Py_t_M_inv.t();
	G.submat(span(IND_Proven_n,IND_Total-1),span(IND_Proven_n,IND_Total-1))=diagmat(M_inv);
	//reorder the position, 
	CharacterVector IND_G_final=union_cpp(IND_Proven,IND_Young);
	IntegerVector reorder_pos1=match(IND_G,IND_G_final)-1;
	uvec reorder_pos=as<arma::uvec>(reorder_pos1);	
	G=G.submat(reorder_pos,reorder_pos);
	return G;
}

// [[Rcpp::export]]
List APY_inverse_cpp(arma::Mat<int> & M,CharacterVector IND_geno, 
								double APY_eigen_threshold=0.95,
								int APY_n_core=0,
								bool re_inverse=false){

Rcout<<"Start constructing inverse-genomic additive relationship matrix(APY method)......"<<endl;
arma::Mat<double> G=G_matrix_cpp(M,false,false,false,false)["G"];
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
CharacterVector IND_Proven=sample(IND_geno,N_selected,false);

arma::Mat<double> G_inv=APY_cpp(G,IND_Proven,IND_geno);
arma::Mat<double> G1;
Rcout<<"Complete constructing inverse-genomic additive relationship matrix(APY method)!"<<endl;
if(re_inverse==true){
G_inv.diag()=G_inv.diag()+0.0001;	
G1 =arma::inv(G_inv);	
G_inv.diag()=G_inv.diag()-0.0001;	
Rcout<<"Complete constructing inverse-(inverse-genomic additive relationship matrix(APY method))!"<<endl;
}	

return List::create(Named("G") = G1,_["Ginv"] = G_inv);
 
}
	





// [[Rcpp::export]]
arma::Mat<double> makeAinv_cpp(arma::Mat<int> Pedigree, bool inbreeding=false){

	Rcout<<"Start constructing inverse-pedigree additive relationship matrix......"<<endl;  
	int IND=Pedigree.n_rows;	
	arma::Col<int> Animal=Pedigree.col(0);
	arma::Col<int> Sire=Pedigree.col(1);
	arma::Col<int> Dam=Pedigree.col(2);
	arma::Mat<double> L(IND,IND);
	arma::Mat<double> A(IND,IND);	
	arma::Mat<double> A_inv(IND,IND);
	A_inv.fill(0);
		
	double tmp,d;
if(inbreeding==true){	

	for(int i=0; i < IND; i++){
	
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

	for(int i=0; i < IND; i++){

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
	Rcout<<"Complete constructing inverse-pedigree additive relationship matrix!"<<endl;   
	 return A_inv; 
}	


// [[Rcpp::export]]
DataFrame cal_homo_inbred_cpp(arma::Mat<int> & data_genumeric, CharacterVector IND_geno){
	
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
arma::Mat<double> makeInbreeding_cpp(arma::Mat<int> Pedigree){
	int IND=Pedigree.n_rows;	
	arma::Col<int> Animal=Pedigree.col(0);
	arma::Col<int> Sire=Pedigree.col(1);
	arma::Col<int> Dam=Pedigree.col(2);
	arma::Mat<double> L(IND,IND);
	arma::Mat<double> A(IND,IND);	

	double tmp;

	for(int i=0; i < IND; i++){
	
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
List makeD_cpp(arma::Mat<int> Pedigree, bool inverse=false){   //稀疏矩阵运行速度会变慢

Rcout<<"Start constructing pedigree dominance relationship matrix......"<<endl;	
	arma::Mat<double> A_matrix(Pedigree.n_rows,Pedigree.n_rows);	
	arma::Mat<double> D_matrix(Pedigree.n_rows,Pedigree.n_rows);	
    D_matrix.fill(0);
	A_matrix.fill(0);	
	arma::Col<int> Animal=Pedigree.col(0);
	arma::Col<int> Sire=Pedigree.col(1);
	arma::Col<int> Dam=Pedigree.col(2);
	int Sire_i_pos,Sire_j_pos,Dam_i_pos,Dam_j_pos,IND=Pedigree.n_rows;
	for(int i=0; i < IND; i++){

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
Rcout<<"Complete constructing pedigree dominance relationship matrix!"<<endl;		
	arma::Mat<double> D_inv;
	if(inverse==true){
		D_matrix.diag()=D_matrix.diag()+0.0001;
		D_inv=arma::inv(D_matrix);
		D_matrix.diag()=D_matrix.diag()-0.0001;	
Rcout<<"Complete constructing inverse-pedigree dominance relationship matrix!"<<endl;			
	}
	return List::create(Named("D") = D_matrix,_["Dinv"] = D_inv);
	}



// openMP, not very accurate
// [[Rcpp::export]]
List gene_dropping_D(arma::Mat<int> Pedigree, bool inverse=false, int iteration=1000, bool diagnoal=true, int cpu_cores=1){
	
	Rcout<<"Start constructing pedigree dominance relationship matrix(gene dropping)......"<<endl; 
     omp_set_num_threads(cpu_cores);
//父母只要有任何一方缺失，显性亲缘关系就为0
	Rcout<<"The max iteration is "<<iteration<<endl;
	
	int n_ind=Pedigree.n_rows;	
	arma::Mat<int> allele(n_ind,3);
	arma::Mat<double> D_sim(n_ind,n_ind);
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
	arma::Mat<double> A=makeA_cpp(Pedigree);
	D_sim.diag()=2-A.diag();
	}
	arma::Mat<double> D_inv;
	Rcout<<"Complete constructing pedigree dominance relationship matrix(gene dropping)!"<<endl;  	
	if(inverse==true){
		D_sim.diag()=D_sim.diag()+0.0001;
		D_inv=arma::inv(D_sim);
		D_sim.diag()=D_sim.diag()-0.0001;
	Rcout<<"Complete constructing inverse-pedigree dominance relationship matrix(gene dropping)!"<<endl;; 		
	}
	return List::create(Named("D") = D_sim,_["Dinv"] = D_inv);
}
	


// [[Rcpp::export]]
List makeHA_cpp(arma::Mat<int> & Pedigree, arma::Mat<int> & M, 
						   CharacterVector IND_geno,
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
				
				arma::Mat<double> P_A=makeA_cpp(Pedigree);
				arma::Mat<double> P_Ainv=makeAinv_cpp(Pedigree);
				arma::Mat<int> M_new=M.rows(pos_geno);
				arma::Mat<double> H_A,H_Ainv;
				arma::Col<double> inbred;

				Rcout<<"Start constructing H_Additive relationship matrix......"<<endl; 						   
				
				arma::Mat<double> A22=P_A.submat(pos_A22,pos_A22);
				arma::Mat<double> A22_inv=arma::inv(A22);		

				if(APY_algorithm==false){		
				arma::Mat<double> A11=P_A.submat(pos_A11,pos_A11);
				arma::Mat<double> A12=P_A.submat(pos_A11,pos_A22);				
				arma::Mat<double> G_A=G_matrix_cpp(M_new)["G"];
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
				int ind_n1=pos_A11.size(),ind_n2=pos_A22.size(),ind_n=pos_A11.size()+pos_A22.size();
				if(direct==true){
				arma::Mat<double> tmp=A12*A22_inv*G_A;	
				H_A=P_A.submat(pos_A,pos_A);
				H_A.submat(span(0,ind_n1-1),span(0,ind_n1-1))=A11+A12*A22_inv*(G_A-A22)*A22_inv*(A12.t());
				H_A.submat(span(0,ind_n1-1),span(ind_n1,ind_n-1))=tmp;
				H_A.submat(span(ind_n1,ind_n-1),span(0,ind_n1-1))=tmp.t();
				H_A.submat(span(ind_n1,ind_n-1),span(ind_n1,ind_n-1))=G_A;
				inbred=H_A.diag()-1;
				}
				if(inverse==true){
				
				arma::Mat<double> G_Ainv;
				G_Ainv=arma::inv(G_A);
				//arma::Mat<double> P_Ainv=makeAinv_cpp(Pedigree);
				H_Ainv=P_Ainv.submat(pos_A,pos_A);
				H_Ainv.submat(pos_H22,pos_H22)=H_Ainv.submat(pos_H22,pos_H22)+G_Ainv-A22_inv;				
				}			
				}else{
				
				arma::Mat<double> G_Ainv=APY_inverse_cpp(M_new,IND_geno,APY_eigen_threshold,APY_n_core,false)["Ginv"];	
		
				arma::Mat<double> P_Ainv=makeAinv_cpp(Pedigree);
				H_Ainv=P_Ainv.submat(pos_A,pos_A);
				H_Ainv.submat(pos_H22,pos_H22)=H_Ainv.submat(pos_H22,pos_H22)+G_Ainv-A22_inv;					
					
					
				}
				
				return List::create(Named("H") =H_A,_["Hinv"] = H_Ainv,Named("inbred") =inbred);
}						   
	
	

// [[Rcpp::export]]
List makeHA_metafounder_cpp(arma::Mat<int> & Pedigree, arma::Mat<int> & M, 
						   std::vector<std::string> meta_ind,
							arma::mat provided_gamma,
						   arma::uvec pos_A11,
						   arma::uvec pos_A22,
						   arma::uvec pos_geno,
						   arma::uvec pos_A,
						   arma::uvec pos_H22,
						   bool direct=false,
						   bool inverse=true,
						   double omega=0.05,
						   bool scaled=false,
						   std::string gamma_method="gls"
						   ){	
				
				arma::Mat<double> P_A=makeA_cpp(Pedigree);
				arma::Mat<double> P_Ainv=makeAinv_cpp(Pedigree);
				arma::Mat<int> M_new=M.rows(pos_geno);
				arma::Mat<double> G_A=G_matrix_cpp(M_new,true,false,true,false)["G"];
				int ind_n1=pos_A11.size(),ind_n2=pos_A22.size(),ind_n=pos_A11.size()+pos_A22.size();
				Rcout<<"Start constructing metafounder-H_Additive relationship matrix......"<<endl; 
				
				arma::Mat<double> A22=P_A.submat(pos_A22,pos_A22);
				arma::Mat<double> A22_inv=inv(A22);
				//arma::Mat<double> identity2(ind_n,ind_n,fill::ones);			
				
				int n_meta=meta_ind.size();
				arma::mat gamma(n_meta,n_meta),Q(ind_n,n_meta,fill::ones),miu,M1;
				M1=conv_to<arma::Mat<double>>::from(M_new);
			    
				if(n_meta>=2){
					
					Q=cal_breed_pro(Pedigree,n_meta);				
				}
				Q=Q.rows(pos_A22);
				if(gamma_method=="naive"){	
					Rcout<<"Estimating gamma via naive method......"<<endl;
					
					
					if(n_meta==1){
						gamma=2*var(mean(M1));
					}else{
						miu=inv(Q.t() *Q)*(Q.t()*M1);
						
						for(int i=0;i<n_meta;i++){
						
							gamma(i,i)=2*sum((miu.row(i)-1)%(miu.row(i)-1))/M1.n_cols;
						
							for(int j=0;j<i;j++){
							
								gamma(i,j)=gamma(j,i)=2*sum((miu.row(i)-1)%(miu.row(j)-1))/M1.n_cols;
							
							}
						}
					}

				}else if(gamma_method=="gls"){
					Rcout<<"Estimating gamma via gls method......"<<endl;
					//arma::Mat<double> identity1(ind_n2,1,fill::ones);
					//arma::Row<double> miu=1/accu(A22_inv)*identity1.t()*A22_inv*M_new;	
					
					if(n_meta==1){
						miu=1/accu(A22_inv)*Q.t()*A22_inv*M1;
						gamma=2*var(miu.row(0));
					}else{
						
						miu=inv(Q.t()*A22_inv*Q)*(Q.t()*A22_inv*M1);
						
						for(int i=0;i<n_meta;i++){
							
							//arma::Row<double> tmpi=miu.row(i);
							gamma(i,i)=2*sum((miu.row(i)-1)%(miu.row(i)-1))/M1.n_cols;
						
							for(int j=0;j<i;j++){
								
								//arma::Row<double> tmpj=miu.row(j);							
								gamma(i,j)=gamma(j,i)=2*sum((miu.row(i)-1)%(miu.row(j)-1))/M1.n_cols;
							
							}
						}		
					}					
				
				}else if(gamma_method=="provided"){	
					Rcout<<"Using provided gamma to perform calculating!"<<endl;
					gamma=provided_gamma;
				}
				
				Rcout<<"Estimated gamma="<<endl;
				Rcout<<gamma<<endl;
				
				P_A=makeA_meta_cpp(Pedigree,gamma); 			 
				//P_Ainv=inv(P_A);
				double k=1+mean(gamma.diag())/2-mean(mean(gamma));
				//P_A=P_A*(1-gamma/2)+gamma*identity2;
				A22=P_A.submat(pos_A22,pos_A22);
				A22_inv=arma::inv(A22);
				//P_A(0,0)=gamma;			
				
				arma::Mat<double> A11=P_A.submat(pos_A11,pos_A11);
				arma::Mat<double> A12=P_A.submat(pos_A11,pos_A22);
				
				Rcout<<"Adjusting G_A , the adjusting parameter omega is: omega= "<<omega<<endl;
				G_A=G_A*(1-omega)+A22*omega;
				
				arma::Mat<double> H_A,H_Ainv;
				arma::Col<double> inbred1=A11.diag(),inbred2=A22.diag(),inbred;
				inbred=join_cols(inbred1,inbred2);
				
				if(direct==true){
				arma::Mat<double> tmp=A12*A22_inv*G_A;	
				H_A=P_A.submat(pos_A,pos_A);
				H_A.submat(span(0,ind_n1-1),span(0,ind_n1-1))=A11+A12*A22_inv*(G_A-A22)*A22_inv*(A12.t());
				H_A.submat(span(0,ind_n1-1),span(ind_n1,ind_n-1))=tmp;
				H_A.submat(span(ind_n1,ind_n-1),span(0,ind_n1-1))=tmp.t();
				H_A.submat(span(ind_n1,ind_n-1),span(ind_n1,ind_n-1))=G_A;
				H_A=H_A;   
				H_A.clean(matrix_threshold);
				
					if(scaled==true){
						
						H_A=H_A/k;	 //使得metafounder的方差组分是 compatable
							
					}
				
				}
				
				if(inverse==true){
				arma::Mat<double> G_Ainv=arma::inv(G_A);
				arma::Mat<double> P_Ainv=arma::inv(P_A);
				H_Ainv=P_Ainv.submat(pos_A,pos_A);
				H_Ainv.submat(pos_H22,pos_H22)=H_Ainv.submat(pos_H22,pos_H22)+G_Ainv-A22_inv;
				H_Ainv.clean(matrix_threshold);
				
					if(scaled==true){
						
						H_Ainv=H_Ainv*k;	 //使得metafounder的方差组分是 compatable
							
					}
				
				}
				
								
				return List::create(Named("H") =H_A,_["Hinv"] = H_Ainv,Named("inbred") =inbred,
									Named("miu") =miu);
}		



// [[Rcpp::export]]
List makeA_metafounder_cpp(arma::Mat<int> & Pedigree, arma::Mat<int> & M, 
						   std::vector<std::string> meta_ind,
							arma::mat provided_gamma,
						   arma::uvec pos_A11,
						   arma::uvec pos_A22,
						   arma::uvec pos_geno,
						   arma::uvec pos_A,
						   arma::uvec pos_H22,
						   bool direct=false,
						   bool inverse=true,
						   double omega=0.05,
						   bool scaled=false,
						   std::string gamma_method="gls"
						   ){	
				
				arma::Mat<double> P_A=makeA_cpp(Pedigree);
				arma::Mat<double> P_Ainv=makeAinv_cpp(Pedigree);
				arma::Mat<int> M_new=M.rows(pos_geno);
				arma::Mat<double> G_A=G_matrix_cpp(M_new,true,false,true,false)["G"];
				int ind_n1=pos_A11.size(),ind_n2=pos_A22.size(),ind_n=pos_A11.size()+pos_A22.size();
				Rcout<<"Start constructing metafounder-H_Additive relationship matrix......"<<endl; 
				
				arma::Mat<double> A22=P_A.submat(pos_A22,pos_A22);
				arma::Mat<double> A22_inv=P_Ainv.submat(pos_A22,pos_A22);
				//arma::Mat<double> identity2(ind_n,ind_n,fill::ones);			
				
				int n_meta=meta_ind.size();
				arma::mat gamma(n_meta,n_meta),Q(ind_n,n_meta,fill::ones),miu,M1;
				M1=conv_to<arma::Mat<double>>::from(M_new);
			    
				if(n_meta>=2){
					
					Q=cal_breed_pro(Pedigree,n_meta);				
				}
				Q=Q.rows(pos_A22);
				if(gamma_method=="naive"){	
					Rcout<<"Estimating gamma via naive method......"<<endl;
					
					
					if(n_meta==1){
						gamma=2*var(mean(M1));
					}else{
						miu=inv(Q.t() *Q)*(Q.t()*M1);
						
						for(int i=0;i<n_meta;i++){
						
							gamma(i,i)=2*sum((miu.row(i)-1)%(miu.row(i)-1))/M1.n_cols;
						
							for(int j=0;j<i;j++){
							
								gamma(i,j)=gamma(j,i)=2*sum((miu.row(i)-1)%(miu.row(j)-1))/M1.n_cols;
							
							}
						}
					}
					
				}else if(gamma_method=="gls"){
					Rcout<<"Estimating gamma via gls method......"<<endl;
					//arma::Mat<double> identity1(ind_n2,1,fill::ones);
					//arma::Row<double> miu=1/accu(A22_inv)*identity1.t()*A22_inv*M_new;	
					
					if(n_meta==1){
						miu=1/accu(A22_inv)*Q.t()*A22_inv*M1;
						gamma=2*var(miu.row(0));
					}else{
						
						miu=inv(Q.t()*A22_inv*Q)*(Q.t()*A22_inv*M1);
						
						for(int i=0;i<n_meta;i++){
							
							//arma::Row<double> tmpi=miu.row(i);
							gamma(i,i)=2*sum((miu.row(i)-1)%(miu.row(i)-1))/M1.n_cols;
						
							for(int j=0;j<i;j++){
								
								//arma::Row<double> tmpj=miu.row(j);							
								gamma(i,j)=gamma(j,i)=2*sum((miu.row(i)-1)%(miu.row(j)-1))/M1.n_cols;
							
							}
						}		
					}					
				
				}else if(gamma_method=="provided"){	
					Rcout<<"Using provided gamma to perform calculating!"<<endl;
					gamma=provided_gamma;
				}
				
				Rcout<<"Estimated gamma="<<endl;
				Rcout<<gamma<<endl;
				
				P_A=makeA_meta_cpp(Pedigree,gamma); 	
				double k=1+mean(gamma.diag())/2-mean(mean(gamma));				
				if(direct==true){
				P_A.clean(matrix_threshold);			
					if(scaled==true){						
						P_A=P_A/k;	 //使得metafounder的方差组分是 compatable							
					}				
				}
				
				if(inverse==true){
				P_Ainv=arma::inv(P_A);
				P_Ainv.clean(matrix_threshold);				
					if(scaled==true){						
						P_Ainv=P_Ainv*k;	 //使得metafounder的方差组分是 compatable
							
					}
				
				}
				
								
				return List::create(Named("A") =P_A,_["Ainv"] = P_Ainv);
}




// [[Rcpp::export]]
List makeHD_cpp(arma::Mat<int> & Pedigree, arma::Mat<int> & M, 
						   CharacterVector IND_geno,
						   arma::uvec pos_A11,
						   arma::uvec pos_A22,
						   arma::uvec pos_geno,
						   arma::uvec pos_A,
						   arma::uvec pos_H22,
						   bool gene_dropping=false,
						   int cpu_cores=1,
						   int gene_dropping_iteration=1000,
						   bool direct=false,
						   bool inverse=true,
						   double omega=0.05){					
				
				Rcout<<"Start constructing H_Dominance relationship matrix......"<<endl; 	
				arma::Mat<double> P_D=makeD_cpp(Pedigree)[0];				

				arma::Mat<int> M_new=M.rows(pos_geno);
				arma::Mat<double> H_D,H_Dinv;
				
				arma::Mat<double> D22=P_D.submat(pos_A22,pos_A22);
				
				arma::Mat<double> D22_inv=arma::inv(D22);			
		
				arma::Mat<double> D11=P_D.submat(pos_A11,pos_A11);
				
				arma::Mat<double> D12=P_D.submat(pos_A11,pos_A22);		
				
				arma::Mat<double> G_D=Deviation_matrix_cpp(M_new)["D"];
                double diag_GD=mean(G_D.diag());
				double offdiag_GD=mean(G_D.elem(trimatl_ind(size(G_D),-1)));
                Rcout<<"The mean of diagnal and off-diagnal of G_D is : "<<diag_GD<<" and "<<offdiag_GD<<endl;
				
				double diag_D22=mean(D22.diag());
				double offdiag_D22=mean(D22.elem(trimatl_ind(size(D22),-1)));
				Rcout<<"The mean of diagnal and off-diagnal of P_D22 is : "<<diag_D22<<" and "<<offdiag_D22<<endl;

				double beta=(diag_D22-offdiag_D22)/(diag_GD-offdiag_GD);
				double alpha=(offdiag_D22-beta*offdiag_GD);
				
                Rcout<<"Blending G_D to P_D , the adjusting parameter is: alpha= "<<alpha<<" , beta= "<<beta<<endl;
                G_D=alpha+beta*G_D;		
				Rcout<<"Adjusting G_D , the adjusting parameter omega is: omega= "<<omega<<endl;
				G_D=G_D*(1-omega)+D22*omega;
				int ind_n1=pos_A11.size(),ind_n2=pos_A22.size(),ind_n=pos_A11.size()+pos_A22.size();
				if(direct==true){
				arma::Mat<double> tmp=D12*D22_inv*G_D;	
				H_D=P_D.submat(pos_A,pos_A);
				H_D.submat(span(0,ind_n1-1),span(0,ind_n1-1))=D11+D12*D22_inv*(G_D-D22)*D22_inv*(D12.t());
				H_D.submat(span(0,ind_n1-1),span(ind_n1,ind_n-1))=tmp;
				H_D.submat(span(ind_n1,ind_n-1),span(0,ind_n1-1))=tmp.t();
				H_D.submat(span(ind_n1,ind_n-1),span(ind_n1,ind_n-1))=G_D;
				
				}
				if(inverse==true){
				
				arma::Mat<double> G_Dinv;
				G_Dinv=arma::inv(G_D);
				arma::Mat<double> P_Dinv=arma::inv(P_D);
				H_Dinv=P_Dinv.submat(pos_A,pos_A);
				H_Dinv.submat(pos_H22,pos_H22)=H_Dinv.submat(pos_H22,pos_H22)+G_Dinv-D22_inv;				
				}			
				Rcout<<"Complete constructing H_Dominance relationship matrix!"<<endl; 	
				return List::create(Named("H") =H_D,_["Hinv"] = H_Dinv);
}


// [[Rcpp::export]]
List makeHD_gene_dropping_cpp(arma::Mat<int> & Pedigree, arma::Mat<int> & M, 
						   CharacterVector IND_geno,
						   arma::uvec pos_A11,
						   arma::uvec pos_A22,
						   arma::uvec pos_geno,
						   arma::uvec pos_A,
						   arma::uvec pos_H22,
						   bool gene_dropping=false,
						   int cpu_cores=1,
						   int gene_dropping_iteration=1000,
						   bool direct=false,
						   bool inverse=true,
						   double omega=0.05){					
				
				Rcout<<"Start constructing H_Dominance relationship matrix......"<<endl; 	

				arma::Mat<double> P_D=gene_dropping_D(Pedigree,false,gene_dropping_iteration,cpu_cores)[0];


				arma::Mat<int> M_new=M.rows(pos_geno);
				arma::Mat<double> H_D,H_Dinv;
				
				arma::Mat<double> D22=P_D.submat(pos_A22,pos_A22);
				
				arma::Mat<double> D22_inv=arma::inv(D22);			
		
				arma::Mat<double> D11=P_D.submat(pos_A11,pos_A11);
				
				arma::Mat<double> D12=P_D.submat(pos_A11,pos_A22);		
				
				arma::Mat<double> G_D=Deviation_matrix_cpp(M_new)["D"];
                double diag_GD=mean(G_D.diag());
				double offdiag_GD=mean(G_D.elem(trimatl_ind(size(G_D),-1)));
                Rcout<<"The mean of diagnal and off-diagnal of G_D is : "<<diag_GD<<" and "<<offdiag_GD<<endl;
				
				double diag_D22=mean(D22.diag());
				double offdiag_D22=mean(D22.elem(trimatl_ind(size(D22),-1)));
				Rcout<<"The mean of diagnal and off-diagnal of P_D22 is : "<<diag_D22<<" and "<<offdiag_D22<<endl;

				double beta=(diag_D22-offdiag_D22)/(diag_GD-offdiag_GD);
				double alpha=(offdiag_D22-beta*offdiag_GD);
				
                Rcout<<"Blending G_D to P_D , the adjusting parameter is: alpha= "<<alpha<<" , beta= "<<beta<<endl;
                G_D=alpha+beta*G_D;		
				Rcout<<"Adjusting G_D , the adjusting parameter omega is: omega= "<<omega<<endl;
				G_D=G_D*(1-omega)+D22*omega;
				int ind_n1=pos_A11.size(),ind_n2=pos_A22.size(),ind_n=pos_A11.size()+pos_A22.size();
				if(direct==true){
				arma::Mat<double> tmp=D12*D22_inv*G_D;	
				H_D=P_D.submat(pos_A,pos_A);
				H_D.submat(span(0,ind_n1-1),span(0,ind_n1-1))=D11+D12*D22_inv*(G_D-D22)*D22_inv*(D12.t());
				H_D.submat(span(0,ind_n1-1),span(ind_n1,ind_n-1))=tmp;
				H_D.submat(span(ind_n1,ind_n-1),span(0,ind_n1-1))=tmp.t();
				H_D.submat(span(ind_n1,ind_n-1),span(ind_n1,ind_n-1))=G_D;
				
				}
				if(inverse==true){
				
				arma::Mat<double> G_Dinv;
				G_Dinv=arma::inv(G_D);
				arma::Mat<double> P_Dinv=arma::inv(P_D);
				H_Dinv=P_Dinv.submat(pos_A,pos_A);
				H_Dinv.submat(pos_H22,pos_H22)=H_Dinv.submat(pos_H22,pos_H22)+G_Dinv-D22_inv;				
				}			
				Rcout<<"Complete constructing H_Dominance relationship matrix!"<<endl; 
				return List::create(Named("H") =H_D,_["Hinv"] = H_Dinv);
}