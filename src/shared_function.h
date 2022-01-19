// [[Rcpp::depends(BH, bigmemory)]]
//[[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppProgress)]]


#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <omp.h>
#include <stdio.h>
#include <numeric>


using namespace Rcpp;
using namespace std;
using namespace arma;

void  cumulativeSum(const std::vector<int> input, std::vector<int>& result);

std::string paste_vec_short( arma::Col<int> & data);

std::vector<int> get_haplotype_set_short( arma::Mat<int> & data_haplotype);

void allele_get_haplotype_set_short( arma::Mat<int> & data_haplotype, Rcpp::List & haplotype_allele, int i_pos);

List define_block_window_kb_cpp(std::vector<int> block1,std::vector<int> block2,std::vector<int> tmp_data_map_3);

std::string single_base_factor_cpp(CharacterVector data, std::string miss_base);

std::string pair_base_factor_cpp(CharacterVector data, std::string miss_base);


std::string	get_blupf90_allele_string_phased_vcf(arma::Col<int> allele_string_col1, arma::Col<int> allele_string_col2,
												 int max_length,std::vector<int> cumsum_haplo_type_num, int n_snp);
												 
std::string	get_blupf90_allele_string_numeric(arma::Row<int> allele_string_row,int max_length);

SEXP make_bigmemory_object_cpp(int nrow, int ncol,std::string file_name,std::string file_path, std::string type);

SEXP make_bigmemory_object_address_cpp(int nrow, int ncol,std::string file_name,std::string file_path, std::string type);	

arma::Mat<int> NumericMatrix_to_arma(NumericMatrix & data_numeric);												 

arma::Mat<int> DataFrame_to_arma(DataFrame & data_numeric);	

IntegerVector get_offspring_generation_cpp(DataFrame ped, CharacterVector IND_base);

List single_pedigree_cpp(CharacterMatrix ped);

												 
std::string	get_blupf90_allele_string_unphased_haplotype(arma::Col<int> allele_string_col1, arma::Col<int> allele_string_col2,
												         int max_length);												 

CharacterVector union_cpp(CharacterVector X, CharacterVector Y);

void delete_bigmemory_file_cpp(std::string  matrix_type,
								std::string bigmemory_data_name,
						       std::string bigmemory_data_path,
							   bool message);


arma::Mat<int> get_allele(arma::Mat<int> Pedigree);

arma::Mat<double> matrix_col3_old(arma::Mat<double> & G,arma::Col<int> IND_geno, bool det, int cpu_cores, double threshold);

SEXP matrix_col3_memory_old(SEXP pBigMat,std::string bigmemory_data_name,std::string bigmemory_data_path,
						arma::Col<int> IND_geno, bool det, int cpu_cores, double threshold);

SEXP matrix_col3_memory_alt(arma::Mat<double> G,std::string bigmemory_data_name,std::string bigmemory_data_path,
						arma::Col<int> IND_geno, bool det, int cpu_cores, double threshold);

arma::Mat<double> matrix_col3(arma::Mat<double> & G,arma::Col<int> IND_geno, bool det, int cpu_cores, double threshold);

SEXP matrix_col3_memory(SEXP pBigMat,std::string bigmemory_data_name,std::string bigmemory_data_path,
						arma::Col<int> IND_geno, bool det, int cpu_cores, double threshold);
