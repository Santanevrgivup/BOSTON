#include <fstream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <ctime>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm> 
// 11
#ifndef _CLOCK_T_DEFINED
typedef long clock_t;
#define _CLOCK_T_DEFINED
#endif
//#define CLOCKS_PER_SEC ((clock_t)1000)

using namespace std;
using namespace arma;
using namespace Rcpp;

//' define a class
//' 
class ComputeQuantiesWithoutCVT {
private:
  // input data
  arma::mat Y;
  arma::mat cell_locs;
  arma::mat D;
  arma::mat tau;
  std::string kernel_type;
  double kernel_param;
  
  // tmp variables
  arma::mat Hinv;
  arma::mat Py;
  arma::mat Hinv_sq;
  arma::vec XtHinvX;
  arma::vec XtHinvX_sq;
  
  // final variables
  arma::vec yPKPy;
  arma::vec trPK;
  arma::vec trPKP;
  arma::vec trPKPK;
  arma::vec trPP;
  arma::vec newInfoM;
  arma::vec ee;
  arma::vec kk;
  arma::vec df;
  
  // parameters
  size_t num_cell;
  size_t num_gene;
  size_t block_size;
  
public:
  // construct function
  ComputeQuantiesWithoutCVT(const arma::mat& Y, const arma::mat& cell_locs, 
                   const arma::mat& D, const arma::mat& tau,
                   const std::string& kernel_type, double kernel_param,
                   size_t block_size = 50)
    : Y(Y), cell_locs(cell_locs), D(D), tau(tau),
      kernel_type(kernel_type), kernel_param(kernel_param),
      block_size(block_size) {
    
    num_cell = Y.n_rows;
    num_gene = Y.n_cols;
    
    // block size
    if (block_size > num_cell) {
      this->block_size = num_cell;
    }
    if (this->block_size < 10) {
      this->block_size = 10;
    }
  }
  
  // main 
  void compute() {
    // initial 
    yPKPy.set_size(num_gene);
    trPK.set_size(num_gene);
    trPKP.set_size(num_gene);
    trPKPK.set_size(num_gene);
    trPP.set_size(num_gene);
    
    yPKPy.zeros();
    trPK.zeros();
    trPKP.zeros();
    trPKPK.zeros();
    trPP.zeros();
    
    //Hinv related quantities
    ComputeHinvQuantities();
    
    // compute
    arma::vec x_sq = pow(cell_locs.col(0), 2.0);
    arma::vec y_sq = pow(cell_locs.col(1), 2.0);
    
    // kernel types
    if (kernel_type == "gaussian") {
      ComputeGaussianKernelQantities(x_sq, y_sq);
    } else if (kernel_type == "linear") {
      ComputeLinearKernelQantities(x_sq, y_sq);
    } else if (kernel_type == "cosine") {
      ComputeCosineKernelQantities(x_sq, y_sq);
    } else {
      std::cerr << "Error: Unknown kernel type: " << kernel_type << std::endl;
      return;
    }
    
    // final results
    ComputeFinalQuantities();
  }
  
  // get results
  const arma::vec& getYPKPy() const { return yPKPy; }
  const arma::vec& getTrPK() const { return trPK; }
  const arma::vec& getTrPKP() const { return trPKP; }
  const arma::vec& getTrPKPK() const { return trPKPK; }
  const arma::vec& getTrPP() const { return trPP; }
  const arma::vec& getNewInfoM() const { return newInfoM; }
  const arma::vec& getEe() const { return ee; }
  const arma::vec& getKk() const { return kk; }
  const arma::vec& getDf() const { return df; }
  
  // block size
  void setBlockSize(size_t new_block_size) {
    block_size = new_block_size;
    if (block_size > num_cell) {
      block_size = num_cell;
    }
    if (block_size < 10) {
      block_size = 10;
    }
  }
  
private:
  // Hinv-related quantities
  void ComputeHinvQuantities() {
    Hinv.set_size(num_cell, num_gene);
    Py.set_size(num_cell, num_gene);
    Hinv_sq.set_size(num_cell, num_gene);
    XtHinvX.set_size(num_gene);
    XtHinvX_sq.set_size(num_gene);
    
    #pragma omp parallel for
    for(size_t ig = 0; ig < num_gene; ig++) {
      arma::vec Hinv_tmp = tau(0, ig) * (1.0 / (D.col(ig) + 1e-5));
      Hinv_tmp += tau(1, ig);
      Hinv.col(ig) = 1.0 / (Hinv_tmp + 1e-5);
      XtHinvX(ig) = sum(Hinv.col(ig));
      
      double dot_Hinv_Y = dot(Hinv.col(ig), Y.col(ig));
      Py.col(ig) = Hinv.col(ig) % Y.col(ig) - Hinv.col(ig) * (dot_Hinv_Y / XtHinvX(ig));
      
      Hinv_sq.col(ig) = Hinv.col(ig) % Hinv.col(ig);
      XtHinvX_sq(ig) = XtHinvX(ig) * XtHinvX(ig);
    }
  }
  
  // shared computation across kernels
  void ComputeStatistics(const arma::mat& HiKHiK, const arma::mat& HinvK, const arma::mat& KPy) {
    clock_t Start = clock();
    arma::vec trHiKHiK(num_gene, arma::fill::zeros);
    arma::vec HinvHinvK(num_gene, arma::fill::zeros);
    
    arma::vec sum_Hinv = sum(Hinv, 0).t();
    arma::vec sum_Hinv_sq_all = sum(Hinv_sq, 0).t();
    
#pragma omp parallel for
    for(size_t ig = 0; ig < num_gene; ig++) {
      const arma::vec& Hinv_col = Hinv.col(ig);
      const arma::vec& Py_col = Py.col(ig);
      const arma::vec& Hinv_sq_col = Hinv_sq.col(ig);
      const arma::vec& HiKHiK_col = HiKHiK.col(ig);
      const arma::vec& HinvK_col = HinvK.col(ig);
      const arma::vec& KPy_col = KPy.col(ig);
      
      double XtHinvX_ig = XtHinvX(ig);
      double XtHinvX_sq_ig = XtHinvX_sq(ig);
      double sum_Hinv_sq = sum_Hinv_sq_all(ig);
      
      yPKPy(ig) = 0.5 * dot(KPy_col, Py_col);
      trHiKHiK(ig) = dot(HiKHiK_col, Hinv_col);
      HinvHinvK(ig) = dot(HinvK_col, Hinv_col);
      trPK(ig) = sum_Hinv(ig) - HinvHinvK(ig) / XtHinvX_ig;
      
      arma::vec tmp_vec = Hinv_sq_col % HinvK_col;
      double sum_tmp_vec = sum(tmp_vec);
      
      trPKP(ig) = sum_Hinv_sq - 2 * sum_tmp_vec / XtHinvX_ig + 
        sum_Hinv_sq * HinvHinvK(ig) / XtHinvX_sq_ig;
      
      arma::vec tmp_vec2 = HinvK_col % Hinv_col;
      double dot_tmp_vec2_HinvK = dot(tmp_vec2, HinvK_col);
      
      trPKPK(ig) = trHiKHiK(ig) - 2.0 * dot_tmp_vec2_HinvK / XtHinvX_ig + 
        HinvHinvK(ig) * HinvHinvK(ig) / XtHinvX_sq_ig;
      
      double dot_Hinv_Hinv_sq = dot(Hinv_col, Hinv_sq_col);
      trPP(ig) = sum_Hinv_sq - 2.0 * dot_Hinv_Hinv_sq / XtHinvX_ig + 
        sum_Hinv_sq * sum_Hinv_sq / XtHinvX_sq_ig;
    }
    
    //std::cout << "Statistical quantities computation time: " 
    //          << (double)(clock() - Start)/CLOCKS_PER_SEC << " (sec.)" << std::endl;
  }
  // Gaussian kernel
  void ComputeGaussianKernelQantities(const arma::vec& x_sq, const arma::vec& y_sq) {
    arma::mat HiKHiK(num_cell, num_gene, arma::fill::zeros);
    arma::mat HinvK(num_cell, num_gene, arma::fill::zeros);
    arma::mat KPy(num_cell, num_gene, arma::fill::zeros);
    
    double sigma_sq = 2.0 * kernel_param * kernel_param;
    
    #pragma omp parallel for schedule(dynamic)
    for(size_t ic_start = 0; ic_start < num_cell; ic_start += block_size) {
      size_t ic_end = ic_start + block_size;
      if (ic_end > num_cell) {
        ic_end = num_cell;
      }
      
      for(size_t ic = ic_start; ic < ic_end; ic++) {
        double x_ic = cell_locs(ic, 0);
        double y_ic = cell_locs(ic, 1);
        double x_ic_sq = x_ic * x_ic;
        double y_ic_sq = y_ic * y_ic;
        
        arma::vec out = x_sq - 2.0 * x_ic * cell_locs.col(0) + x_ic_sq +
          y_sq - 2.0 * y_ic * cell_locs.col(1) + y_ic_sq;
        
        arma::vec tmp_out = exp(-out / sigma_sq);
        arma::vec tmp_out_sq = tmp_out % tmp_out;
        
        HiKHiK.row(ic) = tmp_out_sq.t() * Hinv;
        HinvK.row(ic) = tmp_out.t() * Hinv;
        KPy.row(ic) = tmp_out.t() * Py;
      }
    }// end for
    
    // stats
    ComputeStatistics(HiKHiK, HinvK, KPy);
  }// end gaussian kernel
  
  // linear kernel
  void ComputeLinearKernelQantities(const arma::vec& x_sq, const arma::vec& y_sq) {
    arma::mat HiKHiK(num_cell, num_gene, arma::fill::zeros);
    arma::mat HinvK(num_cell, num_gene, arma::fill::zeros);
    arma::mat KPy(num_cell, num_gene, arma::fill::zeros);
    
    #pragma omp parallel for
    for(size_t ic = 0; ic < num_cell; ic++) {
      double x_ic = cell_locs(ic, 0);
      double y_ic = cell_locs(ic, 1);
      
      arma::vec tmp_out = (x_ic * cell_locs.col(0) + y_ic * cell_locs.col(1)) * 0.5;
      arma::vec tmp_out_sq = tmp_out % tmp_out;
      
      HiKHiK.row(ic) = tmp_out_sq.t() * Hinv;
      HinvK.row(ic) = tmp_out.t() * Hinv;
      KPy.row(ic) = tmp_out.t() * Py;
    }
    
    // stats
    ComputeStatistics(HiKHiK, HinvK, KPy);
  }
  
  // cosine kernel
  void ComputeCosineKernelQantities(const arma::vec& x_sq, const arma::vec& y_sq) {
    arma::mat HiKHiK(num_cell, num_gene, arma::fill::zeros);
    arma::mat HinvK(num_cell, num_gene, arma::fill::zeros);
    arma::mat KPy(num_cell, num_gene, arma::fill::zeros);
    
    double two_pi_over_param = 2.0 * arma::datum::pi / kernel_param;
    
    #pragma omp parallel for
    for(size_t ic = 0; ic < num_cell; ic++) {
      //double x_ic = cell_locs(ic, 0);
      //double y_ic = cell_locs(ic, 1);
      
      //arma::vec out_sqrt = arma::sqrt(x_sq - 2.0 * x_ic * cell_locs.col(0) + x_ic*x_ic +
      //  y_sq - 2.0 * y_ic * cell_locs.col(1) + y_ic*y_ic);
      //cout << out_sqrt <<endl;
      arma::rowvec v1 = cell_locs.row(ic);
      arma::vec out_sqrt = arma::sqrt(pow(v1(0) - cell_locs.col(0), 2.0) + pow(v1(1) - cell_locs.col(1), 2.0));
      arma::vec tmp_out = arma::cos(two_pi_over_param * out_sqrt);
      //if(ic==0){cout<<"K (cosine) = "<< tmp_out.head(20) <<endl;
      arma::vec tmp_out_sq = tmp_out % tmp_out;
      
      HiKHiK.row(ic) = tmp_out_sq.t() * Hinv;
      HinvK.row(ic) = tmp_out.t() * Hinv;
      KPy.row(ic) = tmp_out.t() * Py;
    }
    
    // stats
    ComputeStatistics(HiKHiK, HinvK, KPy);
  }
  
  // final quantities
  void ComputeFinalQuantities() {
    newInfoM.set_size(num_gene);
    ee.set_size(num_gene);
    kk.set_size(num_gene);
    df.set_size(num_gene);
    
    arma::vec newInfoM_p1 = 0.5 * trPKPK;
    newInfoM = newInfoM_p1 - 0.5 * (trPKP % trPKP) / trPP;
    ee = 0.5 * trPK;
    kk = 0.5 * newInfoM / ee;
    df = 2.0 * (ee % ee) / newInfoM;
  }
};

//' Fast compute the testing quantities without covariates, double format
//' input all genes at a time
//' include the loop number of gene in Cpp file
//' @param Yin Working vector
//' @param Din Weight for each gene
//' @param locationin Location for each cell
//' @param kernel_paramin Gaussian kernel parameters
//' @param tauin Initial value for variance component
//' @param kernel_typein The kernel type, including 'gaussian', 'linear' and 'cosine'
//' 
//' @return A list
//' 
//' 
//' @export
// [[Rcpp::export]]
SEXP fastComputeQuantitiesWithoutCVT(SEXP Yin, SEXP locationin, SEXP kernel_paramin, SEXP Din, SEXP tauin, SEXP kernel_typein) {
  try	{
    // variable format convertion
    arma::mat Y = as<arma::mat>(Yin); // a matrix for all genes
    arma::mat cell_locs = as<arma::mat>(locationin);
    //arma::vec kernel_param = as<arma::vec>(kernel_paramin);
    double kernel_param = Rcpp::as<double>(kernel_paramin);
    arma::mat D = as<arma::mat>(Din); // a matrix for all genes
    arma::mat tau = as<arma::mat>(tauin); // a matrix for all genes
    const string kernel_type = Rcpp::as<string>(kernel_typein);
    
    size_t block_size = 50;
    
    
    // definite variables
    const size_t num_cell = Y.n_rows;
    const size_t num_gene = Y.n_cols;
    //const int num_kernel = kernel_param.n_elem;
    
    
    // define object
    ComputeQuantiesWithoutCVT calculator(Y, cell_locs, D, tau, kernel_type, kernel_param, block_size);
    
    // compute
    calculator.compute();
    
    const arma::vec& yPKPy = calculator.getYPKPy();
    //const arma::vec& trPK = calculator.getTrPK();
    //const arma::vec& trPKP = calculator.getTrPKP();
    const arma::vec& trPKPK = calculator.getTrPKPK();
    //const arma::vec& trPP = calculator.getTrPP();
    //const arma::vec& newInfoM = calculator.getNewInfoM();
    const arma::vec& ee = calculator.getEe();
    const arma::vec& kk = calculator.getKk();
    const arma::vec& df = calculator.getDf();
    // return values
    return List::create(Named("S0") = yPKPy, Named("ee") = ee, Named("infoMp1") = 0.5*trPKPK, Named("df") = df, Named("kk") = kk);
  } // end try
  catch (std::exception &ex)
  {
    forward_exception_to_r(ex);
  }
  catch (...)
  {
    ::Rf_error("C++ exception (unknown reason)...");
  }
  return R_NilValue;
} // end funcs

//' Compute the testing quantities without covariates, float format
//' @param yin Working vector
//' @param Pyin The vector P*y
//' @param cov_matin Kernel matrix to be tested
//' @param Din Weight for each gene
//' @param tauin Initial value for variance component
//' 
//' @return A list
//' 
//' 
//' @export
// [[Rcpp::export]]
SEXP ComputeTestQuantRcpp_nocov(SEXP yin, SEXP Pyin, SEXP cov_matin, SEXP Din, SEXP tauin){
	try
	{
		arma::vec y = as<arma::vec>(yin);
		arma::vec Py = as<arma::vec>(Pyin);
		arma::mat cov_mat = as<arma::mat>(cov_matin);
		arma::vec D = as<arma::vec>(Din);
		arma::vec tau = as<arma::vec>(tauin);

		const int num_cell = y.n_elem;
		arma::vec Hinv(num_cell);
		arma::vec one_vec = ones<arma::vec>(num_cell);

		Hinv = tau(0) * (1.0 / (D + 1e-5));
		Hinv += tau(1) * one_vec;
		Hinv = 1.0 / (Hinv + 1e-5); // Hinv is a diagonal matrix
		//arma::vec Hinvy = Hinv % y;

		arma::vec HinvX = Hinv;
		double XtHinvX = sum(HinvX);
	
		arma::mat P = - arma::kron(HinvX, HinvX.t())/XtHinvX;
		P.diag() = P.diag() + Hinv;

		arma::rowvec PKp2 = HinvX.t()*cov_mat;

		arma::mat PK = cov_mat.each_col() % HinvX - arma::kron(HinvX, PKp2)/XtHinvX;
		cout<<"old trace(PK_p1)="<<trace(cov_mat.each_col() % HinvX)<<endl;
		cout<<"old trace(PK_p2)="<<trace(arma::kron(HinvX, PKp2))<<endl;
		cout<<"old trace(PP)="<< accu(P%P)<<endl;
		cout<<"old trace(PK)="<< trace(PK)<<endl;
		double trace_PKP = accu(PK % P);
		cout<<"old trace(PKP)="<< trace_PKP<<endl;
		cout<<"old trace(PKPK)="<< trace(PK * PK)<<endl;
		double newInfoM_p1 = 0.5 * trace(PK * PK);
		cout<<"newInfoM_p1 ="<< newInfoM_p1<<endl;
		
		double newInfoM = newInfoM_p1 - 0.5 * trace_PKP*trace_PKP/accu(P % P);
		double ee = trace(PK) / 2.0;
		double kk = newInfoM / (2.0 * ee);
		double df = 2.0 * ee * ee / newInfoM;

		arma::vec PKPy = PK * Py;

		double S0 = 0.5 * dot(y, PKPy);
		cout<<"S0 = " << S0 <<endl;
		// return values
		return List::create(Named("S0") = S0, Named("ee") = ee, Named("infoMp1") = newInfoM_p1, Named("df") = df, Named("kk") = kk);
	} // end try
	catch (std::exception &ex)
	{
		forward_exception_to_r(ex);
	}
	catch (...)
	{
		::Rf_error("C++ exception (unknown reason)...");
	}
	return R_NilValue;
} // end funcs


//' Variance component estimation with covariates using Average Information algorithm
//' @param Yin Working vector
//' @param Xin Covariate matrix
//' @param Din Weight for each gene
//' @param tauin Initial value for variance component
//' @param fixtauin Variance component to be optimized
//' @param tolin Tolerance
//' 
//' @return A list
//' 
//' @export
// [[Rcpp::export]]
SEXP CovariatesAI(SEXP Yin, SEXP Xin, SEXP Din, SEXP tauin, SEXP fixtauin, SEXP tolin) { /*Average Information*/
  try	{
    arma::vec Y = as<arma::vec>(Yin);
    arma::mat X = as<arma::mat>(Xin);
    arma::vec D = as<arma::vec>(Din);
    arma::vec tau = as<arma::vec>(tauin);
    const uvec fixtau = as<uvec>(fixtauin);
    const int num_cov_mat2 = sum(fixtau == 0);
    const double tol = Rcpp::as<double>(tolin);
    uvec ZERO = (tau < tol);
    
    const int num_cell = X.n_rows;
    const int num_cvt = X.n_cols; // if number of column X isnot equal to 1
    arma::vec Hinv(num_cell);
    arma::vec one_vec = ones<arma::vec>(num_cell);
    
    Hinv = tau(0) * (1.0 / (D + 1e-5));
    Hinv += tau(1) * one_vec;
    Hinv = 1.0 / (Hinv + 1e-5);
    arma::vec HinvY = Hinv % Y;
    arma::mat HinvX = X.each_col() % Hinv;
    arma::mat XtHinvX = X.t() * HinvX;
    arma::mat XtHinvX_inv = inv_sympd(XtHinvX);
    
    arma::mat P = diagmat(Hinv) - HinvX * XtHinvX_inv * HinvX.t();
    
    arma::vec alpha = XtHinvX_inv * HinvX.t() * Y;
    arma::vec eta = Y - tau(0) * (HinvY - HinvX * alpha) / D;
    arma::vec PY = P * Y;
    
    if (num_cov_mat2 > 0) {
      const uvec idxtau = find(fixtau == 0);
      arma::mat AImat(num_cov_mat2, num_cov_mat2); //average information matrix
      //arma::vec PY = P * Y;
      arma::vec score(num_cov_mat2), PAPY;
      for (size_t i = 0; i < num_cov_mat2; i++) {
        PAPY = P * PY;
        score(i) = dot(Y, PAPY) - sum(P.diag());
        for (size_t j = 0; j <= i; j++)	{
          AImat(i, j) = dot(PY, PAPY);
          if (j != i)	{
            AImat(j, i) = AImat(i, j);
          } // end fi
        }	 //end for j
      }		  // end for i
      
      arma::vec Dtau = solve(AImat, score);
      arma::vec tau0 = tau;
      
      tau.elem(idxtau) = tau0.elem(idxtau) + Dtau;
      
      tau.elem(find(ZERO % (tau < tol))).zeros();
      double step = 1.0;
      while (any(tau < 0.0)) {
        step *= 0.5;
        tau.elem(idxtau) = tau0.elem(idxtau) + step * Dtau;
        tau.elem(find(ZERO % (tau < tol))).zeros();
      }
      tau.elem(find(tau < tol)).zeros();
    } // end fi
    // boundary tau 0<= tau <=10
    // tau.elem(find(tau >10.0)).ones();
    // return values
    return List::create(Named("tau") = tau, Named("P") = P, Named("cov") = XtHinvX_inv,	Named("alpha") = alpha, Named("Py") = PY, Named("eta") = eta);
  }
  catch (std::exception &ex)
  {
    forward_exception_to_r(ex);
  }
  catch (...)
  {
    ::Rf_error("C++ exception (unknown reason)...");
  }
  return R_NilValue;
} // end funcs


//' Variance component estimation without covariates using Average Information algorithm, float format
//' @param Yin Working vector
//' @param Xin Covariate matrix
//' @param Din Weight for each gene
//' @param tauin Initial value for variance component
//' @param fixtauin Variance component to be optimized
//' @param tolin Tolerance
//' 
//' @return A list
//' 
//' 
//' @export
// [[Rcpp::export]]
SEXP noCovariatesAI(SEXP Yin, SEXP Xin, SEXP Din, SEXP tauin, SEXP fixtauin, SEXP tolin) { /*Average Information*/
  try	{
    arma::vec Y = as<arma::vec>(Yin);
    arma::mat X = as<arma::mat>(Xin);
    arma::vec D = as<arma::vec>(Din);
    arma::vec tau = as<arma::vec>(tauin);
    const uvec fixtau = as<uvec>(fixtauin);
    const int num_cov_mat2 = sum(fixtau == 0);
    const double tol = Rcpp::as<double>(tolin);
    uvec ZERO = (tau < tol);
    
    const int num_cell = X.n_rows;
    const int num_cvt = X.n_cols; // only suitable for intercept case
    
    arma::vec Hinv(num_cell);
    arma::vec one_vec = ones<arma::vec>(num_cell);
    
    
    Hinv = tau(0) * (1.0 / (D + 1e-5));
    Hinv += tau(1) * one_vec;
    Hinv = 1.0 / (Hinv + 1e-5);
    
    arma::vec HinvY = Hinv % Y;
    
    arma::vec HinvX = Hinv;
    double XtHinvX = sum(HinvX);
    double XtHinvX_inv = 1.0 / XtHinvX;
    arma::vec P_diag = Hinv - (HinvX % HinvX) * XtHinvX_inv;
    double alpha = XtHinvX_inv * dot(HinvX, Y);
    arma::vec eta = Y - tau(0) * (HinvY - HinvX * alpha) / D;
    
    arma::vec PY = HinvY - HinvX * XtHinvX_inv * (HinvX.t() * Y);
    
    
    if (num_cov_mat2 > 0) {
      const uvec idxtau = find(fixtau == 0);
      arma::mat AImat(num_cov_mat2, num_cov_mat2); //average information matrix
      //arma::vec PY = P * Y;
      arma::vec score(num_cov_mat2);
      for (size_t i = 0; i < num_cov_mat2; i++) {
        
        arma::vec PAPY = Hinv % PY - HinvX * XtHinvX_inv * (HinvX.t() * PY);
        
        score(i) = dot(Y, PAPY) - sum(P_diag);
        for (size_t j = 0; j <= i; j++)	{
          AImat(i, j) = dot(PY, PAPY);
          if (j != i)	{
            AImat(j, i) = AImat(i, j);
          } // end fi
        }	 //end for j
      }		  // end for i
      
      arma::vec Dtau = solve(AImat, score);
      arma::vec tau0 = tau;
      
      tau.elem(idxtau) = tau0.elem(idxtau) + Dtau;
      tau.elem(find(ZERO % (tau < tol))).zeros();
      double step = 1.0;
      while (any(tau < 0.0)) {
        step *= 0.5;
        tau.elem(idxtau) = tau0.elem(idxtau) + step * Dtau;
        tau.elem(find(ZERO % (tau < tol))).zeros();
      } // end while
      tau.elem(find(tau < tol)).zeros();
    } // end fi
    // boundary tau 0<= tau <=10
    //tau.elem(find(tau >10.0)).ones();
    // return values
    
    return List::create(Named("tau") = tau, Named("Py") = PY, Named("cov") = XtHinvX_inv,	Named("alpha") = alpha, Named("eta") = eta);
  } // end try
  catch (std::exception &ex)
  {
    forward_exception_to_r(ex);
  }
  catch (...)
  {
    ::Rf_error("C++ exception (unknown reason)...");
  }
  return R_NilValue;
} // end funcs


//' Compute the testing quantities with covariates, fast trace term
//' @param Xin covariates
//' @param Pyin The vector P*y
//' @param Kin Kernel matrix to be tested
//' @param Din Weight for each gene
//' @param tauin Initial value for variance component
//' 
//' @return A list
//' 
//' 
//' @export
// [[Rcpp::export]]
SEXP FastTraceComputeTestQuantRcpp_cov(SEXP Pyin, SEXP Xin, SEXP Kin, SEXP Din, SEXP tauin){
  try
  {
    arma::vec Py = as<arma::vec>(Pyin); // n-dim vector
    arma::mat K = as<arma::mat>(Kin); // nxn matrix
    arma::vec D = as<arma::vec>(Din);// n-dim vector
    arma::mat X = as<arma::mat>(Xin);
    arma::vec tau = as<arma::vec>(tauin); // 2-dim vector
    arma::vec Hinv(D.n_elem);
    arma::vec one_vec = ones<arma::vec>(D.n_elem);
    
    Hinv = tau(0) * (1.0 / (D + 1e-5));
    Hinv += tau(1) * one_vec;
    Hinv = 1.0 / (Hinv + 1e-5); // Hinv is a diagonal matrix
    //arma::vec Hinvy = Hinv % y;
    
    //arma::vec HinvX = Hinv;
    //double XtHinvX = arma::sum(HinvX);
    arma::mat HinvX = X.each_col() % Hinv; // nxd
    arma::mat Hinv2X = X.each_col() % (Hinv%Hinv); // nxd
    arma::mat XtHinvX = X.t() * HinvX;
    arma::mat XtHinv2X = X.t() * Hinv2X;
    arma::mat XtHinvX_inv = inv_sympd(XtHinvX);
    arma::mat HinvXXtHinvX_inv = HinvX*XtHinvX_inv; // nxd
    arma::mat HinvXXtHinvX_invXtHinv2X = HinvXXtHinvX_inv*XtHinv2X;
    // we can reduce the computational burden due to only consider the trace of matrix term, edited by sun, 2021-1-18 11:54:42
    // required quantities
    //double trK = arma::trace(K); // time-cost step, using pre-defined value
    //double trKK = arma::accu(K%K); // time-cost step, using pre-defined value
    //double trHinv = arma::sum(Hinv);
    
    arma::mat KHinvX = K*HinvX; // time-cost step, nxd
    arma::vec KPy = K*Py; // time-cost step/ nx1
    arma::mat HinvK = K.each_col()%Hinv; // nxd
    arma::mat HinvKHinvX = KHinvX.each_col()%Hinv; // nxd
    arma::mat XHinvKHinvX = HinvX.t()*KHinvX; // dxd
    arma::mat HinvHinvXXtHinvX_inv = HinvXXtHinvX_inv.each_col() % Hinv;
    arma::mat HinvXXtHinvX_invXHinvKHinvX = HinvXXtHinvX_inv*XHinvKHinvX;
    arma::mat HinvKHinvXXtHinvX_inv = HinvKHinvX*XtHinvX_inv;
    arma::mat KHinvXXtHinvX_inv = KHinvX*XtHinvX_inv;
    // compute trPK
    double trPK = arma::dot(Hinv,K.diag()) - arma::accu(HinvXXtHinvX_inv % KHinvX);
    //cout<<"new trace(PK_p1)="<< arma::dot(Hinv,K.diag())<<endl;
    //cout<<"new trace(PK_p2)="<< trHinvXKHinvX<<endl;
    // compute trPP
    double trPP = arma::dot(Hinv, Hinv) - 2*arma::accu(HinvHinvXXtHinvX_inv % HinvX) + arma::accu(HinvXXtHinvX_invXtHinv2X%HinvXXtHinvX_inv);
    //cout<<"new trace(PP)="<< trPP<<endl;
    // compute trPKP
    double trPKP = arma::dot(Hinv%Hinv, K.diag()) - 2*arma::accu(HinvXXtHinvX_inv%HinvKHinvX) + arma::accu(HinvXXtHinvX_invXHinvKHinvX%HinvXXtHinvX_inv);
    //cout<<"new trace(PKP)="<< trPKP<<endl;
    // compute trPKPK
    double trPKPK = arma::accu(HinvK%HinvK.t()) - 2*arma::accu(HinvKHinvXXtHinvX_inv%KHinvX) + arma::accu(HinvXXtHinvX_invXHinvKHinvX%KHinvXXtHinvX_inv);
    //cout<<"new trace(PKPK)="<< trPKPK<<endl;
    
    // testing quantities
    double im = 0.5*trPKPK - 0.5*trPKP*trPKP/trPP;
    double ee = trPK / 2.0;
    double kk = im / (2.0 * ee);
    double df = 2.0 * ee * ee / im;
    double ss = 0.5*arma::dot(Py, KPy);
    //cout<<"S0 = " << ss <<endl;
    // return values
    return List::create(Named("S0") = ss, Named("ee") = ee, Named("infoMp1") = 0.5*trPKPK, Named("df") = df, Named("kk") = kk);
  } // end try
  catch (std::exception &ex)
  {
    forward_exception_to_r(ex);
  }
  catch (...)
  {
    ::Rf_error("C++ exception (unknown reason)...");
  }
  return R_NilValue;
} // end funcs

//' Compute the testing quantities with covariates, float format
//' @param yin Working vector
//' @param Pyin The vector P*y
//' @param Xin Covariate matrix, including the intercept
//' @param cov_matin Kernel matrix to be tested
//' @param Din Weight for each gene
//' @param tauin Initial value for variance component
//' 
//' @return A list
//' 
//' 
//' @export
// [[Rcpp::export]]
SEXP ComputeTestQuantRcpp_cov(SEXP yin, SEXP Pyin, SEXP Xin, SEXP cov_matin, SEXP Din, SEXP tauin) {
  try
  {
    arma::vec y = as<arma::vec>(yin);
    arma::vec Py = as<arma::vec>(Pyin);
    arma::mat cov_mat = as<arma::mat>(cov_matin);
    arma::mat X = as<arma::mat>(Xin);
    arma::vec D = as<arma::vec>(Din);
    arma::vec tau = as<arma::vec>(tauin);
    
    const int num_cell = y.n_elem;
    arma::vec Hinv(num_cell);
    arma::vec one_vec = ones<arma::vec>(num_cell);
    
    Hinv = tau(0) * (1.0 / (D + 1e-5));
    Hinv += tau(1) * one_vec;
    Hinv = 1.0 / (Hinv + 1e-5); // Hinv is a diagonal matrix
    //arma::vec Hinvy = Hinv % y;
    arma::mat HinvX = X.each_col() % Hinv;
    arma::mat XtHinvX = X.t() * HinvX;
    arma::mat XtHinvX_inv = inv_sympd(XtHinvX);
    arma::mat P = diagmat(Hinv) - HinvX * XtHinvX_inv * HinvX.t();
    
    // modified by sun, 2019-4-13 16:25:06
    arma::mat PK = P*cov_mat;
    double trace_PKP = accu(PK % P);
    
    // modified by sun, 2019-4-9 12:26:03
    double newInfoM_p1 = 0.5 * trace(PK * PK);
    double newInfoM = newInfoM_p1 - 0.5 * trace_PKP*trace_PKP/accu(P % P);
    double ee = trace(PK) / 2.0;
    double kk = newInfoM / (2.0 * ee);
    double df = 2.0 * ee * ee / newInfoM;
    arma::vec PKPy = PK * Py;
    
    double S0 = 0.5 * dot(y, PKPy);
    double ll = 0.0;
    
    // return values
    return List::create(Named("S0") = S0, Named("ee") = ee, Named("infoMp1") = newInfoM_p1, Named("df") = df, Named("kk") = kk);
  } // end try
  catch (std::exception &ex)
  {
    forward_exception_to_r(ex);
  }
  catch (...)
  {
    ::Rf_error("C++ exception (unknown reason)...");
  }
  return R_NilValue;
} // end funcs


//' Compute Gaussian kernel parameters
//' @param xin spatial location
//' 
//' @return A list
//' 
//' 
//' @export
// [[Rcpp::export]]
SEXP ComputeKernelParamLessMem(SEXP xin) {
  try
  {
    arma::mat x = as<arma::mat>(xin);
    
    unsigned int outrows = x.n_rows;
    
    double out = 0.0, lmin = 9999.0, lmax = 0.0;
    for (size_t i = 0; i < outrows - 1; i++) {
      arma::rowvec v1 = x.row(i);
      for (size_t j = i + 1; j < outrows; j++) {
        out = sqrt(sum(pow(v1 - x.row(j), 2.0)));
        //out = sum(pow(v1 - x.row(j), 2.0));
        if(out > 1e-8){
          if(out < lmin){
            lmin = out;
          }else if( out > lmax){
            lmax = out;
          }
        }
      }// end for j
    }// end for i
    //arma::vec sub_out = out.elem( find(out>1e-6) );
    //double lmin = sub_out.min()/2;
    //double lmax = sub_out.max()*2;
    lmin = lmin/2;
    lmax = 2*lmax;
    arma::vec kparam = exp10( linspace(log10(lmin), log10(lmax), 10) );
    // return values
    return List::create(Named("kparam") = kparam);
  } // end try
  catch (std::exception &ex)
  {
    forward_exception_to_r(ex);
  }
  catch (...)
  {
    ::Rf_error("C++ exception (unknown reason)...");
  }
  return R_NilValue;
} // end funcs



//' Compute the testing quantities for linear mixed model, float format
//' @param covariatesin Covariates matrix
//' @param Xdaggerin The inverse matrix of covariates
//' @param norm_countsin normalized count matrix
//' @param cov_matin Kernel matrix to be testedt
//' 
//' @return A list
//' 
//' 
//' @export
// [[Rcpp::export]]
SEXP ComputeTestQuantRcpp_RLSKAT(SEXP cov_matin, SEXP Xdaggerin, SEXP covariatesin, SEXP norm_countsin){
  try {
    arma::vec Xdagger = as<arma::vec>(Xdaggerin);
    arma::vec covariates = as<arma::vec>(covariatesin);
    arma::mat cov_mat = as<arma::mat>(cov_matin);
    arma::mat norm_counts = as<arma::mat>(norm_countsin);
    
	// calculate SKS
	arma::rowvec SKSp1 = Xdagger.t()*cov_mat;
	arma::mat SKSp2 = cov_mat - arma::kron(covariates, SKSp1);
	arma::rowvec SKSp3 = Xdagger.t()*SKSp2;
	arma::mat SKS = SKSp2 - arma::kron(covariates, SKSp3);

	arma::vec eigval;
	arma::mat eigvec;
	//arma::eig_sym(eigval, eigvec, SKS);
	arma::eig_sym(eigval, eigvec, SKS);

	// svd for B
	arma::mat B = join_rows(cov_mat,covariates);
	arma::mat U;
	arma::vec d;
	arma::mat V;
	arma::svd(U,d,V,B);

	// calculate nominators
	arma::mat SKSc = (SKS*norm_counts) %norm_counts;
	arma::rowvec nominators = arma::sum(SKSc, 0);

	//calculate denonimators
	arma::rowvec norm_countsp1 = Xdagger.t()*norm_counts;
	arma::mat tmpmat = arma::kron(covariates, norm_countsp1);
	arma::mat cde = norm_counts - norm_counts%tmpmat;
    arma::rowvec denonimators = arma::sum(cde, 0);

    // return values
    return List::create(Named("eigval")=eigval, Named("d")=d, Named("nominators")=nominators, Named("denonimators")=denonimators);
  } // end try
  catch (std::exception &ex)
  {
    forward_exception_to_r(ex);
  }
  catch (...)
  {
    ::Rf_error("C++ exception (unknown reason)...");
  }
  return R_NilValue;
} // end funcs


///////////////////////////////////////////////////////////////////////////////////////////
////                             CODE END HERE                                           //
///////////////////////////////////////////////////////////////////////////////////////////
