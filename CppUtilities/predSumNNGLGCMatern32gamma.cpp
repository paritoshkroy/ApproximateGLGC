#include <RcppArmadillo.h>
#include <ComputeDistMat.cpp>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export()]]

arma::mat predSumNNGLGCMatern32gamma(arma::colvec obsY, arma::mat obsX, arma::mat prdX, 
                          arma::mat obsCoords, arma::mat prdNeiDist, 
                          arma::umat prdNeiID, arma::mat z,
                          arma::mat beta, arma::vec sigmasq1, arma::vec lscale1,
                          arma::vec sigmasq2, arma::vec lscale2,
                          arma::vec tausq, arma::vec gamma,
                          unsigned int iterprint)
{
  unsigned int L = beta.n_cols;     // number of posterior samples
  unsigned int N0 = prdNeiID.n_rows;   // number of prediction locations
  int K = prdNeiID.n_cols;   // number of nearest neighbor
  arma::mat ypredsummary = zeros(N0,5);
  arma::vec ypred = zeros(L);
  arma::vec quant = zeros(3);
  arma::vec prob = {0.025,0.5,0.975};
  arma::mat thisCoords = zeros(K,2);
  arma::mat thisDist, thisNeiDist, C1Nei, C2Nei, thisNeiXb, thisz;
  arma::rowvec Weight1, Weight2;
  arma::rowvec c10, c20, thisXdist;
  arma::urowvec thisNeiID;
  arma::colvec thisObsY;
  double m1, m2, v1, v2, z0;
  arma::umat NeiID = prdNeiID - 1;  // index in C is equal to (index in R - 1)
  
  for(unsigned int i=0; i<N0; i++) { // i index for the number of prediction locations
    thisNeiID = NeiID.row(i);
    thisCoords = obsCoords.rows(thisNeiID);
    thisNeiDist = ComputeDistMat(thisCoords);
    thisXdist = prdNeiDist.row(i);
    thisObsY = obsY(thisNeiID);
    thisNeiXb = obsX.rows(thisNeiID) * beta; //obsXb.rows(thisNeiID);
    thisz = z.rows(thisNeiID);
    
    for(unsigned int l=0; l<L; l++) { // l index for the number of MCMC samples
      
      c10 = sigmasq1(l)* ((1+1.732051*thisXdist/lscale1(l)) % exp(-1.732051*thisXdist/lscale1(l)));
      C1Nei = sigmasq1(l)*((1+1.732051*thisNeiDist/lscale1(l)) % exp(-1.732051*thisNeiDist/lscale1(l)));
      
      Weight1 = c10 * inv(C1Nei);
      m1 = as_scalar(Weight1 * thisz.col(l));
      v1 = as_scalar(sigmasq1(l) - Weight1 * c10.t());
      z0 = m1 + pow(v1,0.5) * randn();
      
      // for memmory allocation c0, Sig, Weight, m, v are kept the same
      // these are the same dimensional
      
      c20 = sigmasq2(l)* ((1+1.732051*thisXdist/lscale2(l)) % exp(-1.732051*thisXdist/lscale2(l)));
      C2Nei = sigmasq2(l)*((1+1.732051*thisNeiDist/lscale2(l)) % exp(-1.732051*thisNeiDist/lscale2(l))) + tausq(l)*eye(K,K);
      
      Weight2 = c20 * inv(C2Nei);
      m2 = as_scalar(Weight2 * (thisObsY - thisNeiXb.col(l) - gamma(l) * exp(thisz.col(l))));
      v2 = as_scalar(sigmasq2(l) + tausq(l) - Weight2 * c20.t());
      ypred(l) = as_scalar(prdX.row(i)*beta.col(l)) + gamma(l) * exp(z0) + m2 + pow(v2, 0.5) * randn();
    }
    quant = quantile(ypred,prob);
    ypredsummary(i,0) = mean(ypred);
    ypredsummary(i,1) = stddev(ypred);
    ypredsummary(i,2) = quant(0);
    ypredsummary(i,3) = quant(1);
    ypredsummary(i,4) = quant(2);
    
    // Monitoring the processes
    if(i%iterprint==0){
      Rcout << "Prediction upto the " << i <<"th locations is completed" << std::endl;
    }
    
  }
  Rcout << std::endl;
  Rcout << std::endl;
  Rcout << "Output is a matrix with columns mean, sd, q2.5, q50, q97.5, respectively." << std::endl;
  return ypredsummary;
}
