#include <TMB.hpp>
#include <string>
using namespace tmbutils;

template<class Type>
struct list_SparseMatrix_from_R : vector<SparseMatrix<Type> > {
  list_SparseMatrix_from_R(SEXP x){
    int len_x = LENGTH(x);
    (*this).resize(len_x);
    for (int index = 0; index < len_x; index++){
      SEXP cur = VECTOR_ELT(x, index);
      (*this)(index) = asSparseMatrix<Type>(cur);
    }
  }
};

template<class Type>
struct list_Scalar_from_R : vector<Type> {
  list_Scalar_from_R(SEXP x){
    int len_x = LENGTH(x);
    (*this).resize(len_x);
    for (int index = 0; index < len_x; index++){
      SEXP cur = VECTOR_ELT(x, index);
      (*this)(index) = *REAL(cur);
    }
  }
};

template<class Type>
Type objective_function<Type>::operator() ()
{
  // For random effects
  DATA_STRUCT(X, list_SparseMatrix_from_R); // Design matrix
  DATA_STRUCT(B, list_SparseMatrix_from_R); // Design matrix
  DATA_STRUCT(P, list_SparseMatrix_from_R); // Penalty matrix
  DATA_STRUCT(logPdet, list_Scalar_from_R); // Determinant of (fixed) penalty matrix

  // The last element of u and the last element of alpha are for the
  // Gaussian Family. (An extra variance parameter)
  DATA_STRUCT(u, list_Scalar_from_R);       // Pc prior, u param
  DATA_STRUCT(alpha, list_Scalar_from_R);   // Pc prior, alpha param
  DATA_STRUCT(betaprec, list_Scalar_from_R);// For boundary, beta ~iid N(0,1/betaprec)

  // For fixed effects
  DATA_STRUCT(beta_fixed_prec, list_Scalar_from_R);
  DATA_STRUCT(Xf, list_SparseMatrix_from_R);

  DATA_VECTOR(y); //response variable
  
  int betadim1 = X(0).cols(); // Number of boundary conditions in RE1
  int d1 = P(0).cols(); // Number of spline in RE1

  int betadim2 = X(1).cols(); // Number of boundary conditions in RE2
  int d2 = P(1).cols(); // Number of spline in RE2

  int beta_fixed_dim0 = Xf(0).cols(); // Number of fixed effect0, should be ONE!! (Since it is the intercept)
  int beta_fixed_dim1 = Xf(1).cols(); // Number of fixed effect 1
  int beta_fixed_dim2 = Xf(2).cols(); // Number of fixed effect 2

  // Parameter
  PARAMETER_VECTOR(W); // W = c(U1,U2, beta1, beta2, beta_fixed0, beta_fixed1, beta_fixed2)

  vector<Type> U1(d1);
  vector<Type> U2(d2);

  vector<Type> beta1(betadim1);
  vector<Type> beta2(betadim2);
  vector<Type> beta_fixed0(beta_fixed_dim0);
  vector<Type> beta_fixed1(beta_fixed_dim1);
  vector<Type> beta_fixed2(beta_fixed_dim2);

  for (int i=0;i<d1;i++) U1(i) = W(i);
  for (int i=0;i<d2;i++) U2(i) = W(i + d1);

  for (int i=0;i<betadim1;i++) beta1(i) = W(i + d1 + d2);
  for (int i=0;i<betadim2;i++) beta2(i) = W(i + d1 + d2 + betadim1);

  for (int i=0;i<beta_fixed_dim0;i++) beta_fixed0(i) = W(i + d1 + d2 + betadim1 + betadim2);
  for (int i=0;i<beta_fixed_dim1;i++) beta_fixed1(i) = W(i + d1 + d2 + betadim1 + betadim2 + beta_fixed_dim0);
  for (int i=0;i<beta_fixed_dim2;i++) beta_fixed2(i) = W(i + d1 + d2 + betadim1 + betadim2 + beta_fixed_dim0 + beta_fixed_dim1);

  // The last element is for the variance of Gaussian family 
  PARAMETER_VECTOR(theta);

  // Transformations
  vector<Type> eta = X(0) * beta1 + X(1) * beta2 + B(0) * U1 + B(1) * U2 + Xf(0) * beta_fixed0 + Xf(1) * beta_fixed1 + Xf(2) * beta_fixed2;
  Type sigma1 = exp(-0.5*theta(0));
  REPORT(sigma1);
  Type sigma2 = exp(-0.5*theta(1));
  REPORT(sigma2);
  Type sigma3 = exp(-0.5*theta(2));
  REPORT(sigma3);



  // Log likelihood
  Type ll = 0;
  ll = sum(dnorm(y, eta, sigma3, TRUE));
  REPORT(ll);
  

  // Log prior on W
  Type lpW = 0;
  // Cross product (for each RE and its boundary, and for fixed effect)
  // RE1:
  vector<Type> P1U1 = P(0) * U1;
  Type U1P1U1 = (U1 * P1U1).sum();
  lpW += -0.5 * exp(theta(0)) * U1P1U1; // U part
  Type bb1 = (beta1 * beta1).sum();
  lpW += -0.5 * betaprec(0) * bb1; // Beta part
  // Log determinant
  Type logdet1 = d1 * theta(0) + logPdet(0);
  lpW += 0.5 * logdet1; // P part

  // RE2:
  vector<Type> P2U2 = P(1) * U2;
  Type U2P2U2 = (U2 * P2U2).sum();
  lpW += -0.5 * exp(theta(1)) * U2P2U2; // U part
  Type bb2 = (beta2 * beta2).sum();
  lpW += -0.5 * betaprec(1) * bb2; // Beta part
  // Log determinant
  Type logdet2 = d2 * theta(1) + logPdet(1);
  lpW += 0.5 * logdet2; // P part

  // Fixed Effects;
  Type bbf0 = (beta_fixed0 * beta_fixed0).sum();
  lpW += -0.5 * beta_fixed_prec(0) * bbf0; //
  Type bbf1 = (beta_fixed1 * beta_fixed1).sum();
  lpW += -0.5 * beta_fixed_prec(1) * bbf1; //
  Type bbf2 = (beta_fixed2 * beta_fixed2).sum();
  lpW += -0.5 * beta_fixed_prec(2) * bbf2; //


  // Log prior for theta
  Type lpT = 0;
  // Variance of RE1
  Type phi1 = -log(alpha(0)) / u(0);
  lpT += log(0.5 * phi1) - phi1*exp(-0.5*theta(0)) - 0.5*theta(0);

  // Variance of RE2
  Type phi2 = -log(alpha(1)) / u(1);
  lpT += log(0.5 * phi2) - phi2*exp(-0.5*theta(1)) - 0.5*theta(1);

  // Variance of the Gaussian family
  Type phi3 = -log(alpha(2)) / u(2);
  lpT += log(0.5 * phi3) - phi3*exp(-0.5*theta(2)) - 0.5*theta(2);
  
  // Final result!
  Type logpost = -1 * (ll + lpW + lpT);
  
  return logpost;
}
