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

  // The last element of u and the last element of alpha are for the family
  // (For example: the Gaussian noise variance when family = "Gaussian")
  DATA_STRUCT(u, list_Scalar_from_R);       // Pc prior, u param
  DATA_STRUCT(alpha, list_Scalar_from_R);   // Pc prior, alpha param
  DATA_STRUCT(betaprec, list_Scalar_from_R);// For boundary, beta ~iid N(0,1/betaprec)

  // For fixed effects
  DATA_STRUCT(beta_fixed_prec, list_Scalar_from_R);
  DATA_STRUCT(Xf, list_SparseMatrix_from_R);

  DATA_VECTOR(y); //response variable
  
  DATA_SCALAR(family_type); // Family types: Gaussian - 0, Poisson - 1, Binomial - 2

  vector<int> betadim(X.size());
  int sum_betadim = 0;
  for (int i = 0; i < X.size(); i++){
    betadim(i) = X(i).cols(); // Number of boundary conditions in RE #i
    sum_betadim += betadim(i);
  }

  vector<int> d(P.size());
  int sum_d = 0;
  for (int i = 0; i < P.size(); i++){
    d(i) = P(i).cols(); // Number of spline in RE #i
    sum_d += d(i);
  }

  vector<int> beta_fixed_dim(Xf.size());
  for (int i = 0; i < Xf.size(); i++){
    // beta_fixed_dim(0) is for the intercept 
    // beta_fixed_dim(1) is for the fixed effect 1, etc
    beta_fixed_dim(i) = Xf(i).cols(); 
  }

  // Parameter
  PARAMETER_VECTOR(W); // W = c(U1,U2, beta1, beta2, beta_fixed0, beta_fixed1, beta_fixed2)

  vector<vector<Type>> U(d.size());
  for (int i = 0; i < d.size(); i++){
    vector<Type> cur_U(d(i));
    U(i) = cur_U;
  }

  vector<vector<Type>> beta(betadim.size());
  for (int i = 0; i < betadim.size(); i++){
    vector<Type> cur_beta(betadim(i));
    beta(i) = cur_beta;
  }

  vector<vector<Type>> beta_fixed(beta_fixed_dim.size());
  for (int i = 0; i < beta_fixed_dim.size(); i++){
    vector<Type> cur_beta_fixed(beta_fixed_dim(i));
    beta_fixed(i) = cur_beta_fixed;
  }

  for (int j = 0; j < d.size(); j++){
    if (j == 0){
      for (int i=0;i<d(j);i++) U(j)(i) = W(i);
    }
    else{
      for (int i=0;i<d(j);i++) U(j)(i) = W(i + d(j-1));
    }
  }

  for (int j=0;j<betadim.size();j++){
    if (j == 0){
      for (int i=0;i<betadim(j);i++) beta(j)(i) = W(i + sum_d);
    }
    else{
      int prev_betadim_sum = 0;
      for (int i=0;i<j;i++) prev_betadim_sum += betadim(i);
      for (int i=0;i<betadim(j);i++) beta(j)(i) = W(i + sum_d + prev_betadim_sum);
    }
  }

  for (int j=0;j<beta_fixed_dim.size();j++){
    if (j == 0){
      for (int i=0;i<beta_fixed_dim(j);i++) beta_fixed(j)(i) = W(i + sum_d + sum_betadim);
    }
    else{
      int prev_beta_fixed_dim_sum = 0;
      for (int i=0;i<j;i++) prev_beta_fixed_dim_sum += beta_fixed_dim(i);
      for (int i=0;i<beta_fixed_dim(j);i++) beta_fixed(j)(i) = W(i + sum_d + sum_betadim + prev_beta_fixed_dim_sum);
    }
  }

  // The last element is for the variance of Gaussian family 
  PARAMETER_VECTOR(theta);

  // Transformations
  vector<Type> eta = X(0) * beta(0); // adding intercept

  for (int i = 1; i < X.size(); i++){
    eta += X(i) * beta(i);  // adding each boundary condition for RE
  }

  for (int i = 0; i < B.size(); i++){
    eta += B(i) * U(i); // adding each local spline effect for RE
  }

  for (int i = 0; i < Xf.size(); i++){
    eta += Xf(i) * beta_fixed(i); // adding each fixed effect
  }

  vector<Type> sigma(theta.size());
  for (int i = 0; i < theta.size(); i++){
    sigma(i) = exp(-0.5*theta(i)); // transforming each variance parameter (from log precision to standard deviation)
    REPORT(sigma(i));
  }


  // Log likelihood
  Type ll = 0;
  // Family types: Gaussian - 0, Poisson - 1, Binomial - 2
  if (family_type == 0){
    ll = sum(dnorm(y, eta, sigma(theta.size() - 1), TRUE));
  } 
  else if (family_type == 1){
    ll = sum(dpois(y, exp(eta), TRUE));
  } 
  // else if (family_type == 2){
  //   Initialize a vector of 1s with the size of length(y)
  //   ll = sum(dbinom_robust(y, size, eta, TRUE));
  // } 
  REPORT(ll);
  

  // Log prior on W
  Type lpW = 0;
  // Cross product (for each RE and its boundary, and for fixed effect)
  // For Random Effects:
  for (int i = 0; i < X.size(); i++){
    lpW += -0.5 * exp(theta(i)) * ((U(i) * (P(i) * U(i))).sum()); // U part (spline effect)
    Type bb = (beta(i) * beta(i)).sum();
    lpW += -0.5 * betaprec(i) * bb; // Beta part (boundary condition)
    // Log determinant
    Type logdet = d(i) * theta(i) + logPdet(i);
    lpW += 0.5 * logdet; // P part
  }

  // For Fixed Effects;
  for (int i = 0; i < beta_fixed.size(); i++){
    Type bbf = (beta_fixed(i) * beta_fixed(i)).sum(); // Fixed effect
    lpW += -0.5 * beta_fixed_prec(i) * bbf; //
  }

  // Log prior for theta
  Type lpT = 0;
  // Variance of each random effect (and the family)
  for (int i = 0; i < alpha.size(); i++){
    Type phi = -log(alpha(i)) / u(i);
    lpT += log(0.5 * phi) - phi*exp(-0.5*theta(i)) - 0.5*theta(i);
  }
  
  // Final result!
  Type logpost = -1 * (ll + lpW + lpT);
  
  return logpost;
}
