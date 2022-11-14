#include <TMB.hpp>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template <class Type>
Type nllBioProcess(array<Type> P, vector<Type> meanVec, vector<Type> logPhi, Type logSdP){
    int nrow=P.dim[0];
    int ncol=P.dim[1];
    int n=nrow*ncol;
    vector<int> r(n); r.setZero();   
    vector<int> c(n); c.setZero();
    int idx=0;
    for(int j=0; j<ncol; ++j){
      for(int i=0; i<nrow; ++i){
        r(idx)=i;
	c(idx)=j;
	++idx;
      }
    }
    matrix<Type> Wc(n,n); Wc.setZero();
    matrix<Type> Wd(n,n); Wd.setZero();
    for(int i=0; i<n; ++i){
      for(int j=0; j<n; ++j){
	if((c(i)==c(j))&&(abs(r(i)-r(j))==1)){
	  Wc(i,j)=1;
	  Wc(i,i)-=1;
        }   
 	if( (((r(i)-r(j))==1)&&((c(i)-c(j))==1))||(((r(i)-r(j)==(-1)))&&((c(i)-c(j))==(-1))) ){
       	  Wd(i,j)=1;
	  Wd(i,i)-=1;
	}
      }
    }
    
    array<Type> mP(nrow,ncol);
    for(int i=0; i<nrow; ++i){
      for(int j=0; j<ncol; ++j){
	mP(i,j)=meanVec(j);
      }
    }
    
    vector<Type> phi=exp(logPhi);
    matrix<Type> I(Wc.rows(),Wc.cols());
    I.setIdentity();
    matrix<Type> Q=I-phi(0)*Wc-phi(1)*Wd;
    
    using namespace density;
    return SCALE(GMRF(asSparseMatrix(Q)),exp(logSdP))((P-mP).vec());
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_ARRAY(Y)
    
  PARAMETER_VECTOR(logPhi)
  PARAMETER_VECTOR(mu)
  PARAMETER(logSdProc)   
  PARAMETER(logSdObs)   
  PARAMETER_ARRAY(P)

  Type jnll = nllBioProcess(P,mu,logPhi,logSdProc);  

  for(int i=0; i<Y.dim(0); ++i){
    for(int j=0; j<Y.dim(1); ++j){
      if(!isNA(Y(i,j))){
        jnll += -dnorm(Y(i,j),P(i,j),exp(logSdObs),true);
      }
    }
  }
  return(jnll);
}
