#include <TMB.hpp>

#ifndef REPORT_F
#define REPORT_F(name,F)					\
  if(isDouble<Type>::value && F->current_parallel_region<0) {		\
    Rf_defineVar(Rf_install(#name),					\
	      PROTECT(asSEXP(name)),F->report);			\
    UNPROTECT(1);						\
  }
#endif


template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template <class Type>
Type nllBioProcess(array<Type> P, vector<Type> meanVec, vector<Type> logPhi, Type logSdP, objective_function <Type > *of){
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
    matrix<Type> Wp(n,n); Wp.setZero();    
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
	if(logPhi.size()==3){
          if ((c(i)==(ncol-1)) && (c(j)==(ncol-1)) && (abs(r(i)-r(j))==1) ){
            Wp(i,j)=1;
	    Wp(i,i)-=1;
	  }
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
    if(logPhi.size()==3){
      Q-=phi(2)*Wp;
    }
    REPORT_F(Wp,of)  
    using namespace density;
    return SCALE(GMRF(asSparseMatrix(Q)),exp(logSdP))((P-mP).vec());
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_ARRAY(Y);
  DATA_ARRAY(w);
    
  PARAMETER_VECTOR(logPhi);
  PARAMETER_VECTOR(mu);
  PARAMETER(logSdProc);  
  PARAMETER(logSdObs);   
  PARAMETER_ARRAY(P);

  Type jnll = nllBioProcess(P,mu,logPhi,logSdProc, this);  

  for(int i=0; i<Y.dim(0); ++i){
    for(int j=0; j<Y.dim(1); ++j){
      if(!isNA(Y(i,j))){
	Type yobs = log(Y(i,j)) - log(1.0 - Y(i,j));
        jnll += -dnorm(yobs,P(i,j),w(i,j) * (1.0 + exp(logSdObs)),true) + Y(i,j) * (1 - Y(i,j));
      }
    }
  }

  return(jnll);
}
