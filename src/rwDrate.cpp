#include <TMB.hpp>

template <class Type>
Type dt1(Type x){
  return Type(1.0)/M_PI/(Type(1.0)+x*x);
}

template <class Type>
Type ilogit(Type x){
  return Type(1.0)/(Type(1.0)+exp(-x));
}

template <class Type>
Type nldens(Type x, Type mu, Type sd, Type p){
  Type z=(x-mu)/sd;
  return -log(Type(1.0)/sd*((Type(1.0)-p)*dnorm(z,Type(0.0),Type(1.0),false)+p*dt1(z)));
}


template<class Type>
bool isNA(Type x){
	return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(slope);
  DATA_VECTOR(idx);
  DATA_VECTOR(j);
  DATA_VECTOR(w);
//  DATA_SCALAR(u0);

  PARAMETER(logSdProc);
  PARAMETER(logSdObs);
  PARAMETER(u0);
  PARAMETER_VECTOR(u);
  PARAMETER(logitp);

  int procTime=u.size();
  int obsTime=slope.size();

  //*********************************************************//
  // Parameter transformations:

  Type p = ilogit(logitp); // logit transform mixture proportions
  Type sdProc = exp(logSdProc); // Exponential, to ensure it's positive
  Type sdObs = exp(logSdObs);  // Exponential, to ensure it's positive


  //*********************************************************//
  // Process:

  Type nll=-dnorm(u(0),u0,sdProc,true); // Initialize starting point
  //  Type ans=0;

  for(int i=1;i<procTime;i++){
    nll += -dnorm(u(i),u(i-1),sdProc,true);
  }


  //*********************************************************//
  // Observations:

  for(int i=0;i<obsTime;i++){
    if(!isNA(slope(i))){
	    Type idx1 = CppAD::Integer(idx(i));
	    Type idx0 = idx1-1;
	    Type u1 = u(CppAD::Integer(idx1));
	    Type u00 = u(CppAD::Integer(idx0));
	    Type jj = j(i);
	    Type umu = (u00*jj)+(u1*(1-(jj)));
	    Type ww = 1-w(i);
	    Type pp = p*ww;
	// Mixture of normal and t:
	//    nll += nldens(slope(i), umu, sdObs, p);
	//    nll += nldens(slope(i), umu, sdObs, p*ww);
	// Mixture of normals:
	    nll -= log((1-pp)*dnorm(slope(i), umu, sdObs, false) + pp*dnorm(slope(i), Type(0), Type(0.5), false));
	  }
  }
  return nll;
}

// See examples for mixture model and continuous time at:
// https://github.com/cavios/tshydro/blob/master/tsHydro/src/track.cpp
