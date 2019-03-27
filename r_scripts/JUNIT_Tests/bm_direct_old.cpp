#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
RcppExport SEXP bm_direct (SEXP dat, SEXP pars) 
{
  /* 
  * all objects ordered from 1:(Nnode(phy)+Ntip(phy)) unless noted otherwise
  * dat: a list of elements 
  * len: edge lengths (vector)
  * root: root ID
  * y: tip data 
  * order: order of internal IDs for pruning algorithm
  * pars: rates associated with each branch 
  */
  
  try {
    /* call in parameters associated with 'dat' object */
    Rcpp::List cache(dat);
    
    int root = Rcpp::as<int>(cache["root"]);
    int n = Rcpp::as<int>(cache["n"]);
    
    double drift = Rcpp::as<double>(cache["drift"]);
    
    std::vector<double> len = Rcpp::as<std::vector<double> >(cache["len"]);
    std::vector<double> y = Rcpp::as<std::vector<double> >(cache["y"]);
    std::vector<double> var = Rcpp::as<std::vector<double> >(cache["var"]);
    std::vector<double> known = Rcpp::as<std::vector<double> >(cache["given"]); //
    std::vector<int> intorder = Rcpp::as<std::vector<int> >(cache["intorder"]);
    std::vector<int> tiporder = Rcpp::as<std::vector<int> >(cache["tiporder"]);
    std::vector<int> descR = Rcpp::as<std::vector<int> >(cache["descRight"]);
    std::vector<int> descL = Rcpp::as<std::vector<int> >(cache["descLeft"]);
    
    std::vector<double> rates = Rcpp::as<std::vector<double> > (pars);
    
    std::vector<double> lq;
    lq.assign(n,0.0);
    
    double yi, ri, li, m1, m2, v1, v2, v12, m12, n12, nm, nv, m, mm, v, k;
    
    double const PIx = 4.0*atan(1.0);
    
    std::vector<double> branchinitM;
    std::vector<double> branchinitV;
    std::vector<double> branchbaseM;
    std::vector<double> branchbaseV;
    
    branchinitM.assign(n,0.0);
    branchinitV.assign(n,0.0);
    branchbaseM.assign(n,0.0);
    branchbaseV.assign(n,0.0);
    
    int i, z, cur, d1, d2;
    
    /* mean and variance for leaves */
    z=tiporder.size();
    for(i=0; i<z; i++){ 
      cur=tiporder[i]-1;
      yi=y[cur];
      li=len[cur];
      ri=rates[cur];
      
      branchinitM[cur] = yi;
      branchbaseM[cur] = yi + drift*li;
      branchbaseV[cur] = var[cur] + li*ri;
    }
    
    /* mean, variance, and density for edges */
    z=intorder.size();
    for(i=0; i<z; i++){ 
      cur=intorder[i]-1;
      d1=descR[cur]-1;
      d2=descL[cur]-1;
      m1=branchbaseM[d1];
      m2=branchbaseM[d2];
      v1=branchbaseV[d1];
      v2=branchbaseV[d2];
      
      v12=v1+v2;
      
      m = (((m1*v2) + (m2*v1))/v12);                  // phylogenetic mean expectation
      branchinitM[cur] = m;
      
      v = ((v1*v2)/v12);
      branchinitV[cur] = v;
      
      m12=pow((m1-m2),2);
      lq[cur] = ((-m12/(2*v12)) - (log(2*PIx*v12)/2));
      
      k=known[cur];                                   
      
      if( k == (signed)1 )                            // resolve whether node state is given (k==1)
      {
        nm=y[cur];
        nv=var[cur];
        mm=m;
        
        v12=v+nv;
        m = ((mm*nv) + (nm*v))/v12;
        branchinitM[cur] = m;
        
        v = (v*nv)/v12;
        branchinitV[cur] = v;
        
        m12=pow((mm-nm),2);
        lq[cur]+=((-m12/(2*v12)) - (log(2*PIx*v12)/2));
      }
      
      li=len[cur];            
      branchbaseM[cur] = m + drift*li;
      branchbaseV[cur] = v + rates[cur]*li;
    }
    
    /* compute root */ 
    cur=root-1;
    d1=descR[cur]-1;
    d2=descL[cur]-1;
    m1=branchbaseM[d1];
    m2=branchbaseM[d2];
    v1=branchbaseV[d1];
    v2=branchbaseV[d2];
    v12=v1+v2;
    
    m=(((m1*v2) + (m2*v1))/v12); //
    branchinitM[cur] = m; //
    v=((v1*v2)/v12); //
    branchinitV[cur] = v;
    m12=pow((m1-m2),2);
    lq[cur] = ((-m12/(2*v12)) - (log(2*PIx*v12)/2));
    
    /* compute root lnL (either ML or given) */
    k=known[cur];
    if(k == (signed)1 )                                
    { // given state
      nm=y[cur];                                     
      nv=var[cur] + v;
      n12=pow((nm-m),2);                             
      lq[cur] += ((-n12/(2*nv)) - (log(2*PIx*nv)/2));
    }
    else
    { // ML state
      lq[cur] += (- (log(2*PIx*v)/2));
      
    }
    
    
    /* PREPARE OUTPUT FOR R */
    return Rcpp::List::create(
      Rcpp::Named("initM",branchinitM),
      Rcpp::Named("initV",branchinitV),
      Rcpp::Named("baseM",branchbaseM),
      Rcpp::Named("baseV",branchbaseV),
      Rcpp::Named("lq",lq)
    );
    
    
  } catch( std::exception &ex ) {		
    forward_exception_to_r( ex );
  } catch(...) { 
    ::Rf_error( "C++ exception: unknown reason" ); 
  }
  return R_NilValue; 
}