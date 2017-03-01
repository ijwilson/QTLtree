/** @file */
// @file      Time-stamp: <2013-12-10 14:48:58 nijw>
#ifndef GSL_RAND_H__
#define GSL_RAND_H__

#include "utilityfunctionals.h"
#include "newio.h"

#include <cassert>
#include <stdexcept>  // for underflow_error

// include the header files -- try to use R also
extern "C" {
#ifdef USE_R
#ifndef USE_RMATH
#include <R.h>
#endif
#include <Rmath.h>
#else
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#endif
}

#undef sexp

/** The global random number generator             */
class rng;

/** wrapper for the gsl rng class             */
class rng {
public:
#ifndef USE_R
  // constructor and destructor
  rng(unsigned long int seed) {
    _r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(_r,seed);
  }
  ~rng() {
    gsl_rng_free(_r);
  }
  /** Reset the random number seed    */
  void reset(unsigned long int s) {
    gsl_rng_set(_r,s);
  }
  double next() {
    return gsl_rng_uniform(_r);
  }
  unsigned long int rint(unsigned long int mx) {
    return gsl_rng_uniform_int(_r,mx);
  }
  double normal() {
    return gsl_ran_ugaussian(_r);
  }
  double rgamma(double a, double scale) {
    return gsl_ran_gamma(_r,a,scale);
  }
  double rGamma(double a, double scale) {
    return gsl_ran_gamma(_r,a,scale);
  }
  double rlnorm(double zeta, double sigma) {
    return gsl_ran_lognormal(_r, zeta, sigma);
  }
  unsigned int rBinom(double  pp, unsigned int n) {
    return gsl_ran_binomial(_r,pp,n);
  }
  unsigned int rhypergeometric(unsigned int n1, unsigned int n2, unsigned int t) {
    return gsl_ran_hypergeometric(_r,n1,n2,t);
  }
  /*Function: unsigned int gsl_ran_hypergeometric (const gsl_rng * r, unsigned int n1, unsigned int n2, unsigned int t)

    This function returns a random integer from the hypergeometric distribution. The probability distribution for hypergeometric random variates is,

              p(k) =  C(n_1, k) C(n_2, t - k) / C(n_1 + n_2, t)

    where C(a,b) = a!/(b!(a-b)!) and t <= n_1 + n_2. The domain of k is max(0,t-n_2), ..., min(t,n_1).

    If a population contains n_1 elements of “type 1” and n_2 elements of “type 2” then the hypergeometric distribution gives the probability of 
    obtaining k elements of “type 1” in t samples from the population without replacement. */


  /** return a random poisson deviate with mean mu  */
  int rpoisson(double  mu) {
    return int(gsl_ran_poisson(_r,mu));
  }
  /** return a geometric rng
      p(k) =  p (1-p)^(k-1)
      for k >= 1.
  */
  int rgeometric(double p) {
    return gsl_ran_geometric(_r,p);
  }
  /** Generate exponential random variables         */
  double sexp(void) {
    return gsl_ran_exponential(_r,1.0);
  }
  double sexp(double mu) {
    return gsl_ran_exponential(_r,mu);
  }
  /** Generate a dirichlet random variable          */
  std::vector<double> rdirichlet(const std::vector<double> &a) {
    size_t d=a.size();
    std::vector<double> bb(d);
    //    for (size_t ii=0;ii<d;ii++) bb[ii] = sexp(a[ii]);
    //  double sm=std::accumulate(bb.begin(),bb.end(),0.0);
    //  for (size_t ii=0;ii<d;ii++) bb[ii] /= sm;
     gsl_ran_dirichlet(_r,d,&a[0],&bb[0]);
     return bb;
  }
#ifdef USETNT
#include "tnt/tnt.h"
 /** Generate a dirichlet random variable          */
  void rdirichlet(double *x,const TNT::Array1D<double> &a) {
    size_t d=a.dim();
    gsl_ran_dirichlet(_r,d,&a[0],x);
  }
#endif

#else
  /** this is set for the R random number generator.  Note that
   * there should only be a single member of this class, and that
   * this should be enforced (using a singleton template)        */
  rng() {
#ifdef LOUD
    std::cerr << "reading random seed" << std::endl;
#endif
    GetRNGstate();
  }
  ~rng() {
#ifdef LOUD
    std::cerr << "writing random seed" << std::endl;
#endif
    PutRNGstate();
  }
  double next() {
    return unif_rand();
  }
  unsigned long int rint(unsigned long int mx) {
    return (unsigned long)(next()*mx);
  }
  double normal() {
    return norm_rand();
  }
  double rGamma(double a, double b) {
    return rgamma(a,b);
  }
  int rBinom(double pp, int n) {
    return static_cast<int>(rbinom(static_cast<double>(n),pp));
  }
  int rpoisson(double mu) {
    return int(rpois(mu));
  }
  double sexp() {
    return exp_rand();
  }
  double sexp(double mu) {
   return mu*exp_rand();
  }
  int rgeometric(double p) {
    return int(rgeom(p));
  }
  double rlnorm(double zeta, double sigma) {
    return rlnorm(zeta, sigma);
  }
  std::vector<double> rdirichlet(const std::vector<double> &a) {
    size_t d=a.size();
    std::vector<double> bb(d);
    for (size_t ii=0;ii<d;ii++) bb[ii] = rgamma(a[ii],1.0);
    double sm=std::accumulate(bb.begin(),bb.end(),0.0);
    for (size_t ii=0;ii<d;ii++) bb[ii] /= sm;
    return bb;
  }
  unsigned int rhypergeometric(unsigned int n1, unsigned int n2, unsigned int t) {
    return static_cast<unsigned int>(rhyper(static_cast<double>(n1),static_cast<double>(n2),static_cast<double>(t)));
  }
#endif
  //
  double operator()(double mn, double mx) {
    return mn+(mx-mn)*next();
  }
  double operator()() {
    return next();
  }
  unsigned long operator()(unsigned long N) {
    return rint(N);
  }
  // std::vector<double> operator()(unsigned long n) {
//     std::vector<double> x(n);
//     for (unsigned long i=0;i<n;i++) {
//       x[i]=next();
//     }
//     return x;
//   };


  std::vector<double> normal(unsigned long n) {
    std::vector<double> a(n);
    for (unsigned int i=0;i<n;i++) a[i]= normal();
    return a;
  }

  /** Sample a sorted pair of integers in [from,to] with first>second */
  std::pair<int,int> sample2intsorted(int from, int to) {
    std::pair<int,int> p;
    double len=(double)(to-from+1);
    p.second = from + (int)(next()*len);
    p.first = from + (int)(next()*(len-1.));
    if (p.first>=p.second) {
      p.first++;
    } else {
      std::swap(p.first,p.second);
    }
    return p;
  }
  /** Sample a sorted pair of integers in [0,to) with first>second */
  std::pair<int,int> sample2intsorted(int to) {
    std::pair<int,int> p;
    double len=(double)(to);
    p.second = static_cast<int>(next()*len);
    p.first =  static_cast<int>(next()*(len-1.));
    if (p.first>=p.second) {
      p.first++;
    } else {
      std::swap(p.first,p.second);
    }
    return p;
  }
  /** Sample a pair of (different) integers from  [0,to) */
  std::pair<int,int> sample2int(int to) {
    std::pair<int,int> p;
    p.first = rint(to-1);
    p.second=rint(to-1);
    if (p.first>=p.second) {
      p.first++;
    }
    return p;
  }
#ifndef USE_R
  /** produce a sorted random selection k integers from [0:n) */
  std::vector<int> integer_choose(int k, int n) {
    std::vector<int> a(n);
    for (int i=0;i<n;i++) a[i]=i;
    std::vector<int> dest(k);
    gsl_ran_choose(_r,&dest[0],k,&a[0],n,sizeof(int));
    return dest;
  }
 /** produce a random permutation of k integers  from [0:n) (k=1..n) */
  std::vector<int> integer_permutations(int k, int n) {
    std::vector<int> a(n);
    for (int i=0;i<n;i++) a[i]=i;
    std::vector<int> dest(k);
    gsl_ran_choose(_r,&dest[0],k,&a[0],n,sizeof(int));
    gsl_ran_shuffle(_r,&dest[0], k,sizeof(int));
    return dest;
  }
#else
 /** produce a sorted random selection k integers from [0:n) */
  std::vector<int> integer_permutations(int k, int n) {
    std::vector<int> a(n);
    for (int i=0;i<n;i++) a[i]=i;
    for (int i=0;i<k;i++) {
      int swp=rint(n);
      int tmp=a[swp];
      a[swp]=a[i];
      a[i]=tmp;
    }
    a.resize(k);
    return a;
  }
 /** produce a random permutation of k integers  from [0:n) (k=1..n) */
  std::vector<int> integer_choose(int k, int n) {
    std::vector<int> a = integer_permutations(k,n);
    sort(a.begin(),a.end());
    return a;
  }
#endif


#ifndef USE_R
  gsl_rng *gslrng() {
    return _r;
  }
private:
  gsl_rng *_r;
  rng();  // this is only used with USE_R
#endif
private:
  // functions that are not defined for safety
  rng(const rng &a);
  rng& operator=(const rng &a);
};



template <typename T> size_t gen_from_cump(const T &cp, rng &r);
template <typename T> int gen_from_pr(const T &p, int lst, rng &r);
template <typename T>
int gen_from_pb(const T *const p, int lst, rng &r);

/** Generate from cumulative probabilities in p  */
template <typename T>
size_t gen_from_cump(const T &cp, rng &r)
// a routine to pick a number from 0 to n-1 based on cumulative
// "probabilities" in p
// this only takes deques and vectors...
{
  if (cp.back()<=0.0)
    throw std::runtime_error("sum of probabilities is "
			       "equal to zero in  gen_from_cump");
  double relprob = r();
  size_t where = static_cast<int>(relprob*(cp.size())); // approximate start
  relprob *= cp.back();
  for (;;) {
    if  (relprob <= cp[where]) {
      if (where==0) return 0;
      if (relprob > cp[where-1]) return where;
      else where--;
    } else {
      where++;
      assert(where < cp.size());
      if (where==cp.size()-1) return where;
    }
  }
}
/** Generate from a vector of probabilities.
 * psample gives the probability of the sample selected      */
template <typename T>
int gen_from_p(const T &p, rng &r)
{
  std::vector<double> cprob(p.size());
  std::transform(p.begin(),p.end(),cprob.begin(),cumsum<double>());
  //std::partial_sum(p.begin(),p.end(),cprob.begin());
  return gen_from_cump(cprob,r);
}
/** Generate n samples from a vector of probabilities.
 * psample gives the probability of the sample selected      */
template <typename T>
std::vector<int> gen_from_p(int n,const T &p, rng &r)
{
  std::vector<int> ret;
  ret.reserve(n);
  std::vector<double> cprob(p.size());
  std::transform(p.begin(),p.end(),cprob.begin(),cumsum<double>());
  //std::partial_sum(p.begin(),p.end(),cprob.begin());
  for (int i=0;i<n;i++) ret.push_back(gen_from_cump(cprob,r));
  return ret;
}
/** Generate n samples from a vector of probabilities.
 * psample gives the probability of the sample selected      */
template <typename T>
void gen_from_pn(int *res,int n,T *first, int len, rng &r)
{
  std::vector<double> cprob(len);
  std::transform(first,first+len,cprob.begin(),cumsum<double>());
  for (int i=0;i<n;i++) res[i] = gen_from_cump(cprob,r);
}
#ifdef USETNT
#include "tnt/tnt.h"
/** Generate n samples from a vector of probabilities.
 * psample gives the probability of the sample selected      */
template <typename T>
int gen_from_pTNT(const TNT::Array1D<T> &p, rng &r,double &psample)
{  
  // std::vector<int> ret;
  //ret.reserve(n);
  std::vector<double> cprob(p.dim());
  std::transform(&p[0],&p[0]+p.dim(),cprob.begin(),cumsum<double>());
  int i= gen_from_cump(cprob,r);
  if (i>0) psample=(cprob[i]-cprob[i-1])/cprob.back();
  else psample=cprob[i]/cprob.back();
  assert(psample>=0.0&&psample<=1.0);
  return i;
}
#endif
/** Generate n samples from a vector of probabilities.
 * psample gives the probability of the sample selected      */
template <typename T>
std::pair<int,int> gen_2_from_p(const T &p, rng &r)
{
  std::pair<int,int> ret;
  std::vector<double> cprob(p.size());
  std::transform(p.begin(),p.end(),cprob.begin(),cumsum<double>());
  //std::partial_sum(p.begin(),p.end(),cprob.begin());
  ret.first=gen_from_cump(cprob,r);
  ret.second=gen_from_cump(cprob,r);
  return ret;
}
/** Generate n samples from a vector of probabilities.
 * psample gives the probability of the sample selected      */
template <typename INTTYPE, typename T>
void gen_from_p(INTTYPE *a,int n,const T &p, rng &r)
{
  std::vector<double> cprob(p.size());
  std::transform(p.begin(),p.end(),cprob.begin(),cumsum<double>());
  //std::partial_sum(p.begin(),p.end(),cprob.begin());
  for (int i=0;i<n;i++) a[i]=gen_from_cump(cprob,r);
}
/** Generate from a vector of probabilities.
 * psample gives the probability of the sample selected      */
template <typename T>
int gen_from_p(const T &p, rng &r,double &psample)
{
  std::vector<double> cprob(p.size());
  std::transform(p.begin(),p.end(),cprob.begin(),cumsum<double>());
  int i= gen_from_cump(cprob,r);
  if (i>0) psample=(cprob[i]-cprob[i-1])/cprob.back();
  else psample=cprob[i]/cprob.back();
  assert(psample>=0.0&&psample<=1.0);
  return i;
}
/** Generate from a vector of probabilities.
 * psample gives the probability of the sample selected      */
template <typename T>
int gen_from_p(int n,const T &p, rng &r,double &psample)
{
  std::vector<double> cprob(p.size());
    std::transform(p.begin(),p.end(),cprob.begin(),cumsum<double>());
  int i= gen_from_cump(cprob,r);
  if (i>0) psample=(cprob[i]-cprob[i-1])/cprob.back();
  else psample=cprob[i]/cprob.back();
  assert(psample>=0.0&&psample<=1.0);
  return i;
}
/** this assumes that T has a random access iterator */
/** gives values in [frst,last]                      */
template <typename T>
int gen_from_p(const T &p,int frst, int lst, rng &r)
{
  std::vector<double> cprob(lst-frst+1);
    std::transform(p.begin()+frst,p.begin()+lst+1,cprob.begin(),cumsum<double>());
    //std::partial_sum(p.begin()+frst,p.begin()+lst+1,cprob.begin());
  return frst+gen_from_cump(cprob,r);
}
/** this assumes that T has a random access iterator */
/** gives values in [0,last)                      */
template <typename T>
int gen_from_pr(const T &p, int lst, rng &r)
{
  std::vector<double> cprob(lst);
    std::transform(p.begin(),p.begin()+lst,cprob.begin(),cumsum<double>());
    // std::partial_sum(p.begin(),p.begin()+lst,cprob.begin());
  return gen_from_cump(cprob,r);
}
/** this assumes that T has a random access iterator */
/** gives values in [0,last)                      */
template <typename T>
int gen_from_pr(const T &p, int lst, rng &r, double &psample)
{
  std::vector<double> cprob(lst);
  std::transform(p.begin(),p.begin()+lst,cprob.begin(),cumsum<double>());
   //std::partial_sum(p.begin(),p.begin()+lst,cprob.begin());
  int i = gen_from_cump(cprob,r);
  if (i>0) psample=(cprob[i]-cprob[i-1])/cprob.back();
  else psample=cprob[i]/cprob.back();
  assert(psample>=0.0&&psample<=1.0);
  return i;
}

/** generate proportional to numbers in [first,end)
 * don't need to be random access                           */
template <typename ITOR>
int gen_from_p(ITOR first, ITOR it_end, rng &r, double &psample)
{
  std::vector<double> cprob(std::distance(first,it_end));
  //std::partial_sum(first,it_end,cprob.begin());
  std::transform(first,it_end,cprob.begin(),cumsum<double>());
  int i= gen_from_cump(cprob,r);
  if (i>0) psample=(cprob[i]-cprob[i-1])/cprob.back();
  else psample=cprob[i]/cprob.back();
  assert(psample>=0.0&&psample<=1.0);
  return i;
}
/** this assumes that T has a random access iterator */
/** gives values in [0,last]                         */
/** I would like to get rid .....                    */
template <typename T>
int gen_from_p(const T *p, int lst, rng &r)
{
  std::vector<double> cprob(lst+1);
  std::transform(p,p+lst+1,cprob.begin(),cumsum<double>());
   //std::partial_sum(p,p+lst+1,cprob.begin());
  return gen_from_cump(cprob,r);
}
/** this assumes that T has a random access iterator */
/** gives values in [0,last)                         */
/** I would like to get rid .....                    */
template <typename T>
int gen_from_pb(const T *const p, int lst, rng &r)
{
  std::vector<double> cprob(lst);
   std::transform(p,p+lst,cprob.begin(),cumsum<double>());
   //std::partial_sum(p,p+lst,cprob.begin());
  return gen_from_cump(cprob,r);
}
/** generate proportional to numbers in [first,end)
 * don't need to be random access                           */
template <typename ITOR>
int gen_from_p(ITOR first, ITOR it_end, rng &r)
{
  std::vector<double> cprob(std::distance(first,it_end));
   std::transform(first,it_end,cprob.begin(),cumsum<double>());
   //std::partial_sum(first,it_end,cprob.begin());
  return gen_from_cump(cprob,r);
}
/** generate proportional to numbers in [first,end)
 * don't need to be random access                           */
// template <typename ITOR>
// int gen_from_p(ITOR first, ITOR it_end, rng &r, double &prob)
// {
//   std::vector<double> cprob(std::distance(first,it_end));
//   std::transform(first,it_end,cprob.begin(),cumsum<double>());
//   int i=gen_from_cump(cprob,r);
//   prob=cprob[i]/cprob.back();
//   return i;
// }
/** template class for permuting a set of T's              */
template <class T>
void permute(std::vector<T> &x,rng &r)
{
  for (size_t i=0;i<x.size();i++) {
    unsigned long which = r.rint(x.size());
    std::swap(x[which],x[i]);
  }
}

#ifndef USE_R
template <class T>
void permute_gsl(std::vector<T> &x,rng &r)
{
  gsl_ran_shuffle(r.gslrng(),&x[0],x.size(),sizeof(T));
}

#endif
/** template class for permuting a set of T's              */
template <class T>
void permute(T *x,int size,rng &r)
{
  for (int i=0;i<size;i++) {
    int which = r.rint(size);
    std::swap(x[which],x[i]);
  }
}
// /** template class for permuting a set of T's              */
// template <class T>
// void trypermute(T *x,int size,rng &r,bool goback=false)
// {
//   static *T lastorder;
//   static int lastsize=0;

//   if (goback) {
//     if (lastsize!=size)
//       throw
// 	std::domain_error("can only get back last order when one has been stored");
//     for (int i=0;i<size;i++) x[i]=lastorder[i];
//     return;
//   }
//   if (lastsize!=size) {
//     delete[] lastorder;
//     lastorder=new T[size];
//     lastsize=size;
//     for (int i=0;i<size;i++) lastorder[i]=x[i];
//   }
//   for (int i=0;i<size;i++) {
//     int which = r.rint(size);
//     T tmp=x[which];
//     x[which]=x[i];
//     x[i]=tmp;
//   }
//   return;
// }

/** template class for permuting a set of T's              */
template <class T,class Y>
void permute2(T *x1,Y *x2,int size,rng &r)
{
  for (int i=0;i<size;i++) {
    int which = r.rint(size);
    std::swap(x1[which],x1[i]);
    std::swap(x2[which],x2[i]);
  }
}

#endif
