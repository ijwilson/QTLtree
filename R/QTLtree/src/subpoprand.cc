#include "R.h"
#include "gsl_rand.H"
#include "tnt/tnt.h"
#include <vector>
#include <algorithm>

extern "C" {
 /** get the test statistic - assumes that the vector of cases is sorted  */
  double teststat(const std::vector<int> &cases, int *g,int n) 
  {    
    TNT::Array2D<int> m(2,2,0);
    int CaseIndex=0;
      
    for (size_t i=0;i<n;i++) {
      if (cases[CaseIndex]==i) {
        if (g[2*i]!=NA_INTEGER) {
          m[g[2*i]][0]+=1;
        }
        if (g[2*i+1]!=NA_INTEGER) {
          m[g[2*i+1]][0]+=1;
        }
        CaseIndex++;
      } else {
        if (g[2*i]!=NA_INTEGER) {
          m[g[2*i]][1]+=1;
        }
        if (g[2*i+1]!=NA_INTEGER) {
          m[g[2*i+1]][1]+=1;
        }
      }
    }
  
    double colsum0=m[0][0]+m[1][0];
    double colsum1=m[0][1]+m[1][1];

    if (colsum0==0||colsum1==0) return 0.0;
    double tot=colsum0+colsum1;
  
    double rowsum0=m[0][0]+m[0][1];
    double rowsum1=m[1][0]+m[1][1];
      
    double e00 = rowsum0*colsum0/tot;
    double e10 = rowsum1*colsum0/tot;
    double e01 = rowsum0*colsum1/tot;
    double e11 = rowsum1*colsum1/tot;
    return (m[0][0]-e00)*(m[0][0]-e00)/e00+(m[1][0]-e10)*(m[1][0]-e10)/e10
      +(m[0][1]-e01)*(m[0][1]-e01)/e01+(m[1][1]-e11)*(m[1][1]-e11)/e11;
  } 
  //
  //
  //
  void subpoprand(int *genotype, int *n, int *nSNPs, int *cc, int *region
                  , int *reps,double *obsts,double *p) 
  {    
    rng r;
    int nregions=*std::max_element(region,region+*n)+1;
    std::vector<std::vector<int> > regionSamples(nregions);
    std::vector<int> regionCases(nregions,0);
    std::vector<int> cases;
    // first count the regions and get the case index numbers
    for (int i=0;i<*n;i++) {
      regionSamples[region[i]].push_back(i);
      if (cc[i]==1) {
        regionCases[region[i]]+=1;
        cases.push_back(i);
      }
    }
        
    std::vector<std::vector<int> > replicateSamples(*reps);
    for (int i=0;i<*reps;i++) {
      for (int jj=0;jj<nregions;jj++) {
        permute(regionSamples[jj],r);
        std::copy(regionSamples[jj].begin(),regionSamples[jj].begin()+regionCases[jj]
                  ,std::back_inserter(replicateSamples[i]));
        
      }
      if (cases.size()!=replicateSamples[i].size()) {
        return;
      }
      std::sort(replicateSamples[i].begin(),replicateSamples[i].end());
    }
    for (int i=0;i<*nSNPs;i++) {
      int *g=genotype+(*n)*i*2;
      obsts[i]=teststat(cases,g,*n);
      if (fabs(obsts[i])>0.01) {
        p[i]=0.0;
        for (int jj=0;jj<*reps;jj++) {
          double repts= teststat(replicateSamples[jj],g,*n);
          if (obsts[i]<=repts) p[i]+=1.0; 
        }
      } else p[i]=static_cast<double>(*reps);
    } 
  }
}
                             
