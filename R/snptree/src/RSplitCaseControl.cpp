#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <stdexcept>
#include "splitter.h"
#include "binode.h"
#include "gsl_rand.h"
#include "tnt/tnt.h"


extern "C" {
 
  /** Get the split in a form that is suitable for using within ape                    
   * The nodepos and edgepos give the left hand and right hand 
   * positions of the minimum and maximum range of the haplotypes at each 
   * node and tip (note both are -1 for tips with only a single haplotype 
   * i.e. those that have been successfully separated)                        
   * The positions are calculated for all data, not just those positions that have
   * been separated out, but they always include the starting position (which should guarantee
   * that all positions that the haplotypes have been split on are included).   
   * 
   * data:          carries the set of genetic data
   * samplesize:    the number of samples
   * nSNP:          the number of SNPs
   * positions:     The positions to split at 
   * npos:          The number of positions to split at
   * cases:         A vector of cases
   * ncases:        The number of cases
   * edge:          a vector of edges 
   * ccleaf:        Number of cases and controls at the leaves
   * ccnode:        Number of cases and control at the nodes
   * len:           The number of leaves
   * labels:        The labels at the leaves (in order)
  */
  void GetSplitCC( int *data, int *samplesize, int *nSNP, int *positions, int *npos, int *cases
		   , int *ncases, int *edge, int *ccleaf, int *ccnode,int *len, int *comblabels
		 , int *nodepos) 
  {  
    TNT::Array2D<int> d(*samplesize,*nSNP,data);
    splitter<int> s(d,d.dim1(),d.dim2());
    for (int i=0;i<*npos;i++) s.split(positions[i]);
    *len=s.nleaves();
    
    std::vector<std::pair<int,int> > edges;
    std::vector<std::vector<int> > labels;
    
    s.apesplit(edges,labels);

    int nedges=static_cast<int>(edges.size());
    if (nedges!=2*(*len-1))
      throw std::range_error("problem in GetSplit\n");

    for (int i=0;i<nedges;i++) {
      edge[i]=edges[i].first;
      edge[nedges+i]=edges[i].second;
    }
    // extract cases and controls for the leaves
    int count=0;
    for (int i=0;i<*len;i++) {
      ccleaf[i]=count_intersection(labels[i],cases,*ncases);
      ccleaf[*len+i] = labels[i].size()-ccleaf[i];
      for (size_t j=0;j<labels[i].size();j++) comblabels[count++] = labels[i][j]+1;
    }   
    s.getCaseControlNodes(ccnode,cases,*ncases,nedges/2);
    // now try to get the lengths.  Note that the centre of this split is 
    // positions[0].
    std::vector<int> nodePos;
    s.getNodesPositions(nodePos);
	  for (size_t ii=0;ii<nodePos.size();ii++) nodepos[ii]=nodePos[ii];
  }
  
  /** C code to perform the bulk of the test work                             */
  void splitTestCC(int *data, 
                 int *samplesize, 
                 int *nSNP, 
                 int *positions, 
                 int *npos, 
                 int *cases, int *ncases, 
                 int *reps, int *maxk, 
                 double *teststat, 
                 double *randteststats,
                 int *nterm, 
                 char **statPick) {
    
    TNT::Array2D<int> d(*samplesize, *nSNP, data);
    splitter<int> s(d,d.dim1(),d.dim2());
    for (int i=0;i<*npos;i++) s.split(positions[i]);
    
    std::vector<int> caseVec(cases,cases+*ncases);
    
    std::vector<double> stat=s.getStat(caseVec, *maxk, *statPick);   // use default test statistic
    for (size_t jj=0;jj< *maxk;jj++) teststat[jj] = stat[jj];

    rng r;
    std::vector<int> cc(*samplesize);
    for (int i=0;i<*samplesize;i++) cc[i]=i;
    
    for (int i=0;i<*reps;i++) {
      permute(cc,r);
      caseVec.assign(cc.begin(),cc.begin()+ *ncases);
      std::sort(caseVec.begin(),caseVec.begin()+ *ncases);
      std::vector<double> rstat =  s.getStat(caseVec,*maxk,*statPick);
      for (size_t jj=0;jj< *maxk;jj++) {
        randteststats[i*(*maxk) + jj] = rstat[jj];
      }
    }
    *nterm=s.nleaves();
  }  
}
