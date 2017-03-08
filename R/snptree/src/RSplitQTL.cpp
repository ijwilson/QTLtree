#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <stdexcept>

#include "splitter.h"
#include "gsl_rand.h"
#include "tnt/tnt.h"


extern "C" {
  
  /**
   *
   *  Split something based on QTL values at the nodes of a tree
   *
   */
  
  void splittestQTL( int *data, 
                     int *samplesize, 
                     int *nSNP, 
                     int *positions, 
                     int *npos,
                     double *qtl,
                     int *reps, 
                     int *maxk,
                     double *teststat,
                     double *randteststats,
                     int *nterm, 
                     char **statPick)
  {

    TNT::Array2D<int> d(*samplesize, *nSNP, data);
    splitter<int> s(d,d.dim1(), d.dim2());
    for (int i=0;i<*npos;i++) s.split(positions[i]);
    
    std::vector<double> myqtl(qtl, qtl + *samplesize);
    std::vector<double> stat=s.qtlStat(myqtl, *maxk, *statPick);
    for (size_t jj=0;jj< *maxk;jj++) teststat[jj] = stat[jj];
    
    rng r;
    for (int i=0; i< *reps; i++) {
        permute(myqtl, r);
        std::vector<double> rstat =  s.qtlStat(myqtl, *maxk, *statPick);
        for (size_t jj=0;jj< *maxk;jj++) {
          randteststats[i*(*maxk) + jj] = rstat[jj];
        }
      }

      *nterm=s.nleaves();
  }

  /** Get the split in a form that is suitable for using within ape                    
   * The nodepos and edgepos give the left hand and right hand 
   * positions of the minimum and maximum range of the haplotypes at each 
   * node and tip (note both are -1 for tips with only a single haplotype 
   * i.e. those that have been successfully separated)                        
   * The positions are calculated for all data, not just those positions that have
   * been separated out, but they always include the starting position (which should guarantee
   * that all positions that the haplotypes have been split on are included).   
  */
  void GetQTLSplit( int *data, 
                    int *samplesize, 
                    int *nSNP, 
                    int *positions, 
                    int *npos,
                    double *qtl,
                    int *edge
                    ,int *len
                    ,int *leafcount
                    ,int *comblabels
                    ,int *nodepos
                    ,int *tippos
                    ,double *nodestat
                    ,double *leafstat
                    ) 
  {  
    TNT::Array2D<int> d(*samplesize,*nSNP,data);
    splitter<int> s(d,d.dim1(),d.dim2());
    for (int i=0;i<*npos;i++) s.split(positions[i]);
    *len=s.nleaves();
   
    // get the edges and labels in a form for ape
    std::vector<std::pair<int,int> > edges;
    std::vector<std::vector<int> > labels;
    s.apesplit(edges,labels);     
    int nedges=static_cast<int>(edges.size());
    if (nedges!=2*(*len-1)) {
      throw std::range_error("problem in GetQTLSplit\n");
    }
   
    for (int i=0;i<nedges;i++) {
      edge[i]=edges[i].first;
      edge[nedges+i]=edges[i].second;
    }
    
    s.getLengths(positions[0],nodepos,tippos);
    // get statistics for the tree
    double xbar,s2;
    std::vector<double> qv(qtl,qtl+*samplesize);
    avevar(qv,xbar,s2);
    int count=0;
    for (size_t ii=0;ii<labels.size();ii++) {
      double sum=0.0;
      leafcount[ii] = static_cast<int>(labels[ii].size());
      for (int jj=0;jj<leafcount[ii];jj++) {
        sum += qtl[labels[ii][jj]];
        comblabels[count++] = labels[ii][jj]+1;
      }
      leafstat[ii] =  (sum/leafcount[ii]-xbar)/sqrt(s2/leafcount[ii]);
    }
    std::vector<std::vector<int> > nodelabels;
    s.getNodesLabels(nodelabels);
    for (size_t ii=0;ii<nodelabels.size();ii++) {
      double sum=0.0;
      size_t n= nodelabels[ii].size();
      for (size_t jj=0;jj<n;jj++) sum += qtl[nodelabels[ii][jj]];
      nodestat[ii] = (sum/n - xbar)/sqrt(s2/n);
      leafcount[labels.size()+ii] = nodelabels[ii].size();
    }
  }
}
