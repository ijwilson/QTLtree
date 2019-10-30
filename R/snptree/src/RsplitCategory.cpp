#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <stdexcept>

#include "splitter.h"
#include "gsl_rand.h"
#include "tnt/tnt.h"


extern "C" {
  

  void GetCategorySplit( int *data
                           ,int *samplesize
                           ,int *nSNP
                           ,int *positions
                           ,int *npos
                           ,int  *category
                           ,int *edge
                           ,int *len
                           ,int *leafcount
                           ,int *comblabels
                           ,int *nodepos
                           ,int *tippos
                           ,int *nodecat
                           ,int *leafcat
                           ,int *ncategories) 
  {  
    // Rprintf("Start of function\n");
    TNT::Array2D<int> d(*samplesize,*nSNP,data);
    splitter<int> s(d,d.dim1(),d.dim2());
    for (int i=0;i<*npos;i++) s.split(positions[i]);
    *len=s.nleaves();
    
    // get the edges and labels in a form for ape
    std::vector<std::pair<int,int> > edges;
    std::vector<std::vector<int> > labels;
    s.apesplit(edges,labels);  
    
    int nedges=static_cast<int>(edges.size());
    if (nedges!=2*(*len-1))
      throw std::range_error("problem in GetCategorySplit\n");
    
    for (int i=0;i<nedges;i++) {
      edge[i]=edges[i].first;
      edge[nedges+i]=edges[i].second;
    }
    //   Rprintf("geting lengths\n");
    s.getLengths(positions[0],nodepos,tippos);
    
    for (size_t ii=0;ii<labels.size();ii++) {
      leafcount[ii] = static_cast<int>(labels[ii].size());
      int *counts = leafcat+*ncategories*ii;
      //      Rprintf("=========== leaf %d ================\n",ii);
      for (int jj=0;jj<leafcount[ii];jj++) 
        counts[category[labels[ii][jj]-1]-1]++;
    }
    //    Rprintf("getting node categories\n");
    
    std::vector<std::vector<int> > nodelabels;
    s.getNodesLabels(nodelabels);
    for (size_t ii=0;ii<nodelabels.size();ii++) {
      int *counts = nodecat+*ncategories*ii;
      for (int jj=0;jj<nodelabels[ii].size();jj++) {
        counts[category[nodelabels[ii][jj]-1]-1]++;
      }
    }
  }

 
}
