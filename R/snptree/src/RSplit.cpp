#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <stdexcept>
#include "splitter.h"
#include "binode.h"
#include "gsl_rand.h"
#include "tnt/tnt.h"


extern "C" {
 
 /** Get a split that splits everything (if no information then the 
  * splits are at random                                               */
 void simplesplit( int *data, int *samplesize, int *nSNP, 
                   int *edge,int *comblabels) 
 {  
   TNT::Array2D<int> d(*samplesize,*nSNP,data);
   splitter<int> s(d,d.dim1(),d.dim2());
   for (int i=0;i<*nSNP;i++) s.split(i);
   
   while(s.fullsplit());
   
   std::vector<std::pair<int,int> > edges;
   std::vector<std::vector<int> > labels;
   
   s.apesplit(edges,labels);
   
   
   int nedges=static_cast<int>(edges.size());
   if (nedges!=2*(*samplesize-1))
     throw std::runtime_error("problem in C++ code simplesplit\n");
   
   for (size_t j=0;j<labels.size();j++) comblabels[j] = labels[j][0];
   
   for (int i=0;i<nedges;i++) {
     edge[i]=edges[i].first;
     edge[nedges+i]=edges[i].second;
   }
 }
 
}
