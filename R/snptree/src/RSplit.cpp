#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <stdexcept>
#include "splitter.h"
#include "binode.h"
#include "gsl_rand.h"
#include "tnt/tnt.h"


extern "C" {
  
  void split_simple( int *data, 
                     int *samplesize, 
                     int *nSNP, 
                     int *positions, 
                     int *npos, 
                     int *edge,
                     int *len,
                     int *leafcount, 
                     int *comblabels, 
                     int *nodepos) 
  {  
    TNT::Array2D<int> d(*samplesize,*nSNP,data);
    splitter<int> s(d,d.dim1(),d.dim2());
    for (int i=0;i<*npos;i++) s.split(positions[i]);
    *len=s.nleaves();
    
    std::vector<std::pair<int,int> > edges;
    std::vector<std::vector<int> > labels;
    
    s.apesplit(edges, labels);
    
    // get the lengths of the labels to allow us to pass this information back tp R
    size_t count=0;
    for (size_t ii=0;ii<labels.size();ii++) {
      for (size_t jj=0; jj<labels[ii].size(); jj++) 
        comblabels[count++] = labels[ii][jj]+1;  // get 1 offset not 0
      leafcount[ii] = static_cast<int>(labels[ii].size());
    }
    
    int nedges=static_cast<int>(edges.size());
    if (nedges!=2*(*len-1))
      throw std::range_error("problem in split_simple\n");
    
    for (int i=0;i<nedges;i++) {
      edge[i] = edges[i].first;
      edge[nedges+i] = edges[i].second;
    }
    // now try to get the lengths.  Note that the centre of this split is 
    // positions[0].
    std::vector<int> nodePos;
    s.getNodesPositions(nodePos);
    for (size_t ii=0;ii<nodePos.size();ii++) nodepos[ii]=nodePos[ii];
  }
  
  
  
  
 
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
