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
    TNT::Array2D<int> d(*samplesize, *nSNP, data);     // convert the data into a 2-d tnt array
    splitter<int> s(d, d.dim1(), d.dim2());            // define the splitter object s
    for (int i=0;i<*npos;i++) s.split(positions[i]);   // split at positions
    *len=s.nleaves();                                  // how many leaves do we have on the tree
    
    std::vector<std::pair<int,int> > edges;            // set up date structure for edges
    std::vector<std::vector<int> > labels;             // and for labels  
     
    s.apesplit(edges, labels);
    
    // get the lengths of the labels to allow us to pass this information back to R
    // the information is returned in the correct order for ape which is root->left->right
    size_t count=0;
    for (size_t ii=0; ii < labels.size(); ii++) {
      for (size_t jj=0; jj < labels[ii].size(); jj++) 
        comblabels[count++] = labels[ii][jj]+1;  // get 1 offset not 0
      leafcount[ii] = static_cast<int>(labels[ii].size());
    }
    
    int nedges=static_cast<int>(edges.size());
    if (nedges != 2*(*len-1))
      throw std::range_error("problem in C++ code split_simple\n");
    
    for (int i=0; i<nedges; i++) {
      edge[i] = edges[i].first;
      edge[nedges+i] = edges[i].second;
    }
    // now try to get the lengths.  Note that the centre of this split is positions[0].
    std::vector<int> nodePos;
    s.getNodesPositions(nodePos);
    for (size_t ii=0;ii<nodePos.size();ii++) nodepos[ii]=nodePos[ii];
}
  
  
  
  
 
 /** Get a split that splits everything (if no information then the 
  * splits are at random                                               */
 void bsplit( int *data, int *samplesize, int *nSNP, 
                   int *xxx, int *nleaves) 
 {  
   TNT::Array2D<int> d(*samplesize, *nSNP, data);
   splitter<int> s(d, *samplesize, *nSNP);
   for (int i=0;i<*nSNP;i++) s.split(i);
   
   *nleaves=s.nleaves();
   std::cerr << "we have " << *nleaves << " leaves" << std::endl;
   int nedges = 2*(*nleaves-1);
   
   std::cerr << "we have " << nedges << " edges" << std::endl;
   
   s.edges_positions_counts(xxx, xxx+nedges, xxx+2*nedges, xxx+3*nedges);

 }
  
}
