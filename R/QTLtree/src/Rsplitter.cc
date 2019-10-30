#include "splitter.h"
#include "gsl_rand.h"
#include "tnt/tnt.h"
#include <R.h>
#include <Rmath.h>
#include <iostream>
//#include <sstream>
//#include <iterator>
//#include <algorithm>

extern "C" {
  void splittest(int *data, int *samplesize, int *nSNP, int *positions, int *npos, int *cases
		 , int *ncases, int *reps, double *teststat, double *pval,int *nterm) {
    *pval=0.0;
    TNT::Array2D<int> d(*samplesize,*nSNP,data);
    splitter<int> s(d,d.dim1(),d.dim2());
    for (int i=0;i<*npos;i++) s.split(positions[i]);
    
    int TestS=s.testStat1(cases,*ncases);
    
    rng r;
    std::vector<int> cc(*samplesize);
    for (int i=0;i<*samplesize;i++) cc[i]=i;
    
    for (int i=0;i<*reps;i++) {
      permute(cc,r);
      int stat = s.testStat1(&cc[0],*ncases);
      if (stat<=TestS) *pval+=1.0; 
    }
    *teststat=TestS;
    *pval /= static_cast<double>(*reps);
    *nterm=s.nleaves();
  }  
  /**
   *
   *  Split something based on QTL values at the nodes of a tree
   *
   */
  void splittestQTL( int *data, int *samplesize, int *nSNP, int *positions, int *npos,double *qtl,
                     int *edge
                     ,int *len
                     ,int *leafcount
                     ,int *comblabels
                     ,int *nodepos
                     ,int *tippos
                     ,double *teststat
                     ,int *reps
                     , double *pval,int *nterm)
  {
    		 
    *pval=0.0;
    TNT::Array2D<int> d(*samplesize,*nSNP,data);
    splitter<int> s(d,d.dim1(),d.dim2());
    for (int i=0;i<*npos;i++) s.split(positions[i]);
        
    rng r;
    
    for (int i=0;i<*reps;i++) {
      //      permute(cc,r);
      //int stat = s.testStat1(&cc[0],*ncases);
      double stat=1.0;
      if (stat<= *teststat) *pval+=1.0; 
    }
  
    *pval /= static_cast<double>(*reps);
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
   * labels:        The labels at the leaves (in order
  */
  void GetSplit( int *data, int *samplesize, int *nSNP, int *positions, int *npos, int *cases
		   , int *ncases, int *edge, int *ccleaf, int *ccnode,int *len,int *comblabels
		 , int *nodepos, int *tippos) 
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
      std::cerr << "problem in GetSplit\n";

    for (int i=0;i<nedges;i++) {
      edge[i]=edges[i].first;
      edge[nedges+i]=edges[i].second;
    }
    // extract cases and controls for the leaves
    int count=0;
    for (int i=0;i<*len;i++) {
      ccleaf[i]=count_intersection(labels[i],cases,*ncases);
      ccleaf[*len+i] = labels[i].size()-ccleaf[i];
      for (size_t j=0;j<labels[i].size();j++) comblabels[count++] = labels[i][j];
    }   
    s.getCaseControlNodes(ccnode,cases,*ncases,nedges/2);
    // now try to get the lengths.  Note that the centre of this split is 
    // positions[0].
    s.getLengths(positions[0],nodepos,tippos);
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
  void GetQTLSplit( int *data, int *samplesize, int *nSNP, int *positions, int *npos,double *qtl,
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
    if (nedges!=2*(*len-1))
      std::cerr << "problem in GetSplit\n";
   
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
        comblabels[count++] = labels[ii][jj];
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
			 ,int *ncategories
                    ) 
  {  

    std::cerr << "Start of function" << std::endl;
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
      std::cerr << "problem in GetSplit\n";

    for (int i=0;i<nedges;i++) {
      edge[i]=edges[i].first;
      edge[nedges+i]=edges[i].second;
    }
    std::cerr << "geting lenghts" << std::endl;
    s.getLengths(positions[0],nodepos,tippos);
   
    for (size_t ii=0;ii<labels.size();ii++) {
      leafcount[ii] = static_cast<int>(labels[ii].size());
      int *counts = leafcat+*ncategories*ii;
      std::cerr << "=========== leaf " << ii << " ================\n";
      for (int jj=0;jj<leafcount[ii];jj++) 
	counts[category[labels[ii][jj]-1]-1]++;
    }
    std::cerr << "getting node categories" << std::endl;

    std::vector<std::vector<int> > nodelabels;
    s.getNodesLabels(nodelabels);
    for (size_t ii=0;ii<nodelabels.size();ii++) {
      int *counts = nodecat+*ncategories*ii;
      for (int jj=0;jj<nodelabels[ii].size();jj++) {
	//std::cerr << nodelabels[ii][jj] << std::endl;
	counts[category[nodelabels[ii][jj]-1]-1]++;
      }
    }
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
      std::cerr << "problem in C++ code simplesplit\n";

    for (size_t j=0;j<labels.size();j++) comblabels[j] = labels[j][0];

    for (int i=0;i<nedges;i++) {
      edge[i]=edges[i].first;
      edge[nedges+i]=edges[i].second;
    }
  }
}
