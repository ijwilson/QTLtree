/** @file          *
 * Fastsplit - split a dataset into bifurcating trees and then test the significance of
 * the branch splits when compared to cases and controls  
 */
 
#include "splitter.H"
#include "BeagleData.H"
#include "tnt/tnt.h"
#include "options.H"
#include "newio.H"
#include "gsl_rand.H"
#include "closest.H"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>   // for lower_bound and upper_bound
#include "targets.H"

int main(int argc, char *argv[]) 
{
  int reps,maxk,seed,qtlpos,dataset;
  std::string infile,positionFile,regionFile;
 
  bool regional=true;
  char direction,statisticUsed;
  std::string targetfile;

  // get the version number from SVN
  const char *version=SVN_REV;
  options o("Options used",version);
 
  try {
    o.add(&targetfile,"t","File containing list of targets","genes.target");

    o.add(&regionFile,"R","regions filename (blank for <infile>.region)\n"
          "\t\t\tnot needed unless you use regional randomisation","");

    o.add(&reps,"r","replicates for the randomisation",1000);
    o.add(&maxk,"k","k - the number of disjoint node statistics to collect",10);
    o.add(&seed,"seed","Random Number Seed",1);
    o.add(&direction,"direction","The direction of the tree, <L>eft, <R>ight or <C>entral",'C');
    o.add(&qtlpos,"qtlpos","Which QTL to use (1 offset)",1);
    o.add(&dataset,"d","datasets to use",1);
    //  o.add(&regional,"regional","Use regional randomisation?");
    o.add(&statisticUsed,"s","Statistic to use <A>bsolute, <S>quared, <P>ositive or <N>egative",'A');

    o.readcommandline(argc,argv);
  }  
  catch(std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }
  catch(...) {
    std::cerr << "Exception of unknown type!\n";
  }
  rng r(seed);
  targets tar(targetfile);
  std::cerr << "read " << tar.name.size() << " genes to analyse\n";
 regionFile = "/home/nijw/GAW/fastsplit/data/regions.beagle";

 std::vector<CCData *> d(23);
  for (size_t ii=1;ii<=22;ii++) {
    std::ostringstream oss;
    oss << "/home/nijw/GAW/fastsplit/data/c"<<ii<<".phased";
    d[ii] = new CCData(oss.str().c_str());
    oss.str("");
    oss << "/home/nijw/GAW/phased/c"<<ii<<"_beagle.markers";
    d[ii]->addPositions(oss.str().c_str());
    oss.str("");
    oss << "/home/nijw/GAW/fastsplit/data/qtls." << dataset;
    d[ii]->addQuantitativeTraits(oss.str().c_str());
    d[ii]->addRegions(regionFile.c_str());
  }
    

  std::cout << "# command line: ";
  for (int jj=0;jj<argc;jj++) std::cout << argv[jj] << " ";
  std::cout << "  version " << version << std::endl;
 
  if (qtlpos<1 or qtlpos>d[1]->QuantitativeTraits.dim1()) {
    std::cerr << "the qtl at position " << qtlpos << "does not exist - maximum is " 
              << d[1]->QuantitativeTraits.dim1() << std::endl;
    exit(EXIT_FAILURE);
  }

  // Get the set of randomisations that you are going to use for all trees
  std::vector<double> myqtl(d[1]->QuantitativeTraits[qtlpos-1],d[1]->QuantitativeTraits[qtlpos-1]+d[1]->nsamples); 
  std::vector<std::vector<double> > randqtl(reps,myqtl);

  if (!regional) {
    for (int ii=0;ii<reps;ii++) permute(randqtl[ii],r);
  } else {  // really this should be diploid
    int MaxRegion = *std::max_element(d[1]->region.begin(),d[1]->region.end());
    std::vector<std::vector<int> > RegionalVector(MaxRegion);
    for (size_t i=0;i<d[1]->region.size();i+=2) RegionalVector[d[1]->region[i]-1].push_back(i);
    std::vector<std::vector<int> > RandRegionalVector=RegionalVector;
    for (int ii=0;ii<reps;ii++) {
      for (size_t jj=0;jj<RegionalVector.size();jj++)  {
        permute(RandRegionalVector[jj],r);
        for (size_t kk=0;kk<RegionalVector[jj].size();kk++) {
          randqtl[ii][RegionalVector[jj][kk]] = myqtl[RandRegionalVector[jj][kk]];
          randqtl[ii][RegionalVector[jj][kk]+1] = myqtl[RandRegionalVector[jj][kk]+1];
        }
      }
    }
  }

  for (size_t ii=0;ii<tar.name.size();ii++) {
    int minpos=100000;
    int maxpos=-1;
    int chr = tar.chromosome[ii];
    int startpos=tar.position[ii];
    // find positions to split for this gene
    splitter<int> s(d[chr]->haplotypes,d[chr]->haplotypes.dim1(),d[chr]->haplotypes.dim2());
    int count=0;
    if (direction=='C') {
      bool up=false;
      bool down=false;
     
      closest cl(d[chr]->position,startpos);
      // find positions to use
      int splits=0;
      for (int j=0;;j++) {
        if (j>=cl.size()) break;
        if (d[chr]->position[cl[j]]<tar.gene_start[ii]) {
          down=true;
          continue;
        }  else if (d[chr]->position[cl[j]]>tar.gene_end[ii]) {
          up=true;
          continue;
        }
        if (up and down) break;
        count++;
        s.split(cl[j]);
        splits++;
        if (cl[j]<minpos) minpos=cl[j];
        if (cl[j]>maxpos) maxpos=cl[j];
      } 
    }
  
 else if (direction=='R') {
    std::cerr << "note defined";
  } else if (direction=='L') {
    std::cerr << "note defined";
  }
    if (minpos<=maxpos) {  // was anything done
      // if (nleaves<maxk) maxk=nleaves;
      std::vector<double> stat=s.qtlStat(myqtl,maxk,statisticUsed);
      
      TNT::Array2D<double> RandStatistics(maxk,reps);
      for (int replicate=0;replicate<reps;replicate++) {
        std::vector<double> stmp=s.qtlStat(randqtl[replicate],maxk,statisticUsed);
        for (int k=0;k<maxk;k++) RandStatistics[k][replicate]=stmp[k];
      }
      TNT::Array1D<double>  pval(maxk,0.0);
      for (int k=0;k<maxk;k++) {
        for (int j=0;j<reps;j++)  if (stat[k]>RandStatistics[k][j]) pval[k] += 1.0;
      }
      std::cout << tar.name[ii] << " " << tar.position[ii] << " " << chr << " "  << count << " ";
      for (int k=0;k<maxk;k++) std::cout << stat[k] << " " << pval[k] << " ";
      
      std::cout << minpos << " " << d[chr]->position[minpos] << " " << d[chr]->position[maxpos] 
                << " " << s.nleaves() << std::endl;

    } 
  }
  exit(EXIT_SUCCESS);
}
