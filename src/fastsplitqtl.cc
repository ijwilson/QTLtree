/** @file          *
 * Fastsplit - split a dataset into bifurcating trees and then test the significance of
 * the branch splits when compared to cases and controls  
 */
 
#include "splitter.h"
#include "BeagleData.h"
#include "tnt/tnt.h"
#include "options.h"
#include "newio.h"
#include "gsl_rand.h"
#include "closest.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>   // for lower_bound and upper_bound


int main(int argc, char *argv[]) 
{
  int starting,nSplits,reps,maxk,seed,qtlpos;
  std::string infile,positionFile,regionFile,excludeFile,stem,QTLfile;
  double FirstPos,LastPos,maxdistance;
  bool removeCentre=false,loud=false,regional=false;
  char direction;

  // get the version number from GIT
  const char *version=GIT_REV;
  options o("Options used",version);
 
  try {
    o.add(&stem,"input","Stem of input filenames  (.phased, .position and .trait added to create file names)\n"
          "\t\t\t if this is left blank then the individual file names are used","");
    o.add(&infile,"f","Beagle input filename","beagleOut.phased");
    o.add(&positionFile,"P","positions filename (blank for <infile>.position)","");
    o.add(&QTLfile,"q","Quantitative trait file","");
    o.add(&regionFile,"R","regions filename (blank for <infile>.region)\n"
          "\t\t\tnot needed unless you use regional randomisation","");
    o.add(&starting,"b","first SNP to split",1);
    o.add(&excludeFile,"E","A file that contains SNPs positions to exclude","");
    o.add(&nSplits,"n","Number of SNPs to split",20);
    o.add(&reps,"r","replicates for the randomisation",1000);
    o.add(&maxk,"k","k - the number of disjoint node statistics to collect",10);
    o.add(&seed,"seed","Random Number Seed",1);

    o.add(&FirstPos,"start","The first position to analyse (in MB) - negative for the position of the first SNP",-1.0);
    o.add(&LastPos,"end","The last position to analyse (in MB) - negative for the position of the last SNP",-1.0);
    o.add(&direction,"direction","The direction of the tree, <L>eft, <R>ight or <C>entral",'C');
    o.add(&removeCentre,"removeCentre","Remove the First SNP");
    o.add(&qtlpos,"qtlpos","Which QTL to use (1 offset)",1);
    o.add(&maxdistance,"maxd","Maximum distance to consider (in MB)",10000.0);
    //  o.add(&regional,"regional","Use regional randomisation?");
   
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

  if (stem!="") {
    infile=stem+".phased";
    QTLfile=stem+".trait";
    positionFile=stem+".position";
  }

  std::cout << "# command line: ";
  for (int jj=0;jj<argc;jj++) std::cout << argv[jj] << " ";
  std::cout << "  version " << version << std::endl;

  CCData d(infile);
  d.addPositions(positionFile.c_str()); 
   
  d.addQuantitativeTraits(QTLfile.c_str());
  if (regionFile!="") regional=true;
  if (regional) d.addRegions(regionFile.c_str());

  std::vector<int> ExcludePositions;
  if (excludeFile!="") {
    std::ifstream inEx(excludeFile.c_str());
    if (!inEx) 
      throw std::string("unable to open exclusions file");
    copy( std::istream_iterator<int>(inEx),
          std::istream_iterator<int>(),
          back_inserter(ExcludePositions ) );
  } 

  //get the first and last positions to construct trees at

  int start,final;
  if (FirstPos<0.0) {
    start=0;
  } else {
    std::vector<int>::iterator iii=
      std::lower_bound(d.position.begin(),d.position.end(),FirstPos*1000000.0);
    if (iii==d.position.end()) {
      std::cerr << "No SNPs after first position " << FirstPos << ", Exiting\n";
      exit(EXIT_FAILURE);
    }
    start=std::distance(d.position.begin(),iii); 
  }
  if (LastPos<0.0)  {
    final=d.position.size();
  } else {
    std::vector<int>::iterator iii=
      std::upper_bound(d.position.begin(),d.position.end(),LastPos*1000000.0);
    final=std::distance(d.position.begin(),iii); 
  }

  if (loud)  std::cerr << "Read in " << d.haplotypes.dim2() << " SNPs " 
                       << " from " << d.haplotypes.dim1() << " haplotypes\n";  

  if (qtlpos<1 or qtlpos>d.QuantitativeTraits.dim1()) {
    std::cerr << "the qtl at position " << qtlpos << "does not exist - maximum is " 
              << d.QuantitativeTraits.dim1() << std::endl;
    exit(EXIT_FAILURE);
  }

  // Get the set of randomisations that you are going to use for all trees
  std::vector<double> myqtl(d.QuantitativeTraits[qtlpos-1],d.QuantitativeTraits[qtlpos-1]+d.nsamples); 
  std::vector<std::vector<double> > randqtl(reps,myqtl);

  if (!regional) {
    for (int ii=0;ii<reps;ii++) permute(randqtl[ii],r);
  } else {  // really this should be diploid
    int MaxRegion = *std::max_element(d.region.begin(),d.region.end());
    std::vector<std::vector<int> > RegionalVector(MaxRegion);
    for (size_t i=0;i<d.region.size();i+=2) RegionalVector[d.region[i]-1].push_back(i);
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

  // a little check
  if (regional) {
    int MaxRegion = *std::max_element(d.region.begin(),d.region.end());
    std::vector<double> sm(MaxRegion,0.0);
    for (size_t i=0;i<d.region.size();i++) sm[d.region[i]-1] += d.QuantitativeTraits[qtlpos-1][i];
    for (int ii=0;ii<reps;ii++) {
      std::vector<double> testsm(MaxRegion,0.0);
      for (size_t jj=0;jj<d.region.size();jj++) testsm[d.region[jj]-1] += randqtl[ii][jj];
      for (int jj=0;jj<MaxRegion;jj++) 
        if (fabs(testsm[jj] - sm[jj])>1E-8) std::cerr << "error with randomisation";
    }
  }
  // double sm = std::accumulate(d.QuantitativeTraits[qtlpos-1],d.QuantitativeTraits[qtlpos-1]+d.nsamples,0.0);
  // for (int ii=0;ii<reps;ii++) {
  //   double testsm = std::accumulate(randqtl[ii].begin(),randqtl[ii].end(),0.0);
  //   if (fabs(testsm - sm)>1E-8) std::cerr << "error with randomisation";
  // }

  // for (int ii=0;ii<10;ii++) {
  //   for (int j=0;j<5;j++) std::cerr << randqtl[ii][j] << " ";
  //   std::cerr << std::endl;//yprintvec(randqtl[ii]);
  // }

  for (int i=start;i<final;i++) {
    int minpos=d.position.size(),maxpos=0;
    int splits=0;
    int firstIndex=removeCentre?1:0;
    splitter<int> s(d.haplotypes,d.haplotypes.dim1(),d.haplotypes.dim2());
    if (direction=='C') {
      closest cl(d.position,d.position[i]);
      for (int j=firstIndex;splits<nSplits;j++) {
        if (j>=cl.size()) break;
        if (fabs(d.position[cl[j]] - d.position[i])>maxdistance*1000000.0) break;
        if (std::find(ExcludePositions.begin(),ExcludePositions.end(),d.position[cl[j]])==ExcludePositions.end()) {
          s.split(cl[j]);
          splits++;
          if (cl[j]<minpos) minpos=cl[j];
          if (cl[j]>maxpos) maxpos=cl[j];
        } 
      }
    } else if (direction=='R') {
      for (int j=firstIndex;splits<nSplits;j++) {
        if (static_cast<size_t>(i+j)==d.position.size()) break; // no more positions to split
        if (d.position[i+j]-d.position[i]>maxdistance*1000000.0) break;
        if (std::find(ExcludePositions.begin(),ExcludePositions.end(),d.position[i+j])
            ==ExcludePositions.end()) {
          s.split(i+j);
          splits++;
          if (i+j<minpos) minpos=i+j;
          if (i+j>maxpos) maxpos=i+j;
        } 
      }
    } else if (direction=='L') {
      for (int j=firstIndex;splits<nSplits;j++) {
        if (i-static_cast<int>(j)<0) break; // no more positions to split
        if (d.position[i]-d.position[i-j]>maxdistance*1000000.0) break;
        if (std::find(ExcludePositions.begin(),ExcludePositions.end(),d.position[i-j])
            ==ExcludePositions.end()) {
          s.split(i-j);
          splits++;
          if (i-j<minpos) minpos=i-j;
          if (i-j>maxpos) maxpos=i-j;
        } 
      }
    }
    if (minpos<=maxpos) {  // was anything done

      std::vector<double> stat=s.qtlStat(myqtl,maxk);
      
      TNT::Array2D<double> RandStatistics(maxk,reps);
      for (int replicate=0;replicate<reps;replicate++) {
        std::vector<double> stmp=s.qtlStat(randqtl[replicate],maxk);
        for (int k=0;k<maxk;k++) RandStatistics[k][replicate]=stmp[k];
      }
      TNT::Array1D<double>  pval(maxk,0.0);
      for (int k=0;k<maxk;k++) {
        for (int j=0;j<reps;j++)  if (stat[k]>RandStatistics[k][j]) pval[k] += 1.0;
      }
      std::cout << d.position[i] << " ";
      for (int k=0;k<maxk;k++) std::cout << stat[k] << " " << pval[k] << " ";
      
      std::cout << minpos << " " << d.position[minpos] << " " << d.position[maxpos] 
                << " " << s.nleaves() << std::endl;
    } 
  }
  exit(EXIT_SUCCESS);
}
