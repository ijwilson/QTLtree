#ifndef FSDATA_H__
#define FSDATA_H__

#include <cstdio>
#include <iostream>
#include <map>          // for popmap
#include <algorithm>
#include "read_csv.h"
#include "closest.h"
#include "tnt/tnt.h"

class snp_info;
class ind_info;

void conv(const SimpleCSVReader &csv, const snp_info &snp
                       , const ind_info &id,TNT::Array2D<int> &a,std::vector<int> &positions);
void variantCounts(const SimpleCSVReader &csv, const snp_info &snp
                   , const ind_info &id,TNT::Array2D<int> &a,std::vector<int> &positions);
enum SNP_Type  {Nonsynonymous, Synonymous, unknown};

SNP_Type getsnptype(const std::string &a);
//**************************************************
//    Keeps SNP information
//***************************************************
class snp_info {
public:
  snp_info(const char *filename);
  size_t findIndex(const std::string &name) const {
    return std::distance(snp_name.begin(),
                         std::find(
                                   snp_name.begin(),
                                   snp_name.end(),
                                   name)
                         );
  }
  size_t size() const {
    return snp_name.size();
  }
  std::vector<std::string> snp_name;      
  std::vector<int> chr;
  std::vector<unsigned long> position;
  std::vector<std::string> gene;
  std::vector<SNP_Type> snp_type;
  std::vector<char> variant;
  std::vector<double> maf;
};
//**************************************************
//    Keeps individuals information
//**************************************************
class ind_info {
public:
  ind_info(const char *filename);
  size_t findIndex(const std::string &name) const {
    return std::distance(id.begin(),
                         std::find(
                                   id.begin(),
                                   id.end(),
                                   name)
                         );
  }
  size_t size() const {
    return id.size();
  }
  void makepopmap(int poplength) {
    int count=0;
    for (size_t i=0;i<size();i++) {
      std::string sp=population[i].substr(0,poplength);
      if (popmap.find(sp)==popmap.end()) popmap[sp]=count++;
      popcode.push_back(popmap[sp]);
    }
  }
  std::ostream &printpopmap(std::ostream &o) {
    o << "Population Labels\n-----------------\n";
    std::map<std::string,int>::iterator i=popmap.begin();
    while (i!=popmap.end()) {
      o << "    " << i->first << ": " << i->second << std::endl;
      i++;
    }
    return o;
  }
  void add_pheno(const char *filename);
  size_t nregion() const {return popmap.size();}
  std::vector<double> qtl(int col) const {
    std::vector<double> a(size());
    for (size_t i=0;i<size();i++) a[i]=phenotype[i][col];
    return a;
  }
  std::vector<std::string> id;
  std::vector<int> sex;
  //std::vector<int> smoke;
  //std::vector<int> Affected;
  std::vector<std::vector<double> > phenotype;
  std::vector<int> age;
  std::vector<std::string> population;   // the population name for a sample
  std::vector<int> popcode;              // the populattion code: 0 offset
  std::map<std::string,int> popmap;      // map of population name to popcode
};
//****************************************************
//   Class to hold all the information
//****************************************************
class FSData {
public:
  FSData(const char *snpfilename,const char *idfilename):snps(snpfilename),id(idfilename) {} 
  void add_data(const char *filename,bool anyvariant=true) {  // add genotype data
    SimpleCSVReader geno(filename);
    size_t index = snps.findIndex(geno.colnames[1]);
    //std::cerr << geno.colnames[1] << " " << index << std::endl;
    assert(snps.chr[index]>=1 && snps.chr[index] <=22);
    positions[snps.chr[index]] = std::vector<int>(0);
    gdata[snps.chr[index]] = TNT::Array2D<int>(geno.nrow(),geno.ncol()-1,0);// all zeros to start
    if (anyvariant) 
      conv(geno,snps,id, gdata[snps.chr[index]],positions[snps.chr[index]]);
    else 
      variantCounts(geno,snps,id, gdata[snps.chr[index]],positions[snps.chr[index]]);

  }
  size_t n() const {return id.size();};   // number of individuals
  int **data(int chromosome)  {
    return gdata[chromosome];
  }
  size_t nsnps(int chromosome)  {
    return gdata[chromosome].dim2();
  }
  bool contains(int chr) const {
    if (gdata.find(chr)==gdata.end()) return false;
    return true;
  }
  std::ostream &print(std::ostream &o) const {
    std::map<int, TNT::Array2D<int> >::const_iterator dat=gdata.begin();
    while (dat!=gdata.end()) {
      for (int ii=0;ii<dat->second.dim2();ii++) {
        // o << dat->first << " " << positions[dat->first].at(ii); 
        for (int jj=0;jj<dat->second.dim1();jj++) {
          o <<  dat->second[jj][ii] << " ";
        }
        o << std::endl;
      }
      dat++;
    }
    return o;
  }
  
  std::vector<int> nvar(int chromosome,int posmin, int posmax,int &count);
  ///////////////////////////////////////////////////  
  // DATA
  ///////////////////////////////////////////////////
  snp_info snps;
  ind_info id;
  std::map<int, TNT::Array2D<int> > gdata;    // map of genetic data for chromosome
  std::map<int,std::vector<int> > positions;  // map of positions for chromosome
};


//****************************************************
// A Class to hold the target regions for the 
// QTLtree analysis
//****************************************************
class targets {
public:
  targets(const std::string &filename);
  size_t size() const {
    return name.size();
  }
  ///////////////////////////////////////////////////  
  // DATA
  ///////////////////////////////////////////////////
  std::vector<std::string> name;
  std::vector<char> orientation;
  std::vector<int> chromosome,gene_id,gene_length;
  std::vector<int> position,gene_start,gene_end;
};

#endif
