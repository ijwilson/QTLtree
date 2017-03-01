#include "FSData.h"


SNP_Type getsnptype(const std::string &a) {
  if (a=="Nonsynonymous") return Nonsynonymous;
  else if (a=="Synonymous") return Synonymous;
  else if (a=="Unknown") return unknown;
  else if (a=="None") return unknown;
  else if (a==""||a=="NA") return unknown;
  else {
    std::cerr << "warning, type " << a << " not known\n";
  }
 return unknown;
}
/***************************************************************************/
snp_info::snp_info(const char *filename) {
    SimpleCSVReader snp(filename);  // read from snp file
    for (size_t ii=0;ii<snp.nrow();ii++) {
      snp_name.push_back(snp[ii][0]);
      chr.push_back(atoi(snp[ii][1].c_str()));
      position.push_back(atol(snp[ii][2].c_str()));
      gene.push_back( snp[ii][3]);
      snp_type.push_back( getsnptype(snp[ii][4]));
      variant.push_back( snp[ii][5][0]);
      maf.push_back(atof(snp[ii][5].c_str()));
    }
  }
/***************************************************************************/
ind_info::ind_info(const char *filename) {
    SimpleCSVReader idcsv(filename);  // read from snp file
    for (size_t ii=0;ii<idcsv.nrow();ii++) {
      id.push_back(idcsv[ii][0]);
      sex.push_back(atoi(idcsv[ii][1].c_str()));
      age.push_back(atoi(idcsv[ii][2].c_str()));
      population.push_back( idcsv[ii][3]);
    }
  }
void ind_info::add_pheno(const char *filename)
{
  SimpleCSVReader pheno(filename);  // read from pheno file
  size_t len = pheno[0].size();
  size_t nQTL = len-1;
  std::cerr <<"Reading " << nQTL << " quantitative traits from file " << filename << std::endl;

  for (size_t ii=0;ii<pheno.nrow();ii++) {
    if (pheno[ii][0] != id[ii]) {  // this should not have to be the case but is here
      std::cerr << "Error, names in phenotype and ped file are not in the same order\n";
      exit(EXIT_FAILURE);
    }
    std::vector<double> q(nQTL);
    for (size_t jj=0;jj<nQTL;jj++) {
      assert(pheno[ii].size()==len);
      q[jj] = atof(pheno[ii][1+jj].c_str());
    }
    phenotype.push_back(q);
  }
}
/*************************************************************************************************************
 *  FSData                                    
 *************************************************************************************************************/
std::vector<int> FSData::nvar(int chromosome,int posmin, int posmax,int &count) 
 {
    std::vector<int> x(n(),0);
    count=0;
    for (size_t ii=0;ii<nsnps(chromosome);ii++) {
      //   std::cerr << gdata[chromosome].dim1() << " " << gdata[chromosome].dim2() << std::endl;
      if (posmin<=positions[chromosome][ii] and positions[chromosome][ii]<=posmax) {
        count++;
        for (size_t jj=0;jj<n();jj++) if (gdata[chromosome][jj][ii]!=0) x[jj]++;
      }
    }
    return x;
  }
/***************************************************************************/
void conv(const SimpleCSVReader &csv, const snp_info &snp
                       , const ind_info &id,TNT::Array2D<int> &a,std::vector<int> &positions) {

  positions.resize(0);
  for   (size_t ii=0;ii<csv.nrow();ii++) {  
    if (id.id[ii] != csv[ii][0]) {
      std::cerr << "wrong order for individuals, expected " 
                << id.id[ii] << " and got " << csv[ii][0]  << " exiting\n";
      exit(EXIT_FAILURE);
    }
  }

  for (size_t jj=1;jj<csv.ncol();jj++) {  // loop over SNPs
    size_t index = snp.findIndex(csv.colnames[jj]);
    if (index ==snp.size()) {
      std::cerr << "snp " << csv.colnames[jj] 
                << " not found in snp info, exiting";
      exit(EXIT_FAILURE);
    }
    positions.push_back(snp.position[index]);
    char variant = snp.variant[index];

    for   (size_t ii=0;ii<csv.nrow();ii++) {  // loop over individuals
      int count=0;
      assert(csv[ii][jj][1]=='/');
      if (csv[ii][jj][0]==variant 
          || csv[ii][jj][2]==variant) {
        a[ii][jj-1]=1;
        count++;
      }
      if (count>static_cast<int>(csv.nrow())/2) {
        std::cerr << "potential problem with snp at position " << ii << std::endl;
      }
    }
  }
} 
/***************************************************************************/
void variantCounts(const SimpleCSVReader &csv, const snp_info &snp
                   , const ind_info &id,TNT::Array2D<int> &a,std::vector<int> &positions) {
  positions.resize(0);
  for   (size_t ii=0;ii<csv.nrow();ii++) {  
    if (id.id[ii] != csv[ii][0]) {
      std::cerr << "wrong order for individuals, expected " 
                << id.id[ii] << " and got " << csv[ii][0] 
                << " exiting\n";
      exit(EXIT_FAILURE);
    }
  }

  for (size_t jj=1;jj<csv.ncol();jj++) {  // loop over SNPs
    size_t index = snp.findIndex(csv.colnames[jj]);
    if (index ==snp.size()) {
      std::cerr << "snp " << csv.colnames[jj] 
                << " not found in snp info, exiting";
      exit(EXIT_FAILURE);
    }
    positions.push_back(snp.position[index]);
    char variant = snp.variant[index];

    for   (size_t ii=0;ii<csv.nrow();ii++) {  // loop over individuals
      assert(csv[ii][jj][1]=='/');
      a[ii][jj-1] = (csv[ii][jj][0]==variant)+(csv[ii][jj][2]==variant);
    }      
  }
} 
/*******************************************************************/
targets::targets(const std::string &filename) {
  SimpleCSVReader csv(filename.c_str());
   for (size_t ii=0;ii<csv.nrow();ii++) {
      name.push_back(csv[ii][0]);
      chromosome.push_back(atoi(csv[ii][1].c_str()));
      position.push_back(atoi(csv[ii][2].c_str()));
      gene_start.push_back( atoi(csv[ii][3].c_str()));
      gene_end.push_back( atoi(csv[ii][4].c_str()));
      gene_length.push_back( atoi(csv[ii][5].c_str()));
      orientation.push_back(csv[ii][6][0]);
      gene_id.push_back(atoi(csv[ii][7].c_str()));
    }
}
