#ifndef FSDATA_H__
#define FSDATA_H__


class FSData {
public:

  FSData(const SimpleCSVReader &snps);    // first thing to load for checks

  std::vector<TNT::Array2D<int> > snps; // vector of chromosomes
  std::vector<std::vector<std::string> > snpid;
  std::vector<std::string> id;
  std::vector<std::vector<double> > QTL;  // column major, so QTL[0] is the first QTL
  std::vector<std::string> region;
};

TNT::Array2D<int> InputData(const std::vector<std::vector<std::string> > &xx);
