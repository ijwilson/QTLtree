/**  Quick program to extract the number of variants by gene and by individual, 
can also (by adding the correct flag) output a file with all allele counts      */
#include "options.H"
#include "read_csv.H"
#include "FSData.H"
#include "progressBar.H"
#include "newio.H"

#include <iostream>

int main(int argc, char *argv[]) 
{
  std::string targetFile,pedFile,snpFile,basedir;
  int dataset;
 
  const char *version=SVN_REV;
  
  options o("Options used",version);

  try {
    o.add(&basedir,"b","Base Directory","/users/nijw/GAW/");
    o.add(&targetFile,"t","File containing list of targets","Extras/gene_info.4ago2010");
    o.add(&pedFile,"R","ped filename, gives the list of samples and their regions"
          ,"CD/unrelateds.ped");
    o.add(&snpFile,"snp","list of snps","CD/snp_info");
    o.add(&dataset,"d","datasets to use (when have multiple QTL datasets)",1);
    //    o.add(&phenoFile,"pheno","pheno filename.","CD/unr_phen");
    o.readcommandline(argc,argv);
  }  
  catch(std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }
  catch(...) {
    std::cerr << "Exception of unknown type!\n";
  }
  // Now read the data in raw
  pedFile=basedir+pedFile;
  snpFile=basedir+snpFile;
  FSData fsd(snpFile.c_str(),pedFile.c_str());
  //  fsd.id.add_pheno(pheno.str().c_str());
  
  std::cerr << "Adding chromosomes: ";
  int sumsnps=0;
  for (int i=1;i<=22;i++) {
    std::cerr << i << " ";
    std::ostringstream oss;
    oss << basedir << "CD/c" << i << "_snps.unr";
    fsd.add_data(oss.str().c_str(),false);
    std::cerr << "read " << fsd.nsnps(i) << " snps\n";
    sumsnps+=fsd.nsnps(i);
  }
  std::cerr << "Total of " << sumsnps << " read." << std::endl;
  fsd.print(std::cout);

  exit(EXIT_SUCCESS);
}
