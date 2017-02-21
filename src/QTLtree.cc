/** Program to read GAW data, and GAW phenotype files 
 * and to do my sets of analyses on them    */
 //./QTLtree --b=../python/ --snp=markers.csv --R=pops.ped --pheno=qtl.csv --t=targets.csv --qtlpos=1
#include "options.H"
#include "read_csv.H"
#include "FSData.H"
#include "gsl_rand.H"
#include "progressBar.H"
#include "splitter.H"

#include <iostream>

int main(int argc, char *argv[]) 
{
  std::string targetFile,pedFile,snpFile,phenoFile,basedir,genotypeFile;
  bool keepRandomisations; 
  int reps,maxk,seed,qtlpos,abbrevLength;//,dataset
  char direction,statisticUsed; 
  const char *version=SVN_REV;
  
  options o("Options used",version);

  try {
    o.add(&basedir,"b","Base Directory - the directory the data files","/users/nijw/GAW/CD/");
    o.add(&targetFile,"t","File containing list of targets","../Extras/gene_info.4ago2010");
    o.add(&pedFile,"R","ped filename, gives the list of samples and their regions"
          ,"unrelateds.ped");
    o.add(&snpFile,"snp","snp info filename, gives the snp codes"
          ,"snp_info");
    o.add(&genotypeFile,"g","Genotype files, replace * with chromosome numbers","c*_snps.unr");
    o.add(&phenoFile,"pheno","pheno filename.","unr_phen.1");
    o.add(&reps,"r","replicates for the randomisation",1000);
    o.add(&keepRandomisations,"keep","Keep all the randomisations - useful to look at the distribution of "
          "randomisation statistics but this produces lots of data, use very carefully");
    o.add(&maxk,"k","k - the number of disjoint node statistics to collect",10);
    o.add(&seed,"seed","Random Number Seed",1);
    o.add(&direction,"direction","The direction of the tree, <L>eft, <R>ight or <C>entral",'C');
    o.add(&qtlpos,"qtlpos","Which QTL to use (1 offset but ignore label)",4);
    o.add(&statisticUsed,"s","Statistic to use <A>bsolute, <Z>squared, <P>ositive or <N>egative",'A');
    o.add(&abbrevLength,"al","Length of population abbreviation",4);
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
 
  ostringstream pheno;
  pheno << basedir << phenoFile;
  pedFile=basedir+pedFile;
  snpFile=basedir+snpFile;
  FSData fsd(snpFile.c_str(),pedFile.c_str());
  fsd.id.add_pheno(pheno.str().c_str());

  std::cerr << "Adding chromosomes: ";
  size_t replacePosition=genotypeFile.find('*');
  for (int i=1;i<=22;i++) {
    std::cerr << i << " ";
    std::ostringstream oss;
    oss << basedir << genotypeFile.substr(0,replacePosition) << i 
    << genotypeFile.substr(replacePosition+1);
    fsd.add_data(oss.str().c_str());
  }
    
  std::cerr << std::endl;
  fsd.id.makepopmap(abbrevLength);
  fsd.id.printpopmap(std::cerr);

  targets tar(basedir+targetFile);
  std::cerr << tar.size() << " target genes read from file " << basedir+targetFile << std::endl;
  rng r(seed);
  // now get the different randomisations
  
  // Get the set of randomisations that you are going to use for all trees
  std::vector<double> myqtl = fsd.id.qtl(qtlpos-1);
  std::vector<std::vector<double> > randqtl(reps,myqtl);
  progressBar pbperm(std::cerr,"Getting permutations",reps);  
  size_t MaxRegion = fsd.id.nregion();
  std::vector<std::vector<int> > RegionalVector(MaxRegion);
  for (size_t i=0;i<fsd.n();i++) RegionalVector[fsd.id.popcode[i]].push_back(i);
  std::vector<std::vector<int> > RandRegionalVector=RegionalVector;
 
  for (int ii=0;ii<reps;ii++) {
    pbperm.update(ii);
    for (size_t jj=0;jj<RegionalVector.size();jj++)  {
      permute(RandRegionalVector[jj],r);
      for (size_t kk=0;kk<RegionalVector[jj].size();kk++) {
        randqtl[ii][RegionalVector[jj][kk]] = myqtl[RandRegionalVector[jj][kk]];
      }
    }
  }
  pbperm.finish();
  /*
  std::ofstream qcheck("qtl.check");
  printvector(qcheck,myqtl,"\n");
  for (int ii=0;ii<reps;ii++) printvector(qcheck,randqtl[ii],'\n');
  qcheck.close();
  */
 
  for (size_t ii=0;ii<tar.size();ii++) {
    int chr = tar.chromosome[ii];
    // find positions to split for this gene
    if (fsd.contains(chr)) {   
      int minpos=100000;
      int maxpos=-1;
      int startpos=tar.position[ii];

      splitter<int> s(fsd.data(chr),fsd.n(),fsd.nsnps(chr));
      //    std::cerr << "creating tree for gene " <<  tar.name[ii] <<
      //     " with " << fsd.nsnps(chr) << " snps for " << fsd.n() << " individuals" << std::endl;
        
      int count=0;
      if (direction=='C') {
        bool up=false;
        bool down=false;
        
        closest cl(fsd.positions[chr],startpos);
        //    std::cerr << "positions are " << cl.print(std::cerr);
        // find positions to use
        int splits=0;
        for (int j=0;;j++) {
          if (j>=cl.size()) break;
          if (fsd.positions[chr][cl[j]]<tar.gene_start[ii]) {
            down=true;
            continue;
          }  else if (fsd.positions[chr][cl[j]]>tar.gene_end[ii]) {
            up=true;
            continue;
          }
          if (up and down) break;
          count++;
          //      std::cerr << "splitting at position " << cl[j] << std::endl;
          s.split(cl[j]);
          splits++;
          if (cl[j]<minpos) minpos=cl[j];
          if (cl[j]>maxpos) maxpos=cl[j];
        } 
      }
      
      else if (direction=='R') {
        throw ijwerror("No trees defined for direction ",direction);
      } else if (direction=='L') {
          throw ijwerror("No trees defined for direction ",direction);
      } else {
         throw ijwerror("No trees defined for direction ",direction);
      }
      if (minpos<=maxpos) {  // was anything done
        std::vector<double> stat=s.qtlStat(myqtl,maxk,statisticUsed);
        std::vector<std::vector<double> > RandStatistics;
      
        for (int replicate=0;replicate<reps;replicate++) {
          RandStatistics.push_back(s.qtlStat(randqtl[replicate],maxk,statisticUsed));
        }

        if (keepRandomisations) {      // keep the randomised data
          ostringstream oss;
          oss <<  tar.name[ii] <<".rand";
          std::ofstream o(oss.str().c_str());
          std::copy(stat.begin(),stat.end(),std::ostream_iterator<double>(o," "));
          o << std::endl;
          for (int kkk=0;kkk<reps;kkk++) {
            std::copy(RandStatistics[kkk].begin(),RandStatistics[kkk].end()
                      ,std::ostream_iterator<double>(o," "));
            o << std::endl;
          }
          o.close();
        }

        TNT::Array1D<double>  pval(maxk,0.0);
        for (int j=0;j<reps;j++)  {
          for (int k=0;k<maxk;k++) {
            if (stat[k]>RandStatistics[j][k]) pval[k] += 1.0;
          }
        }
        std::cout << tar.name[ii] << " " << tar.position[ii] << " " << chr << " "  << count << " ";
        for (int k=0;k<maxk;k++) std::cout << stat[k] << " " << pval[k] << " ";
        std::cout << minpos << " " << fsd.positions[chr][minpos]  << " " 
                  << fsd.positions[chr][maxpos] << " " << s.nleaves() << std::endl;
        
      } 
    }
  }
  exit(EXIT_SUCCESS);
}
