#include <iostream>

#include "Naomini/Forward.hpp"
#include "Naomini/Molecule.hpp"
#include "Naomini/MoleculeFactory.hpp"
#include "Naomini/Helpers.hpp"

#include "ExtendedConnectivityFingerPrint.hpp"
#include "SimilarityMeasure.hpp"

namespace Naomini{
void analyseMoleculeSimilarity(MoleculeVector molecules,
    unsigned nofIterations){


/**** calculate an ECFP for each molecule *************************************/
  std::vector<ECFP> fingerprints;
  try{
    for (MoleculePtr mol : molecules) {
      ECFP fp = getECFP(mol, nofIterations);
      fingerprints.push_back(fp);
    }
  }
  catch (const char *err){
    std::cerr << err << std::endl;
  }

/****  calculate the tanimoto coefficient for each pair of molecules  *********/

/****  calculate the cosine coefficient for each pair of molecules  ***********/

/****  calculate the hamming distance for each pair of molecules  *************/

/****  analyze similarities / distances ***************************************/

}
} // end namespace Naomini

using namespace Naomini;
int main(int argc, char *argv[]) {
  // check the number of arguments
  if (argc != 3){
    std::cout << "Usage: " << argv[0] << " <molecule-file> <number of iterations>" << std::endl;
    return 1;
  }

  // read the number of ecfp-iterations from argument 2
  unsigned nofIterations = atoi(argv[2]);
  if (nofIterations > 10){
    std::cout << "Use a sensible number of iterations." << std::endl
        << "Usage: " << argv[0] << " <molecule-file> <number of iterations>" << std::endl;
    return 1;
  }

  // get a molecule for each molecule entry in the provided file
  MoleculeVector mols = MoleculeFactory::getAllMolecules(std::string(argv[1]));
  if ( mols.size() < 1){
    std::cout << "No molecules could be read from \"" << argv[1] <<"\"."  <<   std::endl
        << "Usage: " << argv[0] << " <molecule-file> <number of iterations>" << std::endl;
    return 1;
  }

  analyseMoleculeSimilarity(mols, nofIterations);

  // delete all molecules
  for (MoleculePtr molecule : mols) {
    delete molecule;
  }
  return 0;
}


