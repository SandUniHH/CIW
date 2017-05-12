#include <iomanip>

#include "boost/numeric/ublas/matrix.hpp"

#include "Naomini/MoleculeFactory.hpp"
#include "Naomini/Helpers.hpp"
#include "Naomini/MoleculeDrawer.hpp"
#include "Naomini/Molecule.hpp"

#include "BronKerbosch.hpp"

namespace Naomini {

void computeMCSForAllMolecules(MoleculeVector &mols) {

  for (unsigned i = 0; i < mols.size(); i++) {
    for (unsigned j = i+1; j < mols.size(); j++)
    {
      std::cout << "Molecules: " << mols.at(i) << " " << mols.at(j) << " : \n";
      calculateBronKerbosch(mols.at(i), mols.at(j));
    }
  }
}

}

int main(int argc, char *argv[])
{
  // check the number of arguments
  if (argc != 2){
    std::cout << "Usage: " << argv[0] << " <molecule-file>" << std::endl;
    return 1;
  }

  // get a molecule for each molecule entry in the provided file
  Naomini::MoleculeVector mols = Naomini::MoleculeFactory::getAllMolecules(std::string(argv[1]));
  if ( mols.size() < 1){
    std::cout << "No molecules could be read from \"" << argv[1] <<"\"."  <<   std::endl
        << "Usage: " << argv[0] << " <molecule-file>" << std::endl;
    return 1;
  }

  Naomini::computeMCSForAllMolecules(mols);
  return 0;
}

