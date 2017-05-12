#include <iostream>

#include "Naomini/Atom.hpp"
#include "Naomini/Bond.hpp"
#include "Naomini/Molecule.hpp"
#include "Naomini/MoleculeFactory.hpp"
#include "Naomini/MoleculeDrawer.hpp"
#include "Naomini/Helpers.hpp"

#include "Ullmann.hpp"

namespace Naomini{

int main(int argc, char *argv[]) {
  if (argc != 2){
    std::cout << "Usage: " << argv[0]
              << " <molecule-file> " << std::endl;
    return 1;
  }
  // get a molecule for each molecule entry in the provided file
  MoleculeVector mols = MoleculeFactory::getAllMolecules(std::string(argv[1]));
  Ullmann ullmannAlgorithm(mols.at(0), mols.at(1));
  ullmannAlgorithm.run();

  std::cout << "Substructure matching finished." << std::endl;
  return 0;
}
}
