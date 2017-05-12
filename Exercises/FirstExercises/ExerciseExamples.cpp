#include <iostream>

#include "Naomini/Atom.hpp"
#include "Naomini/Bond.hpp"
#include "Naomini/Molecule.hpp"
#include "Naomini/MoleculeFactory.hpp"
#include "Naomini/MoleculeDrawer.hpp"
#include "Naomini/Helpers.hpp"


namespace Naomini{

BondVector findAllDoubleBonds(MoleculePtr mol){
  BondVector doubleBonds;
  // insert your code here!
  return doubleBonds;
}

AtomVector findAllAmineNitrogens(MoleculePtr mol){
  AtomVector amineNitrogens;
  // insert your code here!
  return amineNitrogens;
}

int main(int argc, char *argv[]) {
  // check number of arguments
  if (argc != 2){
    std::cout << "Usage: " << argv[0] << " <molecule-file> " << std::endl;
    return 1;
  }

  // get a molecule for each molecule entry in the provided file
  MoleculeVector mols = MoleculeFactory::getAllMolecules(std::string(argv[1]));

  // insert your code here!

  // delete all molecules
  for (MoleculePtr molecule : mols) {
    delete molecule;
  }
  return 0;
} // end main

} // end namespace Naomini






