#include <iostream>

#include "Naomini/Atom.hpp"
#include "Naomini/Bond.hpp"
#include "Naomini/Molecule.hpp"
#include "Naomini/MoleculeFactory.hpp"
#include "Naomini/MoleculeDrawer.hpp"
#include "Naomini/Helpers.hpp"

namespace Naomini {

BondPtr* findAllDoubleBonds(MoleculePtr mol){

  BondPtr *bondArray;
  bondArray = (BondPtr*) malloc (100 * sizeof(BondPtr));

  int counter = 0;
  for (unsigned i = 0; i < mol->getBonds().size(); i++) {
    if(mol->getBonds().at(i)->getType() == DOUBLE) {
      bondArray[i] = mol->getBonds().at(i);
      counter++;
    }
  }

  // print double bond information
  if (counter == 0){
    std::cout << "No double bonds found." << std::endl;
  } else {
    std::cout << "Double bonds: " << std::endl;
    for (int i = 0; i < counter; i++) {
      std::cout <<  "    " << mol->getBonds().at(i) << std::endl;
    }
  }
  return bondArray;
}

AtomVector findAllAmineNitrogens(MoleculePtr mol){
  return AtomVector();
}

int main(int argc, char *argv[]) {

  if(argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <Molecule file>" << std::endl;
    return 1;
  }

  MoleculeVector inMols = MoleculeFactory::getAllMolecules(std::string(argv[1]));

  for (MoleculePtr molecule : inMols) {
    BondPtr* bondArray = findAllDoubleBonds(molecule);
    free(bondArray);
  }


  return 0;
}

} //namespace
