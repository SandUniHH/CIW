#include <iostream>

#include "Naomini/Molecule.hpp"
#include "Naomini/Atom.hpp"
#include "Naomini/MoleculeFactory.hpp"

using namespace Naomini;

int main(int argc, char *argv[]) {

  if(argc != 2) {
    // insert error message
    return 1;
  }

  // loads the first molecule of a given file
  MoleculePtr mol = MoleculeFactory::getFirstMolecule(std::string(argv[1]));

  // insert your code here!
  std::cout << "Name: " << mol->getName() << std::endl;

  for (unsigned i = 0; i < mol->getNofAtoms(); i++ ) {
    std::cout << "Name: " << mol->getAtoms().at(i)->getName();
    std::cout << std::endl;
  }

  return 0;
}



