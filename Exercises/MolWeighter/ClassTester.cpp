#include <iostream>

#include "MolWeighter.hpp"
#include "Naomini/Molecule.hpp"
#include "Naomini/MoleculeFactory.hpp"

using namespace Naomini;

int main(int argc, char *argv[]) {


  if(argc != 2) {
    std::cerr << "Usage: <Molecule file>" << std::endl;
    return 1;
  }

  MoleculeVector mols = MoleculeFactory::getAllMolecules(std::string(argv[1]));
  MolWeighter weighter(mols);

  if (weighter.getMoleculesWithWeight(Equal, 350).size() == 0) {
    std::cout << "  SUCCESS: Equal weighter works!" << std::endl;
  } else {
    std::cout << "  ERROR: Equal weighter fails!" << std::endl;
  }
  if (weighter.getMoleculesWithWeight(MoreThan, 350).size() == 701) {
    std::cout << "  SUCCESS: Less than weighter works!" << std::endl;
  } else {
    std::cout << "  ERROR: Less than weighter fails!" << std::endl;
  }
  if (weighter.getMoleculesWithWeight(LessThan, 350).size() == 733) {
    std::cout << "  SUCCESS: More than weighter works!" << std::endl;
  } else {
    std::cout << "  ERROR: More than weighter fails!" << std::endl;
  }

  return 0;
}

